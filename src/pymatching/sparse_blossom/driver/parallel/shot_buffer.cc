// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/driver/parallel/shot_buffer.h"

#include <iostream>
#include <utility>

#include <omp.h>

#ifdef USE_SHMEM
#include <shmem.h>
#endif

namespace pm {

ShotContainer::ShotContainer(
    int num_partitions, int num_virtual_boundaries, int num_observables_in)
    : partition_hits(num_partitions),
      virtual_boundary_hits(num_virtual_boundaries),
      num_observables(num_observables_in),
      res(num_observables_in)
{
    // Reservation must be a genuine worst-case upper bound, not a heuristic: the claim side reads
    // extraction_jobs[i] without a lock, which is only safe if push_back never reallocates. A
    // per-L formula (e.g. ceil(num_partitions/L)) is NOT safe here -- with small L in non-preemptive/
    // post-hoc-chunking mode, every internal fusion node can end up posting its own standalone job,
    // so the true worst case is bounded by the whole tree's node count, not num_partitions/L. Reuse
    // the same bound tasks.reserve() already uses (2*num_partitions-1), plus slack for the (at most 2,
    // for ROUND) cross-rank jobs. L-independent, safe regardless of extract_preemptively/L's value.
    extraction_jobs.reserve(static_cast<size_t>(num_partitions) * 2 + 1);
}

ShotContainer::ShotContainer(ShotContainer&& other) noexcept
    : sparse_shot(std::move(other.sparse_shot)),
      partition_hits(std::move(other.partition_hits)),
      virtual_boundary_hits(std::move(other.virtual_boundary_hits)),
      num_observables(other.num_observables),
      res(std::move(other.res)),
      num_task_roots(other.num_task_roots),
      thread_results(std::move(other.thread_results)),
      extraction_jobs(std::move(other.extraction_jobs)),
      tasks(std::move(other.tasks)),
      local_seam_tasks(std::move(other.local_seam_tasks))
{
    num_roots_done.store(other.num_roots_done.load());
    extraction_posted_count.store(other.extraction_posted_count.load());
    extraction_claim_cursor.store(other.extraction_claim_cursor.load());
    pending_extraction_jobs.store(other.pending_extraction_jobs.load());
#ifdef USE_SHMEM
    cross_rank_tasks = std::move(other.cross_rank_tasks);
#endif
}

ShotContainer& ShotContainer::operator=(ShotContainer&& other) noexcept {
    if (this != &other) {
        sparse_shot = std::move(other.sparse_shot);
        partition_hits = std::move(other.partition_hits);
        virtual_boundary_hits = std::move(other.virtual_boundary_hits);
        num_observables = other.num_observables;
        res = std::move(other.res);
        tasks = std::move(other.tasks);
        local_seam_tasks = std::move(other.local_seam_tasks);
        num_task_roots = other.num_task_roots;
        num_roots_done.store(other.num_roots_done.load());
        thread_results = std::move(other.thread_results);
        extraction_jobs = std::move(other.extraction_jobs);
        extraction_posted_count.store(other.extraction_posted_count.load());
        extraction_claim_cursor.store(other.extraction_claim_cursor.load());
        pending_extraction_jobs.store(other.pending_extraction_jobs.load());
#ifdef USE_SHMEM
        cross_rank_tasks = std::move(other.cross_rank_tasks);
#endif
    }
    return *this;
}

void ShotContainer::clear() {
    sparse_shot.clear();
    for (auto& hits : partition_hits) {
        hits.clear();
    }
    for (auto& hits : virtual_boundary_hits) {
        hits.clear();
    }
    res.reset();
    // Per-shot queue state, reset here to mirror partition_hits/virtual_boundary_hits/res above --
    // called exactly when this container is about to be reused for the next shot, at which point
    // the caller has already established pending_extraction_jobs == 0.
    extraction_jobs.clear();
    extraction_posted_count.store(0, std::memory_order_relaxed);
    extraction_claim_cursor.store(0, std::memory_order_relaxed);
    pending_extraction_jobs.store(0, std::memory_order_relaxed);
}

void ShotContainer::post_extraction_job(const ExtractionJob& job, std::ofstream* t_out) {
    if (DEBUG && t_out) {
        *t_out << "  POST_JOB " << (job.subtree_root->is_fusion ? "vb=" : "p=") << job.subtree_root->part
               << " task=" << job.subtree_root << std::endl << std::flush;
    }
    // Incremented before the job is published (posted_count's release store below) so no thread
    // can ever observe a claimable job that pending_extraction_jobs hasn't already counted.
    pending_extraction_jobs.fetch_add(1, std::memory_order_acq_rel);
    {
        std::lock_guard<std::mutex> lock(extraction_post_mutex);
        extraction_jobs.push_back(job);
    }
    extraction_posted_count.fetch_add(1, std::memory_order_release);
}

ExtractionJob* ShotContainer::try_claim_extraction_job() {
    size_t claimed = extraction_claim_cursor.load(std::memory_order_relaxed);
    for (;;) {
        size_t posted = extraction_posted_count.load(std::memory_order_acquire);
        if (claimed >= posted) {
            return nullptr;
        }
        if (extraction_claim_cursor.compare_exchange_weak(
                claimed, claimed + 1, std::memory_order_acq_rel, std::memory_order_relaxed)) {
            return &extraction_jobs[claimed];
        }
        // claimed was updated to the current cursor value by the failed CAS; retry.
    }
}

void ShotContainer::reset() {
    current_buffer_round.store(-1, std::memory_order_release);
    for (auto& task : tasks) {
        task.reset();
    }
    for (auto& seam : local_seam_tasks) {
        seam.reset();
    }
#ifdef USE_SHMEM
    // status_shm and signal_shm are reset at end of each shot; only done_shm needs resetting
    for (auto& crt : cross_rank_tasks) {
        crt.reset();
    }
#endif
}

ShotBuffer::ShotBuffer(
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader_in,
    std::unique_ptr<stim::MeasureRecordWriter> writer_in,
    int num_partitions,
    int num_virtual_boundaries,
    int num_observables)
    : reader(std::move(reader_in)), writer(std::move(writer_in))
{
#ifdef ENABLE_SHOT_BUFFERS
    const int num_containers = NUM_BUFFERS_PER_UNIT;
#else
    const int num_containers = 1;
#endif
    buffer.reserve(static_cast<size_t>(num_containers));
    for (size_t i = 0; i < num_containers; ++i) {
        buffer.emplace_back(
            num_partitions, num_virtual_boundaries, num_observables);
    }
}

void ShotBuffer::write_result_and_get_next_shot(
    int shot_container_id,
    std::vector<int>& node_part_id
#ifdef USE_SHMEM
    , bool write_results
#endif
) {
    std::unique_lock<std::mutex> lock(m);
    cv.wait(lock, [&] {
        return next_shot_container_id == shot_container_id;
    });
    auto& shot = buffer[shot_container_id];
#ifdef USE_SHMEM
    if (write_results) {
#endif
        if (DEBUG) {
            std::cout << "  T" << omp_get_thread_num() << " writing results for shot buffer " << shot_container_id << std::endl;
        }
        for (size_t k = 0; k < shot.num_observables; k++) {
            writer->write_bit(shot.res.obs_crossed[k]);
        }
        writer->write_end();
#ifdef USE_SHMEM
    }
#endif
    read_shot(shot_container_id, node_part_id);
#ifdef ENABLE_SHOT_BUFFERS
    if (++next_shot_container_id >= NUM_BUFFERS_PER_UNIT) {
        next_shot_container_id = 0;
    }
#endif
    lock.unlock();
    cv.notify_all();
}

// Reads a shot
// Should only before multi-threading or when a thread hold the lock
void ShotBuffer::read_shot(int shot_container_id, std::vector<int>& node_part_id) {
    auto& shot = buffer[shot_container_id];
    shot.clear();
#ifdef ENABLE_SHOT_BUFFERS
    if (last_shot_container_id < 0) {
#endif
        bool shot_read = pm::start_and_read_entire_record_buffered(*reader, shot.sparse_shot);
        if (shot_read) {
            for (auto det : shot.sparse_shot.hits) {
                int part_id = node_part_id[det];
                if (part_id >= 0) {
                    shot.partition_hits[part_id].push_back(det);
                } else {
                    shot.virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot.current_buffer_round++;
            shot.current_buffer_round.notify_all();
        } else {
#ifdef ENABLE_SHOT_BUFFERS
            last_shot_container_id = shot_container_id - 1;
            if (last_shot_container_id < 0) {
                last_shot_container_id = NUM_BUFFERS_PER_UNIT - 1;
            }
#else
            shot.current_buffer_round.store(-1);
            shot.current_buffer_round.notify_all();
#endif
        }
#ifdef ENABLE_SHOT_BUFFERS
    } else if (last_shot_container_id == shot_container_id) {
        for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
            buffer[i].current_buffer_round.store(-1);
            shot.current_buffer_round.notify_all();
        }
    }
#endif
}

void ShotBuffer::reset() {
    // Reset ShotBuffer counters
    next_shot_container_id = 0;
    last_shot_container_id = -1;
    // Reset per-ShotContainer state
    for (int i=0; i < buffer.size(); ++i) {
        buffer[i].reset();
    }
}

}  // namespace pm
