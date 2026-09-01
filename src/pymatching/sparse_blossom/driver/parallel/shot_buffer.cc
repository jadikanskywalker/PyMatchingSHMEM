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

#include <immintrin.h>
#include <iostream>
#include <utility>

#include <omp.h>

#ifdef USE_SHMEM
#include <shmem.h>
#endif

#ifdef SCOREP_USER_ENABLE
#include <scorep/SCOREP_User.h>
#endif

namespace pm {

ShotIOResource::ShotIOResource(
    int num_partitions,
    int num_virtual_boundaries,
    int num_observables,
    int num_threads,
    int num_extraction_units)
    : partition_hits(num_partitions),
      virtual_boundary_hits(num_virtual_boundaries),
      res(num_observables),
      thread_results(num_threads),
      thread_extended_results(num_threads, pm::ExtendedMatchingResult(num_observables))
{
    // Pre-size to exactly num_extraction_units slots, up front, and never grow afterward --
    // reserve() first means emplace_back() below never reallocates (the claim side reads
    // extraction_jobs[i] without a lock, only safe if the backing storage never moves), and
    // emplace_back() constructs each slot in place, so ExtractionJob's embedded atomic is never
    // moved or copied.
    extraction_jobs.reserve(static_cast<size_t>(num_extraction_units));
    for (int i = 0; i < num_extraction_units; ++i) {
        extraction_jobs.emplace_back();
    }
}

ShotIOResource::ShotIOResource(ShotIOResource&& other) noexcept
    : sparse_shot(std::move(other.sparse_shot)),
      partition_hits(std::move(other.partition_hits)),
      virtual_boundary_hits(std::move(other.virtual_boundary_hits)),
      obs_mask(other.obs_mask),
      res(std::move(other.res)),
      thread_results(std::move(other.thread_results)),
      thread_extended_results(std::move(other.thread_extended_results)),
      extraction_jobs(std::move(other.extraction_jobs))
      // extraction_jobs moves as a plain std::vector -- always fine regardless of ExtractionJob
      // containing an atomic, since vector move never touches individual elements.
{
    num_roots_done.store(other.num_roots_done.load());
    extraction_posted_count.store(other.extraction_posted_count.load());
    extraction_claim_cursor.store(other.extraction_claim_cursor.load());
    pending_extraction_jobs.store(other.pending_extraction_jobs.load());
}

ShotIOResource& ShotIOResource::operator=(ShotIOResource&& other) noexcept {
    if (this != &other) {
        sparse_shot = std::move(other.sparse_shot);
        partition_hits = std::move(other.partition_hits);
        virtual_boundary_hits = std::move(other.virtual_boundary_hits);
        obs_mask = other.obs_mask;
        res = std::move(other.res);
        thread_results = std::move(other.thread_results);
        thread_extended_results = std::move(other.thread_extended_results);
        extraction_jobs = std::move(other.extraction_jobs);
        num_roots_done.store(other.num_roots_done.load());
        extraction_posted_count.store(other.extraction_posted_count.load());
        extraction_claim_cursor.store(other.extraction_claim_cursor.load());
        pending_extraction_jobs.store(other.pending_extraction_jobs.load());
    }
    return *this;
}

void ShotIOResource::post_extraction_job(int shot_container_id, Task* subtree_root, std::ofstream* t_out) {
    if (DEBUG && t_out) {
        *t_out << "  POST_JOB " << (subtree_root->is_fusion ? "vb=" : "p=") << subtree_root->part
               << " task=" << subtree_root << std::endl << std::flush;
    }
    // Gates ShotContainer reuse (see pending_extraction_jobs's own comment); ordering relative to
    // the slot write below doesn't matter -- only ready's release/acquire pair below governs
    // visibility of the job's own data.
    pending_extraction_jobs.fetch_add(1, std::memory_order_acq_rel);
    // fetch_add hands this poster a unique, exclusively-owned slot index -- no other poster can
    // ever write the same index, so no lock is needed to guard the write into extraction_jobs[idx].
    size_t idx = extraction_posted_count.fetch_add(1, std::memory_order_relaxed);
    extraction_jobs[idx].shot_container_id = shot_container_id;
    extraction_jobs[idx].subtree_root = subtree_root;
    extraction_jobs[idx].ready.store(1, std::memory_order_release);
}

ExtractionJob* ShotIOResource::try_claim_extraction_job(size_t num_extraction_units) {
    // fetch_add unconditionally hands this claimer a unique slot index -- no CAS/retry needed,
    // unlike posting there's no shared data write to protect here.
    size_t idx = extraction_claim_cursor.fetch_add(1, std::memory_order_relaxed);
    if (idx >= num_extraction_units) {
        // Every slot this shot will ever have has already been claimed by someone (possibly
        // still in-flight) -- nothing left to wait for.
        return nullptr;
    }
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_DEFINE(extraction_slot_wait);
    SCOREP_USER_REGION_BEGIN(extraction_slot_wait, "Extraction Slot Wait", SCOREP_USER_REGION_TYPE_COMMON);
#endif
    int spin_count = 1;
    constexpr int max_spin_count = 1024;
    while (extraction_jobs[idx].ready.load(std::memory_order_acquire) == 0) {
        for (int i = 0; i < spin_count; ++i) {
            _mm_pause();
        }
        if (spin_count < max_spin_count) {
            spin_count *= 2;
        }
    }
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_END(extraction_slot_wait);
#endif
    return &extraction_jobs[idx];
}

void ShotIOResource::reset() {
    for (auto& tr : thread_results) {
        tr = pm::MatchingResult{};
    }
    for (auto& ter : thread_extended_results) {
        ter.reset();
    }
    res.reset();
    obs_mask = pm::MatchingResult{};
    num_roots_done.store(0, std::memory_order_relaxed);
    // Slots are reused across --num_repeats repeats, never regrown -- reset each in place rather
    // than clearing the vector.
    for (auto& job : extraction_jobs) {
        job.subtree_root = nullptr;
        job.ready.store(0, std::memory_order_relaxed);
    }
    extraction_posted_count.store(0, std::memory_order_relaxed);
    extraction_claim_cursor.store(0, std::memory_order_relaxed);
    pending_extraction_jobs.store(0, std::memory_order_relaxed);
}

ShotContainer::ShotContainer(ShotContainer&& other) noexcept
    : tasks(std::move(other.tasks)),
      local_seam_tasks(std::move(other.local_seam_tasks))
{
#ifdef USE_SHMEM
    cross_rank_tasks = std::move(other.cross_rank_tasks);
#endif
}

ShotContainer& ShotContainer::operator=(ShotContainer&& other) noexcept {
    if (this != &other) {
        tasks = std::move(other.tasks);
        local_seam_tasks = std::move(other.local_seam_tasks);
#ifdef USE_SHMEM
        cross_rank_tasks = std::move(other.cross_rank_tasks);
#endif
    }
    return *this;
}

void ShotContainer::reset() {
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
    int num_partitions_in,
    int num_virtual_boundaries_in,
    int num_observables_in)
    : reader(std::move(reader_in)), writer(std::move(writer_in)),
      num_partitions(num_partitions_in),
      num_virtual_boundaries(num_virtual_boundaries_in),
      num_observables(num_observables_in)
{
#ifdef ENABLE_SHOT_BUFFERS
    const int num_containers = NUM_BUFFERS_PER_UNIT;
#else
    const int num_containers = 1;
#endif
    buffer.reserve(static_cast<size_t>(num_containers));
    for (size_t i = 0; i < num_containers; ++i) {
        buffer.emplace_back();
    }
}

void ShotBuffer::write_shot_result(
    int shot_container_id,
    size_t shot_id
#ifdef USE_SHMEM
    , bool write_results
#endif
) {
    std::unique_lock<std::mutex> lock(m);
    cv.wait(lock, [&] {
        return next_shot_container_id == shot_container_id;
    });
#ifdef USE_SHMEM
    if (write_results) {
#endif
        if (DEBUG) {
            std::cout << "  T" << omp_get_thread_num() << " writing results for shot " << shot_id << std::endl;
        }
        auto& r = io_resources[shot_id].res;
        for (int k = 0; k < num_observables; k++) {
            writer->write_bit(r.obs_crossed[k]);
        }
        writer->write_end();
#ifdef USE_SHMEM
    }
#endif
#ifdef ENABLE_SHOT_BUFFERS
    if (++next_shot_container_id >= NUM_BUFFERS_PER_UNIT) {
        next_shot_container_id = 0;
    }
#endif
    lock.unlock();
    cv.notify_all();
}

void ShotBuffer::read_all_shots_and_create_IO_resources(
    std::vector<int>& node_part_id, int num_threads, int num_task_roots_in,
    int num_extraction_units_in) {
    num_task_roots = num_task_roots_in;
    num_extraction_units = num_extraction_units_in;
    io_resources.clear();
    stim::SparseShot sparse_shot;
    while (pm::start_and_read_entire_record_buffered(*reader, sparse_shot)) {
        io_resources.emplace_back(
            num_partitions, num_virtual_boundaries, num_observables, num_threads, num_extraction_units);
        auto& io = io_resources.back();
        for (auto det : sparse_shot.hits) {
            int part_id = node_part_id[det];
            if (part_id >= 0) {
                io.partition_hits[part_id].push_back(det);
            } else {
                io.virtual_boundary_hits[-1 * part_id - 1].push_back(det);
            }
        }
        sparse_shot.clear();
    }
}

size_t ShotBuffer::num_shots() const {
    return io_resources.size();
}

std::vector<uint64_t>& ShotBuffer::hits(size_t shot_id, bool is_vb, int id) {
    return is_vb ? io_resources[shot_id].virtual_boundary_hits[id]
                 : io_resources[shot_id].partition_hits[id];
}

pm::ExtendedMatchingResult& ShotBuffer::result(size_t shot_id) {
    return io_resources[shot_id].res;
}

pm::MatchingResult& ShotBuffer::obs_mask(size_t shot_id) {
    return io_resources[shot_id].obs_mask;
}

std::vector<pm::MatchingResult>& ShotBuffer::thread_results(size_t shot_id) {
    return io_resources[shot_id].thread_results;
}

std::vector<pm::ExtendedMatchingResult>& ShotBuffer::thread_extended_results(size_t shot_id) {
    return io_resources[shot_id].thread_extended_results;
}

int64_t ShotBuffer::fetch_add_roots_done(size_t shot_id, int64_t n) {
    return io_resources[shot_id].num_roots_done.fetch_add(n, std::memory_order_acq_rel);
}

bool ShotBuffer::all_roots_done(size_t shot_id) const {
    return io_resources[shot_id].num_roots_done.load(std::memory_order_acquire) == num_task_roots;
}

void ShotBuffer::post_extraction_job(size_t shot_id, int shot_container_id, Task* subtree_root, std::ofstream* t_out) {
    io_resources[shot_id].post_extraction_job(shot_container_id, subtree_root, t_out);
}

ExtractionJob* ShotBuffer::try_claim_extraction_job(size_t shot_id) {
    return io_resources[shot_id].try_claim_extraction_job((size_t)num_extraction_units);
}

int ShotBuffer::pending_extraction_jobs(size_t shot_id) const {
    return io_resources[shot_id].pending_extraction_jobs.load(std::memory_order_acquire);
}

void ShotBuffer::mark_extraction_job_finished(size_t shot_id) {
    io_resources[shot_id].pending_extraction_jobs.fetch_sub(1, std::memory_order_acq_rel);
}

void ShotBuffer::reset() {
    next_shot_container_id = 0;
    last_shot_container_id = -1;
    for (auto& io : io_resources) {
        io.reset();
    }
    for (auto& sc : buffer) {
        sc.reset();
    }
}

}  // namespace pm
