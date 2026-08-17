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

#ifndef PYMATCHING2_SHOT_BUFFER_H
#define PYMATCHING2_SHOT_BUFFER_H

#include <atomic>
#include <condition_variable>
#include <fstream>
#include <mutex>
#include <vector>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"
#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"

namespace pm {

#ifdef USE_SHMEM
enum ShotStatus : uint64_t { READY, PUT_SUMMARY, PUT_RESULT, WROTE_RESULT };
#endif

// Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md Design §3): a
// self-describing extraction job, unconditional across build configs since checkpoint chains and
// post-hoc chunking both run under plain USE_THREADS too, not just USE_SHMEM. Any thread can process
// a job regardless of which shot/root it's nominally working on. A job always names a fully-closed
// unit subtree -- closed meaning both its boundaries (if any) have already been "divided" (see
// DecodingUnit::divide_vb) before the job is posted, so it's safe to extract independently, by any
// thread, at any time. Cross-rank fusion extraction is handled inline by the resolving thread
// instead of going through this queue (posting would only add overhead there).
struct ExtractionJob {
    int shot_container_id;
    Task* subtree_root;
};

struct LocalSeamContinuationJob {
    size_t shot_container_id;
    size_t special_task_id; // id of completed seam in child's SpecialTask vector
    Task* child_task;
};

// ShotContainer isolates everything needed to solve
// a single shot in parallel.
struct ShotContainer {
    alignas(64) std::atomic<int64_t> current_buffer_round{-1};  // Used for idle thread spin-wait until new shot read

    stim::SparseShot sparse_shot;

    std::vector<std::vector<uint64_t>> partition_hits;
    std::vector<std::vector<uint64_t>> virtual_boundary_hits;

    int num_observables;
    pm::MatchingResult obs_mask;
    pm::ExtendedMatchingResult res;

    // Set during build_tasks_*; count of chain tops (parent==nullptr) across tasks + CRTs. Always 1
    // for round-partitioning (a single true root regardless of extract_preemptively/L -- checkpoint
    // chain-links are never independent roots, see this-is-a-broader-purrfect-crystal.md Design §2).
    int num_task_roots{0};
    // Incremented by each thread when it finishes all its roots for a shot.
    // Last thread (cumulative count reaches num_task_roots) resets to 0 and writes result.
    alignas(64) std::atomic<int64_t> num_roots_done{0};
    // Per-thread partial MatchingResult accumulator (indexed by omp thread id).
    // Sized to num_threads during build_tasks_*.
    std::vector<pm::MatchingResult> thread_results;

    // Unit-checkpointed extraction queue (Design §3-6). Posting (rare, ~num_partitions/L per shot) is
    // mutex-guarded; claiming (hot, up to ~2*num_partitions steal attempts per shot) is non-locking
    // via a fetch_add cursor. extraction_jobs' capacity is reserved once in the constructor and never
    // grown, so a claimer reading extraction_jobs[i] without the lock never races a reallocation.
    std::mutex extraction_post_mutex;
    std::vector<ExtractionJob> extraction_jobs;
    alignas(64) std::atomic<size_t> extraction_posted_count{0};  // published (release) job count
    alignas(64) std::atomic<size_t> extraction_claim_cursor{0};  // next unclaimed job index
    // Gates ShotContainer reuse until every posted job has actually been executed (not just
    // claimed) -- incremented on post, decremented once a claimed job finishes processing.
    alignas(64) std::atomic<int> pending_extraction_jobs{0};

    std::mutex local_seam_continuation_post_mutex;
    std::vector<LocalSeamContinuationJob> local_seam_continuation_jobs;
    alignas(64) std::atomic<size_t> local_seam_continuation_posted_count{0};  // published (release) job count
    alignas(64) std::atomic<size_t> local_seam_continuation_claim_cursor{0};  // next unclaimed job index
    // Gates ShotContainer reuse until every posted job has actually been executed (not just
    // claimed) -- incremented on post, decremented once a claimed job finishes processing.
    alignas(64) std::atomic<int> local_seam_continuation_extraction_jobs{0};

    // t_out: optional per-thread debug stream; when DEBUG and non-null, logs the vb/task address of
    // the job being posted (see this-is-a-broader-purrfect-crystal.md) so posting/draining can be
    // traced end-to-end alongside decode_shots()'s existing per-thread traces.
    void post_extraction_job(const ExtractionJob& job, std::ofstream* t_out = nullptr);
    ExtractionJob* try_claim_extraction_job();

    std::vector<Task> tasks;  // Tasks handle dynamic fusion tree synchonization
    // Local-observable-boundary seams (OBS partitioning only) -- never SHMEM-specific, unlike
    // cross_rank_tasks, so unconditional even though only build_tasks_for_obs_patch_partitioning
    // populates it today.
    std::vector<LocalSeamTask> local_seam_tasks;
#ifdef USE_SHMEM
    std::vector<CrossRankTask> cross_rank_tasks;
#endif

    ShotContainer(int num_partitions, int num_virtual_boundaries, int num_observables);

    ShotContainer(const ShotContainer&) = delete;
    ShotContainer& operator=(const ShotContainer&) = delete;

    ShotContainer(ShotContainer&& other) noexcept;
    ShotContainer& operator=(ShotContainer&& other) noexcept;

    void clear();
    void reset();
};

// Rotating buffer of shot containers with ordered shot
// reading and result writing.
struct ShotBuffer {
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader;
    std::unique_ptr<stim::MeasureRecordWriter> writer;

    std::vector<ShotContainer> buffer;

    alignas(64) int64_t next_shot_container_id{ 0 };
    alignas(64) int64_t last_shot_container_id{ -1 };

    // Lock for result writing & shot reading
    std::mutex m;
    std::condition_variable cv;

    ShotBuffer(
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        int num_partitions,
        int num_virtual_boundaries,
        int num_observables);

    // Reads next shot into shot_container_id
    void write_result_and_get_next_shot(int shot_container_id, std::vector<int> &node_part_id
#ifdef USE_SHMEM
        , bool write_results=true
#endif
    );

    void read_shot(int shot_container_id, std::vector<int> &node_part_id);

    void reset();
};

}

#endif // PYMATCHING2_SHOT_BUFFER_H