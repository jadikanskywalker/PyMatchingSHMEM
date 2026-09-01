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
// enum ShotStatus : uint64_t { READY, PUT_SUMMARY, PUT_RESULT, WROTE_RESULT };
#endif

// Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md Design §3): a
// self-describing extraction job. Any thread can process a job; a job always names a fully-closed
// unit subtree -- closed meaning its outer boundaries have been "divided" (see DecodingUnit::
// divide_vb). Cross-rank fusion extraction of received window is handled inline avoiding posting.
// ready gates a claimer reading subtree_root/shot_container_id: a slot exists (default-constructed,
// unposted) as soon as ShotIOResource itself is built, so a thread can claim slot i before job i is
// actually posted -- see ShotIOResource::try_claim_extraction_job. alignas(64) bumps the whole
// struct's alignment so adjacent slots in ShotIOResource::extraction_jobs never share a cache line
// (each slot is spun on by a different thread while waiting).
struct ExtractionJob {
    int shot_container_id;
    Task* subtree_root;
    pm::MatchingResult obs_mask;
    pm::ExtendedMatchingResult res;
    alignas(64) std::atomic<int> ready{0};

    ExtractionJob() = default;
    ExtractionJob(const ExtractionJob&) = delete;
    ExtractionJob& operator=(const ExtractionJob&) = delete;

    // ready makes the implicit move ctor/assign deleted; std::vector<ExtractionJob>::reserve()
    // still needs this to at least compile (even though, in practice, reserve() is only ever
    // called on an empty vector here -- see ShotIOResource's ctor -- so there's nothing to
    // actually move at runtime). Mirrors ShotIOResource's own atomic-skipping move pattern.
    ExtractionJob(ExtractionJob&& other) noexcept
        : shot_container_id(other.shot_container_id),
          subtree_root(other.subtree_root),
          obs_mask(other.obs_mask),
          res(std::move(other.res)) {
        ready.store(other.ready.load(std::memory_order_relaxed), std::memory_order_relaxed);
    }
    ExtractionJob& operator=(ExtractionJob&& other) noexcept {
        if (this != &other) {
            shot_container_id = other.shot_container_id;
            subtree_root = other.subtree_root;
            obs_mask = other.obs_mask;
            res = std::move(other.res);
            ready.store(other.ready.load(std::memory_order_relaxed), std::memory_order_relaxed);
        }
        return *this;
    }
};

// ShotIOResource isolates everything needs to store hits (input for a shot), manage extraction jobs,
// and build results (output for a shot). This is decoupled from ShotContainer to allow multiple shots
// to be active in the same buffer, avoiding a global current_buffer_round atomic -> inter-shot barrier.
// Instead, partition tasks are marked ready when their extraction unit is extracted; a thread cannot 
// only steal the partition task for the next shot when ready. This effectively moves inter-shot
// synchronization to the task level (granular) instead of all threads racing on a single atomic.
struct ShotIOResource {
    // -- Input --
    // Hits
    stim::SparseShot sparse_shot;
    // Partitioned hits
    std::vector<std::vector<uint64_t>> partition_hits;
    std::vector<std::vector<uint64_t>> virtual_boundary_hits;
    // -- Output --
    // num_observables deliberately NOT stored here -- it's a fixed per-run constant, already
    // cached once on ShotBuffer::num_observables; no need for every shot to carry its own copy.
    pm::MatchingResult obs_mask;
    pm::ExtendedMatchingResult res;
    // Per-thread partial MatchingResult accumulator
    std::vector<pm::MatchingResult> thread_results;
    // Per-thread partial ExtendedMatchingResult accumulator (>64-obs path only, but always sized
    // so callers don't need to branch). extract_paths_from_match_edges (mwpm.cc) writes directly
    // into obs_crossed/weight with plain, non-atomic ops -- multiple threads can be inside
    // process_extraction_job/divide_vb/extract_crt_received_window simultaneously for the same
    // shot (try_claim_extraction_job's lock-free cursor lets any idle thread claim any posted
    // job), so writing straight into a single shared `res` races. Mirrors thread_results exactly:
    // each thread accumulates into its own slot, combined once single-threaded at the end of the
    // shot under the same pending_extraction_jobs==0 gate that already exists there.
    std::vector<pm::ExtendedMatchingResult> thread_extended_results;
    // num_task_roots also deliberately NOT stored here -- same reasoning as num_observables
    // above: a fixed per-run constant (always 1 for ROUND, my_obs_count for OBS, never varies
    // per shot), cached once on ShotBuffer::num_task_roots instead.
    // Incremented by each thread when it finishes all its roots for a shot.
    // Thread to reach num_task_roots writes result
    alignas(64) std::atomic<int64_t> num_roots_done{0};
    // Unit-checkpointed extraction queue (Design §3-6). extraction_jobs is pre-sized to exactly
    // num_extraction_units slots at construction (a static, verified-exact per-shot constant, see
    // ShotBuffer::num_extraction_units) and never grown afterward -- posting and claiming both
    // reduce to a single fetch_add each (see post_extraction_job/try_claim_extraction_job), no
    // mutex needed: each poster's fetch_add hands it a unique, exclusively-owned slot index.
    std::vector<ExtractionJob> extraction_jobs;
    alignas(64) std::atomic<size_t> extraction_posted_count{0};  // next slot index to post into
    alignas(64) std::atomic<size_t> extraction_claim_cursor{0};  // next slot index to claim/wait on
    // Gates ShotContainer reuse until every posted job has actually been executed (not just
    // claimed) -- incremented on post, decremented once a claimed job finishes processing.
    alignas(64) std::atomic<int> pending_extraction_jobs{0};

    ShotIOResource(
        int num_partitions,
        int num_virtual_boundaries,
        int num_observables,  // sizes res/obs_mask/thread_extended_results only, not stored
        int num_threads,
        int num_extraction_units);

    // Custom move ctor/assign needed so std::vector<ShotIOResource> can grow via emplace_back
    // during read_all_shots_and_create_IO_resources (io_resources is sized incrementally, one shot
    // at a time, since the total shot count isn't known upfront) -- mirrors ShotContainer's
    // existing pattern: manually .load()/.store() every atomic. extraction_jobs itself moves as a
    // plain std::vector (always fine regardless of ExtractionJob containing an atomic -- vector
    // move never touches individual elements).
    ShotIOResource(const ShotIOResource&) = delete;
    ShotIOResource& operator=(const ShotIOResource&) = delete;
    ShotIOResource(ShotIOResource&& other) noexcept;
    ShotIOResource& operator=(ShotIOResource&& other) noexcept;

    // Per-value parameters, not a whole ExtractionJob&: ExtractionJob now owns an atomic (ready),
    // so there's nothing sensible to construct-and-pass-by-reference at a call site -- callers
    // only ever have a (shot_container_id, subtree_root) pair to publish. t_out: optional
    // per-thread debug stream; when DEBUG and non-null.
    void post_extraction_job(int shot_container_id, Task* subtree_root, std::ofstream* t_out = nullptr);
    // num_extraction_units: the exact, static total this shot will ever post (see ShotBuffer::
    // num_extraction_units) -- once claim_cursor reaches it, every slot has been claimed and no
    // more will ever exist, so this returns nullptr immediately rather than waiting. Otherwise
    // spins on the claimed slot's own ready flag until that specific job is posted.
    ExtractionJob* try_claim_extraction_job(size_t num_extraction_units);

    // Per-repeat (--num_repeats) reset only -- never called per-shot. Resets thread_results,
    // thread_extended_results, res, obs_mask, num_roots_done, and the whole extraction queue.
    // Deliberately does NOT touch sparse_shot/partition_hits/virtual_boundary_hits (read once
    // upfront and reused verbatim across repeats).
    void reset();
};

// ShotContainer isolates everything needed to decode a single shot in parallel,
// all the way down to reserved ephemeral field slots on DetectorNodes.
struct ShotContainer {
    std::vector<Task> tasks;  // Tasks handle dynamic fusion tree synchonization
    // Local-observable-boundary seams (OBS partitioning only) -- never SHMEM-specific, unlike
    // cross_rank_tasks, so unconditional even though only build_tasks_for_obs_patch_partitioning
    // populates it today.
    std::vector<LocalSeamTask> local_seam_tasks;
#ifdef USE_SHMEM
    std::vector<CrossRankTask> cross_rank_tasks;
#endif

    // No constructor params -- tasks/local_seam_tasks/cross_rank_tasks are sized/reserved later
    // by build_tasks_for_round_partitioning/build_tasks_for_obs_patch_partitioning, which have
    // the actual partition/vb counts needed for a correct worst-case reservation. Nothing here
    // depends on partition/vb/observable counts now that Tier B construction has moved to
    // ShotIOResource.
    ShotContainer() = default;

    ShotContainer(const ShotContainer&) = delete;
    ShotContainer& operator=(const ShotContainer&) = delete;

    ShotContainer(ShotContainer&& other) noexcept;
    ShotContainer& operator=(ShotContainer&& other) noexcept;

    void reset();
};

// Rotating buffer of shot containers with ordered shot
// reading and result writing.
struct ShotBuffer {
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader;
    std::unique_ptr<stim::MeasureRecordWriter> writer;

    std::vector<ShotContainer> buffer;

    std::vector<ShotIOResource> io_resources;

    // Stashed from the constructor so read_all_shots_and_create_IO_resources (called separately,
    // once the total shot count is known) doesn't need them re-passed. num_observables also
    // serves as the single source of truth callers should read directly (ShotIOResource itself
    // no longer stores a per-shot copy).
    int num_partitions{0};
    int num_virtual_boundaries{0};
    int num_observables{0};
    // Fixed for the whole run (always 1 for ROUND, my_obs_count for OBS) -- set once by
    // read_all_shots_and_create_IO_resources, read directly by callers instead of a per-shot
    // ShotIOResource field (it never varies per shot, same reasoning as num_observables above).
    int num_task_roots{0};
    // Exact total extraction jobs a shot will ever post (ceil(my_partition_count/L) for ROUND,
    // summed per-obs-patch for OBS) -- also fixed for the whole run, set once by
    // read_all_shots_and_create_IO_resources, and used both to pre-size each ShotIOResource's
    // extraction_jobs slots and as try_claim_extraction_job's exact termination bound.
    int num_extraction_units{0};

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

    // Writes shot_id's finished result to the output file once it's this container's turn
    // (unchanged cv/next_shot_container_id turn-taking -- serializes *output* order, an
    // orthogonal concern from the per-partition decode desynchronization this design targets;
    // shot_container_id cycles 0..N-1 across shots 0..total-1 in lockstep with real shot order,
    // so per-container turn-taking still correctly serializes writes in global shot order even
    // though containers are now reused many times). No longer reads a new shot afterward --
    // reading is fully upfront now, see read_all_shots_and_create_IO_resources.
    void write_shot_result(int shot_container_id, size_t shot_id
#ifdef USE_SHMEM
        , bool write_results=true
#endif
    );

    // Reads every shot in the input file upfront, partitioning each shot's hits by
    // node_part_id and appending one fully-populated ShotIOResource per shot to io_resources.
    // num_threads sizes each ShotIOResource's thread_results/thread_extended_results at
    // construction, rather than those being set post-hoc by task-tree-building code.
    // num_task_roots_in/num_extraction_units_in are stashed directly onto ShotBuffer::
    // num_task_roots/num_extraction_units (fixed per-run constants, not per-shot -- see the field
    // comments above), and num_extraction_units_in also sizes each new ShotIOResource's
    // extraction_jobs slots.
    void read_all_shots_and_create_IO_resources(
        std::vector<int>& node_part_id, int num_threads, int num_task_roots_in,
        int num_extraction_units_in);

    // -- Helpers below decouple DecodingUnit from ShotIOResource's internal layout: every one of
    // these is keyed by shot_id (Tier B, one slot per real shot) rather than shot_container_id
    // (Tier A, the small reused Task-tree pool) -- see docs/decentralized_shot_sync_design.md. --

    // Total real shot count -- a property of the whole buffer, not of one shot_id, but still
    // routed through a helper rather than a direct io_resources.size() read: a future
    // streaming/circular-buffer redesign (design doc §10) wouldn't have a plain "vector size ==
    // total shots" shape, and callers should never need to know that today's implementation does.
    size_t num_shots() const;
    // Returns partition_hits[id] (is_vb=false) or virtual_boundary_hits[id] (is_vb=true) for shot_id.
    std::vector<uint64_t>& hits(size_t shot_id, bool is_vb, int id);
    pm::ExtendedMatchingResult& result(size_t shot_id);
    pm::MatchingResult& obs_mask(size_t shot_id);
    // No num_observables(shot_id) helper -- it's a fixed per-run constant, read directly off
    // ShotBuffer::num_observables at call sites instead (never per-shot).
    std::vector<pm::MatchingResult>& thread_results(size_t shot_id);
    std::vector<pm::ExtendedMatchingResult>& thread_extended_results(size_t shot_id);
    int64_t fetch_add_roots_done(size_t shot_id, int64_t n);
    // True once every root of this shot's tree has resolved (compares against the cached
    // num_task_roots internally, rather than making every caller re-derive the comparison).
    bool all_roots_done(size_t shot_id) const;
    // No num_task_roots(shot_id) helper -- fixed per-run constant, read directly off
    // ShotBuffer::num_task_roots at call sites instead (never per-shot).
    void post_extraction_job(size_t shot_id, int shot_container_id, Task* subtree_root, std::ofstream* t_out = nullptr);
    // Forwards the cached num_extraction_units bound -- callers never need to query "how many have
    // posted so far" separately; nullptr means every slot for this shot has been claimed.
    ExtractionJob* try_claim_extraction_job(size_t shot_id);
    int pending_extraction_jobs(size_t shot_id) const;
    // Called once a claimed job has actually finished processing (not just claimed) -- decrements
    // the counter post_extraction_job incremented at post time.
    void mark_extraction_job_finished(size_t shot_id);

    void reset();
};

}

#endif // PYMATCHING2_SHOT_BUFFER_H