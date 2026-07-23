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

#ifndef PYMATCHING2_DECODING_UNIT_H
#define PYMATCHING2_DECODING_UNIT_H

#include <condition_variable>
#include <mutex>
#include <vector>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/shot_buffer.h"

#ifdef USE_SHMEM
#include <shmem.h>
#endif

#include "pymatching/sparse_blossom/flooder/blossom_child.h"

namespace pm {

// Forward declared rather than #included: user_graph.h conditionally #includes this
// header (to get SharedMatchingGraph), so a full #include here would cycle back before
// UserGraph is declared. A forward declaration is sufficient since UserGraph only
// appears by value in the DecodingUnit constructor's declaration; decoding_unit.cc
// (the definition) and namespaced_main.cc (the call site) both #include user_graph.h
// directly, giving them the complete type where it's actually needed.
class UserGraph;

#ifdef USE_SHMEM
#define SHMEM_NUM_ATOMICS_PER_CROSS_RANK_FUSION 3

struct FusionSummary {
    GraphFillRegion* regions_ptr_base; // index where solving partition's arena begins
    DetectorNode* static_nodes_base; // Base address of graph.nodes on the sender PE
    size_t blossom_children_size;
    size_t regions_to_unmatch_size;
    GraphFillRegion* regions_to_unmatch[]; 
    // uint64_t shmem_bitmap[]; // implicitly follows
};
#endif

// A decoding unit is responsible for decoding a connected decoding graph in parallel
//   Connected decoding graphs are partitioned for parallel solving
//   A decoding task involves solving a partition or fusing two solved partition along a virtual boundary
//   Decoding units implement multi-threading and OpenSHMEM behavoior
struct DecodingUnit {
    // Graph
    SharedMatchingGraph graph;

    // Shots
    std::shared_ptr<ShotBuffer> shot_buffer;

    // Params
    bool ensure_search_flooder_included;
    bool enable_correlations;
#ifdef ENABLE_DRAW_FLAGS
    bool draw_frames;
#endif

    // Solvers
    size_t num_threads;
    // size_t num_partition_units;
    size_t num_solvers_per_buffer;
    std::vector<std::shared_ptr<Mwpm>> solvers;

    std::vector<size_t> my_partition_task_ids;

#ifdef USE_SHMEM
    
    int pid;
    int n_pes;

    // Symmetric memory
    DetectorNodeEphemeralFields* node_ephemeral_fields_ptr{ nullptr }; // ephemeral field buffers for DetectorNodes
    GraphFillRegion* regions_ptr{ nullptr }; // buffers for region SHMEMArenas
    BlossomChild* child_edges_ptr{ nullptr }; // blossom parent-child relationship buffer
    // uint64_t* atomics_ptr{ nullptr }; // shot counter for each shot container
    uint64_t* task_status_ptr{ nullptr }; // status and signal for each fusion task
    FusionSummary* task_fusion_summary_ptr { nullptr }; // fusion summary for each fusion task
    // Info used for symmetric memory accesses
    size_t nodes_nelems_per_buffer;
    size_t regions_nelems_per_solver;
    size_t child_edges_nelems_per_solver;
    size_t task_fusion_summary_size_per_task;
    size_t regions_matched_to_vb_nelems;
    int num_cross_rank_fusions = 0;  // ROUND=2, OBS=num_seams; set in constructor before alloc

    inline FusionSummary* get_fusion_summary_ptr(size_t shot_container_id, bool iamleft) {
        bool even = !(pid % 2);
        bool take_second_slot = (even && iamleft) || (!even && !iamleft);
        return reinterpret_cast<FusionSummary*>(
            reinterpret_cast<char*>(task_fusion_summary_ptr)
            + (shot_container_id * num_cross_rank_fusions + take_second_slot) * task_fusion_summary_size_per_task
        );
    }

    inline uint64_t* get_task_status_ptr(size_t shot_container_id, bool iamleft) {
        // There are two uint64s per slot: status and signal
        //   Simple left/right case: num_cross_rank_fusions=2 slots per buffer
        //   Even PEs map slots (to left, to right), odd PEs (to right, to left)
        //   This ensures tasks correspond on each PE
        bool even = !(pid % 2);
        bool take_second_slot = (even && iamleft) || (!even && !iamleft);
        return task_status_ptr + SHMEM_NUM_ATOMICS_PER_CROSS_RANK_FUSION*(num_cross_rank_fusions*shot_container_id + take_second_slot);
    }

    // OBS strategy: index by seam index s directly (both PEs compute seam_infos identically)
    inline uint64_t* get_task_status_ptr_for_seam(size_t shot_id, int seam_idx) {
        return task_status_ptr
            + SHMEM_NUM_ATOMICS_PER_CROSS_RANK_FUSION
              * (num_cross_rank_fusions * shot_id + seam_idx);
    }
    inline FusionSummary* get_fusion_summary_ptr_for_seam(size_t shot_id, int seam_idx) {
        return reinterpret_cast<FusionSummary*>(
            reinterpret_cast<char*>(task_fusion_summary_ptr)
            + (num_cross_rank_fusions * shot_id + seam_idx)
              * task_fusion_summary_size_per_task);
    }

    inline DetectorNodeEphemeralFields* get_node_fields_ptr(size_t shot_container_id, int partition_id) {
        auto& bounds = graph.partition_bounds[partition_id];
        return node_ephemeral_fields_ptr + shot_container_id * nodes_nelems_per_buffer + bounds.first;
    }

    inline GraphFillRegion* get_regions_ptr(size_t shot_container_id, int partition_id) {
        return regions_ptr + (shot_container_id * graph.num_partitions + partition_id) * regions_nelems_per_solver;
    }

    inline BlossomChild* get_child_edges_ptr(size_t shot_container_id, int partition_id) {
        return child_edges_ptr + (shot_container_id * graph.num_partitions + partition_id) * child_edges_nelems_per_solver;
    }

    inline size_t get_cross_rank_fusion_idx(size_t shot_container_id, size_t index) {
        return num_cross_rank_fusions * shot_container_id + index;
    }

#endif

    // Initialization Functions
    // dem_for_drawing is only consulted when draw_frames is true, to extract detector
    // coordinates for visualization; it may be null when the graph was loaded from a
    // binary graph cache (--graph_cache_path) rather than parsed fresh from a DEM, in
    // which case passing draw_frames=true will throw.
    DecodingUnit(
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        pm::UserGraph user_graph,
        weight_int num_distinct_weights,
        bool ensure_search_flooder_included,
        bool enable_correlations
#ifdef ENABLE_DRAW_FLAGS
        , bool draw_frames
        , const stim::DetectorErrorModel* dem_for_drawing = nullptr
#endif
        );

    ~DecodingUnit();

    void build_tasks_for_round_partitioning();
#ifdef USE_SHMEM
    void build_tasks_for_obs_patch_partitioning();
#endif

    void build_solvers();

    inline int get_solver_id(int shot_container_id, int partition)
    {
        return num_solvers_per_buffer * shot_container_id + partition;
    }


    // Decoding Functions
#ifdef USE_SHMEM
    void send_solution_to_remote_pe(size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &task, std::ofstream &t_out);
    bool get_solution_from_remote_pe(
        size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &task, std::ofstream &t_out,
        std::vector<uint64_t>& hitsref); // returns whether solving is necessary
    // Extracts a cross-rank fusion's received window immediately, inline on the resolving thread --
    // not queued (queuing would only add overhead here, since this thread is already doing the
    // work). Computes the received partition/vb range directly from crt's own fields
    // (division-strategy-aware), then extracts it the same way an ordinary leaf/vb would be.
    // Assumes crt->part's own vb has already been divided (see divide_vb) by the caller.
    void extract_crt_received_window(ShotContainer& shot, size_t shot_container_id, CrossRankTask& crt, int tid, std::ofstream* t_out = nullptr);
#endif

    // Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md). Universal
    // across build configs -- both plain USE_THREADS and USE_SHMEM+USE_THREADS builds post to and
    // drain the same job queue.

    // Separates the two subgraphs joined at vb_id: loops every node in graph.vb_bounds[vb_id]
    // (not just ones with "hits" -- a blossom can span the vb without either side registering a hit
    // exactly there) and shatters/extracts any non-null region_that_arrived_top. Synchronous/inline,
    // never queued -- this is what makes it safe to later post an extraction job for either side
    // independently: once divided, neither side's regions reference across the boundary anymore.
    // anchor_partition must be a genuine descendant leaf id (caller's own vb_solver_offset-resolved
    // partition, not an arbitrary fixed choice -- see definition for why).
    // prune_target: the fusion task whose vb this is (Task* for a checkpoint, CrossRankTask* for a
    // cross-rank fusion -- both derive from TaskBase). A blossom shattered here can span farther than
    // this vb and destroy a region still referenced in prune_target->regions_matched_to_virtual_
    // boundary (kept there for a LATER setup() call to consume) -- passing prune_target makes divide_vb
    // prune any now-stale pointers out of that list afterward. Pass nullptr only when no future
    // setup() call will ever read that list again (post-hoc chunking, run after the whole tree is
    // already fully solved).
    // t_out: optional per-thread debug stream; when DEBUG and non-null, logs the vb being divided and
    // the anchor/prune_target so divide timing can be traced end-to-end (see this-is-a-broader-
    // purrfect-crystal.md).
    void divide_vb(ShotContainer& shot, int shot_container_id, int vb_id, int anchor_partition, int tid, TaskBase* prune_target, std::ofstream* t_out = nullptr);

    // Shatter+extract an entire unit subtree (job.subtree_root), accumulating into
    // shot.thread_results[tid] (bit-packed) or shot.res directly (extended observables, under omp
    // critical). Any solver from this job's shot_container_id block works as a scratch accumulator
    // (region ownership is globally node-indexed, not solver-private) -- tid must be the
    // *processing* thread's own, since get_solver_id depends on tid under plain USE_THREADS.
    // Decrements pending_extraction_jobs when done. Precondition: subtree_root is a fully-closed
    // unit (see divide_vb) -- every vb this subtree touches on its way to being posted has already
    // been divided by the poster, so no vb-boundary handling is needed inside the walk itself...
    // except each internal fusion's *own* vb, still extracted here via the simpler (pre-existing,
    // not yet fixed -- tracked separately) virtual_boundary_hits-based approach.
    // t_out: optional per-thread debug stream; when DEBUG and non-null, logs the vb/task address of
    // the job being processed.
    void process_extraction_job(ShotContainer& shot, const ExtractionJob& job, int tid, std::ofstream* t_out = nullptr);

    // Non-preemptive (extract_preemptively == false) post-hoc chunking: a single top-down walk of the
    // balanced-over-units tree, driven entirely by the extraction-role tags (see
    // decoding_task.h/build_tasks_for_round_partitioning) -- no dynamic leaf-counting needed. At a
    // unit root, post it and return; at a connector, divide its own vb inline (making both children
    // independently safe to post/recurse into) and recurse into both.
    void post_hoc_chunk_and_post(Task* node, ShotContainer& shot, int shot_container_id, int tid, std::ofstream* t_out = nullptr);

    void decode_shots();
    void reset();
};

}  // namespace pm

#endif  // PYMATCHING2_DECODING_UNIT_H