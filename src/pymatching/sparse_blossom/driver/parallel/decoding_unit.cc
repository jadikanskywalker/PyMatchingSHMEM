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

#include "pymatching/sparse_blossom/driver/parallel/decoding_unit.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <omp.h>
#include <set>
#include <unordered_map>
#include <vector>
#include <format>

// #include "profiling/profiling_json.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"

#ifdef SCOREP_USER_ENABLE
#include <scorep/SCOREP_User.h>
#endif

// #define SHMEM_MAX_SYNC_STEPS 5
// #define SHMEM_SYNC_SUMMARY 0
// #define SHMEM_SYNC_RES 1

pm::DecodingUnit::DecodingUnit(
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
    std::unique_ptr<stim::MeasureRecordWriter> writer,
    pm::UserGraph user_graph,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included,
    bool enable_correlations
#ifdef ENABLE_DRAW_FLAGS
    ,
    bool draw_frames,
    const stim::DetectorErrorModel* dem_for_drawing
#endif
    )
    : ensure_search_flooder_included(ensure_search_flooder_included),
      enable_correlations(enable_correlations)
#ifdef ENABLE_DRAW_FLAGS
      , draw_frames(draw_frames)
#endif
{
#ifdef USE_SHMEM
    n_pes = shmem_n_pes();
    pid = shmem_my_pe();
#endif
    // --- Create shared matching graph ---
#ifdef ENABLE_SHOT_BUFFERS
    const int num_shot_containers = NUM_BUFFERS_PER_UNIT;
#else
    const int num_shot_containers = 1;
#endif
#ifdef USE_SHMEM
    nodes_nelems_per_buffer = user_graph.nodes.size();
    node_ephemeral_fields_ptr = static_cast<DetectorNodeEphemeralFields*>(
        shmem_malloc(num_shot_containers * nodes_nelems_per_buffer * sizeof(DetectorNodeEphemeralFields)));
    if (node_ephemeral_fields_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric detector node buffer.");
    }
#endif
    graph = user_graph.to_shared_matching_graph(
        num_distinct_weights
#ifdef USE_SHMEM
        ,
        node_ephemeral_fields_ptr
#endif
    );
    if (graph.num_partitions <= 0) {
        throw std::invalid_argument("Graph partitioning produced no partitions. Check --rounds_per_partition.");
    }
#ifdef USE_SHMEM
    // Bound k
    int p_per_pe = graph.num_partitions / n_pes;
    if (n_pes == 2 && config_parallel::k > graph.num_partitions/2) {
        config_parallel::k = graph.num_partitions/2;
        std::cout << "NOTE: k set to " << config_parallel::k << " for 2 PEs" << std::endl << std::flush;
    } else {
        if (p_per_pe < 2) {
            throw std::invalid_argument("The number of partitions per PE should be >= 2 for more than 2 ranks.");
        } else if (config_parallel::division_strategy == config_parallel::ROUND && config_parallel::k > p_per_pe / 2) {
            // Bounding k for correctness
            config_parallel::k = p_per_pe / 2;
            std::cout << "NOTE: k bounded to " << config_parallel::k << std::endl << std::flush;
        } else if (config_parallel::division_strategy == config_parallel::OBS && config_parallel::k > graph.p_per_obs_patch) {
            // Bounding k for correctness
            config_parallel::k = graph.p_per_obs_patch;
            std::cout << "NOTE: k bounded to " << config_parallel::k << std::endl << std::flush;
        }
    }
    // --- Compute num_cross_rank_fusions before allocations ---
    if (config_parallel::division_strategy == config_parallel::OBS) {
        if (DEBUG) std::cout << "K_p: " << graph.p_per_obs_patch << "; K_vb: " << graph.vb_per_obs_patch << std::endl;
        num_cross_rank_fusions = (int)graph.num_virtual_boundaries - graph.vb_per_obs_patch * (int)graph.num_obs_patches;
    } else {
        num_cross_rank_fusions = 2;  // 2
    }
    // --- Allocate sychronization & summary memory ---
    task_status_ptr = static_cast<uint64_t*>(shmem_align(sizeof(uint64_t), num_cross_rank_fusions * SHMEM_NUM_ATOMICS_PER_CROSS_RANK_FUSION * num_shot_containers * sizeof(uint64_t)));
    regions_nelems_per_solver = graph.node_part_id.size() / graph.num_partitions * SHMEM_ARENA_BUFFER_FACTOR;
    if (regions_nelems_per_solver % 64 > 0) {
        regions_nelems_per_solver =
            (regions_nelems_per_solver / 64 + 1) * 64;  // make multiple of 64
    }
    // regions_nelems_per_solver / 64 gives number of uint64_t for bitmap
    child_edges_nelems_per_solver = regions_nelems_per_solver; // Could reduce this
    regions_matched_to_vb_nelems = graph.node_part_id.size() / graph.num_rounds * SHMEM_INTERSECTION_BUFFER_FACTOR;
    child_edges_ptr = static_cast<BlossomChild*>(shmem_malloc(child_edges_nelems_per_solver * graph.num_partitions * num_shot_containers * sizeof(BlossomChild)));
    if (child_edges_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric blossom child buffer.");
    }
    task_fusion_summary_size_per_task = sizeof(FusionSummary) +
                                        regions_matched_to_vb_nelems * sizeof(GraphFillRegion*) +
                                        (std::max(2, config_parallel::k) * regions_nelems_per_solver / 8); /* bit map in bytes (== nelems/64 * 8) */
    task_fusion_summary_ptr = static_cast<FusionSummary*>(shmem_malloc(task_fusion_summary_size_per_task * num_cross_rank_fusions * num_shot_containers));
    if (task_status_ptr == nullptr || task_fusion_summary_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric atomics buffer.");
    }
#endif
    // --- Create shot buffer ---
    shot_buffer = std::make_shared<pm::ShotBuffer>(
        std::move(reader),
        std::move(writer),
        graph.num_partitions,
        graph.num_virtual_boundaries,
        graph.graph_ptr->num_observables);
    // --- Fill buffer with shots ---
    if (DEBUG) {
        std::cout << "DEBUG: Reading shots \n" << std::flush;
    }
    for (int i = 0; i < shot_buffer->buffer.size(); ++i) {
        shot_buffer->read_shot(i, graph.node_part_id);
    }
    // --- Set num_threads ---
    if (DEBUG) {
        std::cout << "DEBUG: Setting num threads\n" << std::flush;
    }
    int max_threads = omp_get_max_threads();
    // --- Populate my_partition_task_ids (ROUND only -- OBS populates this itself, from inside
    // build_tasks_for_obs_patch_partitioning, recording each leaf's *actual* tasks[] index as it's
    // created. A fixed p_base formula here would assume every observable occupies a fixed-size block
    // in tasks[], which no longer holds once seams/CRTs/continuation fusions get interleaved into
    // tasks[] by the preemptive-extraction worklist construction -- see now-its-time-to-jazzy-
    // turing.md.) ---
    if (config_parallel::division_strategy == config_parallel::OBS) {
#ifdef USE_SHMEM
        const int base         = (int)graph.num_obs_patches / n_pes;
        const int rem          = (int)graph.num_obs_patches % n_pes;
        const int my_obs_start = base * pid + std::min(pid, rem);
        const int my_obs_count = base + (pid < rem ? 1 : 0);
#else
        const int my_obs_start = 0;
        const int my_obs_count =  graph.num_obs_patches;
#endif
        // Global partition id of my own first local partition, and how many I own -- contiguous by
        // construction (obs-major numbering, per-PE contiguous observable ownership). Distinct from
        // my_partition_task_ids, which holds indices into this PE's own tasks[] array, not global
        // partition ids -- populated later, see above.
        my_partitions_start = my_obs_start * (int)graph.p_per_obs_patch;
        my_partition_count  = my_obs_count * (int)graph.p_per_obs_patch;
        my_partition_task_ids.reserve((size_t)my_partition_count);
    } else {
#ifdef USE_SHMEM
        const int p_base   = (int)graph.num_partitions / n_pes;
        const int p_leftover = (int)graph.num_partitions % n_pes;
        my_partitions_start = p_base * pid + std::min(pid, p_leftover);
        my_partition_count = p_base + (pid < p_leftover);
#else
        my_partitions_start = 0;
        my_partition_count = graph.num_partitions;
#endif
        my_partition_task_ids.reserve(my_partition_count);
    }
    // my_partition_task_ids isn't populated yet at this point for either strategy (OBS fills it in
    // from build_tasks_for_obs_patch_partitioning, ROUND from build_tasks_for_round_partitioning, both
    // called after this constructor) -- size off my_partition_count instead, which is always equal to
    // the eventual my_partition_task_ids.size() once construction finishes.
    num_threads = std::min(max_threads, my_partition_count);
    // solvers only exist for this PE's own partitions -- see build_solvers() and the plan for this
    // refactor (a full Mwpm is unnecessary for a remote partition, only its SHMEMArena is).
    num_solvers_per_buffer = (size_t)my_partition_count;
    std::cout << "graph.num_partitions = " << graph.num_partitions << "; num_threads: " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
#ifdef USE_SHMEM
    // --- Allocate buffer for GraphFillRegion arenas ---
    // Sized globally (graph.num_partitions, not num_solvers_per_buffer/my_partition_count): the
    // underlying SHMEM region buffer must stay addressable by every PE's own global partition ids for
    // RMA puts, independent of how many local Mwpm solver objects this PE actually constructs.
    if (DEBUG) {
        std::cout << "DEBUG: Allocating Regions\n" << std::flush;
    }
    regions_ptr = static_cast<GraphFillRegion*>(
        shmem_malloc(regions_nelems_per_solver * graph.num_partitions * num_shot_containers * sizeof(GraphFillRegion)));
    if (regions_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric region buffer.");
    }
    if (DEBUG) {
        std::cout << "PE" << pid << " symmetric allocations:" << std::endl
                  << "  node_ephemeral_fields_ptr: " << node_ephemeral_fields_ptr << " (" << nodes_nelems_per_buffer * sizeof(DetectorNodeEphemeralFields) << " bytes)" << std::endl
                  << "  task_status_ptr: " << task_status_ptr << " (" << (graph.num_partitions-1) * 2 * sizeof(uint64_t) << " bytes)" << std::endl
                  << "  task_fusion_summary_ptr: " << task_fusion_summary_ptr << " (" << task_fusion_summary_size_per_task * (graph.num_partitions-1) << " bytes)" << std::endl
                  << "  regions_ptr: " << regions_ptr << " (" << regions_nelems_per_solver * (num_threads) * 2 * num_shot_containers * sizeof(GraphFillRegion) << " bytes)" << std::endl
                  << "  regions_nelems_per_solver: " << regions_nelems_per_solver << std::endl
                  << "  child_edges_ptr: " << child_edges_ptr << " (" << child_edges_nelems_per_solver * graph.num_partitions * num_shot_containers * sizeof(BlossomChild) << " bytes)" << std::endl
                  << std::flush;
    }
#endif
    // --- Build solvers (before tasks: task construction assigns each task's own ->solver, which
    // needs solvers[]/remote_arenas already built) ---
    build_solvers();
    // --- Build local merge-tree tasks and cross-rank fusions ---
    if (config_parallel::division_strategy == config_parallel::OBS) {
        build_tasks_for_obs_patch_partitioning();
    } else {
        build_tasks_for_round_partitioning();
    }
#ifdef ENABLE_DRAW_FLAGS
    if (draw_frames) {
        if (dem_for_drawing == nullptr) {
            throw std::invalid_argument(
                "--draw_frames requires the DEM to be available, but the graph was loaded from "
                "--graph_cache_path instead of parsed fresh. Remove --draw_frames or drop the graph cache.");
        }
        auto coords = pm::pick_coords_for_drawing_from_dem(*dem_for_drawing, 20);
        for (auto& s : solvers)
            s->coords = coords;
#ifdef USE_SHMEM
        if (pid == 0) {
            std::filesystem::create_directories("out_parallel/frames/0");
        }
        shmem_barrier_all();
#endif
    }
#endif
}

pm::DecodingUnit::~DecodingUnit() {
#ifdef USE_SHMEM
    shmem_free(node_ephemeral_fields_ptr);
    shmem_free(regions_ptr);
    shmem_free(child_edges_ptr);
    shmem_free(task_status_ptr);
    shmem_free(task_fusion_summary_ptr);
#endif
}

namespace {
// Shared by build_tasks_for_round_partitioning() and build_tasks_for_obs_patch_partitioning() --
// see the extraction-unit design comment at the top of build_tasks_for_round_partitioning() for the
// full rationale (unit-checkpointed extraction, is_extraction_unit_root/connector tagging).
struct RangeInfo { Task* task; size_t leftmost; size_t rightmost; };

// Pairs adjacent entries of `level`, carrying an odd one forward to the next level, until one root
// remains. Used both for each unit's own internal subtree (raw leaves in `level`,
// tag_connector=false) and for combining units into a balanced tree (unit roots in `level`,
// tag_connector=true). `on_fusion`, if given, is called right after each new fusion Task is
// constructed with (fusion, vb_id) -- vb_id is the *local* position (relative to whatever range
// `level` covers), letting a caller set fields the Task constructor doesn't default correctly on its
// own (OBS partitioning's vb_marker; round partitioning needs no callback, since its constructor
// default -- vb_marker=part -- is already correct).
//
// `solver_base`: every fusion's own solver is "the partition directly left of the fusion," i.e.
// solvers[some_p_offset + vb_id] for whatever global-partition offset is correct in the caller's own
// scope (round partitioning's own p_offset; OBS's per-observable obs_p_offset, which is *not* the
// same value as the vb-id offset OBS passes as `p_offset` above, since K_vb = K_p - 1). Rather than
// passing that offset in and making pair_up reach into DecodingUnit state to resolve it, each caller
// precomputes the starting pointer into DecodingUnit::solvers once and passes it in directly --
// solver_base[vb_id] is then plain pointer arithmetic, no lookup call needed.
RangeInfo pair_up(
    std::vector<Task>& tasks, int p_offset, std::vector<RangeInfo> level, bool tag_connector,
    const std::shared_ptr<pm::Mwpm>* solver_base,
    const std::function<void(Task*, size_t)>& on_fusion = {}) {
    while (level.size() > 1) {
        std::vector<RangeInfo> next_level;
        for (size_t k = 0; k + 1 < level.size(); k += 2) {
            const RangeInfo& L_op = level[k];
            const RangeInfo& R_op = level[k + 1];
            size_t vb_id = L_op.rightmost;
            tasks.emplace_back((int)vb_id + p_offset, L_op.task, R_op.task, solver_base[vb_id].get());
            Task* fusion = &tasks.back();
            if (tag_connector) fusion->is_extraction_unit_connector = true;
            if (on_fusion) on_fusion(fusion, vb_id);
            next_level.push_back({fusion, L_op.leftmost, R_op.rightmost});
        }
        if (level.size() % 2) next_level.push_back(level.back());
        level = std::move(next_level);
    }
    return level[0];
}

// Groups `n_leaves` consecutive leaves starting at tasks[leaf_base_idx] into L-sized units, builds
// each as its own balanced subtree (via pair_up, tag_connector=false), and tags each unit's root
// is_extraction_unit_root. Returns the unit roots, left for the caller to combine (a balanced tree
// for non-preemptive extraction, or a sequential chain for preemptive -- see
// build_tasks_for_round_partitioning for both).
std::vector<RangeInfo> build_extraction_units(
    std::vector<Task>& tasks, int p_offset, size_t leaf_base_idx, size_t n_leaves, int L,
    const std::shared_ptr<pm::Mwpm>* solver_base,
    const std::function<void(Task*, size_t)>& on_fusion = {}) {
    std::vector<RangeInfo> unit_roots;
    for (size_t ustart = 0; ustart < n_leaves; ustart += (size_t)L) {
        size_t count = std::min((size_t)L, n_leaves - ustart);
        std::vector<RangeInfo> level;
        for (size_t k = 0; k < count; ++k)
            level.push_back({&tasks[leaf_base_idx + ustart + k], ustart + k, ustart + k});
        RangeInfo unit_root = pair_up(tasks, p_offset, std::move(level), /*tag_connector=*/false,
                                       solver_base, on_fusion);
        unit_root.task->is_extraction_unit_root = true;
        unit_roots.push_back(unit_root);
    }
    return unit_roots;
}
}  // namespace

// Builds balanced fusion tree assuming round-based partitioning
void pm::DecodingUnit::build_tasks_for_round_partitioning() {
    // n_pes > graph.num_partitions can leave a PE with nothing to do -- bail before touching
    // build_extraction_units/pair_up with zero leaves.
    if (my_partition_count == 0) return;
#ifdef ENABLE_SHOT_BUFFERS
    const int num_shot_containers = NUM_BUFFERS_PER_UNIT;
#else
    const int num_shot_containers = 1;
#endif
    for (int i = 0; i < my_partition_count; ++i)
        my_partition_task_ids.push_back(i);
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // p_offset (the vb-id-to-part-space offset) and my_partitions_start (the global-partition-to-
    // solvers[]-index offset) are the same value for round partitioning -- reuse the member directly.
    const int p_offset = my_partitions_start;
    // Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md). Every unit (L
    // consecutive leaves, fewer for a ragged last unit) is built as its own small balanced subtree,
    // tagging its root is_extraction_unit_root -- this step is unconditional, regardless of
    // extract_preemptively. The units are then combined either as a sequential chain
    // (extract_preemptively: F1=fuse(U0,U1), F2=fuse(F1,U2), ...) or as a balanced tree over the units
    // (otherwise) -- every inter-unit fusion built either way is tagged is_extraction_unit_connector.
    // Task::parent (set by the Task(vb,left,right) constructor) makes each fusion's natural parent its
    // combiner, so the existing climb-to-parent logic in decode_shots() needs no special-casing.
    //
    // vb id (part, for shot.virtual_boundary_hits[] lookup) is the position-based "rightmost leaf
    // index of the left operand" -- unique across the whole structure, matching how the original
    // balanced-tree code numbers vbs by position.
    //
    // Every task's own solver is assigned once at construction (TaskBase::solver, see decoding_task.h)
    // rather than recomputed later. For ROUND partitioning, p_offset (the vb-id-to-part-space offset)
    // and my_partitions_start (the global-partition-to-solvers[]-index offset) are the same value, so
    // solver_base -- the pointer into solvers[] for this shot_container_id's own local block --
    // already points at "solver for global partition p_offset + 0", and solver_base[vb_id] is exactly
    // solvers[]'s entry for the fusion's own vb id. (OBS partitioning's per-observable tree needs its
    // own, separately-derived solver_base since its p_offset there is a *vb*-space offset, not a
    // partition-space one -- see build_tasks_for_obs_patch_partitioning.)
    //
    // This also matters for correctness, not just indexing: since a fusion's own solver must be *some
    // partition within its own subtree* for Arena safety (see GraphFlooder::create_blossom), always
    // resolving to the vb's own position (rather than, e.g., always collapsing back to the chain's
    // original leftmost leaf) keeps a preemptive chain's later links from reusing an earlier link's
    // solver after that earlier unit has already been divided off and posted for concurrent
    // extraction.
    // Build a full fusion tree: less than 2*N tasks. Reserve to keep element addresses stable.
    for (int shot_container_id=0; shot_container_id < num_shot_containers; shot_container_id++) {  // Lazy repeat per shot container -- CAN BE IMPROVED!
        auto& shot_container = shot_buffer->buffer[shot_container_id];
        auto& tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * graph.num_partitions - 1));
        const std::shared_ptr<pm::Mwpm>* solver_base = &solvers[num_solvers_per_buffer * shot_container_id];
        for (int task_id : my_partition_task_ids) {
            tasks.emplace_back(task_id + p_offset, solver_base[task_id].get());
        }

        int L = config_parallel::L;
        size_t n_leaves = my_partition_task_ids.size();
        std::vector<RangeInfo> unit_roots =
            build_extraction_units(tasks, p_offset, /*leaf_base_idx=*/0, n_leaves, L, solver_base);

        if (config_parallel::extract_preemptively) {
            // Sequential chain over unit_roots.
            RangeInfo chain = unit_roots[0];
            for (size_t k = 1; k < unit_roots.size(); ++k) {
                const RangeInfo& next_unit = unit_roots[k];
                size_t vb_id = chain.rightmost;
                tasks.emplace_back((int)vb_id + p_offset, chain.task, next_unit.task, solver_base[vb_id].get());
                Task* fusion = &tasks.back();
                fusion->is_extraction_unit_connector = true;
                if (k > 1) {
                    // chain.task is the previous chain-link fusion, not a raw leaf: its own vb
                    // (chain.task->part) was already divided and its left/older unit posted for
                    // (possibly concurrent) extraction by the time THIS fusion resolves. Restrict
                    // vb_left so this fusion's flooding can never grow back into that now-separate
                    // unit -- the Task(vb, left, right) constructor otherwise inherits
                    // left_child->vb_left unconditionally, which for a chain link still spans all
                    // the way back to the very first unit's own original vb_left.
                    fusion->vb_left = chain.task->part;
                }
                chain = {fusion, chain.leftmost, next_unit.rightmost};
            }
            // If only one unit exists, no connectors were built -- unit_roots[0].task is already
            // tagged is_extraction_unit_root, and root detection's degenerate-case branch handles it.
        } else {
            // Balanced tree over unit_roots -- same pair_up algorithm, one level higher. If only one
            // unit exists, pair_up is a no-op (returns unit_roots[0] unchanged, no task created).
            pair_up(tasks, p_offset, unit_roots, /*tag_connector=*/true, solver_base);
        }
    }
    if (DEBUG) {
        std::cout << "DEBUG:" << std::endl;
        for (auto& buffer : shot_buffer->buffer) {
            std::cout << "Buffer" << std::endl;
            for (Task& t : buffer.tasks) {
                std::cout << "--addr: " << &t << std::endl
                        << "  part: " << t.part << std::endl
                        << "  vb_left: " << t.vb_left << std::endl
                        << "  vb_right: " << t.vb_right << std::endl
                        << "  is_fusion: " << t.is_fusion << std::endl
                        << "  is_extraction_unit_connector: " << t.is_extraction_unit_connector << std::endl
                        << "  is_extraction_unit_root: " << t.is_extraction_unit_root << std::endl
                        << "  solver: " << t.solver << std::endl
                        << "  left_child: " << t.left_child << std::endl
                        << "  right_child: " << t.right_child << std::endl
                        << "  parent: " << t.parent << std::endl;
            }
        }
    }
    // Completion-counting init: universal across build configs. ROUND partitioning always has
    // exactly one true tree root regardless of extract_preemptively/L -- checkpoint chain-links are
    // never independent roots (they're never pushed as a root, never set i_solved_root), so
    // num_roots_done (only incremented at genuine parent==nullptr detection) can only ever reach 1,
    // not 1+num_checkpoints. Checkpoint/chunk completion is tracked entirely via
    // pending_extraction_jobs instead (see shot_buffer.h).
    for (int i = 0; i < num_shot_containers; ++i) {
        shot_buffer->buffer[i].num_task_roots = 1;
        shot_buffer->buffer[i].thread_results.assign(num_threads, pm::MatchingResult{});
        shot_buffer->buffer[i].num_roots_done.store(0, std::memory_order_relaxed);
    }
#ifdef USE_SHMEM
    // --- Build cross-rank fusions (ROUND topology) ---
    //   THIS IS ALL BASED ON SIMPLE ROUND BASED FUSION ACROSS PEs
    int my_partitions_end = my_partitions_start + my_partition_count;
    for (int i=0; i < num_shot_containers; ++i) {
        shot_buffer->buffer[i].cross_rank_tasks.reserve(2);
        if (pid > 0) { // Add cross-rank fusion on left
            uint64_t* task_status_p = get_task_status_ptr(i, false);
            FusionSummary* fusion_summary_p = get_fusion_summary_ptr(i, false);
            int vb = my_partitions_start-1;
            // Anchor on my own first local partition, right at this boundary -- same "a partition
            // this CRT touches, on my own local side" pattern as OBS's CRT anchor.
            shot_buffer->buffer[i].cross_rank_tasks.emplace_back(
                vb,
                false,
                (vb-config_parallel::k < -1) ? -1 : vb-config_parallel::k,
                (vb+config_parallel::k < (int)graph.num_virtual_boundaries) ? vb+config_parallel::k : (int)graph.num_virtual_boundaries,
                pid-1,
                task_status_p,     // status_shm
                task_status_p + 1, // signal_shm
                task_status_p + 2, // done_shm
                fusion_summary_p,
                solver_for(i, my_partitions_start)
            );
            shot_buffer->buffer[i].tasks.back().add_special_task(&shot_buffer->buffer[i].cross_rank_tasks.back());
            if (DEBUG) {
                auto& t = shot_buffer->buffer[i].cross_rank_tasks.back();
                std::cout << "PE" << pid << " Left cross task:" << std::endl
                        << "    vb: " << t.part << std::endl
                        << "    child: " << t.child << std::endl
                        << "    iamleft: " << t.iamleft << std::endl
                        << "    vb_left: " << t.vb_left << std::endl
                        << "    vb_right: " << t.vb_right << std::endl
                        << "    other_pid: " << t.other_pid << std::endl
                        << "    status_shm: " << t.status_shm << std::endl
                        << "    signal_shm: " << t.signal_shm << std::endl
                        << "    fusion_summary_shm: " << t.fusion_summary_shm << std::endl << std::flush;
            }
        }
        if (pid < n_pes-1) { // Add cross-rank fusion on right
            uint64_t* task_status_p = get_task_status_ptr(i, true);
            FusionSummary* fusion_summary_p = get_fusion_summary_ptr(i, true);
            int vb = my_partitions_end-1;
            // Anchor on my own last local partition, right at this boundary.
            shot_buffer->buffer[i].cross_rank_tasks.emplace_back(
                vb,
                true,
                (vb-config_parallel::k < -1) ? -1 : vb-config_parallel::k,
                (vb+config_parallel::k < (int)graph.num_virtual_boundaries) ? vb+config_parallel::k : (int)graph.num_virtual_boundaries,
                pid+1,
                task_status_p,     // status_shm
                task_status_p + 1, // signal_shm
                task_status_p + 2, // done_shm
                fusion_summary_p,
                solver_for(i, my_partitions_end - 1)
            );
            shot_buffer->buffer[i].tasks.back().add_special_task(&shot_buffer->buffer[i].cross_rank_tasks.back());
            if (DEBUG) {
                auto& t = shot_buffer->buffer[i].cross_rank_tasks.back();
                std::cout << "PE" << pid << "  Right cross task:" << std::endl
                    << "    vb: " << t.part << std::endl
                    << "    child: " << t.child << std::endl
                    << "    iamleft: " << t.iamleft << std::endl
                    << "    vb_left: " << t.vb_left << std::endl
                    << "    vb_right: " << t.vb_right << std::endl
                    << "    other_pid: " << t.other_pid << std::endl
                    << "    status_shm: " << t.status_shm << std::endl
                    << "    signal_shm: " << t.signal_shm << std::endl
                    << "    done_shm: " << t.done_shm << std::endl
                    << "    fusion_summary_shm: " << t.fusion_summary_shm << std::endl << std::flush;
            }
        }
    }
#endif
}

void pm::DecodingUnit::build_tasks_for_obs_patch_partitioning() {
    if (DEBUG) std::cout << "DEBUG: initializing obs-patch tasks" << std::endl;

    const int K_p       = graph.p_per_obs_patch;
    const int K_vb      = graph.vb_per_obs_patch;
    const int num_seams = (int)graph.num_virtual_boundaries - K_vb * (int)graph.num_obs_patches;

    // Determine my obs patch range (same formula as constructor)
#ifdef USE_SHMEM
    const int base         = (int)graph.num_obs_patches / n_pes;
    const int rem          = (int)graph.num_obs_patches % n_pes;
    const int my_obs_start = base * pid + std::min(pid, rem);
    const int my_obs_count = base + (pid < rem ? 1 : 0);
#else
    const int my_obs_start = 0;
    const int my_obs_count = graph.num_obs_patches;
    // Every existing "PE" << pid DEBUG print (both here and in the CRT-only branches below) compiles
    // unchanged in a non-SHMEM build with this shadow -- there's only ever one implicit "PE" (0), which
    // owns every observable.
    const int pid = 0;
#endif
    if (DEBUG) std::cout << "my_obs_start: " << my_obs_start << ", n=" << my_obs_count << "\n" << std::flush;

    // For each seam: determine the two obs patches it connects and the local partition
    // range it touches (used for vb_left/vb_right on cross-PE fusion tasks).
    // unit_lo/unit_hi (populated below, after window conflict resolution) are this seam's extraction-
    // unit span for preemptive extraction -- see now-its-time-to-jazzy-turing.md Design §1.
    struct SeamInfo { int oi, oj, vb_left, vb_right; int unit_lo{-1}, unit_hi{-1}; };
    std::vector<SeamInfo> seam_infos(num_seams);
    for (int s = 0; s < num_seams; ++s) {
        auto [first, last] = graph.vb_bounds[K_vb * graph.num_obs_patches + s];
        std::map<int, std::pair<int,int>> obs_part_range;  // obs_id → (min_local_p, max_local_p)
#ifdef DEBUG
        std::string p_string;
#endif
        for (size_t ni = first; ni <= last; ++ni) {
            for (auto* nbr : graph.graph_ptr->nodes[ni].neighbors) {
                if (nbr == nullptr) continue; // boundary
                int nbr_idx = (int)(nbr - graph.graph_ptr->nodes.data());
                int part_id = graph.node_part_id[nbr_idx];
                int obs_id, low_p, high_p;
                if (part_id >= 0) { // partition node
                    obs_id = part_id / K_p;
                    low_p = part_id % K_p;
                    high_p = low_p;
                    if (DEBUG) p_string += "  obs" + std::to_string(obs_id);
                } else {
                    if (K_vb == 0) continue;  // no intra-obs VBs possible
                    int global_vb = -part_id - 1;
                    int obs_id = global_vb / K_vb;
                    if (obs_id >= (int)graph.num_obs_patches) continue;  // cross-obs seam node
                    int local_vb = global_vb % K_vb;
                    low_p  = local_vb - 1;
                    high_p = local_vb;
                    if (DEBUG) p_string += "  obs" + std::to_string(obs_id);
                }
                auto [it, inserted] = obs_part_range.try_emplace(obs_id, low_p, high_p);
                if (!inserted) {
                    it->second.first  = std::min(it->second.first,  low_p);
                    it->second.second = std::max(it->second.second, high_p);
                }
            }
        }

        auto it = obs_part_range.begin();
        int oi    = it->first;
        int vb_left = std::max(it->second.first - config_parallel::k - 1, -1);
        int vb_right = std::min(it->second.second + config_parallel::k, K_p-1);
        ++it;
        int oj = (it != obs_part_range.end()) ? it->first : oi;
        seam_infos[s] = { oi, oj, vb_left, vb_right };
        if (DEBUG) std::cout << "PE" << pid << " seam " << s << ": obs " << oi << " -- obs " << oj
                             << " node_bounds [" << first << ", " << last << "]\n"
                             << " vb_left=" << vb_left << " vb_right=" << vb_right << "\n    "
                             << p_string << "\n" 
                             << std::flush;
    }

    // --- Global seam window conflict resolution ---
    // All PEs compute this identically (same seam_infos, same iteration order)
    // so every PE ends up with the same adjusted windows.
    for (int si_idx = 0; si_idx < num_seams; ++si_idx) {
        auto& si = seam_infos[si_idx];
        for (int sj_idx = si_idx + 1; sj_idx < num_seams; ++sj_idx) {
            auto& sj = seam_infos[sj_idx];
            if (si.oi != sj.oi && si.oi != sj.oj &&
                si.oj != sj.oi && si.oj != sj.oj)
                continue;
            if (si.vb_left <= sj.vb_left) {
                if (si.vb_right > sj.vb_left) {
                    int mid = (si.vb_right + sj.vb_left + 1) / 2;
                    if (DEBUG) std::cout << "  Global adjust: seam " << si_idx
                        << " vb_right " << si.vb_right << " -> " << mid
                        << ", seam " << sj_idx
                        << " vb_left " << sj.vb_left << " -> " << mid << "\n" << std::flush;
                    si.vb_right = mid;
                    sj.vb_left  = mid;
                }
            } else {
                if (sj.vb_right > si.vb_left) {
                    int mid = (sj.vb_right + si.vb_left + 1) / 2;
                    if (DEBUG) std::cout << "  Global adjust: seam " << sj_idx
                        << " vb_right " << sj.vb_right << " -> " << mid
                        << ", seam " << si_idx
                        << " vb_left " << si.vb_left << " -> " << mid << "\n" << std::flush;
                    sj.vb_right = mid;
                    si.vb_left  = mid;
                }
            }
        }
    }

    // --- Extraction-unit span per seam (preemptive extraction only; harmless to always compute) ---
    // [si.vb_left+1, si.vb_right] (already k-padded, symmetric across both oi's and oj's own local
    // partition numbering -- see now-its-time-to-jazzy-turing.md Design §1) maps directly to the range
    // of L-sized extraction units this seam's deferred span must cover, using the same ustart/L
    // grouping build_extraction_units itself uses.
    const int L = config_parallel::L;
    for (auto& si : seam_infos) {
        si.unit_lo = (si.vb_left + 1) / L;
        si.unit_hi = si.vb_right / L;
    }
    if (DEBUG) {
        for (size_t s = 0; s < seam_infos.size(); ++s) {
            auto& si = seam_infos[s];
            std::cout << "PE" << pid << " seam_infos[" << s << "]: oi=" << si.oi << " oj=" << si.oj
                      << " vb_left=" << si.vb_left << " vb_right=" << si.vb_right
                      << " unit_lo=" << si.unit_lo << " unit_hi=" << si.unit_hi << "\n" << std::flush;
        }
    }
    // Which observable(s) of each seam are local to this PE -- computed once, reused every shot_id
    // (identical graph structure every shot) and by both the preemptive and non-preemptive branches.
    struct SeamOwnership { bool oi_local, oj_local; };
    std::vector<SeamOwnership> seam_ownership(num_seams);
    for (int s = 0; s < num_seams; ++s) {
        const auto& si = seam_infos[s];
        const int loi = si.oi - my_obs_start, loj = si.oj - my_obs_start;
        seam_ownership[s] = { (loi >= 0 && loi < my_obs_count), (loj >= 0 && loj < my_obs_count) };
    }
    // No overlap restriction remains: any number of seams (local or cross-rank, in any combination) may
    // share an extraction unit -- see now-its-time-to-jazzy-turing.md Phase 1.6, which removed the
    // narrower same-unit_hi-both-CRT assertion Phase 1.5 needed (CRT construction no longer commits
    // immediately, so it can stack/defer exactly like a local seam).

#ifdef USE_SHMEM
    const int base_pe = (int)graph.num_obs_patches / n_pes;
    const int rem_pe  = (int)graph.num_obs_patches % n_pes;
    auto other_pid_for = [&](int remote_obs) {
        if (remote_obs < rem_pe * (base_pe + 1)) return remote_obs / (base_pe + 1);
        return rem_pe + (remote_obs - rem_pe * (base_pe + 1)) / base_pe;
    };
#endif

#ifdef ENABLE_SHOT_BUFFERS
    const int num_shot_containers = NUM_BUFFERS_PER_UNIT;
#else
    const int num_shot_containers = 1;
#endif
    for (int shot_id = 0; shot_id < num_shot_containers; ++shot_id) {
        auto& tasks = shot_buffer->buffer[shot_id].tasks;
        auto& local_seams = shot_buffer->buffer[shot_id].local_seam_tasks;
        tasks.reserve(static_cast<size_t>(my_obs_count * (2*K_p - 1) + num_seams + 1));
        local_seams.reserve((size_t)num_seams);
#ifdef USE_SHMEM
        auto& crt = shot_buffer->buffer[shot_id].cross_rank_tasks;
        crt.reserve((size_t)num_seams);
#endif

        if (!config_parallel::extract_preemptively) {
            // ==================== Non-preemptive ====================
            // Each observable becomes one full balanced tree; local seams/CRTs attach to whole-
            // observable roots via add_special_task instead of fusing them into the tree.
            // post_hoc_chunk_and_post (unmodified) does all the actual unit-by-unit division+posting
            // at extraction time, walking this static tree top-down.
            std::vector<Task*> obs_roots(my_obs_count, nullptr);
            for (int lo = 0; lo < my_obs_count; ++lo) {
                const int o = my_obs_start + lo;
                const int obs_p_offset = o * K_p;
                const int obs_vb_offset = o * K_vb;
                const int tree_start = lo * (2*K_p - 1);
                const std::shared_ptr<pm::Mwpm>* solver_base =
                    &solvers[num_solvers_per_buffer * shot_id + (obs_p_offset - my_partitions_start)];
                for (int lp = 0; lp < K_p; ++lp) {
                    if (shot_id == 0) my_partition_task_ids.push_back(tasks.size());
                    tasks.emplace_back(obs_p_offset + lp, lp - 1, lp, solver_base[lp].get());
                }
                auto on_fusion = [&](Task* fusion, size_t vb_id) { fusion->vb_marker = (int)vb_id; };
                std::vector<RangeInfo> unit_roots = build_extraction_units(
                    tasks, obs_vb_offset, /*leaf_base_idx=*/tree_start, (size_t)K_p,
                    config_parallel::L, solver_base, on_fusion);
                RangeInfo obs_root = pair_up(tasks, obs_vb_offset, unit_roots, /*tag_connector=*/true,
                                              solver_base, on_fusion);
                obs_roots[lo] = obs_root.task;
            }
            if (DEBUG) {
                for (int lo = 0; lo < my_obs_count; ++lo) std::cout << "  " << obs_roots[lo];
                std::cout << "\n" << std::flush;
            }

            // Roots never merge under the special-tasks-attach model -- obs_roots[lo] stays this
            // observable's own tree root for the whole function's lifetime, so seams/CRTs attach
            // directly against it (no union-find over a group_root[] indirection needed anymore).
            for (int s = 0; s < num_seams; ++s) {
                const auto& si = seam_infos[s];
                const int loi = si.oi - my_obs_start, loj = si.oj - my_obs_start;
                const bool oi_local = seam_ownership[s].oi_local, oj_local = seam_ownership[s].oj_local;
                if (oi_local && oj_local) {
                    int global_vb = K_vb * (int)graph.num_obs_patches + s;
                    pm::Mwpm* seam_solver = solver_for(shot_id, si.oi * K_p + si.vb_left + 1);
                    auto& seam = local_seams.emplace_back(global_vb, si.vb_left, si.vb_right, seam_solver);
#ifdef ENABLE_DRAW_FLAGS
                    seam.obs_patch_ids = { si.oi, si.oj };
#endif
                    obs_roots[loi]->add_special_task(&seam);
                    obs_roots[loj]->add_special_task(&seam);
                }
#ifdef USE_SHMEM
                else if (oi_local || oj_local) {
                    const int local_lo = oi_local ? loi : loj;
                    Task* local_root = obs_roots[local_lo];
                    const bool iamleft = oi_local;
                    const int remote_obs = oi_local ? si.oj : si.oi;
                    int other_pid = other_pid_for(remote_obs);
                    int global_vb = K_vb * (int)graph.num_obs_patches + s;
                    uint64_t* sp = get_task_status_ptr_for_seam(shot_id, s);
                    FusionSummary* fp = get_fusion_summary_ptr_for_seam(shot_id, s);
                    pm::Mwpm* crt_solver = solver_for(shot_id, (oi_local ? si.oi : si.oj) * K_p + si.vb_left + 1);
                    crt.emplace_back(global_vb, iamleft, si.vb_left, si.vb_right,
                                      other_pid, sp, sp+1, sp+2, fp, crt_solver);
                    crt.back().left_global_offset  = { (size_t)si.oi * K_p, (size_t)si.oi * K_vb };
                    crt.back().right_global_offset = { (size_t)si.oj * K_p, (size_t)si.oj * K_vb };
                    local_root->add_special_task(&crt.back());
                    if (DEBUG) {
                        auto& t = crt.back();
                        std::cout << "PE" << pid << " OBS cross-rank task s=" << s << ":\n"
                                  << "    obs_a: " << si.oi << "  obs_b: " << si.oj << "\n"
                                  << "    vb=" << t.part << " vb_left=" << t.vb_left << " vb_right=" << t.vb_right << "\n"
                                  << "    iamleft=" << t.iamleft << " other_pid=" << t.other_pid << "\n"
                                  << "    child: " << t.child << " (part=" << t.child->part << ")  me: " << &t << "\n"
                                  << std::flush;
                    }
                }
#endif
            }
            shot_buffer->buffer[shot_id].num_task_roots = my_obs_count;
        } else {
            // ==================== Preemptive: build first, attach second ====================
            // See now-its-time-to-jazzy-turing.md's task-building plan. Each observable's own chain is
            // built completely, independently, with no seam awareness at all (Step 1) -- there is no
            // longer any reason one observable's chain construction needs to know about another's, since
            // seams/CRTs attach via add_special_task after the fact rather than being woven into the
            // tree's own pointers as they're built. Step 2a then marks which chain-links fall inside an
            // active seam's span, Step 2b narrows vb_left accordingly, and Step 2c builds and attaches
            // each seam/CRT.
            std::vector<std::vector<Task*>> chain_link_at_unit(my_obs_count);
            std::vector<Task*> obs_roots(my_obs_count, nullptr);

            // ---- Step 1: build every observable's own chain, independently ----
            for (int lo = 0; lo < my_obs_count; ++lo) {
                const int o = my_obs_start + lo;
                const int obs_p_offset = o * K_p;
                const int obs_vb_offset = o * K_vb;
                const std::shared_ptr<pm::Mwpm>* solver_base =
                    &solvers[num_solvers_per_buffer * shot_id + (obs_p_offset - my_partitions_start)];
                for (int lp = 0; lp < K_p; ++lp) {
                    if (shot_id == 0) my_partition_task_ids.push_back(tasks.size());
                    tasks.emplace_back(obs_p_offset + lp, lp - 1, lp, solver_base[lp].get());
                }
                auto on_fusion = [&](Task* fusion, size_t vb_id) { fusion->vb_marker = (int)vb_id; };
                size_t leaf_base_idx = tasks.size() - (size_t)K_p;
                std::vector<RangeInfo> unit_roots = build_extraction_units(
                    tasks, obs_vb_offset, leaf_base_idx, (size_t)K_p, config_parallel::L,
                    solver_base, on_fusion);

                auto& links = chain_link_at_unit[lo];
                links.resize(unit_roots.size());
                links[0] = unit_roots[0].task;
                RangeInfo chain = unit_roots[0];
                for (size_t k = 1; k < unit_roots.size(); ++k) {
                    const RangeInfo& next_unit = unit_roots[k];
                    size_t vb_id = chain.rightmost;
                    tasks.emplace_back((int)vb_id + obs_vb_offset, chain.task, next_unit.task, solver_base[vb_id].get());
                    Task* fusion = &tasks.back();
                    fusion->is_extraction_unit_connector = true;
                    fusion->vb_marker = (int)vb_id;
                    // vb_left narrowing is NOT applied here -- deferred to step 2b, which needs every
                    // seam's defer_division marks in place first (step 2a).
                    links[k] = fusion;
                    chain = {fusion, chain.leftmost, next_unit.rightmost};
                }
                obs_roots[lo] = chain.task;
            }

            // ---- Step 2a: mark defer_division strictly between unit_lo and unit_hi (both excluded) ----
            // Neither the connector closing unit_lo nor the one closing unit_hi is deferred: both
            // divide normally when solved -- unit_hi's own divide just happens after its special_tasks
            // (attached in step 2c) are processed first, a decode-time (phase three) ordering concern,
            // not something this flag itself needs to encode.
            auto mark_deferred = [&](int lo, int unit_lo, int unit_hi) {
                auto& links = chain_link_at_unit[lo];
                for (int u = unit_lo + 1; u <= unit_hi - 1 && u < (int)links.size(); ++u)
                    links[u]->defer_division = true;
            };
            // A CRT's own window (unlike a local seam's) gets *sent to, and destroyed on the sender's
            // own side after sending by, the remote PE* -- so unlike a local seam, it also needs a
            // buffer partition between its own window and its attachment connector's own right edge
            // (that connector's own vb_right, never narrowed, always equals its own unit's natural
            // right edge). Without one, the connector's own right-edge partition (needed intact by the
            // *next* connector, which hasn't resolved yet) is exactly what the CRT would send/destroy.
            // The left edge needs no equivalent check: whatever a connector's own vb_left anchors to
            // was already safely divided by an earlier, already-resolved connector (see step 2b) --
            // only the right edge is still "live," shared with what comes after.
            std::vector<int> crt_attach_unit(num_seams, -1);
            for (int s = 0; s < num_seams; ++s) {
                const auto& si = seam_infos[s];
                if (seam_ownership[s].oi_local) mark_deferred(si.oi - my_obs_start, si.unit_lo, si.unit_hi);
                if (seam_ownership[s].oj_local) mark_deferred(si.oj - my_obs_start, si.unit_lo, si.unit_hi);
#ifdef USE_SHMEM
                if (seam_ownership[s].oi_local != seam_ownership[s].oj_local) {  // exactly one side local
                    const int local_lo = seam_ownership[s].oi_local ? (si.oi - my_obs_start) : (si.oj - my_obs_start);
                    auto& links = chain_link_at_unit[local_lo];
                    int u = si.unit_hi;
                    while (u + 1 < (int)links.size() && links[u]->vb_right <= si.vb_right) {
                        links[u]->defer_division = true;
                        ++u;
                    }
                    crt_attach_unit[s] = u;
                }
#endif
            }

            // ---- Step 2b: vb_left narrowing via a running "last-divided anchor" per observable ----
            // Plain "narrow to my immediate predecessor's .part" is only correct between two connectors
            // that both actually divide. Inside (and at the far end of) a deferred span, nothing has
            // actually been divided since the connector that closed unit_lo, so every connector from
            // unit_lo+1 through unit_hi inclusive -- deferred or not -- anchors back to that same .part,
            // not to whichever connector immediately precedes it. This single forward pass reproduces
            // ordinary one-step-back narrowing exactly wherever there's no deferred span, and handles a
            // deferred span (and the unit_lo==0 edge case, where no anchor is ever established) with no
            // extra case-work.
            for (int lo = 0; lo < my_obs_count; ++lo) {
                auto& links = chain_link_at_unit[lo];
                int anchor = 0;
                bool have_anchor = false;
                for (size_t u = 1; u < links.size(); ++u) {
                    if (have_anchor) links[u]->vb_left = anchor;
                    if (!links[u]->defer_division) {
                        anchor = links[u]->part;
                        have_anchor = true;
                    }
                }
            }

            // ---- Step 2c: build and attach, now that defer_division and vb_left are both final ----
            // The special task's own vb_left/vb_right are read straight off the connector it attaches
            // to (chain_link_at_unit[...][unit_hi]) -- that connector's own (now-final) vb_left is
            // exactly the anchor value step 2b computed, and its vb_right was never touched by
            // narrowing in the first place. Both sides' own chain_link_at_unit arrays hold identical
            // values at the same unit index for a local-local seam (same L/K_p/leaf numbering and
            // identical defer_division patterns on both sides, per Phase 1 Design §1's symmetry), so
            // either side can be read -- use oi's.
            for (int s = 0; s < num_seams; ++s) {
                const auto& si = seam_infos[s];
                const bool oi_local = seam_ownership[s].oi_local, oj_local = seam_ownership[s].oj_local;
                if (oi_local && oj_local) {
                    const int loi = si.oi - my_obs_start, loj = si.oj - my_obs_start;
                    Task* attach_i = chain_link_at_unit[loi][si.unit_hi];
                    Task* attach_j = chain_link_at_unit[loj][si.unit_hi];
                    int global_vb = K_vb * (int)graph.num_obs_patches + s;
                    pm::Mwpm* seam_solver = solver_for(shot_id, si.oi * K_p + si.vb_left + 1);
                    auto& seam = local_seams.emplace_back(
                        global_vb, attach_i->vb_left, attach_i->vb_right, seam_solver);
#ifdef ENABLE_DRAW_FLAGS
                    seam.obs_patch_ids = { si.oi, si.oj };
#endif
                    attach_i->add_special_task(&seam);
                    attach_j->add_special_task(&seam);
                }
#ifdef USE_SHMEM
                else if (oi_local || oj_local) {
                    const int local_lo = oi_local ? (si.oi - my_obs_start) : (si.oj - my_obs_start);
                    // Not necessarily si.unit_hi -- may have been pushed forward in step 2a to keep a
                    // buffer partition between the CRT's own window and its attachment's right edge.
                    Task* attach = chain_link_at_unit[local_lo][crt_attach_unit[s]];
                    const bool iamleft = oi_local;
                    const int remote_obs = oi_local ? si.oj : si.oi;
                    const int other_pid = other_pid_for(remote_obs);
                    const int global_vb = K_vb * (int)graph.num_obs_patches + s;
                    uint64_t* sp = get_task_status_ptr_for_seam(shot_id, s);
                    FusionSummary* fp = get_fusion_summary_ptr_for_seam(shot_id, s);
                    pm::Mwpm* crt_solver = solver_for(shot_id, (oi_local ? si.oi : si.oj) * K_p + si.vb_left + 1);
                    // Unlike a local seam (which inherits its attachment connector's own vb_left/
                    // vb_right -- the maximum range of already-accumulated, undivided graph it's safe
                    // to use), a CRT must stay narrowly windowed to si.vb_left/si.vb_right: this is
                    // what actually gets sent to (and, on the sender's own side, destroyed after
                    // sending by) the remote PE. Inheriting the connector's own broader span would
                    // send/destroy far more than the seam's own boundary window -- fatally so when two
                    // CRTs share the same attachment connector (a coarser --extraction_unit_size than
                    // the seam spacing): the first one to send would destroy the *other* CRT's own
                    // still-needed region too, leaving it nothing to send.
                    crt.emplace_back(global_vb, iamleft, si.vb_left, si.vb_right,
                                      other_pid, sp, sp+1, sp+2, fp, crt_solver);
                    crt.back().left_global_offset  = { (size_t)si.oi * K_p, (size_t)si.oi * K_vb };
                    crt.back().right_global_offset = { (size_t)si.oj * K_p, (size_t)si.oj * K_vb };
                    attach->add_special_task(&crt.back());
                }
#endif
            }
            shot_buffer->buffer[shot_id].num_task_roots = my_obs_count;
        }

        // Every locally-owned observable's own tree/chain has exactly one top under the special-tasks-
        // attach model (seams/CRTs attach beside the tree, never replacing a root), so num_task_roots is
        // set directly (my_obs_count, both branches above) -- no scan needed.
        {
            auto& sc = shot_buffer->buffer[shot_id];
            sc.thread_results.assign(num_threads, pm::MatchingResult{});
            sc.num_roots_done.store(0, std::memory_order_relaxed);
            if (DEBUG) std::cout << "PE" << pid << " num_task_roots=" << sc.num_task_roots << "\n" << std::flush;
        }
    }

    if (DEBUG) {
        std::string tasks = "DEBUG obs-patch tasks:\n";
        for (auto& buffer : shot_buffer->buffer) {
            tasks += "Buffer\n";
            for (Task& t : buffer.tasks) {
                tasks += "--part: " + (std::string)((t.is_fusion) ? "f" : "p") + std::to_string(t.part)
                       + "  vb_left: " + std::to_string(t.vb_left)
                       + "  vb_right: " + std::to_string(t.vb_right)
                       + "  vb_marker: " + std::to_string(t.vb_marker)
                       + "  solver: " + std::format("{:p}", static_cast<void*>(t.solver)) + "\n"
                       + "    is_extraction_unit_root: " + std::to_string(t.is_extraction_unit_root)
                       + "  is_extraction_unit_connector: " + std::to_string(t.is_extraction_unit_connector)
                       + "  defer_division: " + std::to_string(t.defer_division) + "\n"
                       + "    left_child: " + std::format("{:p}",static_cast<void*>(t.left_child)) + ((t.left_child) ? "(" + (std::string)(t.left_child->is_fusion ? "f" : "p") + std::to_string(t.left_child->part) + ")" : "")
                       + "  me: " + std::format("{:p}", static_cast<void*>(&t))
                       + "  right_child: " + std::format("{:p}", static_cast<void*>(t.right_child)) + ((t.right_child) ? "(" + (std::string)(t.right_child->is_fusion ? "f" : "p") + std::to_string(t.right_child->part) + ")" : "") + "\n"
                       + "    parent: f" + std::format("{:p}", static_cast<void*>(t.parent))
                       + "\n";
                tasks += "    special_tasks (" + std::to_string(t.special_tasks.size()) + "):";
                for (SpecialTask* st : t.special_tasks) {
                    tasks += "  " + std::format("{:p}", static_cast<void*>(st))
                           + "(" + (std::string)(st->is_local_seam_fusion() ? "seam" : "crt")
                           + " vb=" + std::to_string(st->part)
                           + " vb_left=" + std::to_string(st->vb_left)
                           + " vb_right=" + std::to_string(st->vb_right) + ")";
                }
                tasks += "\n";
            }
        }
        std::cout << tasks << std::flush;

        // Verification for now-its-time-to-jazzy-turing.md Verification step 3: every entry must be a
        // valid, in-range leaf (never a fusion) index into buffer[0].tasks, covering each leaf exactly
        // once. Only populated once (shot_id == 0), so only buffer 0 is meaningful here.
        std::string pt = "PE" + std::to_string(pid) + " my_partition_task_ids (n="
            + std::to_string(my_partition_task_ids.size()) + ", my_partition_count="
            + std::to_string(my_partition_count) + "):\n";
        auto& buf0 = shot_buffer->buffer[0];
        std::set<size_t> seen;
        bool all_leaves = true, all_unique = true, all_in_range = true;
        for (size_t idx : my_partition_task_ids) {
            pt += " " + std::to_string(idx);
            if (idx >= buf0.tasks.size()) { all_in_range = false; continue; }
            if (buf0.tasks[idx].is_fusion) all_leaves = false;
            if (!seen.insert(idx).second) all_unique = false;
        }
        pt += "\n  all_leaves=" + std::to_string(all_leaves)
            + " all_unique=" + std::to_string(all_unique)
            + " all_in_range=" + std::to_string(all_in_range) + "\n";
        std::cout << pt << std::flush;
    }
}

// Build a full solver for each of this PE's own local partitions, for each shot container. A full
// Mwpm/GraphFlooder is unnecessary for a remote partition -- only its SHMEMArena (bitmap + buffer
// pointer) is, see the remote_arenas block below and region_arena_for().
void pm::DecodingUnit::build_solvers() {
    if (ensure_search_flooder_included || enable_correlations) {
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    }
#ifdef ENABLE_SHOT_BUFFERS
    const int num_shot_containers = NUM_BUFFERS_PER_UNIT;
#else
    const int num_shot_containers = 1;
#endif
    solvers.clear();
    solvers.reserve(static_cast<size_t>(num_solvers_per_buffer * num_shot_containers));
    if (DEBUG) std::cout << "num_solvers_per_buffer: " << num_solvers_per_buffer << std::endl << std::flush;
    for (size_t idx = 0; idx < num_shot_containers; ++idx) {
        for (size_t t = 0; t < num_solvers_per_buffer; ++t) {
            int global_partition = my_partitions_start + (int)t;
            // Each solver shares the same MatchingGraph via shared_ptr.
#ifdef USE_SHMEM
            if (DEBUG) std::cout << "solver: " << idx*num_solvers_per_buffer + t << "  " << get_regions_ptr(idx, global_partition) << std::endl << std::flush;
#endif
            solvers.emplace_back(
                std::make_shared<pm::Mwpm>(pm::GraphFlooder(
                    graph.graph_ptr,
                    idx
#ifdef USE_SHMEM
                    ,
                    get_regions_ptr(idx, global_partition),
                    regions_nelems_per_solver
#endif
#ifdef ENABLE_DRAW_FLAGS
                    , &graph.node_part_id
#endif
                    )));
            solvers.back()->flooder.sync_negative_weight_observables_and_detection_events();
        }
    }
#ifdef USE_SHMEM
    // Bare region-tracking arenas for every partition NOT owned by this PE. get_solution_from_remote_pe
    // copies the sender's bitmap into one of these and stamps owner_arena on every rebased region,
    // purely for later del() bookkeeping -- no other Mwpm/GraphFlooder machinery is ever touched for a
    // remote partition. Ordered to match remote_index(): every global partition outside my own local
    // range, in order. Reserved to the exact final size up front so addresses stay stable -- owner_arena
    // pointers into this vector must never dangle.
    size_t num_remote = (size_t)graph.num_partitions - (size_t)my_partition_count;
    remote_arenas.clear();
    remote_arenas.reserve(num_remote * num_shot_containers);
    for (size_t idx = 0; idx < num_shot_containers; ++idx) {
        for (int p = 0; p < (int)graph.num_partitions; ++p) {
            if (p >= my_partitions_start && p < my_partitions_start + my_partition_count) continue;
            remote_arenas.emplace_back(get_regions_ptr(idx, p), regions_nelems_per_solver);
        }
    }
#endif
}

#ifdef USE_SHMEM
SHMEMArena<pm::GraphFillRegion>& pm::DecodingUnit::region_arena_for(int shot_container_id, int global_partition) {
    if (global_partition >= my_partitions_start && global_partition < my_partitions_start + my_partition_count) {
        return solver_for(shot_container_id, global_partition)->flooder.region_arena;
    }
    size_t num_remote = (size_t)graph.num_partitions - (size_t)my_partition_count;
    return remote_arenas[num_remote * (size_t)shot_container_id + (size_t)remote_index(global_partition)];
}
#endif

namespace {
// Shatter+extract one hits list, handling extended-vs-bit-packed accumulation identically
// everywhere it's needed (process_extraction_job's tree walk, extract_crt_received_window).
// divide_vb's per-region variant is separate since it operates on regions directly, not hit lists.
void accumulate_hits(pm::Mwpm& solver, const std::vector<uint64_t>& hits, bool extended, pm::MatchingResult& local_res) {
    if (extended) {
        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(solver, hits);
    } else {
        local_res += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(solver, hits);
    }
}
}  // namespace

#ifdef USE_SHMEM
// returns true if valid, false o.w.
// Ranges are left inclusive, right exclusive
template <typename T>
inline bool ptr_in_range(const T* p, const std::pair<T*, T*>& range) {
    return p >= range.first && p < range.second;
}

bool check_pointers_for_self_and_all_descendents(
    pm::GraphFillRegion* root,
    std::pair<pm::GraphFillRegion*, pm::GraphFillRegion*> region_range,
    std::pair<pm::DetectorNode*, pm::DetectorNode*> p_node_range,
    std::pair<pm::DetectorNode*, pm::DetectorNode*> vb_node_range,
    std::vector<pm::BlossomChild>* discovered_edges) {
    if (!ptr_in_range(root, region_range)) return false;
    for (pm::DetectorNode* node : root->shell_area) {
        if (!ptr_in_range(node, p_node_range) && !ptr_in_range(node, vb_node_range)) return false;
    }
    for (pm::RegionEdge edge : root->blossom_children) {
        if (!check_pointers_for_self_and_all_descendents(
                edge.region,
                region_range,
                p_node_range,
                vb_node_range,
                discovered_edges)) {
            return false;
        }
        discovered_edges->push_back(pm::BlossomChild{(size_t)(root - region_range.first), edge});
    }
    return true;
}

void pm::DecodingUnit::send_solution_to_remote_pe(size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &t, std::ofstream &t_out) {
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_DEFINE(sender_wait);
    SCOREP_USER_REGION_DEFINE(solution_isolation);
    SCOREP_USER_REGION_DEFINE(putmems);
    SCOREP_USER_REGION_DEFINE(shatter);
#endif
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_BEGIN();
#endif
    size_t other_pid = t.other_pid;
    int p_start, p_end, p_k;
    size_t my_vb_offset = 0;
    if (config_parallel::division_strategy == config_parallel::ROUND) { // ROUND
        p_start = (t.iamleft) ? t.part - config_parallel::k + 1 : t.part + 1;
        p_end = p_start + config_parallel::k - 1;
        if (p_start < 0) {
            p_start = 0;
        }
        if (p_end >= graph.num_partitions) {
            p_end = graph.num_partitions - 1;
        }
        p_k = p_end - p_start + 1;
    } else { // OBS
        auto& my_off    = t.iamleft ? t.left_global_offset  : t.right_global_offset;
        size_t my_p_offset  = my_off.first;
        my_vb_offset        = my_off.second;
        p_start = t.vb_left + 1 + (int)my_p_offset;
        p_end   = t.vb_right + (int)my_p_offset;
        p_k = p_end - p_start + 1;
    }

    if (DEBUG)
        t_out << "    sending (p_start" << p_start << ", p_k=" << p_k << ", p_end=" << p_end << ") to " << other_pid << std::endl;

    // We use the LOCAL slot buffer to construct the payload, then PUT it to the REMOTE slot buffer.
    FusionSummary*& fusion_summary_base = t.fusion_summary_shm;

    // Populate FusionSummary header and regions_to_unmatch first; bitmap area is used as temporary checked map.
    if (t.regions_to_unmatch.size() > regions_matched_to_vb_nelems) {
        throw std::invalid_argument("Rank " + std::to_string(pid) + " Thread " + std::to_string(omp_get_thread_num())
                                    + " Shot " + std::to_string(pid) + " Cross-Rank Fusion " + std::to_string(t.part)
                                    + ": t.regions_to_unmatch (" + std::to_string(t.regions_to_unmatch.size())
                                    + ") exceeds regions_matched_to_vb_nelems (" + std::to_string(regions_matched_to_vb_nelems) + ")"
                                    );
    }
    fusion_summary_base->regions_ptr_base = regions_ptr;
    fusion_summary_base->static_nodes_base = graph.graph_ptr->nodes.data();

    // 1. Isolate Solution
    auto* regions_start_ptr = get_regions_ptr(shot_container_id, p_start);
    size_t regions_total_bytes = p_k * regions_nelems_per_solver * sizeof(GraphFillRegion);
    size_t nodes_nelems_total = graph.partition_bounds[p_end].second - graph.partition_bounds[p_start].first + 1;

    if (DEBUG) {
        t_out << "  isolating solution" << std::endl
              << "    node_start_bound: " << graph.partition_bounds[p_start].first << std::endl
              << "    node_end_bound: " << graph.partition_bounds[p_end].second << " (" << nodes_nelems_total << ")" << std::endl
              << std::flush;
    }

    std::pair<pm::GraphFillRegion *, pm::GraphFillRegion *> region_range = {regions_start_ptr, regions_start_ptr + p_k * regions_nelems_per_solver};
    std::pair<pm::DetectorNode*, pm::DetectorNode*> p_node_range = {
        graph.graph_ptr->nodes.data() + graph.partition_bounds[p_start].first,
        graph.graph_ptr->nodes.data() + graph.partition_bounds[p_end].second + 1
    };
    // Nodes in the fusion VB itself may appear in shell_area of VB-matched regions; accept them as valid.
    auto& fvb_bounds = graph.vb_bounds[t.part];
    std::pair<pm::DetectorNode*, pm::DetectorNode*> vb_node_range = {
        graph.graph_ptr->nodes.data() + fvb_bounds.first,
        graph.graph_ptr->nodes.data() + fvb_bounds.second + 1
    };
    
    BlossomChild* child_edges_buff_base = get_child_edges_ptr(shot_container_id, p_start);
    size_t child_edges_counter = 0;

    std::vector<pm::BlossomChild> discovered_child_edges;
    discovered_child_edges.reserve(64);
    
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BEGIN(solution_isolation, "Sender Solution Isolation", SCOREP_USER_REGION_TYPE_COMMON);
#endif
    // Validation Pass: scan all live regions in the k-partition send window. p_start is always one of
    // my own local partitions (a PE only ever sends its own data), so a direct local lookup works --
    // not t.solver, since the window-scan loop below needs consecutive partitions starting exactly at
    // p_start, which isn't guaranteed to equal t.solver's own independently-chosen anchor.
    auto& solver = *solver_for(shot_container_id, p_start);
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_to_unmatch);
    // reset to 1's for safety
    size_t bitmap_words = solver.flooder.region_arena.shmem_bitmap.size() * p_k;
    for (size_t w = 0; w < bitmap_words; ++w) {
        bitmap_base[w] = ~0ULL;
    }
    // validate all blossom roots
    for (size_t i = 0; i < p_k; ++i) {
        int p_scan = p_start + i;
        auto& solver_scan = *solver_for(shot_container_id, p_scan);
        auto& bitmap = solver_scan.flooder.region_arena.shmem_bitmap;
        GraphFillRegion* p_regions_base = get_regions_ptr(shot_container_id, p_scan);

        for (size_t word_id = 0; word_id < bitmap.size(); ++word_id) {
            uint64_t& word = bitmap[word_id];
            if (word == ~0ULL) {
                continue;
            }
            for (size_t bit = 0; bit < 64; ++bit) {
                if ((word >> bit) & 1ULL) {
                    continue;
                }
                size_t local_idx = word_id * 64 + bit;
                GraphFillRegion* r = p_regions_base + local_idx;
                pm::GraphFillRegion*& blossom_root = (r->blossom_parent_top) ? r->blossom_parent_top :  r;

                bool valid = false;

                if (ptr_in_range(blossom_root, region_range)) {
                    size_t root_idx = (size_t)(blossom_root - region_range.first);
                    size_t root_word = root_idx / 64;
                    uint64_t root_mask = 1ULL << (root_idx % 64);
                    // If root bit is 0, it is already processed.
                    if ((bitmap_base[root_word] & root_mask) == 0ULL) {
                        continue;
                    }
                    // Mark root as seen
                    bitmap_base[root_word] &= ~root_mask;

                    // if (DEBUG) t_out << "    checking blossom_root: " << blossom_root << std::endl << std::flush;
                    discovered_child_edges.clear();

                    valid = check_pointers_for_self_and_all_descendents(
                        blossom_root,
                        region_range,
                        p_node_range,
                        vb_node_range,
                        &discovered_child_edges);
                    if (DEBUG && !valid) t_out << "      blossom_root (" << blossom_root << ") not valid, shell_area.size(): " << blossom_root->shell_area.size() << std::endl
                                               << std::flush;

                    // Also check match partner if present.
                    if (valid && blossom_root->match.region) { // check matched region
                        if (ptr_in_range(blossom_root->match.region, region_range)) {
                            size_t match_idx = (size_t)(blossom_root->match.region - region_range.first);
                            bitmap_base[match_idx / 64] &= ~(1ULL << (match_idx % 64));

                            valid = check_pointers_for_self_and_all_descendents(
                                blossom_root->match.region,
                                region_range,
                                p_node_range,
                                vb_node_range,
                                &discovered_child_edges);
                        } else {
                            valid = false;
                        }
                        if (DEBUG && !valid) t_out << "      blossom_root.match.region (" << blossom_root->match.region << ") not valid" << std::endl << std::flush; 
                    } else if (valid && blossom_root->match.edge.loc_to != nullptr) { // region matched to vb
                        valid = false;
                    }
                }
                if (!valid) {
                    if (DEBUG) t_out << "      SHATTERING blossom_root\n" << std::flush;
                    res += solver.shatter_blossom_and_extract_matches(blossom_root);
                } else {
                    for (const auto& child_edge : discovered_child_edges) {
                        if (child_edges_counter < p_k * child_edges_nelems_per_solver) {
                            child_edges_buff_base[child_edges_counter++] = child_edge;
                        } else {
                            throw std::invalid_argument("Rank " + std::to_string(pid) + " Thread " + std::to_string(omp_get_thread_num())
                                                        + " Shot " + std::to_string(pid) + " Cross-Rank Fusion " + std::to_string(t.part)
                                                        + ": # of blossom child edges exceeds k * child_edges_nelems_per_solver");
                        }
                    }
                }
            }
        }
    }

    // The validation pass above can shatter a region still referenced in t's own (by now, via
    // CrossRankTask::setup(), genuinely populated) regions_matched_to_virtual_boundary -- the one
    // shatter in the whole CRT path with no divide_vb call of its own to fold this into (divide_vb
    // handles it for every other shatter site). See prune_stale_regions_matched_to_vb's own comment.
    prune_stale_regions_matched_to_vb(&t);

#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_END(solution_isolation);
    SCOREP_USER_REGION_BEGIN(putmems, "Sender Putmems", SCOREP_USER_REGION_TYPE_COMMON);
#endif

    if (DEBUG) t_out << "  isolated solution" << std::endl << std::flush;
#ifdef ENABLE_DRAW_FLAGS
    if (draw_frames) {
        // t.part is a vb/fusion id, not a partition -- same class of bug as the dbg_solver fix
        // elsewhere in this file (using it directly here would alias whatever real local task happens
        // to already own that solver index). t.solver is my own already-assigned local anchor.
        draw_frame(*t.solver, pm::MwpmEvent::no_event(), 1001, true, omp_get_thread_num());
    }
#endif

    // Send nodes and regions
    auto* node_fields_base = get_node_fields_ptr(shot_container_id, p_start);
    if (DEBUG) {
        t_out << "    sending DetectorNodeEphemeralFields" << std::endl
              << "      Node Fields Base: " << node_fields_base << std::endl
              << "      Nelems: " << nodes_nelems_total << std::endl << std::flush;
    }
    shmem_ctx_putmem_signal_nbi(t.context_shm, 
                                node_fields_base, 
                                node_fields_base, 
                                nodes_nelems_total * sizeof(DetectorNodeEphemeralFields), 
                                t.signal_shm, 
                                1, SHMEM_SIGNAL_ADD, 
                                other_pid);

    if (DEBUG) {
        t_out << "    sending GraphFillRegions" << std::endl
              << "      Regions Base: " << regions_start_ptr << std::endl
              << "      Nelems: " << p_k * regions_nelems_per_solver << std::endl << std::flush;
    }
    shmem_ctx_putmem_signal_nbi(t.context_shm, 
                                regions_start_ptr, 
                                regions_start_ptr, 
                                regions_total_bytes, 
                                t.signal_shm, 
                                1, SHMEM_SIGNAL_ADD, 
                                other_pid);

    // Send BlossomChild array
    //   Sent as one block from the base of partition_start
    if (DEBUG) {
        t_out << "    sending BlossomChild array" << std::endl
              << "      Child Edges Base: " << child_edges_buff_base << std::endl
              << "      Size: " << child_edges_counter * sizeof(BlossomChild) << std::endl
              << "      Nelems: " << child_edges_counter << std::endl << std::flush;
    }
    shmem_ctx_putmem_signal_nbi(t.context_shm, 
                                child_edges_buff_base, 
                                child_edges_buff_base, 
                                child_edges_counter * sizeof(BlossomChild), 
                                t.signal_shm, 
                                1, SHMEM_SIGNAL_ADD, 
                                other_pid);
    fusion_summary_base->blossom_children_size = child_edges_counter;

    // Validate & copy regions_to_unmatch, skipping any that were shattered during solution isolation.
    size_t valid_rtu_count = 0;
    for (GraphFillRegion* region : t.regions_to_unmatch) {
        // Check that the region is still allocated in its solver's bitmap (bit=0 means taken).
        ptrdiff_t rtu_offset = region - regions_start_ptr;
        bool allocated = false;
        if (rtu_offset >= 0 && rtu_offset < (ptrdiff_t)(p_k * regions_nelems_per_solver)) {
            size_t i_solver  = (size_t)rtu_offset / regions_nelems_per_solver;
            size_t local_idx = (size_t)rtu_offset % regions_nelems_per_solver;
            auto& bm = solver_for(shot_container_id, p_start + (int)i_solver)->flooder.region_arena.shmem_bitmap;
            allocated = !((bm[local_idx / 64] >> (local_idx % 64)) & 1ULL);
        }
        if (!allocated) {
            if (DEBUG) t_out << "      SKIPPING shattered region " << region << std::endl << std::flush;
            continue;
        }
        if (valid_rtu_count < regions_matched_to_vb_nelems) {
            fusion_summary_base->regions_to_unmatch[valid_rtu_count] = region;
            ++valid_rtu_count;
            if (DEBUG) t_out << "      " << region << std::endl << std::flush;

        }
    }
    fusion_summary_base->regions_to_unmatch_size = valid_rtu_count;
    bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_to_unmatch
                               + sizeof(GraphFillRegion*) * valid_rtu_count);

    // Copy Bitmap
    size_t total_bitmap_bytes = 0;
    for (size_t i = 0; i < p_k; ++i) {
        int p = p_start + i;
        auto& solver = *solver_for(shot_container_id, p);
        auto& bitmap = solver.flooder.region_arena.shmem_bitmap;
        size_t bitmap_len = bitmap.size();
        // Copy bitmap to FusiionSummary
        size_t bitmap_size_bytes = bitmap_len * sizeof(uint64_t);
        memcpy(bitmap_base, bitmap.data(), bitmap_size_bytes);
        bitmap_base += bitmap_len; // Pointer arithmetic on uint64_t*
        total_bitmap_bytes += bitmap_size_bytes;
    }

    // Send FusionSummary payload (Header + Regions List + Bitmap)
    size_t summary_payload_size = sizeof(FusionSummary) 
                                + fusion_summary_base->regions_to_unmatch_size * sizeof(GraphFillRegion*) 
                                + total_bitmap_bytes;
    if (DEBUG) {
        t_out << "    FusionSummary size: " << summary_payload_size << std::endl
              << "      FusionSummary base: " << fusion_summary_base << std::endl
              << "      Regions Ptr Base: " << fusion_summary_base->regions_ptr_base << std::endl
              << "      Static Nodes Base: " << fusion_summary_base->static_nodes_base << std::endl
              << "      Regions To Unmatch Size: " << fusion_summary_base->regions_to_unmatch_size << std::endl
              << "      Blossom Children Size: " << fusion_summary_base->blossom_children_size << std::endl
              << "      Bitmap Size Bytes: " << total_bitmap_bytes << std::endl << std::flush;
    }
    shmem_ctx_putmem_signal_nbi(t.context_shm, 
                                fusion_summary_base, 
                                fusion_summary_base, 
                                summary_payload_size, 
                                t.signal_shm, 
                                1, SHMEM_SIGNAL_ADD, 
                                other_pid);
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_END(putmems);
#endif
    
    auto& shot = shot_buffer->buffer[shot_container_id];
    std::vector<std::vector<uint64_t>*> hits;
    for (int p_i = p_start; p_i <= p_end; ++p_i) {
        if (DEBUG) t_out << "  p" << p_i << std::flush;
        hits.emplace_back(&shot.partition_hits[p_i]);
    }
    for (int vb_i = t.vb_left + 1 + (int)my_vb_offset; vb_i < t.vb_right + (int)my_vb_offset; ++vb_i) {
        if (DEBUG) t_out << "  vb" << vb_i << std::flush;
        hits.emplace_back(&shot.virtual_boundary_hits[vb_i]);
    }
    if (DEBUG) t_out << std::endl << std::flush;

    // Ensure completion
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BEGIN(sender_wait, "Sender Ctx Quiet", SCOREP_USER_REGION_TYPE_COMMON);
#endif
    shmem_ctx_quiet(t.context_shm);
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_END(sender_wait);
#endif

    // Cleanup sent regions
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BEGIN(shatter, "Sender Shatter Regions", SCOREP_USER_REGION_TYPE_COMMON);
#endif
    if (DEBUG) t_out << "  shattering sent blossoms" << std::endl << std::flush;
    for (std::vector<uint64_t>* hitsref : hits) {
        for (uint64_t i : *hitsref) {
            auto& node_state = solver.flooder.graph.nodes[i].state(shot_container_id);
            // Only shatter if it hasn't been shattered yet (region_that_arrived is still set)
            if (node_state.region_that_arrived) {
                solver.shatter_blossom_and_extract_matches(node_state.region_that_arrived_top);
            }
        }
    }
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_END(shatter);
#endif
    if (DEBUG) t_out << "  shattered sent blossoms" << std::endl
                     << "  sent all data to " << other_pid << std::endl << std::flush;
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_END();
#endif
}

bool pm::DecodingUnit::get_solution_from_remote_pe(
    size_t shot_container_id, CrossRankTask &t, std::ofstream &t_out,
    std::vector<uint64_t> &hitsref) {
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_DEFINE(receiver_wait);
#endif
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_BEGIN();
#endif
    size_t other_pid = t.other_pid;
    int p_k, p_start, p_end;
    if (config_parallel::division_strategy == config_parallel::ROUND) {
        p_k = config_parallel::k;
        p_start = (t.iamleft) ? t.part + 1 : t.part - p_k + 1;
        if (p_start < 0) {
            p_k += p_start;
            p_start = 0;
        }
        if (p_start + p_k > graph.num_partitions) {
            p_k = graph.num_partitions - p_start;
        }
        p_end = p_start + p_k - 1;
    } else { // OBS
        auto& rem_off        = t.iamleft ? t.right_global_offset : t.left_global_offset;
        size_t remote_p_offset = rem_off.first;
        p_start = t.vb_left + 1 + (int)remote_p_offset;
        p_end   = t.vb_right + (int)remote_p_offset;
        p_k     = p_end - p_start + 1;
    }

    // auto& solver = *t.solver;

    // // Isolating local window from the rest of the graph
    // if (DEBUG) t_out << "    isolating solution" << std::endl << std::flush;
    // int my_vb = (config_parallel::division_strategy == config_parallel::ROUND)
    //                 ? (t.iamleft ? t.vb_left : t.vb_right)
    //                 : t.part;
    // if (my_vb >= 0 && my_vb < graph.num_virtual_boundaries) {
    //     auto& vb_bounds = graph.vb_bounds[my_vb];
    //     size_t num_nodes = vb_bounds.second - vb_bounds.first + 1;
    //     pm::DetectorNodeEphemeralFields* vb_fields_base = node_ephemeral_fields_ptr + shot_container_id * nodes_nelems_per_buffer + vb_bounds.first;
    //     for (int i = 0; i < num_nodes; ++i) {
    //         if ((vb_fields_base+i)->region_that_arrived_top) {
    //             if (DEBUG) t_out << "        SHATTERING: " << (vb_fields_base+i)->region_that_arrived_top << std::endl;
    //             res += solver.shatter_blossom_and_extract_matches((vb_fields_base+i)->region_that_arrived_top);
    //         }
    //     }
    // }

    if (DEBUG) t_out << "    getting (p_start=" << p_start << ", p_k=" << p_k << ", p_end=" << p_end << ") data from " << other_pid << std::endl << std::flush;

#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BEGIN(receiver_wait, "Receiver Wait Until", SCOREP_USER_REGION_TYPE_COMMON);
#endif
    // Wait for signal (4 puts expected)
    shmem_wait_until(t.signal_shm, SHMEM_CMP_EQ, 4);
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_END(receiver_wait);
#endif

    // Read local FusionSummary buffer (which was populated by remote PE)
    FusionSummary*& fusion_summary_base = t.fusion_summary_shm;
    
    // Pointers for rebasing
    GraphFillRegion*& remote_regions_ptr_base = fusion_summary_base->regions_ptr_base;
    DetectorNode*& remote_static_nodes_base = fusion_summary_base->static_nodes_base;
    DetectorNode* local_static_nodes_base = graph.graph_ptr->nodes.data();
    if (DEBUG) t_out << "  local_status_nodes_base: " << local_static_nodes_base << std::endl;

    // Define memory ranges for OOB check
    GraphFillRegion* k_block_regions_base = get_regions_ptr(shot_container_id, p_start);
    GraphFillRegion* k_block_regions_end = k_block_regions_base + p_k * regions_nelems_per_solver;
    // Define GLOBAL memory range for stricter-but-permissive OOB check (allow pointers to ANY partition owned by this PE/process)
    if (DEBUG) t_out << "  k_block_regions_base: " << k_block_regions_base << std::endl
                     << "  k_block_regions_end: " << k_block_regions_end << " (" << p_k*regions_nelems_per_solver << ")" << std::endl << std::flush;

    auto& nodes_start_bounds = graph.partition_bounds[p_start];
    auto& nodes_end_bounds = graph.partition_bounds[p_end];
    size_t nodes_nelems_total = nodes_end_bounds.second - nodes_start_bounds.first + 1;
    ptrdiff_t nodes_start_bound_with_vb, nodes_end_bound_with_vb;
    if (config_parallel::division_strategy == config_parallel::ROUND) {
        nodes_start_bound_with_vb = (t.iamleft) ? (ptrdiff_t)graph.vb_bounds[p_start-1].first : (ptrdiff_t)nodes_start_bounds.first;
        nodes_end_bound_with_vb = (!t.iamleft) ? (ptrdiff_t)graph.vb_bounds[p_end].second : (ptrdiff_t)nodes_end_bounds.second;
    } else { // OBS: partition range only; fusion VB accepted separately below
        nodes_start_bound_with_vb = (ptrdiff_t)nodes_start_bounds.first;
        nodes_end_bound_with_vb   = (ptrdiff_t)nodes_end_bounds.second;
    }
    auto& fvb_bounds = graph.vb_bounds[t.part];
    ptrdiff_t fvb_start = (ptrdiff_t)fvb_bounds.first;
    ptrdiff_t fvb_end   = (ptrdiff_t)fvb_bounds.second;
    if (DEBUG) t_out << "  node_start_bound: " << nodes_start_bounds.first << std::endl
                     << "  node_end_bound: " << nodes_end_bounds.second << " (" << nodes_nelems_total << ")" << std::endl
                     << "  nodes_start_bound_with_vb: " << nodes_start_bound_with_vb << std::endl
                     << "  nodes_end_bound_with_vb: " << nodes_end_bound_with_vb << std::endl << std::flush;

    // Helpers for pointer rebasing with OOB check
    auto rebase_region_ptr = [&](GraphFillRegion* ptr) -> GraphFillRegion* {
        // Rebase to local memory
        GraphFillRegion* rebased = (GraphFillRegion*)((char*)ptr - (char*)remote_regions_ptr_base + (char*)regions_ptr);
        // Check if remote pointer is within the valid range of the transmitted block OR global
        if (rebased < k_block_regions_base || rebased >= k_block_regions_end) {
            if (DEBUG) t_out << "      PRUNED OOB region ptr: " << ptr << " offset to base: " << (char*)ptr - (char*)remote_regions_ptr_base << std::endl << std::flush;
            return nullptr;
        }
        return rebased;
    };
    auto rebase_node_ptr = [&](DetectorNode* ptr) -> DetectorNode* {
        ptrdiff_t index = ptr - remote_static_nodes_base;
        bool in_partition = index >= nodes_start_bound_with_vb && index <= nodes_end_bound_with_vb;
        bool in_fvb       = index >= fvb_start && index <= fvb_end;
        if (!in_partition && !in_fvb) {
            if (DEBUG) t_out << "      PRUNED OOB node ptr: " << ptr << " index: " << index << std::endl << std::flush;
            return nullptr;
        }
        return local_static_nodes_base + index;
    };

    if (DEBUG) {
        t_out << "    FusionSummary recieved" << std::endl
              << "      FusionSummary Base: " << fusion_summary_base << std::endl
              << "      Regions Ptr Base (remote): " << remote_regions_ptr_base << std::endl
              << "      Static Nodes Base (remote): " << remote_static_nodes_base << std::endl
              << "      Regions To Unmatch Size: " << fusion_summary_base->regions_to_unmatch_size << std::endl
              << "      Blossom Children Size: " << fusion_summary_base->blossom_children_size << std::endl << std::flush;
    }

    // 1. Reconstruct regions_to_unmatch using post-prune liveness.
    for (size_t i = 0; i < fusion_summary_base->regions_to_unmatch_size; ++i) {
        GraphFillRegion* rebased = rebase_region_ptr(fusion_summary_base->regions_to_unmatch[i]);
        if (rebased) {
            if (DEBUG) t_out << "      " << rebased << std::endl << std::flush;
            t.regions_to_unmatch.push_back(rebased);
        }
    }

    // 2. Copy Bitmaps & Reconstruct GraphFillRegions (Iterate k partitions)
    if (DEBUG) {
        t_out << "    reconstructing solution state\n"
              << "      Remote Regions To Unmatch Size: " << fusion_summary_base->regions_to_unmatch_size << std::endl << std::flush;
    }
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_to_unmatch
                                        + sizeof(GraphFillRegion*) * fusion_summary_base->regions_to_unmatch_size);
    // p_start..p_start+p_k-1 is the REMOTE side's own partition range -- genuinely needs remote arena
    // access (region_arena_for), the one place in this whole refactor that does: a full Mwpm is never
    // built for a remote partition, only its bare SHMEMArena (see build_solvers()/region_arena_for()),
    // since this loop only needs the bitmap + later del() bookkeeping via owner_arena below, not any
    // other Mwpm/GraphFlooder machinery.
    size_t bitmap_len = region_arena_for(shot_container_id, p_start).shmem_bitmap.size();
    for (size_t i = 0; i < p_k; ++i) {
        int p = p_start + i;
        auto& region_arena = region_arena_for(shot_container_id, p);
        auto& bitmap_dst = region_arena.shmem_bitmap;
        // Copy bitmap
        memcpy(bitmap_dst.data(), bitmap_base + i*bitmap_len, bitmap_dst.size() * sizeof(uint64_t));

        GraphFillRegion* p_regions_base = get_regions_ptr(shot_container_id, p);

        for (size_t word_id = 0; word_id < bitmap_len; ++word_id) {
            uint64_t word = bitmap_dst[word_id];
            if (word != ~0ULL) {
                for (size_t bit = 0; bit < 64; ++bit) {
                    if (!((word >> bit) & 1ULL)) { // taken
                        size_t index = word_id * 64 + bit;
                        GraphFillRegion* r = p_regions_base + index;
                        // Rebase pointers with OOB check
                        if (r->blossom_parent) {
                            r->blossom_parent = rebase_region_ptr(r->blossom_parent);
                            if (DEBUG && !r->blossom_parent) t_out << "        blossom_parent\n" << std::flush;
                        }
                        if (r->blossom_parent_top) {
                            r->blossom_parent_top = rebase_region_ptr(r->blossom_parent_top);
                            if (DEBUG && !r->blossom_parent_top) t_out << "        blossom_parent_top\n" << std::flush;
                        }
                        if (r->match.region) {
                            r->match.region = rebase_region_ptr(r->match.region);
                            if (DEBUG && !r->match.region) t_out << "        match.region\n" << std::flush;
                        }
                        if (r->match.edge.loc_from) {
                            r->match.edge.loc_from = rebase_node_ptr(r->match.edge.loc_from);
                            if (DEBUG && !r->match.edge.loc_from) t_out << "        match.edge.loc_from\n" << std::flush;
                        }
                        if (r->match.edge.loc_to) {
                            r->match.edge.loc_to = rebase_node_ptr(r->match.edge.loc_to);
                            if (DEBUG && !r->match.edge.loc_to) t_out << "        match.edge.loc_to\n" << std::flush;
                        }
                        // Reset vectors (they contain pointers that need manual reconstruction)
                        new (&r->blossom_children) std::vector<RegionEdge>();
                        new (&r->shell_area) std::vector<DetectorNode*>();
                        r->shrink_event_tracker.clear();
                        r->alt_tree_node = nullptr;
                        r->owner_arena = &region_arena;
                    }
                }
            }
        }
    }
    if (DEBUG) {
        t_out << "    rebased regions" << std::endl << std::flush;
    }

    // 4. Reconstruct Blossom Children
    //   Blossom children for all K partitions are packed in the buffer starting at `partition_start`.
    BlossomChild* child_edges_buff_base = get_child_edges_ptr(shot_container_id, p_start);

    for (size_t i=0; i < fusion_summary_base->blossom_children_size; ++i) {
        BlossomChild& child = child_edges_buff_base[i];
        GraphFillRegion* parent = k_block_regions_base + child.blossom_parent;

        RegionEdge local_edge = child.region_edge;
        local_edge.region = rebase_region_ptr(local_edge.region);
        
        if (local_edge.region) {
            if (local_edge.edge.loc_from)
                local_edge.edge.loc_from = rebase_node_ptr(local_edge.edge.loc_from);
            if (local_edge.edge.loc_to)
                local_edge.edge.loc_to = rebase_node_ptr(local_edge.edge.loc_to);
            parent->blossom_children.push_back(local_edge);
        } else {
             if (DEBUG) t_out << "      PRUNED OOB Blossom Child" << std::endl << std::flush;
        }
    }
    if (DEBUG) {
        t_out << "    reconstructed blossom children: " << fusion_summary_base->blossom_children_size << std::endl << std::flush;
    }

    // 5. Rebase DetectorNodeEphemeralFields (Contiguous Block)
    DetectorNodeEphemeralFields* node_fields_base = get_node_fields_ptr(shot_container_id, p_start);
    for (size_t i = 0; i < nodes_nelems_total; ++i) {
        DetectorNodeEphemeralFields& fields = node_fields_base[i];
        // bool clear_node = false;
        bool cleared_region_that_arrived = false;
        if (fields.region_that_arrived) {
            fields.region_that_arrived = rebase_region_ptr(fields.region_that_arrived);
            if (fields.region_that_arrived) {
                fields.region_that_arrived->shell_area.push_back(&graph.graph_ptr->nodes[nodes_start_bounds.first + i]);
            } else {
                cleared_region_that_arrived = true;
            }
        }
        if (fields.region_that_arrived_top) {
            fields.region_that_arrived_top = rebase_region_ptr(fields.region_that_arrived_top);
            if (fields.region_that_arrived_top) {
                if (cleared_region_that_arrived) {
                   throw std::invalid_argument("Rank " + std::to_string(pid) + " Thread " + std::to_string(omp_get_thread_num())
                                                + " Shot " + std::to_string(pid) + " Cross-Rank Fusion " + std::to_string(t.part)
                                                + ": reciever cleared node's region_that_arrived but not region_that_arrived_top");
                }
            } else if (!cleared_region_that_arrived) {
                throw std::invalid_argument("Rank " + std::to_string(pid) + " Thread " + std::to_string(omp_get_thread_num())
                            + " Shot " + std::to_string(pid) + " Cross-Rank Fusion " + std::to_string(t.part)
                            + ": reciever cleared node's region_that_arrived_top but not region_that_arrived");
            }
        }
        if (fields.reached_from_source) {
            fields.reached_from_source = rebase_node_ptr(fields.reached_from_source);
        }
    }

    if (DEBUG) {
        t_out << "    rebased DetectorNodeEphemeralFields" << std::endl
              << "      Node Fields Base: " << node_fields_base << std::endl
              << "      Nelems: " << nodes_nelems_total << std::endl << std::flush;
    }
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_END();
#endif
    return true;
}

void pm::DecodingUnit::extract_crt_received_window(ShotContainer& shot, size_t shot_container_id, CrossRankTask& crt, int tid, std::ofstream* t_out) {
    // Computes the received partition/vb range directly from crt's own fields -- independent of
    // whatever range get_solution_from_remote_pe used internally for rebasing, since that's a
    // different concern (this needs the actual shot.partition_hits[]/virtual_boundary_hits[]
    // extraction range, division-strategy-aware).
    int p_lo, p_hi, vb_lo, vb_hi;
    if (config_parallel::division_strategy == config_parallel::ROUND) {
        int p_k = config_parallel::k;
        p_lo = (crt.iamleft) ? crt.part + 1 : crt.part - p_k + 1;
        if (p_lo < 0) {
            p_k += p_lo;
            p_lo = 0;
        }
        if (p_lo + p_k > graph.num_partitions) {
            p_k = graph.num_partitions - p_lo;
        }
        p_hi = p_lo + p_k - 1;
        vb_lo = crt.vb_left + 1;
        vb_hi = crt.vb_right - 1;
    } else {
        // OBS: known pre-existing gap, not fixed as part of this plan (OBS support is explicitly
        // out of scope -- see this-is-a-broader-purrfect-crystal.md). The old (pre-this-session)
        // code had two extraction-time recomputations of this same range that disagreed with each
        // other on left/right offset selection and off-by-ones (extended vs bit-packed branches) --
        // this mirrors the old bit-packed branch's formula, not independently re-derived/verified.
        auto& rem_off = crt.iamleft ? crt.right_global_offset : crt.left_global_offset;
        p_lo = crt.vb_left + 1 + (int)rem_off.first;
        p_hi = crt.vb_right + (int)rem_off.first;
        vb_lo = crt.vb_left + 1 + (int)rem_off.second;
        vb_hi = crt.vb_right + (int)rem_off.second - 1;
    }

    if (DEBUG && t_out) {
        *t_out << "  EXTRACT_CRT part=" << crt.part
               << " p=[" << p_lo << ", " << p_hi << "]"
               << " vb=[" << vb_lo << ", " << vb_hi << "]"
               << " solver=" << crt.solver << std::endl << std::flush;
    }
    auto& solver = *crt.solver;
    bool extended = shot.num_observables > sizeof(pm::obs_int) * 8;
    pm::MatchingResult local_res{};
    auto do_walk = [&]() {
        for (int p = p_lo; p <= p_hi; ++p) accumulate_hits(solver, shot.partition_hits[p], extended, local_res);
        for (int vb = vb_lo; vb <= vb_hi; ++vb) accumulate_hits(solver, shot.virtual_boundary_hits[vb], extended, local_res);
    };
    if (extended) {
        // No synchronization needed: this thread owns solver exclusively here (only the resolving
        // thread ever touches its own CRT/connector's solver+shot.res at this point -- see plan
        // this-is-a-broader-purrfect-crystal.md Design §10).
        do_walk();
        solver.extract_paths_from_match_edges(solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
        solver.flooder.match_edges.clear();
    } else {
        do_walk();
        shot.thread_results[tid] += local_res;
    }
}
#endif

void pm::DecodingUnit::prune_stale_regions_matched_to_vb(TaskBase* t) {
    auto& list = t->regions_matched_to_virtual_boundary;
    list.erase(
        std::remove_if(list.begin(), list.end(),
            [](pm::GraphFillRegion* r) { return !r->allocated; }),
        list.end());
}

void pm::DecodingUnit::divide_vb(ShotContainer& shot, int shot_container_id, int vb_id, pm::Mwpm* solver_arg, int tid, TaskBase* prune_target, std::ofstream* t_out) {
    if (DEBUG && t_out) {
        *t_out << "  DIVIDE_VB vb=" << vb_id << " solver=" << solver_arg
               << " prune_target=" << static_cast<void*>(prune_target) << std::endl << std::flush;
    }
    // Loop every node in the vb's range (not just ones with "hits" -- a blossom can span the vb
    // without either side registering a hit exactly there) and shatter any region that reached it.
    // Any solver works as scratch in principle (region ownership is globally node-indexed via
    // node.state(...), not solver-private) -- but the caller's own ->solver must still be a genuine
    // descendant leaf's solver, NOT an arbitrary fixed choice: a fixed anchor can collide with
    // whichever solver this thread (or another) is legitimately still using for its own in-flight
    // work at this exact moment -- a real, previously-hit crash.
    auto& vb_bounds = graph.vb_bounds[vb_id];
    auto& solver = *solver_arg;
    bool extended = shot.num_observables > sizeof(pm::obs_int) * 8;
    if (extended) {
        // No synchronization needed: this thread owns solver exclusively here (see Design §10).
        for (size_t i = vb_bounds.first; i <= vb_bounds.second; ++i) {
            auto* region = graph.graph_ptr->nodes[i].state(shot_container_id).region_that_arrived_top;
            if (region) solver.shatter_blossom_and_extract_match_edges(region, solver.flooder.match_edges);
        }
        solver.extract_paths_from_match_edges(solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
        solver.flooder.match_edges.clear();
    } else {
        pm::MatchingResult local_res{};
        for (size_t i = vb_bounds.first; i <= vb_bounds.second; ++i) {
            auto* region = graph.graph_ptr->nodes[i].state(shot_container_id).region_that_arrived_top;
            if (region) local_res += solver.shatter_blossom_and_extract_matches(region);
        }
        shot.thread_results[tid] += local_res;
    }
    // A blossom shattered here can span farther than this vb and destroy a region still referenced
    // in prune_target->regions_matched_to_virtual_boundary (kept there for a LATER setup() call --
    // the next checkpoint up the chain, or a chained special task -- to consume). GraphFillRegion::
    // allocated (set/cleared solely by the owning arena's own alloc/del) makes this a plain liveness
    // check, not an explicitly-collected destroyed-list intersection. prune_target == nullptr (post-
    // hoc chunking, where the whole tree is already fully solved and no future setup() will ever read
    // that list again) skips this entirely.
    if (prune_target) {
        prune_stale_regions_matched_to_vb(prune_target);
    }
}

// Walks start's own left_child chain backward, collecting every deferred connector (is_extraction_
// unit_connector && defer_division) until hitting the first non-deferred boundary -- the connector
// whose own vb was already divided normally, marking the true edge of "already settled" territory (or
// a raw leaf/unit root if the deferred span reaches all the way back to the observable's own start).
// Collected in right-to-left (nearest-to-start-first) order, then divided+posted in reverse (left-to-
// right) order, since each one's own vb divide only needs its own immediate left/right sides to be
// independently addressable, not any relative ordering between different deferred links.
void pm::DecodingUnit::back_divide_walk(
    Task* start, ShotContainer& shot, int shot_container_id, int tid, std::ofstream* t_out) {
    std::vector<Task*> deferred;
    Task* cur = start->left_child;
    while (cur != nullptr && cur->is_extraction_unit_connector && cur->defer_division) {
        deferred.push_back(cur);
        cur = cur->left_child;
    }
    for (auto it = deferred.rbegin(); it != deferred.rend(); ++it) {
        Task* link = *it;
        divide_vb(shot, shot_container_id, link->part, link->solver, tid, link, t_out);
        Task* left = link->left_child;
        Task* closed_unit = left->is_extraction_unit_connector ? left->right_child : left;
        shot.post_extraction_job(ExtractionJob{shot_container_id, closed_unit}, t_out);
    }
}

void pm::DecodingUnit::process_extraction_job(ShotContainer& shot, const ExtractionJob& job, int tid, std::ofstream* t_out) {
    Task* root = job.subtree_root;
    if (DEBUG && t_out) {
        *t_out << "  PROCESS_JOB " << (root->is_fusion ? "vb=" : "p=") << root->part
               << " task=" << static_cast<void*>(root) << std::endl << std::flush;
    }
    // Any solver from this shot_container_id block would work as scratch here (region ownership is
    // globally node-indexed, not solver-private) -- use the subtree root's own already-assigned solver.
    auto& solver = *root->solver;
    bool extended = shot.num_observables > sizeof(pm::obs_int) * 8;
    pm::MatchingResult local_res{};

    auto do_walk = [&]() {
        std::vector<Task*> to_visit = {root};
        while (!to_visit.empty()) {
            Task* curr = to_visit.back(); to_visit.pop_back();
            if (!curr->is_fusion) {
                accumulate_hits(solver, shot.partition_hits[curr->part], extended, local_res);
            } else {
                // Internal (non-checkpoint) fusion vb -- still the simpler hits-list approach for
                // now (tracked separately as a correctness gap, same as divide_vb's approach fixes
                // for checkpoint/split-point boundaries specifically).
                accumulate_hits(solver, shot.virtual_boundary_hits[curr->part], extended, local_res);
                // left_child/right_child are TaskBase* (a preemptive-OBS chain-link's child can be a
                // CrossRankTask -- see decoding_task.h), but this walk only ever runs on non-preemptive/
                // ROUND subtrees, which never attach a CRT as a child (CRTs there sit above a subtree
                // root via CrossRankTask's own parent-chain walk instead) -- always genuinely Task* here.
                if (curr->left_child) to_visit.push_back(static_cast<Task*>(curr->left_child));
                if (curr->right_child && curr->right_child != curr->left_child)
                    to_visit.push_back(static_cast<Task*>(curr->right_child));
            }
        }
    };

    if (extended) {
        do_walk();
        solver.extract_paths_from_match_edges(solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
        solver.flooder.match_edges.clear();
    } else {
        do_walk();
        shot.thread_results[tid] += local_res;  // each thread only ever writes its own slot -- no race
    }
    shot.pending_extraction_jobs.fetch_sub(1, std::memory_order_acq_rel);
}

void pm::DecodingUnit::post_hoc_chunk_and_post(Task* node, ShotContainer& shot, int shot_container_id, int tid, std::ofstream* t_out) {
    if (node->is_extraction_unit_root) {
        shot.post_extraction_job(ExtractionJob{shot_container_id, node}, t_out);
        return;
    }
    // node->is_extraction_unit_connector: divide its own vb inline first -- this is what makes both
    // children independently safe to post/recurse into, regardless of order or which thread
    // eventually drains them -- then recurse into each. Single top-down pass, every node visited
    // exactly once, so there's no risk of double-posting a subtree.
    divide_vb(shot, shot_container_id, node->part, node->solver, tid, /*prune_target=*/nullptr, t_out);
    // Non-preemptive OBS/ROUND trees never attach a CRT as a child (see do_walk's own comment above) --
    // always genuinely Task* here.
    post_hoc_chunk_and_post(static_cast<Task*>(node->left_child), shot, shot_container_id, tid, t_out);
    if (node->right_child != node->left_child) {
        post_hoc_chunk_and_post(static_cast<Task*>(node->right_child), shot, shot_container_id, tid, t_out);
    }
}

// Core parallel decoding loop
void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_BEGIN();
#endif
#ifdef USE_SHMEM
    if (DEBUG) {
        std::string ps = "PE" + std::to_string(pid) + " partition task ids: ";
        for (auto p : my_partition_task_ids)
            ps += std::to_string(p) + " ";
        std::cout << ps << std::endl << std::flush;
    }
#endif
#pragma omp parallel
    {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_DEFINE(local_decoding);
        SCOREP_USER_REGION_DEFINE(cross_rank_fusion);
        SCOREP_USER_REGION_DEFINE(solution_extraction);
        SCOREP_USER_REGION_DEFINE(shot_decode);
        SCOREP_USER_REGION_DEFINE(shot_spin_wait);
        SCOREP_USER_REGION_DEFINE(shot_iteration);
#endif
        const int tid = omp_get_thread_num();
#ifdef ENABLE_SHOT_BUFFERS
        const int num_shot_containers = NUM_BUFFERS_PER_UNIT;
        // Idle-helper participation gate (see this-is-a-broader-purrfect-crystal.md Design §11 and
        // i-see-here-s-the-delegated-pelican.md): bounds how many threads ever spin on the extraction
        // queue at once to roughly the number of extraction units that could ever exist, so idle
        // threads for a small job don't pay any contention cost at all. Strided (not a tid < N
        // prefix) so helper duty doesn't concentrate on the low-tid threads that also own the next
        // shot's earliest static-leaf partitions. Pure function of already-fixed values (num_threads,
        // graph.num_partitions, config_parallel::L) -- identical, deterministic result on every thread.
        const int num_extraction_units =
            ((int)graph.num_partitions + config_parallel::L - 1) / config_parallel::L;  // ceil
        const int helper_stride = std::max(1, (int)num_threads / std::max(1, num_extraction_units));
        if (DEBUG && tid == 0) std::cout << "helper_stride: " << helper_stride << std::endl;
#else
        const int num_shot_containers = 1;
#endif
        std::ofstream t_out;
        if (BARE_DEBUG || DEBUG) {
            std::string t_out_dir = "out_parallel/";
#ifdef USE_SHMEM
            t_out_dir += "p" + std::to_string(pid);
#endif
            std::filesystem::create_directories(t_out_dir);
            std::string t_out_name = t_out_dir + "/t" + std::to_string(tid) + ".out";
            t_out.open(t_out_name);
            std::cout << "T" << tid << " of " << num_threads << std::endl;
        }
        // Start decoding
        size_t shot_container_id = 0;
        int shot_buffer_round =
            shot_buffer->buffer[0].current_buffer_round.load();  // how many times buffer has looped
        int shot_id = shot_buffer_round * num_shot_containers;
        // Thread-local list of local roots solved during the steal loop.
        // Cleared at the start of each shot; a thread may solve multiple roots when
        // it exhausts all its assigned partition leaves (via next_p_id). Universal across build
        // configs -- round-partitioning always has exactly one true root (num_task_roots == 1), so
        // in practice at most one thread ends up with a non-empty list per shot.
        struct RootInfo { Task* task; };
        std::vector<RootInfo> roots_i_solved;
        bool i_solved_last_root = false;
        try {
            while (true) {
                if (BARE_DEBUG) {
                    t_out << std::endl
                          << "Starting shot " << shot_id << ", buffer round " << shot_buffer_round << ", buffer_id "
                          << shot_container_id << std::endl << std::flush;
                }
#ifdef ENABLE_DRAW_FLAGS
                if (draw_frames) {
                    std::string frames_out = "out_parallel/frames/" + std::to_string(shot_id);
#ifndef USE_SHMEM
                    std::filesystem::create_directories(frames_out);
#else
                    frames_out += "/p" + std::to_string(pid);
                    std::filesystem::create_directories(frames_out);
#endif
                    std::filesystem::create_directories(frames_out + "/t" + std::to_string(tid));
                }
#endif
                auto& shot = shot_buffer->buffer[shot_container_id];
                int shot_current_buffer_round = shot.current_buffer_round.load();
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_BEGIN(shot_iteration, "Shot Iteration", SCOREP_USER_REGION_TYPE_COMMON);
                SCOREP_USER_REGION_BEGIN(shot_spin_wait, "Shot Spin Wait", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                while (shot_current_buffer_round < shot_buffer_round) {  // wait
                    if (shot_current_buffer_round < 0) {
                        break;
                    }
                    shot.current_buffer_round.wait(shot_current_buffer_round, std::memory_order_relaxed);
                    shot_current_buffer_round = shot.current_buffer_round.load();
                }
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_END(shot_spin_wait);
#endif
                if (shot_current_buffer_round < 0) {
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_END(shot_iteration);
#endif
                    break;
                }
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_BEGIN(shot_decode, "Shot Decode", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                // we divide partition tasks into sets based on the number of threads available
                Task* t = &shot.tasks[my_partition_task_ids[tid]];
                size_t next_p_inc = num_threads; // cannot inc by two for 1 threads --- does not work for odd
                size_t next_p_id = tid+next_p_inc;
                if (DEBUG) t_out << "solver: " << t->solver << "\n";
                bool stolen = t->try_to_steal(shot_buffer_round);
                roots_i_solved.clear();
                // Declared per-shot (not per-root): a CRT can now be resolved mid-climb, by a thread
                // that never itself reaches a root this shot (a sibling subtree's own solver might win
                // the race up to the shared parent first) -- see the report_done restructuring below.
#ifdef USE_SHMEM
                std::vector<CrossRankTask*> crts_i_handled;
                pm::MatchingResult& my_result = shot.thread_results[tid];
#endif
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_BEGIN(local_decoding, "Local Decoding", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                while (stolen) {  // Got task
                    if (DEBUG) t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl << std::flush;
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    // Solve task
                    t->setup();
                    if (DEBUG) t_out << "  solver: " << t->solver << "\n";
                    pm::Mwpm& solver = *t->solver;
                    solver.prepare_for_task(t, shot_id);
#ifdef USE_SHMEM
                    if (DEBUG) t_out << "  solver bounds: " << solver.flooder.vb_left << "(vb_left) " << solver.flooder.vb_right << " (vb_right)" << std::endl << std::flush;
#ifdef ENABLE_DRAW_FLAGS
                    if (draw_frames && config_parallel::division_strategy == config_parallel::OBS) {
                        // t is never itself a seam/CRT now (those attach via special_tasks instead) --
                        // always the plain "leaf partition" case.
                        size_t K_p = graph.num_partitions / graph.num_obs_patches;
                        size_t K_vb = K_p - 1;
                        size_t patch = (t->is_fusion) ? t->part / K_vb : t->part / K_p;
                        solver.flooder.p_offsets  = { K_p  * patch };
                        solver.flooder.vb_offsets = { K_vb * patch };
                    }
#endif // ENABLE_DRAW_FLAGS
#endif // USE_SHMEM
                    pm::process_timeline_until_completion(
                        solver,
                        hitsref,
#ifdef ENABLE_DRAW_FLAGS
                        draw_frames,
#endif
                        true,
                        tid);
                    if (DEBUG) t_out << "  solved " << (t->is_fusion ? "f" : "p") << t->part << std::endl << std::flush;
                    t->mark_solved();
                    // Process any attached special tasks (seams/CRTs) in order, then t's own deferred
                    // span (if any), then t's own divide+post (if it isn't itself deferred) -- see
                    // now-its-time-to-jazzy-turing.md Phase three §§2-4.
                    for (SpecialTask* st : t->special_tasks) {
                        if (st->is_local_seam_fusion()) {
                            auto* seam = static_cast<LocalSeamTask*>(st);
                            if (seam->try_to_steal(0)) {
                                seam->setup();
                                pm::Mwpm& seam_solver = *seam->solver;
                                seam_solver.prepare_for_task(seam, shot_id);
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames && config_parallel::division_strategy == config_parallel::OBS) {
                                    size_t K_p = graph.num_partitions / graph.num_obs_patches;
                                    size_t K_vb = K_p - 1;
                                    seam_solver.flooder.p_offsets.clear();
                                    seam_solver.flooder.vb_offsets.clear();
                                    for (int op : seam->obs_patch_ids) {
                                        seam_solver.flooder.p_offsets.push_back((size_t)op * K_p);
                                        seam_solver.flooder.vb_offsets.push_back((size_t)op * K_vb);
                                    }
                                }
#endif
                                pm::process_timeline_until_completion(
                                    seam_solver, shot.virtual_boundary_hits[seam->part],
#ifdef ENABLE_DRAW_FLAGS
                                    draw_frames,
#endif
                                    true, tid);
                                seam->mark_solved();
                                divide_vb(shot, (int)shot_container_id, seam->part, seam->solver, tid, seam, &t_out);
                                // Round-tagged, not a plain bool -- see decoding_task.h LocalSeamTask::
                                // ready's own comment for why a stale value would otherwise be possible.
                                seam->ready.store(shot_buffer_round, std::memory_order_release);
                            } else {
                                while (seam->ready.load(std::memory_order_acquire) != shot_buffer_round) {}
                            }
                        }
#ifdef USE_SHMEM
                        else if (st->is_cross_rank_fusion()) {
                            auto* crt = static_cast<CrossRankTask*>(st);
                            if (DEBUG) t_out << "Trying cross-rank fusion vb=" << crt->part
                                             << " iamleft=" << crt->iamleft
                                             << " other_pid=" << crt->other_pid << "\n" << std::flush;
                            if (BARE_DEBUG) t_out << "    waiting until PE done" << std::endl << std::flush;
                            crt->wait_until_done(pid, shot_buffer_round-1);
                            if (crt->try_to_steal(pid)) {
                                if (BARE_DEBUG) t_out << "  Stole CRT with " << crt->other_pid << std::endl << std::flush;
                                crt->setup();
                                auto& crt_solver = *crt->solver;
                                crt_solver.prepare_for_task(crt, shot_id);
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames && config_parallel::division_strategy == config_parallel::OBS) {
                                    crt_solver.flooder.p_offsets  = { crt->left_global_offset.first,  crt->right_global_offset.first  };
                                    crt_solver.flooder.vb_offsets = { crt->left_global_offset.second, crt->right_global_offset.second };
                                }
#endif
                                auto& crt_hitsref = shot.virtual_boundary_hits[crt->part];
                                get_solution_from_remote_pe(shot_container_id, *crt, t_out, crt_hitsref);
                                if (BARE_DEBUG) t_out << "  Solving CRT " << crt->part
                                    << " bounds " << crt_solver.flooder.vb_left << " " << crt_solver.flooder.vb_right << std::endl << std::flush;
                                pm::process_timeline_until_completion(
                                    crt_solver,
                                    crt_hitsref,
#ifdef ENABLE_DRAW_FLAGS
                                    draw_frames,
#endif
                                    true,
                                    tid);
                                crt->mark_solved(pid);
                                // Extract inline, right here on the resolving thread -- not queued
                                // (queuing would only add overhead, since this thread is already
                                // doing the work). Divide crt->part first (now that it's fully
                                // solved) so the received window can be safely extracted alongside
                                // it, in the same breath.
                                divide_vb(shot, (int)shot_container_id, crt->part, crt->solver, tid, /*prune_target=*/crt, &t_out);
                                extract_crt_received_window(shot, shot_container_id, *crt, tid, &t_out);
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames) draw_frame(crt_solver, pm::MwpmEvent::no_event(), 1001, true, tid);
#endif
                            } else {
                                if (BARE_DEBUG) t_out << "  Sending CRT data to " << crt->other_pid << std::endl;
                                crt->setup();
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames) {
                                    // crt->solver is assigned once at construction, anchored on my own
                                    // local side of the boundary (see build_tasks_for_obs_patch_
                                    // partitioning) -- using it directly here is what makes the old bug
                                    // (using crt->part, a vb/fusion id, as if it were a partition index,
                                    // aliasing whatever real local task happened to already own that
                                    // solver) structurally impossible to reintroduce.
                                    auto& dbg_solver = *crt->solver;
                                    dbg_solver.prepare_for_task(crt, shot_id);
                                    if (config_parallel::division_strategy == config_parallel::OBS) {
                                        dbg_solver.flooder.p_offsets  = { crt->left_global_offset.first,  crt->right_global_offset.first  };
                                        dbg_solver.flooder.vb_offsets = { crt->left_global_offset.second, crt->right_global_offset.second };
                                    }
                                    draw_frame(dbg_solver, pm::MwpmEvent::no_event(), 1000, true, tid);
                                }
#endif
                                send_solution_to_remote_pe(shot_container_id, my_result, *crt, t_out);
                            }
                            crts_i_handled.push_back(crt);
                        }
#endif  // USE_SHMEM
                    }
                    if (!t->special_tasks.empty()) {
                        // Each attached special task's own divide_vb call (above) pruned *that*
                        // special task's own regions_matched_to_virtual_boundary -- t's own list,
                        // populated earlier by t->setup() before any special task ran, is now stale.
                        // Write back whatever the last-processed special task ended up with (itself
                        // built by pulling from whichever predecessor add_special_task chained it to
                        // -- see decoding_task.h's own comment), so t's future parent's setup() call
                        // reads the fully-pruned, up-to-date state.
                        t->regions_matched_to_virtual_boundary = std::move(t->special_tasks.back()->regions_matched_to_virtual_boundary);
                        back_divide_walk(t, shot, (int)shot_container_id, tid, &t_out);
                    }
                    // Unit-checkpointed extraction: when a chain checkpoint resolves (and isn't itself
                    // waiting on a deferred span -- see defer_division/back_divide_walk above), divide
                    // its own vb immediately (inline, not queued -- this is what makes it safe to
                    // extract either side independently afterward) and post exactly one job, for
                    // the unit that's now closed on BOTH sides. The *other* child (the new unit
                    // just joined in) stays completely untouched -- it's still needed, intact, by
                    // the *next* fusion's own solve, and must not be extracted yet. This also
                    // covers the true root/last chain-link (itself tagged too) -- its own "closing"
                    // unit gets divided+posted here; its still-reserved right child gets posted
                    // separately once it's recognized as the root below (nothing more will ever
                    // fuse with it).
                    if (config_parallel::extract_preemptively && t->is_extraction_unit_connector && !t->defer_division) {
                        divide_vb(shot, (int)shot_container_id, t->part, t->solver, tid, /*prune_target=*/t, &t_out);
                        Task* left = t->left_child;
                        Task* closed_unit = left->is_extraction_unit_connector
                            ? left->right_child   // Fk, k>1: U(k-1), just closed on its right
                            : left;                // F1: U0, left_child IS the unit itself
                        shot.post_extraction_job(ExtractionJob{(int)shot_container_id, closed_unit}, &t_out);
                    }
                    if (t->parent != nullptr) {
                        Task* local_parent = t->parent;
                        Task* sibling = (local_parent->left_child == t) ? local_parent->right_child : local_parent->left_child;
                        if (DEBUG) t_out << "    t->parent->part=" << local_parent->part << "  t->parent->solver=" << local_parent->solver << "\n" << std::flush;
                        stolen = local_parent->try_to_steal(0);
                        t = local_parent;
                        // Try to steal sibling or descendent of sibling
                        if (!stolen && !sibling->is_fusion) {
                            stolen = sibling->try_to_steal(shot_buffer_round);
                            t = sibling;
                        }
                    } else {
                        // t is a local tree root -- parent is always nullptr here now (a CrossRankTask
                        // never sits at anyone's parent under the special-tasks-attach model).
                        stolen = false;
                        roots_i_solved.push_back({t});
                    }
                    while (!stolen && next_p_id < my_partition_task_ids.size()) {
                        t = &shot.tasks[my_partition_task_ids[next_p_id]];
                        stolen = t->try_to_steal(shot_buffer_round);
                        next_p_id += next_p_inc;
                    }
                    if (DEBUG && t != nullptr) {
                        t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen
                                << std::endl << std::flush;
                    }
                }
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_END(local_decoding);
#endif
                if (!roots_i_solved.empty()) {
                    for (auto& [root_task] : roots_i_solved) {
#ifdef SCOREP_USER_ENABLE
                        SCOREP_USER_REGION_BEGIN(solution_extraction, "Solution Extraction", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                        // Post-solve extraction for this root's own subtree. Non-preemptive: the whole
                        // tree is a balanced-over-units shape that was never divided/posted yet --
                        // post_hoc_chunk_and_post walks it top-down via the extraction-role tags.
                        // Preemptive: if the root is a chain connector, its own divide+"closing unit"
                        // post already happened above at mark_solved time -- only its still-reserved
                        // right child (the final unit, which nothing more will ever fuse with) needs
                        // posting now. Otherwise (degenerate single-unit case, no connectors were ever
                        // built), nothing has been divided/posted yet -- post root_task itself. Neither
                        // branch needs to know about special_tasks -- any seam/CRT attached to root_task
                        // already ran (and divided its own, separate vb id) during the climb above.
                        if (BARE_DEBUG) t_out << "T" << tid << " extracting solution for root part=" << root_task->part << std::endl << std::flush;
                        if (!config_parallel::extract_preemptively) {
                            post_hoc_chunk_and_post(root_task, shot, (int)shot_container_id, tid, &t_out);
                        } else if (root_task->is_extraction_unit_connector) {
                            shot.post_extraction_job(
                                ExtractionJob{(int)shot_container_id, root_task->right_child}, &t_out);
                        } else {
                            shot.post_extraction_job(ExtractionJob{(int)shot_container_id, root_task}, &t_out);
                        }
#ifdef SCOREP_USER_ENABLE
                        SCOREP_USER_REGION_END(solution_extraction);
#endif
                    } // end per-root loop

                    // Last thread (cumulative count == num_task_roots) combines and writes.
                    int n_my = (int)roots_i_solved.size();
                    int prev = shot.num_roots_done.fetch_add(n_my, std::memory_order_acq_rel);
                    if (prev + n_my == shot.num_task_roots) {
                        i_solved_last_root = true;
                    }
                }
#ifdef ENABLE_SHOT_BUFFERS
                // This thread's own static leaf stride is exhausted without ever reaching a local
                // root this shot -- help drain the extraction queue instead of idling until the
                // shot completes (real parallel extraction, not just the single root-reaching
                // thread's own serial drain; also a stress test for extraction-job independence).
                // Not attempted for NUM_BUFFERS_PER_UNIT > 1: there, an idle thread has a genuine
                // choice (keep helping this shot vs. advance into the next shot's data) that this
                // mechanism doesn't address -- see profiling-reveals-...milner.md §7's
                // active_workers/floor mechanism for that, separate future work.
                //
                // Plain busy-spin, not .wait()/.notify() -- deliberately: a spin loop always re-reads
                // the current value, so there's no "asleep" state to miss a wakeup from (see
                // i-see-here-s-the-delegated-pelican.md for the three wait/notify bugs this replaces).
                // Gated to a strided subset of threads (tid % helper_stride == 0) so the number of
                // threads simultaneously spinning on these shared atomics is bounded by how many
                // extraction units could ever exist, not by how many threads happen to be idle --
                // profiling showed naive spinning on a shared atomic scales badly at high thread
                // counts. Strided (not a contiguous tid < N prefix) so helper duty doesn't
                // concentrate on the low-tid threads that also own the next shot's earliest
                // partitions -- most threads never enter this loop at all and are immediately free
                // for the next shot the moment current_buffer_round advances.
                if (tid % helper_stride == 0) {
                    while (true) {
                        if (ExtractionJob* job = shot.try_claim_extraction_job()) {
                            process_extraction_job(shot, *job, tid, &t_out);
                            continue;
                        }
                        // Generation check: with NUM_BUFFERS_PER_UNIT == 1, this shot_container gets
                        // reused in place for the next round the moment it fully completes
                        // (ShotContainer::clear() resets extraction_jobs/extraction_posted_count/
                        // extraction_claim_cursor; num_roots_done is reset separately just before
                        // that). A thread that's lagging (e.g. lost every sibling-steal race this
                        // round) can reach this loop after some *other* thread already finished this
                        // round entirely and moved the container on to a newer one -- in that case
                        // num_roots_done/extraction_posted_count no longer even belong to this
                        // thread's round, so continuing to spin on them is unsafe (races the next
                        // round's own posts/resets). Bail: this round is already done without this
                        // thread's help, exactly as if it had helped and finished.
                        if (shot.num_roots_done.load(std::memory_order_acquire) == shot.num_task_roots
                            || shot.current_buffer_round.load(std::memory_order_acquire) != shot_buffer_round) {
                            break;
                        }
                    }
                }
#else
                // Drain this shot's entire extraction queue before checking completion.
                // num_task_roots is always 1, so exactly one thread reaches this point per shot,
                // and every job-posting site above runs on this same thread's own climb/CRT-
                // handling path strictly before this loop -- nothing can post a new job for this
                // shot after the drain starts.
                while (ExtractionJob* job = shot.try_claim_extraction_job()) {
                    process_extraction_job(shot, *job, tid, &t_out);
                }
#endif
                // report_done tells the other PE "my local window buffer is free for your next shot's
                // send_solution_to_remote_pe to reuse" -- must not fire until every job this shot
                // posted (back_divide_walk, t's own divide+post) has actually been *processed*, not
                // just posted, or the other PE could race ahead into next shot's send while this
                // shot's own extraction is still reading/writing the same regions. Every thread with a
                // non-empty crts_i_handled needs its own explicit wait here -- it can't rely on being
                // the "last root" thread (only one thread ever is) or on the queue-drain loop above
                // (SHMEM's idle-helper exits on root-completion, not job-completion).
#ifdef USE_SHMEM
                if (!crts_i_handled.empty()) {
                    while (shot.pending_extraction_jobs.load(std::memory_order_acquire) != 0) {}
                    if (BARE_DEBUG) t_out << "    reporting done" << std::endl << std::flush;
                    for (auto* crt : crts_i_handled) crt->report_done(shot_buffer_round);
                }
#endif
                if (i_solved_last_root) {
                    i_solved_last_root = false;
                    // Idle-helper threads (above) may still be concurrently draining -- the loop just
                    // above returning empty only proves nothing is left to *claim*, not that every
                    // claimed job has *finished*. Spin for that too before proceeding: write_result_
                    // and_get_next_shot's ShotContainer::clear() must never run concurrently with a
                    // still-in-flight process_extraction_job call. Plain busy-spin (this thread is
                    // never cross-round-stale -- it's the one causing the round to complete -- so no
                    // generation check needed); bounded by however long the last in-flight
                    // process_extraction_job calls, if any, take to finish, typically negligible.
                    while (shot.pending_extraction_jobs.load(std::memory_order_acquire) != 0) {
                    }
                    // Negative-weight correction: a whole-graph constant
                    // (graph.negative_weight_*_set), identical across every solver, folded in
                    // exactly once per shot here -- any of my own local solvers works (see
                    // process_extraction_job). solvers[] only covers my own partitions now, so anchor
                    // on my own first local partition rather than a fixed global index.
                    auto& any_solver = *solver_for((int)shot_container_id, my_partitions_start);
                    if (shot.num_observables > sizeof(pm::obs_int) * 8) {
                        if (!any_solver.flooder.negative_weight_detection_events.empty()) {
                            pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                any_solver, any_solver.flooder.negative_weight_detection_events);
                            any_solver.extract_paths_from_match_edges(
                                any_solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
                            any_solver.flooder.match_edges.clear();
                        }
                        for (auto& obs : any_solver.flooder.negative_weight_observables)
                            *(shot.res.obs_crossed.data() + obs) ^= 1;
                        shot.res.weight += any_solver.flooder.negative_weight_sum;
                    } else {
                        pm::MatchingResult combined{};
                        for (int ti = 0; ti < num_threads; ++ti) {
                            combined += shot.thread_results[ti];
                            shot.thread_results[ti] = {}; // reset
                        }
                        if (!any_solver.flooder.negative_weight_detection_events.empty()) {
                            combined += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                any_solver, any_solver.flooder.negative_weight_detection_events);
                        }
                        combined.obs_mask ^= any_solver.flooder.negative_weight_obs_mask;
                        combined.weight   += any_solver.flooder.negative_weight_sum;
                        if (DEBUG) t_out << "   combined obs_mask: " << combined.obs_mask << std::endl << std::flush;
                        pm::fill_bit_vector_from_obs_mask(
                            combined.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
                        shot.res.weight = combined.weight;
                    }
                    // Reset before unlocking: prevents next-shot threads racing on this counter.
                    shot.num_roots_done.store(0, std::memory_order_release);
                    shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
                }
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_END(shot_decode);
                SCOREP_USER_REGION_END(shot_iteration);
#endif
                // Move on to next shot buffer
                ++shot_id;
#ifdef ENABLE_SHOT_BUFFERS
                ++shot_container_id;
                if (shot_container_id >= NUM_BUFFERS_PER_UNIT) {
                    ++shot_buffer_round;
                    shot_container_id = 0;
                }
#else
                ++shot_buffer_round;
#endif
            }
        } catch (const std::exception& e) {
#pragma omp critical
            {
                std::cerr << "ERROR: Shot " << shot_id << " Thread " << tid << " caught exception: " << e.what()
                          << std::endl;
            }
        } catch (...) {
#pragma omp critical
            {
                std::cerr << "ERROR: Shot " << shot_id << " Thread " << tid << " caught unknown exception."
                          << std::endl;
            }
        }
    }
#ifdef USE_SHMEM
    // ensure all PE's done before exiting
    shmem_barrier_all();
#endif
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_END();
#endif
}

void pm::DecodingUnit::reset() {
    shot_buffer->reset();
    for (int i=0; i < shot_buffer->buffer.size(); ++i) {
        shot_buffer->read_shot(i, graph.node_part_id);
    }
#ifdef USE_SHMEM
    shmem_barrier_all();
#endif
}
