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
#include <omp.h>
#include <set>
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
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included,
    bool enable_correlations
#ifdef ENABLE_DRAW_FLAGS
    ,
    bool draw_frames
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
    auto user_graph =
        pm::detector_error_model_to_user_graph(detector_error_model, enable_correlations, num_distinct_weights);
#ifdef USE_SHMEM
    nodes_nelems_per_buffer = user_graph.nodes.size();
    node_ephemeral_fields_ptr = static_cast<DetectorNodeEphemeralFields*>(
        shmem_malloc(NUM_BUFFERS_PER_UNIT * nodes_nelems_per_buffer * sizeof(DetectorNodeEphemeralFields)));
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
    // atomics_ptr = static_cast<uint64_t*>(shmem_align(sizeof(uint64_t), NUM_BUFFERS_PER_UNIT * sizeof(uint64_t)));
    task_status_ptr = static_cast<uint64_t*>(shmem_align(sizeof(uint64_t), num_cross_rank_fusions * SHMEM_NUM_ATOMICS_PER_CROSS_RANK_FUSION * NUM_BUFFERS_PER_UNIT * sizeof(uint64_t)));
    regions_nelems_per_solver = graph.node_part_id.size() / graph.num_partitions * SHMEM_ARENA_BUFFER_FACTOR;
    if (regions_nelems_per_solver % 64 > 0) {
        regions_nelems_per_solver =
            (regions_nelems_per_solver / 64 + 1) * 64;  // make multiple of 64
    }
    // regions_nelems_per_solver / 64 gives number of uint64_t for bitmap
    child_edges_nelems_per_solver = regions_nelems_per_solver; // Could reduce this
    regions_matched_to_vb_nelems = graph.node_part_id.size() / graph.num_rounds * SHMEM_INTERSECTION_BUFFER_FACTOR;
    child_edges_ptr = static_cast<BlossomChild*>(shmem_malloc(child_edges_nelems_per_solver * graph.num_partitions * NUM_BUFFERS_PER_UNIT * sizeof(BlossomChild)));
    if (child_edges_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric blossom child buffer.");
    }
    task_fusion_summary_size_per_task = sizeof(FusionSummary) + 
                                        regions_matched_to_vb_nelems * sizeof(GraphFillRegion*) + 
                                        (std::max(2, config_parallel::k) * regions_nelems_per_solver / 8); /* bit map in bytes (== nelems/64 * 8) */
    task_fusion_summary_ptr = static_cast<FusionSummary*>(shmem_malloc(task_fusion_summary_size_per_task * num_cross_rank_fusions * NUM_BUFFERS_PER_UNIT));
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
#ifdef USE_SHMEM
    // std::cout << "NOTE: For now, ensure the number of partitions and the number of ranks is a power of 2. This ensures clean divisions in the task tree.\n" << std::flush;
    num_partition_units = 1; // FIX THIS??
    num_solvers_per_buffer = graph.num_partitions;
    // --- Populate my_partition_task_ids ---
    if (config_parallel::division_strategy == config_parallel::OBS) {
        const int base         = (int)graph.num_obs_patches / n_pes;
        const int rem          = (int)graph.num_obs_patches % n_pes;
        const int my_obs_start = base * pid + std::min(pid, rem);
        const int my_obs_count = base + (pid < rem ? 1 : 0);
        my_partition_task_ids.reserve(graph.p_per_obs_patch * my_obs_count);
        for (int lo = 0; lo < my_obs_count; ++lo) {
            size_t p_base = lo * (2 * graph.p_per_obs_patch - 1);
            for (size_t p = p_base; p < p_base + graph.p_per_obs_patch; ++p)
                my_partition_task_ids.push_back(p);
        }
    } else {
        const int p_base   = (int)graph.num_partitions / n_pes;
        const int my_count = p_base + (pid == n_pes - 1 ? (int)graph.num_partitions % n_pes : 0);
        my_partition_task_ids.reserve(my_count);
        for (int i = 0; i < my_count; ++i)
            my_partition_task_ids.push_back(i);
    }
    num_threads = std::min(max_threads, (int)my_partition_task_ids.size());
#else
    num_threads =
        (max_threads > graph.num_partitions) ? graph.num_partitions : max_threads;  // max num_partitions threads
    num_partition_units = graph.num_partitions / num_threads + (graph.num_partitions % num_threads > 0);
    num_solvers_per_buffer = num_threads*num_partition_units;
    for (int i = 0; i < graph.num_partitions; ++i) {
        my_partition_task_ids.push_back(i);
    }
#endif
    std::cout << "graph.num_partitions = " << graph.num_partitions << "; num_threads: " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
    // --- Build local merge-tree tasks and cross-rank fusions ---
#ifdef USE_SHMEM
    if (config_parallel::division_strategy == config_parallel::OBS) {
        build_tasks_for_obs_patch_partitioning();
    } else {
        build_tasks_for_round_partitioning();
    }
#else
    build_tasks_for_round_partitioning();
#endif
#ifdef USE_SHMEM
    // --- Allocate buffer for GraphFillRegion arenas ---
    if (DEBUG) {
        std::cout << "DEBUG: Allocating Regions\n" << std::flush;
    }
    regions_ptr = static_cast<GraphFillRegion*>(
        shmem_malloc(regions_nelems_per_solver * num_solvers_per_buffer * NUM_BUFFERS_PER_UNIT * sizeof(GraphFillRegion)));
    if (regions_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric region buffer.");
    }
    if (DEBUG) {
        std::cout << "PE" << pid << " symmetric allocations:" << std::endl
                  << "  node_ephemeral_fields_ptr: " << node_ephemeral_fields_ptr << " (" << nodes_nelems_per_buffer * sizeof(DetectorNodeEphemeralFields) << " bytes)" << std::endl
                //   << "  atomics_ptr: " << atomics_ptr << " (" << NUM_BUFFERS_PER_UNIT * sizeof(uint64_t) << " bytes)" << std::endl
                  << "  task_status_ptr: " << task_status_ptr << " (" << (graph.num_partitions-1) * 2 * sizeof(uint64_t) << " bytes)" << std::endl
                  << "  task_fusion_summary_ptr: " << task_fusion_summary_ptr << " (" << task_fusion_summary_size_per_task * (graph.num_partitions-1) << " bytes)" << std::endl
                  << "  regions_ptr: " << regions_ptr << " (" << regions_nelems_per_solver * (num_threads) * 2 * NUM_BUFFERS_PER_UNIT * sizeof(GraphFillRegion) << " bytes)" << std::endl
                  << "  regions_nelems_per_solver: " << regions_nelems_per_solver << std::endl
                  << "  child_edges_ptr: " << child_edges_ptr << " (" << child_edges_nelems_per_solver * graph.num_partitions * NUM_BUFFERS_PER_UNIT * sizeof(BlossomChild) << " bytes)" << std::endl
                  << std::flush;
    }
#endif
    // --- Build solvers ---
    build_solvers();
#ifdef ENABLE_DRAW_FLAGS
    if (draw_frames) {
        auto coords = pm::pick_coords_for_drawing_from_dem(detector_error_model, 20);
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
    // shmem_free(atomics_ptr);
    shmem_free(task_status_ptr);
    shmem_free(task_fusion_summary_ptr);
#endif
}

// Builds balanced fusion tree assuming round-based partitioning
void pm::DecodingUnit::build_tasks_for_round_partitioning() {
    if (my_partition_task_ids.empty()) return;
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
#ifdef USE_SHMEM
    const int p_offset = graph.num_partitions / n_pes * pid;
#else
    const int p_offset = 0;
#endif
    // Build a full fusion tree: less than 2*N tasks. Reserve to keep element addresses stable.
    for (int shot_container_id=0; shot_container_id < NUM_BUFFERS_PER_UNIT; shot_container_id++) {  // Lazy repeat per shot container -- CAN BE IMPROVED!
        auto& shot_container = shot_buffer->buffer[shot_container_id];
        auto& tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * graph.num_partitions - 1));
        for (int task_id : my_partition_task_ids) {
            tasks.emplace_back(task_id + p_offset);
        }
        int task_id = my_partition_task_ids.size();
        // Save odd trailing task
        int tail_idx = -1;
        if (my_partition_task_ids.size() % 2) {
            tail_idx = my_partition_task_ids.size() - 1;
        }
        // Add each level of fusions
        int last_step_starts = 0;
        int this_step_starts = task_id;
        int start = 0;
        int step = 2;
        while (start < my_partition_task_ids.size() - 1) {
            int counter = 0;
            int i;
            for (i = start; i < my_partition_task_ids.size() - 1; i += step) {
                Task* left_child = &(tasks)[last_step_starts + 2 * counter];
                int right_child_idx;
                if (last_step_starts + 2 * counter + 1 < this_step_starts) {
                    right_child_idx = last_step_starts + 2 * counter + 1;
                } else {
                    if (tail_idx >= 0) {  // fuse last one with tail
                        right_child_idx = tail_idx;
                        tail_idx = -1;
                    } else {  // save tail
                        tail_idx = last_step_starts + 2 * counter;
                        break;
                    }
                }
                Task* right_child = &(tasks)[right_child_idx];
                tasks.emplace_back(i + p_offset, left_child, right_child);
                ++task_id;
                ++counter;
            }
            if (last_step_starts + 2 * counter < this_step_starts) {  // save tail
                tail_idx = last_step_starts + 2 * counter;
            }
            last_step_starts = this_step_starts;
            this_step_starts = task_id;
            start += step / 2;
            step *= 2;
        }
    }
    if (DEBUG) {
        std::cout << "DEBUG:" << std::endl;
        for (auto& buffer : shot_buffer->buffer) {
            std::cout << "Buffer" << std::endl;
            for (Task& t : buffer.tasks) {
                std::cout << "--part: " << t.part << std::endl
                        << "  vb_left: " << t.vb_left << std::endl
                        << "  vb_right: " << t.vb_right << std::endl
                        << "  is_fusion: " << t.is_fusion << std::endl
                        << "  child_bit: " << t.child_bit << std::endl
                        << "  left_child: " << t.left_child << std::endl
                        << "  right_child: " << t.right_child << std::endl
                        << "  parent: " << t.parent << std::endl;
            }
        }
    }
#ifdef USE_SHMEM
    // --- Build cross-rank fusions (ROUND topology) ---
    //   THIS IS ALL BASED ON SIMPLE ROUND BASED FUSION ACROSS PEs
    int my_partitions_start = my_partition_task_ids.front() + p_offset;
    int my_partitions_end   = (int)my_partition_task_ids.back() + 1 + p_offset;
    for (int i=0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        shot_buffer->buffer[i].cross_rank_tasks.reserve(2);
        // ROUND always has exactly one chain top regardless of CRT count.
        shot_buffer->buffer[i].num_task_roots = 1;
        shot_buffer->buffer[i].thread_results.assign(num_threads, pm::MatchingResult{});
        shot_buffer->buffer[i].num_roots_done.store(0, std::memory_order_relaxed);
        if (pid > 0) { // Add cross-rank fusion on left
            uint64_t* task_status_p = get_task_status_ptr(i, false);
            FusionSummary* fusion_summary_p = get_fusion_summary_ptr(i, false);
            int vb = my_partitions_start-1;
            shot_buffer->buffer[i].cross_rank_tasks.emplace_back(
                vb,
                &shot_buffer->buffer[i].tasks.back(),
                false,
                (vb-config_parallel::k < -1) ? -1 : vb-config_parallel::k,
                (vb+config_parallel::k < (int)graph.num_virtual_boundaries) ? vb+config_parallel::k : (int)graph.num_virtual_boundaries,
                pid-1,
                task_status_p,     // status_shm
                task_status_p + 1, // signal_shm
                task_status_p + 2, // done_shm
                fusion_summary_p
            );
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
            shot_buffer->buffer[i].cross_rank_tasks.emplace_back(
                vb,
                &shot_buffer->buffer[i].tasks.back(),
                true,
                (vb-config_parallel::k < -1) ? -1 : vb-config_parallel::k,
                (vb+config_parallel::k < (int)graph.num_virtual_boundaries) ? vb+config_parallel::k : (int)graph.num_virtual_boundaries,
                pid+1,
                task_status_p,     // status_shm
                task_status_p + 1, // signal_shm
                task_status_p + 2, // done_shm
                fusion_summary_p
            );
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

#ifdef USE_SHMEM
void pm::DecodingUnit::build_tasks_for_obs_patch_partitioning() {
    if (DEBUG) std::cout << "DEBUG: initializing obs-patch tasks" << std::endl;

    const int K_p       = graph.p_per_obs_patch;
    const int K_vb      = graph.vb_per_obs_patch;
    const int num_seams = (int)graph.num_virtual_boundaries - K_vb * (int)graph.num_obs_patches;

    // Determine my obs patch range (same formula as constructor)
    const int base         = (int)graph.num_obs_patches / n_pes;
    const int rem          = (int)graph.num_obs_patches % n_pes;
    const int my_obs_start = base * pid + std::min(pid, rem);
    const int my_obs_count = base + (pid < rem ? 1 : 0);
    if (DEBUG) std::cout << "PE" << pid << " my_obs_start: " << my_obs_start << ", n=" << my_obs_count << "\n" << std::flush;

    // For each seam: determine the two obs patches it connects and the local partition
    // range it touches (used for vb_left/vb_right on cross-PE fusion tasks).
    struct SeamInfo { int oi, oj, vb_left, vb_right; };
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

    for (int shot_id = 0; shot_id < NUM_BUFFERS_PER_UNIT; ++shot_id) {
        auto& tasks = shot_buffer->buffer[shot_id].tasks;
        tasks.reserve(static_cast<size_t>(my_obs_count * (2*K_p - 1) + num_seams + 1));

        std::vector<Task*> obs_roots(my_obs_count, nullptr);

        // --- Balanced trees for my obs patches only ---
        for (int lo = 0; lo < my_obs_count; ++lo) {
            const int o           = my_obs_start + lo;
            const int obs_p_offset  = o * K_p;
            const int obs_vb_offset = o * K_vb;
            const int tree_start    = lo * (2*K_p - 1);  // local offset into tasks[]

            // Leaf tasks: explicit local vb bounds so is_active() uses shared slot indices
            for (int lp = 0; lp < K_p; ++lp)
                tasks.emplace_back(obs_p_offset + lp, lp - 1, lp);

            int task_id = K_p;
            int tail_idx = -1;
            if (K_p % 2) tail_idx = K_p - 1;

            int last_step_starts = 0;
            int this_step_starts = task_id;
            int start = 0;
            int step  = 2;
            while (start < K_p - 1) {
                int counter = 0;
                int i;
                for (i = start; i < K_p - 1; i += step) {
                    Task* left_child = &tasks[tree_start + last_step_starts + 2 * counter];
                    int right_child_idx;
                    if (last_step_starts + 2 * counter + 1 < this_step_starts) {
                        right_child_idx = last_step_starts + 2 * counter + 1;
                    } else {
                        if (tail_idx >= 0) {
                            right_child_idx = tail_idx;
                            tail_idx = -1;
                        } else {
                            tail_idx = last_step_starts + 2 * counter;
                            break;
                        }
                    }
                    Task* right_child = &tasks[tree_start + right_child_idx];
                    tasks.emplace_back(obs_vb_offset + i, left_child, right_child);
                    tasks.back().vb_marker = i;
                    tasks.back().vb_solver_offset = lo + my_obs_start;
                    ++task_id;
                    ++counter;
                }
                if (last_step_starts + 2 * counter < this_step_starts)
                    tail_idx = last_step_starts + 2 * counter;
                last_step_starts = this_step_starts;
                this_step_starts = task_id;
                start += step / 2;
                step  *= 2;
            }

            obs_roots[lo] = &tasks.back();
        }

        // --- Local seam tasks (skip cross-PE seams) ---
        std::vector<Task*> group_root(my_obs_count);
        for (int lo = 0; lo < my_obs_count; ++lo) {
            group_root[lo] = obs_roots[lo];
            if (DEBUG) std::cout << "  " << group_root[lo];
        }
        if (DEBUG) std::cout << "\n" << std::flush;

        std::vector<std::pair<size_t, size_t>> my_remote_seams; // pair (seam_idx, oi_local)
        my_remote_seams.reserve(num_seams);
        for (int s = 0; s < num_seams; ++s) {
            const auto& si = seam_infos[s];
            const int loi = si.oi - my_obs_start;
            const int loj = si.oj - my_obs_start;
            bool oi_local = (loi >= 0 && loi < my_obs_count);
            size_t mine = oi_local + (loj >= 0 && loj < my_obs_count);
            if (mine == 1) {
                my_remote_seams.emplace_back(s, oi_local);  // cross-PE seam
                continue;
            }
            if (mine == 0)
                continue;

            int global_vb = K_vb * (int)graph.num_obs_patches + s;
            Task* ri = group_root[loi];
            Task* rj = group_root[loj];

            if (ri == rj) {
                tasks.emplace_back(global_vb, ri, ri);
                ri->only_child = true;
            } else {
                tasks.emplace_back(global_vb, ri, rj);
            }
            tasks.back().vb_solver_offset = -global_vb + si.oi * K_p + si.vb_left + 1;
            // tasks.back().seam_vb_slot          = K_vb + s;
#ifdef ENABLE_DRAW_FLAGS
            tasks.back().left_obs_patch_id = si.oi;
            tasks.back().right_obs_patch_id = si.oj;
#endif

            Task* new_root = &tasks.back();
            for (int lo2 = 0; lo2 < my_obs_count; ++lo2)
                if (group_root[lo2] == ri || group_root[lo2] == rj)
                    group_root[lo2] = new_root;
        }

        // --- Cross-PE seam tasks ---
        auto& crt = shot_buffer->buffer[shot_id].cross_rank_tasks;
        crt.reserve(my_remote_seams.size());
        const int base_pe = (int)graph.num_obs_patches / n_pes;
        const int rem_pe  = (int)graph.num_obs_patches % n_pes;
        for (size_t seam_i = 0; seam_i < my_remote_seams.size(); ++seam_i) {
            const int s = my_remote_seams[seam_i].first;
            const bool oi_local = my_remote_seams[seam_i].second;
            auto& si = seam_infos[s];
            const int loi = si.oi - my_obs_start;
            const int loj = si.oj - my_obs_start;

            const int local_lo   = oi_local ? loi : loj;
            Task* local_child    = group_root[local_lo];  // top of local computation (post-Union-Find)
            const bool iamleft   = oi_local;              // oi < oj always (map sorted), so oi_local ↔ left
            const int remote_obs = oi_local ? si.oj : si.oi;
            
            // Check remaining seams partition ranges; constrain overlapping pairs at midpoint.
            // seam_j's task hasn't been built yet, and seam_i's emplace_back comes after this
            // loop, so mutating seam_infos in place here is safe for both sides.
            for (size_t seam_j = seam_i+1; seam_j < my_remote_seams.size(); ++seam_j) {
                const bool oi_local_j = my_remote_seams[seam_j].second;
                auto& sj = seam_infos[my_remote_seams[seam_j].first];
                const int local_loj   = oi_local_j ? sj.oi - my_obs_start : sj.oj - my_obs_start;
                const int remote_obs_j = oi_local_j ? sj.oj : sj.oi;
                // Skip unless the two seams share a local obs OR share the same remote obs.
                // The remote-obs check ensures the PE owning the two distinct local obs patches
                // makes the same adjustment as the PE owning the single shared remote obs.
                if (local_lo != local_loj && remote_obs != remote_obs_j)
                    continue;
                if (si.vb_left <= sj.vb_left) {
                    if (si.vb_right > sj.vb_left) {
                        int mid = (si.vb_right + sj.vb_left + 1) / 2;
                        if (DEBUG) std::cout << "  NOTE: PE" << pid << " adjusting seam bounds\n"
                                             << "    s=" << s << ".vb_right " << si.vb_right << " -> " << mid << "\n"
                                             << "    s=" << my_remote_seams[seam_j].first << ".vb_left " << sj.vb_left << " -> " << mid << "\n" << std::flush;
                        si.vb_right = mid;
                        sj.vb_left  = mid;
                    }
                } else {
                    if (sj.vb_right > si.vb_left) {
                        int mid = (sj.vb_right + si.vb_left + 1) / 2;
                        if (DEBUG) std::cout << "  NOTE: PE" << pid << " adjusting seam bounds\n"
                                             << "    s=" << s << ".vb_left " << si.vb_left << " -> " << mid << "\n"
                                             << "    s=" << my_remote_seams[seam_j].first << ".vb_right " << sj.vb_right << " -> " << mid << "\n" << std::flush;
                        sj.vb_right = mid;
                        si.vb_left  = mid;
                    }
                }
            }

            // Determine which PE owns remote_obs
            int other_pid;
            if (remote_obs < rem_pe * (base_pe + 1))
                other_pid = remote_obs / (base_pe + 1);
            else
                other_pid = rem_pe + (remote_obs - rem_pe * (base_pe + 1)) / base_pe;

            const int global_vb = K_vb * (int)graph.num_obs_patches + s;
            uint64_t*      sp = get_task_status_ptr_for_seam(shot_id, s);
            FusionSummary* fp = get_fusion_summary_ptr_for_seam(shot_id, s);

            crt.emplace_back(global_vb, local_child, iamleft,
                             si.vb_left, si.vb_right,
                             other_pid, sp, sp+1, sp+2, fp);
            crt.back().left_global_offset  = { (size_t)si.oi * K_p, (size_t)si.oi * K_vb };
            crt.back().right_global_offset = { (size_t)si.oj * K_p, (size_t)si.oj * K_vb };

            if (DEBUG) {
                auto& t = crt.back();
                std::cout << "PE" << pid << " OBS cross-rank task s=" << s << ":\n"
                          << "    obs_a: " << si.oi << "  obs_b: " << si.oj << "\n"
                          << "    vb=" << t.part << " vb_left=" << t.vb_left << " vb_right=" << t.vb_right << "\n"
                          << "    iamleft=" << t.iamleft << " other_pid=" << t.other_pid << "\n"
                          << "    child: " << t.child << "  me: " << &t << "  child->parent: " << t.child->parent << "\n"
                          << std::flush;
            }
        }
        // Count chain tops (parent==nullptr) as independent roots for this PE.
        // CRT constructor chain-walk already linked CRTs above local roots.
        {
            auto& sc = shot_buffer->buffer[shot_id];
            sc.thread_results.assign(num_threads, pm::MatchingResult{});
            sc.num_roots_done.store(0, std::memory_order_relaxed);
            int roots = 0;
            for (const auto& task : sc.tasks)
                if (task.parent == nullptr) ++roots;
            for (const auto& crt_t : sc.cross_rank_tasks)
                if (crt_t.parent == nullptr) ++roots;
            sc.num_task_roots = roots;
            if (DEBUG) std::cout << "PE" << pid << " num_task_roots=" << roots << "\n" << std::flush;
        }
    }

    if (DEBUG) {
        std::string tasks = "DEBUG obs-patch tasks:\n";
        for (auto& buffer : shot_buffer->buffer) {
            tasks += "Buffer\n";
            for (Task& t : buffer.tasks) {
                tasks += "--part: " + (std::string)((t.is_fusion) ? "f" : "p") + std::to_string(t.part)
                       + "  vb_left: " + std::to_string(t.vb_left)
                       + "  vb_right: " + std::to_string(t.vb_right) + "\n"
                       + "    left_child: " + std::format("{:p}",static_cast<void*>(t.left_child)) + ((t.left_child) ? "(" + (std::string)(t.left_child->is_fusion ? "f" : "p") + std::to_string(t.left_child->part) + ")" : "")
                       + "  me: " + std::format("{:p}", static_cast<void*>(&t))
                       + "  right_child: " + std::format("{:p}", static_cast<void*>(t.right_child)) + ((t.left_child) ? "(" + (std::string)(t.left_child->is_fusion ? "f" : "p") + std::to_string(t.right_child->part) + ")" : "") + "\n"
                       + "    parent: f" + std::format("{:p}", static_cast<void*>(t.parent)) + (t.child_bit == 1 ? "  left" : "  right")
                       + "  only_child: " + std::to_string(t.only_child)
                       + "\n";
            }
        }
        std::cout << tasks << std::flush;
    }
}
#endif

// Build solver for each shot container for each thread
void pm::DecodingUnit::build_solvers() {
    if (ensure_search_flooder_included || enable_correlations) {
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    }
    solvers.clear();
    solvers.reserve(static_cast<size_t>(num_solvers_per_buffer * NUM_BUFFERS_PER_UNIT));
    if (DEBUG) std::cout << "num_solvers_per_buffer: " << num_solvers_per_buffer << std::endl << std::flush;
    for (size_t idx = 0; idx < NUM_BUFFERS_PER_UNIT; ++idx) {
        for (size_t t = 0; t < num_solvers_per_buffer; ++t) {
            // Each solver shares the same MatchingGraph via shared_ptr.
#ifdef USE_SHMEM
            if (DEBUG) std::cout << "solver: " << idx*num_solvers_per_buffer + t << "  " << get_regions_ptr(idx, t) << std::endl << std::flush;
#endif
            solvers.emplace_back(
                std::make_shared<pm::Mwpm>(pm::GraphFlooder(
                    graph.graph_ptr,
                    idx
#ifdef USE_SHMEM
                    ,
                    get_regions_ptr(idx, t),
                    regions_nelems_per_solver
#endif
#ifdef ENABLE_DRAW_FLAGS
                    , &graph.node_part_id
#endif
                    )));
            solvers.back()->flooder.sync_negative_weight_observables_and_detection_events();
        }
    }
}

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
    
    // Validation Pass: scan all live regions in the k-partition send window.
    size_t solver_id = get_solver_id(shot_container_id, p_start);
    auto& solver = *solvers[solver_id];
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_to_unmatch);
    // reset to 1's for safety
    size_t bitmap_words = solvers[solver_id]->flooder.region_arena.shmem_bitmap.size() * p_k;
    for (size_t w = 0; w < bitmap_words; ++w) {
        bitmap_base[w] = ~0ULL;
    }
    // validate all blossom roots
    for (size_t i = 0; i < p_k; ++i) {
        int p_scan = p_start + i;
        auto& solver_scan = *solvers[solver_id + i];
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
                    if (valid && blossom_root->match.region) {
                        if (ptr_in_range(blossom_root->match.region, region_range)) {
                            size_t match_idx = (size_t)(blossom_root->match.region - region_range.first);
                            bitmap_base[match_idx / 64] &= ~(1ULL << (match_idx % 64));
                            if (DEBUG) t_out << "      checking matched region: " << blossom_root->match.region << std::endl << std::flush;

                            valid = check_pointers_for_self_and_all_descendents(
                                blossom_root->match.region,
                                region_range,
                                p_node_range,
                                vb_node_range,
                                &discovered_child_edges);
                        } else
                            valid = false;
                        if (DEBUG && !valid) t_out << "      blossom_root's match (" << blossom_root->match.region << ") not valid" << std::endl << std::flush; 
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

    if (DEBUG) t_out << "  isolated solution" << std::endl << std::flush;
#ifdef ENABLE_DRAW_FLAGS
    if (draw_frames) {
        draw_frame(*solvers[get_solver_id(shot_container_id, t.part)], pm::MwpmEvent::no_event(), 1001, true, omp_get_thread_num());
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
            auto& bm = solvers[solver_id + i_solver]->flooder.region_arena.shmem_bitmap;
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
        auto& solver = *solvers[solver_id + i];
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
    
    auto& shot = shot_buffer->buffer[shot_container_id];
    std::vector<std::vector<uint64_t>*> hits;
    for (int p_i = p_start; p_i <= p_end; ++p_i) {
        if (DEBUG) t_out << "  p" << p_i << std::flush;
        hits.emplace_back(&shot.partition_hits[p_i]);
        shot.i_solved_p[p_i] = false;
    }
    for (int vb_i = t.vb_left + 1 + (int)my_vb_offset; vb_i < t.vb_right + (int)my_vb_offset; ++vb_i) {
        if (DEBUG) t_out << "  vb" << vb_i << std::flush;
        hits.emplace_back(&shot.virtual_boundary_hits[vb_i]);
        shot.i_solved_vb[vb_i] = false;
    }
    shot.i_solved_vb[t.part] = false;
    if (DEBUG) t_out << std::endl << std::flush;

    // Ensure completion
    shmem_ctx_quiet(t.context_shm);

    // Cleanup sent regions
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
    if (DEBUG) t_out << "  shattered sent blossoms" << std::endl
                     << "  sent all data to " << other_pid << std::endl << std::flush;
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_END();
#endif
}

bool pm::DecodingUnit::get_solution_from_remote_pe(size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &t, std::ofstream &t_out, std::vector<uint64_t> &hitsref) {
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_FUNC_BEGIN();
#endif
    size_t other_pid = t.other_pid;
    int p_k, p_start, p_end;
    size_t remote_vb_offset = 0;
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
        remote_vb_offset       = rem_off.second;
        p_start = t.vb_left + 1 + (int)remote_p_offset;
        p_end   = t.vb_right + (int)remote_p_offset;
        p_k     = p_end - p_start + 1;
    }

    size_t solver_id = get_solver_id(shot_container_id, p_start);
    auto& solver = *solvers[solver_id];

    if (DEBUG) t_out << "    isolating solution" << std::endl << std::flush;
    int my_vb = (config_parallel::division_strategy == config_parallel::ROUND)
                    ? (t.iamleft ? t.vb_left : t.vb_right)
                    : t.part;
    if (my_vb >= 0 && my_vb < graph.num_virtual_boundaries) {
        auto& vb_bounds = graph.vb_bounds[my_vb];
        size_t num_nodes = vb_bounds.second - vb_bounds.first + 1;
        pm::DetectorNodeEphemeralFields* vb_fields_base = node_ephemeral_fields_ptr + shot_container_id * nodes_nelems_per_buffer + vb_bounds.first;
        for (int i = 0; i < num_nodes; ++i) {
            if ((vb_fields_base+i)->region_that_arrived_top) {
                if (DEBUG) t_out << "        SHATTERING: " << (vb_fields_base+i)->region_that_arrived_top << std::endl;
                res += solver.shatter_blossom_and_extract_matches((vb_fields_base+i)->region_that_arrived_top);
            }
        }
    }

    if (DEBUG) t_out << "    getting (p_start=" << p_start << ", p_k=" << p_k << ", p_end=" << p_end << ") data from " << other_pid << std::endl << std::flush;

    // Wait for signal (4 puts expected)
    shmem_wait_until(t.signal_shm, SHMEM_CMP_EQ, 4);

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
    size_t bitmap_len = solvers[solver_id]->flooder.region_arena.shmem_bitmap.size();
    for (size_t i = 0; i < p_k; ++i) {
        int p = p_start + i;
        auto& solver = *solvers[solver_id + i];
        auto& bitmap_dst = solver.flooder.region_arena.shmem_bitmap;
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
                        r->owner_arena = &solver.flooder.region_arena;
                    }
                }
            }
        }
        shot_buffer->buffer[shot_container_id].i_solved_p[p] = true;
    }
    for (int vb_i = t.vb_left + 1 + (int)remote_vb_offset; vb_i < t.vb_right + (int)remote_vb_offset; ++vb_i)
        shot_buffer->buffer[shot_container_id].i_solved_vb[vb_i] = true;
    shot_buffer->buffer[shot_container_id].i_solved_vb[t.part] = false;
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
#endif

// Core parallel decoding loop
void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
#ifdef USE_SHMEM
    if (DEBUG) {
        std::string ps = "PE" + std::to_string(pid) + " partition task ids: ";
        for (auto p : my_partition_task_ids)
            ps += std::to_string(p) + " ";
        std::cout << ps << std::endl << std::flush;
    }
    // shmem_barrier_all(); // needed to avoid races on symmetric data (signals)
#endif
#pragma omp parallel
    {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_DEFINE(local_decoding);
        SCOREP_USER_REGION_DEFINE(cross_rank_fusion);
        SCOREP_USER_REGION_DEFINE(solution_extraction);
#endif
        const int tid = omp_get_thread_num();
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
        int shot_id = shot_buffer_round * NUM_BUFFERS_PER_UNIT;
#ifdef USE_SHMEM
        // Thread-local list of local roots solved during the steal loop.
        // Cleared at the start of each shot; a thread may solve multiple roots when
        // it exhausts all its assigned partition leaves (via next_p_id).
        struct RootInfo { Task* task; int solver_id; };
        std::vector<RootInfo> roots_i_solved;
#endif
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
                while (shot_current_buffer_round < shot_buffer_round) {  // wait
                    if (shot_current_buffer_round < 0) {
                        break;
                    }
                    _mm_pause();
                    shot_current_buffer_round = shot.current_buffer_round.load();
                }
                if (shot_current_buffer_round < 0) {
                    break;
                }
                // we divide partition tasks into sets based on the number of threads available
                Task* t = &shot.tasks[my_partition_task_ids[tid]];
                size_t next_p_inc = num_threads; // cannot inc by two for 1 threads --- does not work for odd
                size_t next_p_id = tid+next_p_inc;
                int solver_id = get_solver_id(shot_container_id, t->part, tid);
                if (DEBUG) t_out << "solvers[" << solver_id << "]\n";
                bool stolen = t->try_to_steal_leaf(shot_buffer_round);
#ifdef USE_SHMEM
                roots_i_solved.clear();
#else
                bool i_solved_root = false;
#endif
#ifdef SCOREP_USER_ENABLE
                SCOREP_USER_REGION_BEGIN(local_decoding, "Local Decoding", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                while (stolen) {  // Got task
                    if (DEBUG) t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl << std::flush;
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    // Solve task
                    t->setup();
                    if (DEBUG) t_out << "  solvers[" << solver_id << "]\n";
                    pm::Mwpm& solver = *solvers[solver_id];
                    solver.prepare_for_task(t, shot_id);
#ifdef USE_SHMEM
                    if (DEBUG) t_out << "  solver bounds: " << solver.flooder.vb_left << "(vb_left) " << solver.flooder.vb_right << " (vb_right)" << std::endl << std::flush;
#ifdef ENABLE_DRAW_FLAGS
                    if (draw_frames && config_parallel::division_strategy == config_parallel::OBS) {
                        size_t K_p = graph.num_partitions / graph.num_obs_patches;
                        size_t K_vb = K_p - 1;
                        if (t->is_fusion && t->left_obs_patch_id >= 0 && t->right_obs_patch_id >= 0) {
                            // Local seam task
                            solver.flooder.p_offsets  = { (size_t)t->left_obs_patch_id * K_p,  (size_t)t->right_obs_patch_id * K_p  };
                            solver.flooder.vb_offsets = { (size_t)t->left_obs_patch_id * K_vb, (size_t)t->right_obs_patch_id * K_vb };
                        } else {
                            // Leaf partition
                            size_t patch = (t->is_fusion) ? t->part / K_vb : t->part / K_p;
                            solver.flooder.p_offsets  = { K_p  * patch };
                            solver.flooder.vb_offsets = { K_vb * patch };
                        }
                    }
#endif
#endif
                    pm::process_timeline_until_completion(
                        solver,
                        hitsref,
#ifdef ENABLE_DRAW_FLAGS
                        draw_frames,
#endif
                        true,
                        tid);
#ifdef USE_SHMEM
                    if (DEBUG) t_out << "  solved " << (t->is_fusion ? "f" : "p") << t->part << std::endl << std::flush;
                    (t->is_fusion) ? shot.i_solved_vb[t->part] = true : shot.i_solved_p[t->part] = true;
#endif
                    t->mark_solved();
                    if (t->parent != nullptr
#ifdef USE_SHMEM
                        && !t->parent->is_cross_rank_fusion
#endif
                    ) {
                        Task* local_parent = static_cast<Task*>(t->parent);
                        bool iamleft = t->child_bit == 1;
                        Task* sibling = (iamleft) ? local_parent->right_child : local_parent->left_child;
#ifdef USE_SHMEM
                        solver_id = get_solver_id(shot_container_id, local_parent->part + local_parent->vb_solver_offset, tid);
                        if (DEBUG) t_out << "    t->parent->part=" << local_parent->part << "  t->parent->vb_solver_offset=" << local_parent->vb_solver_offset << "\n" << std::flush;
#else
                        solver_id = get_solver_id(shot_container_id, local_parent->part, tid);
#endif
                        stolen = t->try_to_steal_parent();
                        t = local_parent;
                        // Try to steal sibling or descendent of sibling
                        if (!stolen && !sibling->is_fusion) {
                            solver_id = get_solver_id(shot_container_id, sibling->part, tid);
                            stolen = sibling->try_to_steal_leaf(shot_buffer_round);
                            t = sibling;
                        }
                    } else {
                        // t is a local tree root: parent==nullptr or parent is a CrossRankTask
                        stolen = false;
#ifdef USE_SHMEM
                        roots_i_solved.push_back({t, solver_id});
#else
                        i_solved_root = true;
#endif
                    }
                    while (!stolen && next_p_id < my_partition_task_ids.size()) {
                        t = &shot.tasks[my_partition_task_ids[next_p_id]];
                        solver_id = get_solver_id(shot_container_id, t->part, tid);
                        stolen = t->try_to_steal_leaf(shot_buffer_round);
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
#ifdef USE_SHMEM
                if (!roots_i_solved.empty()) {
                    // Collect all CRTs handled across roots this thread solved this shot
                    std::vector<CrossRankTask*> crts_i_handled;
                    pm::MatchingResult& my_result = shot.thread_results[tid];

                    for (auto& [root_task, root_solver_id] : roots_i_solved) {
                        auto& root_solver = *solvers[root_solver_id];
                        root_solver.flooder.match_edges.clear();
#ifdef SCOREP_USER_ENABLE
                        SCOREP_USER_REGION_BEGIN(cross_rank_fusion, "Cross Rank Fusion", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                        // Walk the CRT chain above this local root
                        std::vector<CrossRankTask*> root_crts;
                        TaskBase* chain_node = root_task->parent;
                        while (chain_node != nullptr) {
                            auto* crt = static_cast<CrossRankTask*>(chain_node);
                            if (DEBUG) t_out << "Trying cross-rank fusion vb=" << crt->part
                                             << " iamleft=" << crt->iamleft
                                             << " other_pid=" << crt->other_pid << "\n" << std::flush;
                            if (crt->try_to_steal(pid, t_out)) {
                                if (BARE_DEBUG) t_out << "  Stole CRT with " << crt->other_pid << std::endl << std::flush;
                                crt->setup();
                                int crt_sid = get_solver_id(shot_container_id,
                                    (int)((crt->iamleft) ? crt->left_global_offset.first
                                                         : crt->right_global_offset.first));
                                auto& crt_solver = *solvers[crt_sid];
                                crt_solver.prepare_for_task(crt, shot_id);
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames && config_parallel::division_strategy == config_parallel::OBS) {
                                    crt_solver.flooder.p_offsets  = { crt->left_global_offset.first,  crt->right_global_offset.first  };
                                    crt_solver.flooder.vb_offsets = { crt->left_global_offset.second, crt->right_global_offset.second };
                                }
                                // if (draw_frames) draw_frame(crt_solver, pm::MwpmEvent::no_event(), 1000, true, tid);
#endif
                                auto& crt_hitsref = shot.virtual_boundary_hits[crt->part];
                                get_solution_from_remote_pe(shot_container_id, my_result, *crt, t_out, crt_hitsref);
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
                                shot.i_solved_vb[crt->part] = true;
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames) draw_frame(crt_solver, pm::MwpmEvent::no_event(), 1001, true, tid);
#endif
                            } else {
                                if (BARE_DEBUG) t_out << "  Sending CRT data to " << crt->other_pid << std::endl;
                                crt->setup();
#ifdef ENABLE_DRAW_FLAGS
                                if (draw_frames) {
                                    auto& dbg_solver = *solvers[get_solver_id(shot_container_id, crt->part)];
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
                            root_crts.push_back(crt);
                            crts_i_handled.push_back(crt);
                            chain_node = crt->parent;
                        }
#ifdef SCOREP_USER_ENABLE
                        SCOREP_USER_REGION_END(cross_rank_fusion);
                        SCOREP_USER_REGION_BEGIN(solution_extraction, "Solution Extraction", SCOREP_USER_REGION_TYPE_COMMON);
#endif
                        if (BARE_DEBUG) t_out << "T" << tid << " extracting solution for root part=" << root_task->part << std::endl << std::flush;

                        // Shatter this root's subtree via task-tree traversal.
                        // Each thread only visits nodes IT solved, avoiding cross-thread races on solvers.
                        if (shot.num_observables > sizeof(pm::obs_int) * 8) {
                            // Extended observables: accumulate under omp critical since shot.res is shared
#pragma omp critical
                            {
                                std::vector<Task*> to_visit = {root_task};
                                while (!to_visit.empty()) {
                                    Task* curr = to_visit.back(); to_visit.pop_back();
                                    if (!curr->is_fusion) {
                                        if (shot.i_solved_p[curr->part]) {
                                            pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                                root_solver, shot.partition_hits[curr->part]);
                                            shot.i_solved_p[curr->part] = false;
                                        }
                                    } else {
                                        if (shot.i_solved_vb[curr->part]) {
                                            pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                                root_solver, shot.virtual_boundary_hits[curr->part]);
                                            shot.i_solved_vb[curr->part] = false;
                                        }
                                        if (curr->left_child) to_visit.push_back(curr->left_child);
                                        if (!curr->only_child && curr->right_child) to_visit.push_back(curr->right_child);
                                    }
                                }
                                for (auto* crt : root_crts) {
                                    // NOT BACKWARD COMPATBILE WITH ROUND
                                    auto& remote_offset = (crt->iamleft) ? crt->left_global_offset : crt->right_global_offset;
                                    size_t p = crt->vb_left + 1 + remote_offset.first; // inclusive
                                    const size_t p_end = crt->vb_right + 1 + remote_offset.first; // exclusive
                                    size_t vb = crt->vb_left + 1 + remote_offset.second; // inclusive
                                    const size_t vb_end = crt->vb_right + remote_offset.second; // exclusive
                                    if (shot.i_solved_vb[crt->part]) {
                                        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                            root_solver, shot.virtual_boundary_hits[crt->part]);
                                        shot.i_solved_vb[crt->part] = false;
                                        for (; p < p_end; ++p) {
                                            pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                                root_solver, shot.partition_hits[p]);
                                        }
                                        for (; vb < vb_end; ++vb) {
                                            pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                                root_solver, shot.virtual_boundary_hits[vb]);
                                        }
                                    }
                                }
                                if (!root_solver.flooder.negative_weight_detection_events.empty())
                                    shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                        root_solver, root_solver.flooder.negative_weight_detection_events);
                                root_solver.extract_paths_from_match_edges(
                                    root_solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
                                for (auto& obs : root_solver.flooder.negative_weight_observables)
                                    *(shot.res.obs_crossed.data() + obs) ^= 1;
                                shot.res.weight += root_solver.flooder.negative_weight_sum;
                            }
                        } else {
                            // Bit-packed case: accumulate thread-locally into thread_results[tid]
                            std::vector<Task*> to_visit = {root_task};
                            while (!to_visit.empty()) {
                                Task* curr = to_visit.back(); to_visit.pop_back();
                                if (!curr->is_fusion) {
                                    if (shot.i_solved_p[curr->part]) {
                                        if (DEBUG) t_out << "  p" << curr->part;
                                        my_result += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                            root_solver, shot.partition_hits[curr->part]);
                                        shot.i_solved_p[curr->part] = false;
                                    }
                                } else {
                                    if (shot.i_solved_vb[curr->part]) {
                                        if (DEBUG) t_out << "  vb" << curr->part;
                                        my_result += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                            root_solver, shot.virtual_boundary_hits[curr->part]);
                                        shot.i_solved_vb[curr->part] = false;
                                    }
                                    if (curr->left_child) to_visit.push_back(curr->left_child);
                                    if (curr->right_child && curr->right_child != curr->left_child) to_visit.push_back(curr->right_child);
                                }
                            }
                            for (auto* crt : root_crts) {
                                if (shot.i_solved_vb[crt->part]) {
                                    if (DEBUG) t_out << "  crt_vb" << crt->part;
                                    my_result += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                        root_solver, shot.virtual_boundary_hits[crt->part]);
                                    shot.i_solved_vb[crt->part] = false;
                                    // NOT BACKWARD COMPATiBLE WITH ROUND
                                    // Need to use different index calculations for p and vb
                                    auto& remote_offset = (crt->iamleft) ? crt->right_global_offset : crt->left_global_offset;
                                    size_t p = crt->vb_left + 1 + remote_offset.first;
                                    const size_t p_end = crt->vb_right + 1 + remote_offset.first; // exclusive
                                    size_t vb = crt->vb_left + 1 + remote_offset.second;
                                    const size_t vb_end = crt->vb_right + remote_offset.second; // exclusive
                                    for (; p < p_end; ++p) {
                                        if (DEBUG) t_out << "  p" << p;
                                        my_result += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                            root_solver, shot.partition_hits[p]);
                                    }
                                    for (; vb < vb_end; ++vb) {
                                        if (DEBUG) t_out << "  vb" << vb;
                                        my_result += pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                            root_solver, shot.virtual_boundary_hits[vb]);
                                    }
                                }
                                
                            }
                            if (!root_solver.flooder.negative_weight_detection_events.empty()) {
                                if (DEBUG) t_out << "  negative detection events\n" << std::flush;
                                my_result += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                    root_solver, root_solver.flooder.negative_weight_detection_events);
                            }
                            my_result.obs_mask ^= root_solver.flooder.negative_weight_obs_mask;
                            my_result.weight   += root_solver.flooder.negative_weight_sum;
                            if (DEBUG) t_out << "\n  obs_mask: " << my_result.obs_mask << std::endl << std::flush;
                        }
#ifdef SCOREP_USER_ENABLE
                        SCOREP_USER_REGION_END(solution_extraction);
#endif
                    } // end per-root loop

                    // Debug: verify all solvers clean after shattering
                    for (int p = 0; p < graph.num_partitions; ++p) {
                        auto& chk = *solvers[get_solver_id(shot_container_id, p)];
                        for (auto word : chk.flooder.region_arena.shmem_bitmap) {
                            if (word != ~0ULL) {
                                if (DEBUG) t_out << "  ERROR: solver p" << p << " not empty\n" << std::flush;
                                else std::cout << "  ERROR: solver p" << p << " not empty\n" << std::flush;
                                break;
                            }
                        }
                    }

                    // Report done before waiting: all-report-then-all-wait avoids deadlock
                    // for same-PE-pair CRTs (e.g., ROUND left+right CRTs).
                    if (BARE_DEBUG) t_out << "    reporting done" << std::endl << std::flush;
                    for (auto* crt : crts_i_handled) crt->report_done(shot_buffer_round);
                    if (BARE_DEBUG) t_out << "    waiting until other PEs done" << std::endl << std::flush;
                    for (auto* crt : crts_i_handled) crt->wait_until_done(pid, shot_buffer_round);
                    if (BARE_DEBUG) t_out << "    all done" << std::endl << std::flush;

                    // Last thread (cumulative count == num_task_roots) combines and writes.
                    int n_my = (int)roots_i_solved.size();
                    int prev = shot.num_roots_done.fetch_add(n_my, std::memory_order_acq_rel);
                    if (prev + n_my == shot.num_task_roots) {
                        if (shot.num_observables <= sizeof(pm::obs_int) * 8) {
                            pm::MatchingResult combined{};
                            for (int ti = 0; ti < num_threads; ++ti) {
                                combined += shot.thread_results[ti];
                                shot.thread_results[ti] = {}; // reset
                            }
                            if (DEBUG) t_out << "   combined obs_mask: " << combined.obs_mask << std::endl << std::flush;
                            pm::fill_bit_vector_from_obs_mask(
                                combined.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
                            shot.res.weight = combined.weight;
                        }
                        // Reset before unlocking: prevents next-shot threads racing on this counter.
                        shot.num_roots_done.store(0, std::memory_order_release);
                        shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
                    }
                }
#else
                if (i_solved_root) {
                    auto& solver = *solvers[solver_id];
                    solver.flooder.match_edges.clear();
                    pm::MatchingResult bit_packed_res;
                    if (shot.num_observables > sizeof(pm::obs_int) * 8) {
                        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                            solver, shot.sparse_shot.hits);
                        if (!solver.flooder.negative_weight_detection_events.empty())
                            shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                solver, solver.flooder.negative_weight_detection_events);
                        solver.extract_paths_from_match_edges(
                            solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
                        for (auto& obs : solver.flooder.negative_weight_observables)
                            *(shot.res.obs_crossed.data() + obs) ^= 1;
                        shot.res.weight += solver.flooder.negative_weight_sum;
                    } else {
                        bit_packed_res =
                            pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                solver, shot.sparse_shot.hits);
                        if (!solver.flooder.negative_weight_detection_events.empty())
                            bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                solver, solver.flooder.negative_weight_detection_events);
                        bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
                        pm::fill_bit_vector_from_obs_mask(
                            bit_packed_res.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
                        shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
                    }
                    shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
                }
#endif
                // Move on to next shot buffer
                ++shot_id;
#if NUM_BUFFERS_PER_UNIT > 1
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
