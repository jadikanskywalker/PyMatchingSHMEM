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
#include <vector>

#include "pymatching/sparse_blossom/driver/user_graph.h"

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
    // --- Allocate sychronization & summary memory ---
    atomics_ptr = static_cast<uint64_t*>(shmem_align(sizeof(uint64_t), NUM_BUFFERS_PER_UNIT * sizeof(uint64_t)));
    task_status_ptr = static_cast<uint64_t*>(shmem_align(sizeof(uint64_t), SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER * 2 * NUM_BUFFERS_PER_UNIT * sizeof(uint64_t))); // Space for two cross-rank fusions
    regions_nelems_per_solver = graph.node_part_id.size() / graph.num_partitions * SHMEM_ARENA_BUFFER_FACTOR;
    if (regions_nelems_per_solver >= 64) {
        regions_nelems_per_solver =
            (regions_nelems_per_solver / 64) * 64;  // make multiple of 64
    } else {
        regions_nelems_per_solver = 64;
    }
    // regions_nelems_per_solver / 64 gives number of uint64_t for bitmap
    child_edges_nelems_per_solver = regions_nelems_per_solver;
    size_t regions_matched_to_vb_nelems = graph.node_part_id.size() / graph.num_rounds * SHMEM_INTERSECTION_BUFFER_FACTOR;
    child_edges_ptr = static_cast<BlossomChild*>(shmem_malloc(child_edges_nelems_per_solver * graph.num_partitions * NUM_BUFFERS_PER_UNIT * sizeof(BlossomChild)));
    if (child_edges_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric blossom child buffer.");
    }
    task_fusion_summary_size_per_task = sizeof(FusionSummary) + 
                                        regions_matched_to_vb_nelems * sizeof(GraphFillRegion*) + 
                                        (regions_nelems_per_solver / 8); /* bit map in bytes (== nelems/64 * 8) */
    task_fusion_summary_ptr = static_cast<FusionSummary*>(shmem_malloc(task_fusion_summary_size_per_task * SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER * NUM_BUFFERS_PER_UNIT));
    if (atomics_ptr == nullptr || task_status_ptr == nullptr || task_fusion_summary_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric atomics buffer.");
    }
    // closest_nb_nelems_per_buffer = graph.node_part_id.size() / graph.num_rounds * SHMEM_INTERSECTION_BUFFER_FACTOR;
    // closest_nb_ptr = static_cast<std::pair<size_t, size_t>*>(
    //     shmem_malloc(closest_nb_nelems_per_buffer * NUM_BUFFERS_PER_UNIT * sizeof(std::pair<size_t, size_t>)));
#endif
    // --- Create shot buffer ---
    shot_buffer = std::make_shared<pm::ShotBuffer>(
#ifdef USE_SHMEM
        atomics_ptr,
#endif
        std::move(reader),
        std::move(writer),
        graph.num_partitions,
        graph.num_virtual_boundaries,
        graph.graph_ptr->num_observables);
    // --- Read first NUM_BUFFERS_PER_UNIT shots ---
    if (DEBUG) {
        std::cout << "DEBUG: Reading shots \n" << std::flush;
    }
    for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        shot_buffer->read_shot(i, graph.node_part_id);
    }
    // --- Set num_threads ---
    if (DEBUG) {
        std::cout << "DEBUG: Setting num threads\n" << std::flush;
    }
    int max_threads = omp_get_max_threads();
#ifdef USE_SHMEM
    int p_per_pe = graph.num_partitions / n_pes;
    if (n_pes == 2 && config_parallel::k > graph.num_partitions/2) {
        config_parallel::k = graph.num_partitions/2;
        std::cout << "NOTE: k set to " << config_parallel::k << " for 2 PEs" << std::endl << std::flush;
    } else {
        if (p_per_pe < 2) {
            throw std::invalid_argument("The number of partitions per PE should be >= 2 for more than 2 ranks.");
        } else if (p_per_pe == 2 && config_parallel::k > 2) {
            // Bounding k for correctness
            config_parallel::k = 2; // (n_pes > 2) ? 1 : 2;
            std::cout << "NOTE: k bounded to " << config_parallel::k << std::endl << std::flush;
        } else if (p_per_pe / 2 < config_parallel::k) {
            // Bounding k for correctness
            config_parallel::k = p_per_pe / 2;
            std::cout << "NOTE: k bounded to " << config_parallel::k << std::endl << std::flush;
        }
    }
    std::cout << "NOTE: For now, ensure the number of partitions and the number of ranks is a power of 2. This ensures clean divisions in the task tree.\n" << std::flush;
    num_threads = (max_threads < p_per_pe) ? max_threads : p_per_pe; 
    num_partition_units = 1; // FIX THIS??
    num_solvers_per_buffer = graph.num_partitions;
    // --- Populate my_partitions (should probably make multiple of power of 2)
    //   THIS IS ALL BASED ON SIMPLE ROUND_BASE FUSION ACROSS PEs
    my_partitions.reserve(p_per_pe);
    int my_partitions_start = p_per_pe * pid;
    int my_partitions_end = my_partitions_start + p_per_pe;
    if (pid == n_pes-1) {
        my_partitions_end += graph.num_partitions % n_pes;
    }
    for (int i = my_partitions_start; i < my_partitions_end; ++i) {
        my_partitions.push_back(i);
    }
#else
    num_threads =
        (max_threads > graph.num_partitions) ? graph.num_partitions : max_threads;  // max num_partitions threads
    num_partition_units = graph.num_partitions / num_threads + (graph.num_partitions % num_threads > 0);
    num_solvers_per_buffer = num_threads*num_partition_units;
    for (int i = 0; i < graph.num_partitions; ++i) {
        my_partitions.push_back(i);
    }
#endif
    std::cout << "graph.num_partitions = " << graph.num_partitions << "; num_threads: " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
    // --- Build my tasks ---
    build_tasks_for_round_partitioning();
#ifdef USE_SHMEM
    // --- Build my cross-rank fusions ---
    //   THIS IS ALL BASED ON SIMPLE ROUND BASED FUSION ACROSS PEs
    for (int i=0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        shot_buffer->buffer[i].cross_rank_tasks.reserve(2);
        if (pid > 0) { // Add cross-rank fusion on left
            uint64_t* task_status_p = get_task_status_ptr(i, false);
            FusionSummary* fusion_summary_p = get_fusion_summary_ptr(i, false);
            int vb = my_partitions_start-1;
            shot_buffer->buffer[i].cross_rank_tasks.emplace_back(
                vb,
                &shot_buffer->buffer[i].tasks.back(),
                false,
                (vb-config_parallel::k < -1) ? -1 : vb-config_parallel::k,
                (vb+config_parallel::k < graph.num_virtual_boundaries) ? vb+config_parallel::k : graph.num_virtual_boundaries,
                pid-1,
                task_status_p,
                task_status_p + 1,
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
                (vb+config_parallel::k < graph.num_virtual_boundaries) ? vb+config_parallel::k : graph.num_virtual_boundaries,
                pid+1,
                task_status_p,
                task_status_p + 1,
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
                    << "    fusion_summary_shm: " << t.fusion_summary_shm << std::endl << std::flush;
            }
        }
    }
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
                  << "  atomics_ptr: " << atomics_ptr << " (" << NUM_BUFFERS_PER_UNIT * sizeof(uint64_t) << " bytes)" << std::endl
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
    shmem_free(atomics_ptr);
    shmem_free(task_status_ptr);
    shmem_free(task_fusion_summary_ptr);
#endif
}

// Builds balanced fusion tree assuming round-based partitioning
void pm::DecodingUnit::build_tasks_for_round_partitioning() {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // Build a full fusion tree: less than 2*N tasks. Reserve to keep element addresses stable.
    for (int shot_container_id=0; shot_container_id < NUM_BUFFERS_PER_UNIT; shot_container_id++) {  // Lazy repeat per shot container -- CAN BE IMPROVED!
        auto& shot_container = shot_buffer->buffer[shot_container_id];
        auto& tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * graph.num_partitions - 1));
        for (int task_id : my_partitions) {
            tasks.emplace_back(task_id);
        }
        int task_id = my_partitions.size();
        // Save odd trailing task
        int tail_idx = -1;
        if (my_partitions.size() % 2) {
            tail_idx = my_partitions.size() - 1;
        }
        // Add each level of fusions
        int last_step_starts = 0;
        int this_step_starts = task_id;
        int start = 0;
        int step = 2;
        while (start < my_partitions.size() - 1) {
            int counter = 0;
            int i;
            for (i = start; i < my_partitions.size() - 1; i += step) {
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
                tasks.emplace_back(i + my_partitions[0], left_child, right_child
// #ifdef USE_SHMEM
//                     , task_status_ptr + (shot_container_id * (graph.num_partitions-1) + i) * 2,
//                     task_status_ptr + (shot_container_id * (graph.num_partitions-1) + i) * 2 + 1,
//                     0, 0 // WRONG NEED TO FIX
// #endif
                );
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
// #ifdef USE_SHMEM
//         task_id = graph.num_partitions;
//         for (size_t i=0; i < graph.num_partitions-1; ++i) {
//             Task* left_child = &(tasks[i]);
//             Task* right_child = &(tasks[i+1]);
//             tasks.emplace_back(
//                 i,
//                 left_child,
//                 right_child,
//                 task_status_ptr + (shot_container_id * (graph.num_partitions-1) + i) * 2,
//                 task_status_ptr + (shot_container_id * (graph.num_partitions-1) + i) * 2 + 1,
//                 i,
//                 i+1
//             );
//         }
// #endif
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
}

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

// bool check_pointers_for_self_and_all_descendents(pm::GraphFillRegion* root, std::pair<pm::GraphFillRegion*, pm::GraphFillRegion*> region_range, std::pair<pm::DetectorNode*, pm::DetectorNode*> node_range) {
//     if (!ptr_in_range(root, region_range)) return false;
//     for (pm::DetectorNode* node : root->shell_area) {
//         if (!ptr_in_range(node, node_range)) return false;
//     }
//     for (pm::RegionEdge edge : root->blossom_children) {
//         if (!check_pointers_for_self_and_all_descendents(edge.region, region_range, node_range)) return false;
//     }
//     return true;
// }

bool check_pointers_for_self_and_all_descendents(
    pm::GraphFillRegion* root,
    std::pair<pm::GraphFillRegion*, pm::GraphFillRegion*> region_range,
    std::pair<pm::DetectorNode*, pm::DetectorNode*> node_range,
    std::vector<pm::BlossomChild>* discovered_edges) {
    if (!ptr_in_range(root, region_range)) return false;
    for (pm::DetectorNode* node : root->shell_area) {
        if (!ptr_in_range(node, node_range)) return false;
    }
    for (pm::RegionEdge edge : root->blossom_children) {
        if (!check_pointers_for_self_and_all_descendents(
                edge.region,
                region_range,
                node_range,
                discovered_edges)) {
            return false;
        }
        discovered_edges->push_back(pm::BlossomChild{(size_t)(root - region_range.first), edge});
    }
    return true;
}

void pm::DecodingUnit::send_solution_to_remote_pe(size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &t, std::ofstream &t_out) {
    size_t other_pid = t.other_pid;
    int k = config_parallel::k;
    int partition_start = (t.iamleft) ? t.part - k + 1 : t.part + 1;
    if (partition_start < 0) {
        k += partition_start; // reduce k
        partition_start = 0;
    }
    if (partition_start + k > graph.num_partitions) {
        k = graph.num_partitions - partition_start;
    }

    if (DEBUG)
        t_out << "    sending (p" << partition_start << ", k=" << k << ") to " << other_pid << std::endl;

    // We use the LOCAL slot buffer to construct the payload, then PUT it to the REMOTE slot buffer.
    FusionSummary*& fusion_summary_base = t.fusion_summary_shm;
    
    // 1. Populate FusionSummary
    fusion_summary_base->regions_to_unmatch_size = t.regions_to_unmatch.size();
    fusion_summary_base->regions_ptr_base = regions_ptr;
    fusion_summary_base->static_nodes_base = graph.graph_ptr->nodes.data();
    
    // Copy regions to unmatch (which are the regions matched to the vb in the child task)
    for (size_t i=0; i<t.regions_to_unmatch.size(); ++i) {
        fusion_summary_base->regions_to_unmatch[i] = t.regions_to_unmatch[i];
        if (DEBUG) t_out << "      " << t.regions_to_unmatch[i] << std::endl << std::flush;
    }

    // 2. Isolate Solution
    auto* regions_start_ptr = get_regions_ptr(shot_container_id, partition_start);
    size_t regions_total_bytes = k * regions_nelems_per_solver * sizeof(GraphFillRegion);

    auto& nodes_start_bounds = graph.partition_bounds[partition_start];
    auto& nodes_end_bounds = graph.partition_bounds[partition_start + k - 1];
    size_t nodes_nelems_total = nodes_end_bounds.second - nodes_start_bounds.first + 1;

    size_t isolate_start_bound = (t.iamleft) ? nodes_start_bounds.first : graph.vb_bounds[partition_start-1].first;
    size_t isolate_end_bound = (t.iamleft) ? graph.vb_bounds[partition_start + k - 1].second + 1 : nodes_end_bounds.second + 1;
    if (DEBUG) {
        t_out << "  isolating solution" << std::endl
              << "    isolate_start_bound: " << isolate_start_bound << std::endl
              << "    node_start_bound: " << nodes_start_bounds.first << std::endl
              << "    node_end_bound: " << nodes_end_bounds.second << " (" << nodes_nelems_total << ")" << std::endl
              << "    isolate_end_bound: " << isolate_end_bound << std::endl << std::flush;
    
    }
    // isolate_solution_before_sending(
    //     res,
    //     solver,
    //     shot_buffer->buffer[shot_container_id], shot_container_id,
    //     (std::pair<size_t, size_t>){partition_start, partition_start+k},
    //     (std::pair<pm::GraphFillRegion *, pm::GraphFillRegion *>){regions_start_ptr, regions_start_ptr + k * regions_nelems_per_solver},
    //     (std::pair<pm::DetectorNode*, pm::DetectorNode*>){graph.graph_ptr->nodes.data() + isolate_start_bound, graph.graph_ptr->nodes.data() + isolate_end_bound},
    //     t_out);
    auto& solver = *solvers[get_solver_id(shot_container_id, partition_start)];
    auto& shot = shot_buffer->buffer[shot_container_id];
    std::vector<std::vector<uint64_t>*> hits;
    size_t p = partition_start;
    for (; p < partition_start + k - 1; ++p) {
        if (DEBUG) t_out << "  p" << p << std::flush;
        hits.emplace_back(&shot.partition_hits[p]);
        if (DEBUG) t_out << "  vb" << p << std::flush;
        hits.emplace_back(&shot.virtual_boundary_hits[p]);
        // relinquish sent solution
        shot.i_solved_p[p] = false;
        shot.i_solved_vb[p] = false;
    }
    if (DEBUG) t_out << "  p" << p << std::endl << std::flush;
    hits.emplace_back(&shot.partition_hits[p]);
    shot.i_solved_p[p] = false;

    std::pair<pm::GraphFillRegion *, pm::GraphFillRegion *> region_range = {regions_start_ptr, regions_start_ptr + k * regions_nelems_per_solver};
    std::pair<pm::DetectorNode*, pm::DetectorNode*> node_range = {graph.graph_ptr->nodes.data() + isolate_start_bound, graph.graph_ptr->nodes.data() + isolate_end_bound};
    BlossomChild* child_edges_buff_base = get_child_edges_ptr(shot_container_id, partition_start);
    size_t child_edges_counter = 0;
    std::vector<pm::BlossomChild> discovered_child_edges;
    discovered_child_edges.reserve(64);
    
    std::vector<pm::GraphFillRegion*> blossom_roots_checked;
    
    // Validation Pass: scan all live regions in the k-partition send window.
    for (size_t i = 0; i < k; ++i) {
        int p_scan = partition_start + i;
        auto& solver_scan = *solvers[get_solver_id(shot_container_id, p_scan)];
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
                pm::GraphFillRegion* blossom_root = r->blossom_parent_top ? r->blossom_parent_top : r;

                // If we've already validated this root, skip re-checking.
                if (std::find(blossom_roots_checked.begin(), blossom_roots_checked.end(), blossom_root) != blossom_roots_checked.end()) {
                    continue;
                }

                if (DEBUG) t_out << "    checking blossom_root: " << blossom_root << std::endl << std::flush;
                discovered_child_edges.clear();

                bool valid = check_pointers_for_self_and_all_descendents(
                    blossom_root,
                    region_range,
                    node_range,
                    &discovered_child_edges);

                // Also check match partner if present.
                if (valid && blossom_root->match.region) {
                    if (DEBUG) t_out << "      checking matched region: " << blossom_root->match.region << std::endl << std::flush;
                    valid = check_pointers_for_self_and_all_descendents(
                        blossom_root->match.region,
                        region_range,
                        node_range,
                        &discovered_child_edges);
                }

                if (!valid) {
                    if (DEBUG) t_out << "      SHATTERING blossom_root\n" << std::flush;
                    res += solver.shatter_blossom_and_extract_matches(blossom_root);
                } else {
                    for (const auto& child_edge : discovered_child_edges) {
                        if (child_edges_counter < child_edges_nelems_per_solver) {
                            child_edges_buff_base[child_edges_counter++] = child_edge;
                        }
                    }
                    blossom_roots_checked.push_back(blossom_root);
                    if (blossom_root->match.region) {
                        blossom_roots_checked.push_back(blossom_root->match.region);
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
    auto* node_fields_base = get_node_fields_ptr(shot_container_id, partition_start);
    if (DEBUG) {
        t_out << "    sending DetectorNodeEphemeralFields" << std::endl
              << "      Node Fields Base: " << node_fields_base << std::endl
              << "      Nelems: " << nodes_nelems_total << std::endl;
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
              << "      Nelems: " << k * regions_nelems_per_solver << std::endl;
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
              << "      Nelems: " << child_edges_counter << std::endl;
    }
    shmem_ctx_putmem_signal_nbi(t.context_shm, 
                                child_edges_buff_base, 
                                child_edges_buff_base, 
                                child_edges_counter * sizeof(BlossomChild), 
                                t.signal_shm, 
                                1, SHMEM_SIGNAL_ADD, 
                                other_pid);

    // Copy Bitmap & Construct BlossomChild array (sparse edge list from bitmap)
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_to_unmatch 
                                        + sizeof(GraphFillRegion*) * fusion_summary_base->regions_to_unmatch_size);
    size_t total_bitmap_bytes = 0;
    for (size_t i = 0; i < k; ++i) {
        int p = partition_start + i;
        auto& solver = *solvers[get_solver_id(shot_container_id, p)];
        auto& bitmap = solver.flooder.region_arena.shmem_bitmap;
        size_t bitmap_len = bitmap.size();
        // Copy bitmap to FusiionSummary
        size_t bitmap_size_bytes = bitmap_len * sizeof(uint64_t);
        memcpy(bitmap_base, bitmap.data(), bitmap_size_bytes);
        bitmap_base += bitmap_len; // Pointer arithmetic on uint64_t*
        total_bitmap_bytes += bitmap_size_bytes;
    }
    fusion_summary_base->blossom_children_size = child_edges_counter;

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
              << "      Bitmap Size Bytes: " << total_bitmap_bytes << std::endl;
    }
    shmem_ctx_putmem_signal_nbi(t.context_shm, 
                                fusion_summary_base, 
                                fusion_summary_base, 
                                summary_payload_size, 
                                t.signal_shm, 
                                1, SHMEM_SIGNAL_ADD, 
                                other_pid);
    
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
    if (DEBUG) t_out << "  shattered sent blossoms" << std::endl << std::flush;
    // for (size_t i = 0; i < k; ++i) {
    //     int p = partition_start + i;
    //     auto& solver = *solvers[get_solver_id(shot_container_id, p)];
    //     auto& bitmap_src = solver.flooder.region_arena.shmem_bitmap;
    //     GraphFillRegion* p_regions = get_regions_ptr(shot_container_id, p);   
    //     for (size_t word_id = 0; word_id < bitmap_src.size(); ++word_id) {
    //         uint64_t word = bitmap_src[word_id];
    //         if (word != ~0ULL) {
    //             for (size_t bit = 0; bit < 64; ++bit) {
    //                 if (!((word >> bit) & 1ULL)) { // taken
    //                     size_t index = word_id * 64 + bit;
    //                     GraphFillRegion* r = p_regions + index;
    //                     r->cleanup_shell_area();
    //                     r->owner_arena->del(r);
    //                     // if (DEBUG) t_out << "    deleted r=" << r << std::endl << std::flush;
    //                 }
    //             }
    //         }
    //     }
    // }
    
    if (DEBUG) t_out << "    sent all data to " << other_pid << std::endl;
}

bool pm::DecodingUnit::get_solution_from_remote_pe(size_t shot_container_id, CrossRankTask &t, std::ofstream &t_out, std::vector<uint64_t> &hitsref) {
    size_t other_pid = t.other_pid;
    int k = config_parallel::k;
    int partition_start = (t.iamleft) ? t.part + 1 : t.part - k + 1;
    if (partition_start < 0) {
        k += partition_start; // reduce k
        partition_start = 0;
    }
    if (partition_start + k > graph.num_partitions) {
        k = graph.num_partitions - partition_start;
    }
    
    if (DEBUG) t_out << "    getting (p" << partition_start << ", k=" << k << ") data from " << other_pid << std::endl << std::flush;

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
    GraphFillRegion* k_block_regions_base = get_regions_ptr(shot_container_id, partition_start);
    GraphFillRegion* k_block_regions_end = k_block_regions_base + k * regions_nelems_per_solver;
    // Define GLOBAL memory range for stricter-but-permissive OOB check (allow pointers to ANY partition owned by this PE/process)
    if (DEBUG) t_out << "  k_block_regions_base: " << k_block_regions_base << std::endl
                     << "  k_block_regions_end: " << k_block_regions_end << " (" << k*regions_nelems_per_solver << ")" << std::endl << std::flush;
    
    auto& nodes_start_bounds = graph.partition_bounds[partition_start];
    auto& nodes_end_bounds = graph.partition_bounds[partition_start + k - 1];
    size_t nodes_nelems_total = nodes_end_bounds.second - nodes_start_bounds.first + 1;
    // DetectorNode* k_block_nodes_base = &graph.graph_ptr->nodes[node_start_bounds.first];
    // DetectorNode* k_block_nodes_end = k_block_nodes_base + nodes_nelems_total;
    ptrdiff_t nodes_start_bound_with_vb = (t.iamleft) ? graph.vb_bounds[partition_start-1].first : nodes_start_bounds.first;
    ptrdiff_t nodes_end_bound_with_vb = (!t.iamleft) ? graph.vb_bounds[partition_start + k - 1].second : nodes_end_bounds.second;
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
        // if (!ptr) return nullptr;
        // Rebase to local memory
        // DetectorNode* rebased = (DetectorNode*)((char*)ptr - (char*)remote_static_nodes_base + (char*)local_static_nodes_base);
        ptrdiff_t index = ptr - remote_static_nodes_base;
        // Check if remote pointer is within the valid range of the transmitted block
        if (index < nodes_start_bound_with_vb || index > nodes_end_bound_with_vb) {
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

    // // 2. Reconstruct regions_to_unmatch
    // for (size_t i = 0; i < fusion_summary_base->regions_to_unmatch_size; ++i) {
    //     GraphFillRegion* rebased = rebase_region_ptr(fusion_summary_base->regions_to_unmatch[i]);
    //     if (rebased) {
    //         if (DEBUG) t_out << "      " << rebased << std::endl << std::flush;
    //         t.regions_to_unmatch.push_back(rebased);
    //     }
    // }

    // 3. Copy Bitmaps & Reconstruct GraphFillRegions (Iterate all k partitions)
    if (DEBUG) {
        t_out << "    reconstructing solution state\n"
              << "      Remote Regions To Unmatch Size: " << fusion_summary_base->regions_to_unmatch_size << std::endl << std::flush;
    }
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_to_unmatch 
                                        + sizeof(GraphFillRegion*) * fusion_summary_base->regions_to_unmatch_size);
    // std::vector<uint64_t*> bitmap_payload_bases(k, nullptr);
    size_t bitmap_len = solvers[partition_start]->flooder.region_arena.shmem_bitmap.size();
    // auto is_region_live = [&](GraphFillRegion* region) -> bool {
    //     if (region == nullptr) {
    //         return false;
    //     }
    //     ptrdiff_t idx = region - k_block_regions_base;
    //     if (idx < 0 || (size_t)idx >= (size_t)k * regions_nelems_per_solver) {
    //         return false;
    //     }
    //     size_t part_offset = (size_t)idx / regions_nelems_per_solver;
    //     size_t local_idx = (size_t)idx % regions_nelems_per_solver;
    //     uint64_t* part_bitmap = bitmap_payload_bases[part_offset];
    //     return ((part_bitmap[local_idx / 64] >> (local_idx % 64)) & 1ULL) == 0ULL;
    // };
        // 6. Copy updated bitmap payloads into solver arenas after node-field rebasing.
    for (size_t i = 0; i < k; ++i) {
        int p = partition_start + i;
        auto& solver = *solvers[get_solver_id(shot_container_id, p)];
    }
    for (size_t i = 0; i < k; ++i) {
        int p = partition_start + i;
        auto& solver = *solvers[get_solver_id(shot_container_id, p)];
        auto& bitmap_dst = solver.flooder.region_arena.shmem_bitmap;
        memcpy(bitmap_dst.data(), bitmap_base + i*bitmap_len, bitmap_dst.size() * sizeof(uint64_t));

        GraphFillRegion* p_regions_base = get_regions_ptr(shot_container_id, p);

        for (size_t word_id = 0; word_id < bitmap_len; ++word_id) {
            uint64_t word = bitmap_dst[word_id];
            if (word != ~0ULL) {
                for (size_t bit = 0; bit < 64; ++bit) {
                    if (!((word >> bit) & 1ULL)) { // taken
                        size_t index = word_id * 64 + bit;
                        GraphFillRegion* r = p_regions_base + index;
                        // bool has_oob_ref = false;
                        // Rebase pointers with OOB check
                        if (r->blossom_parent) {
                            r->blossom_parent = rebase_region_ptr(r->blossom_parent);
                            // if (r->blossom_parent == nullptr) has_oob_ref = true;
                        }
                        if (r->blossom_parent_top) {
                            r->blossom_parent_top = rebase_region_ptr(r->blossom_parent_top);
                            // if (r->blossom_parent_top == nullptr) has_oob_ref = true;
                        }
                        if (r->match.region) {
                            r->match.region = rebase_region_ptr(r->match.region);
                            // if (r->match.region == nullptr) has_oob_ref = true;
                        }
                        if (r->match.edge.loc_from) {
                            r->match.edge.loc_from = rebase_node_ptr(r->match.edge.loc_from);
                            // if (r->match.edge.loc_from == nullptr) has_oob_ref = true;
                        }
                        if (r->match.edge.loc_to) {
                            r->match.edge.loc_to = rebase_node_ptr(r->match.edge.loc_to);
                            // if (r->match.edge.loc_to == nullptr) has_oob_ref = true;
                        }
                        // if (has_oob_ref) {
                        //     bitmap_payload[word_id] |= (1ULL << bit);  // free/delete region slot
                        //     r->blossom_parent = nullptr;
                        //     r->blossom_parent_top = nullptr;
                        //     r->match.region = nullptr;
                        //     r->match.edge.loc_from = nullptr;
                        //     r->match.edge.loc_to = nullptr;
                        // }
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
        // For now claim partition/vb
        shot_buffer->buffer[shot_container_id].i_solved_p[p] = true;
        if (p < graph.num_virtual_boundaries) shot_buffer->buffer[shot_container_id].i_solved_vb[p] = true;
    }
    if (partition_start + k - 1 < graph.num_virtual_boundaries) {
        shot_buffer->buffer[shot_container_id].i_solved_vb[partition_start + k - 1] = false;
    }
    if (DEBUG) {
        t_out << "    rebased regions" << std::endl << std::flush;
    }

    // 3. Reconstruct regions_to_unmatch using post-prune liveness.
    for (size_t i = 0; i < fusion_summary_base->regions_to_unmatch_size; ++i) {
        GraphFillRegion* rebased = rebase_region_ptr(fusion_summary_base->regions_to_unmatch[i]);
        if (rebased) {
            if (DEBUG) t_out << "      " << rebased << std::endl << std::flush;
            t.regions_to_unmatch.push_back(rebased);
        }
    }

    // 4. Reconstruct Blossom Children
    //   Blossom children for all K partitions are packed in the buffer starting at `partition_start`.
    BlossomChild* child_edges_buff_base = get_child_edges_ptr(shot_container_id, partition_start);

    for (size_t i=0; i < fusion_summary_base->blossom_children_size; ++i) {
        BlossomChild& child = child_edges_buff_base[i];
        GraphFillRegion* parent = k_block_regions_base + child.blossom_parent;
        // if (!is_region_live(parent)) {
        //     continue;
        // }

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
    DetectorNodeEphemeralFields* node_fields_base = get_node_fields_ptr(shot_container_id, partition_start);
    for (size_t i = 0; i < nodes_nelems_total; ++i) {
        DetectorNodeEphemeralFields& fields = node_fields_base[i];
        // bool clear_node = false;
        if (fields.region_that_arrived) {
            fields.region_that_arrived = rebase_region_ptr(fields.region_that_arrived);
            if (fields.region_that_arrived) {
                fields.region_that_arrived->shell_area.push_back(&graph.graph_ptr->nodes[nodes_start_bounds.first + i]);
            }
        }
        if (fields.region_that_arrived_top) {
            fields.region_that_arrived_top = rebase_region_ptr(fields.region_that_arrived_top);
            // if (fields.region_that_arrived_top == nullptr || !is_region_live(fields.region_that_arrived_top)) {
            //     clear_node = true;
            // }
        }
        if (fields.reached_from_source) {
            fields.reached_from_source = rebase_node_ptr(fields.reached_from_source);
            // if (fields.reached_from_source == nullptr) {
            //     clear_node = true;
            // }
        }
        // if (clear_node) {
        //     fields.region_that_arrived = nullptr;
        //     fields.region_that_arrived_top = nullptr;
        //     fields.reached_from_source = nullptr;
        // }
    }

    if (DEBUG) {
        t_out << "    rebased DetectorNodeEphemeralFields" << std::endl
              << "      Node Fields Base: " << node_fields_base << std::endl
              << "      Nelems: " << nodes_nelems_total << std::endl << std::flush;
    }
    
    return true;
}
#endif

inline void pm::DecodingUnit::extract_match_edges(pm::Mwpm& solver, pm::ShotContainer& shot, std::vector<uint64_t>& hitsref, size_t tid, std::ostream& t_out) {
#ifdef USE_SHMEM
    pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
        solver, hitsref);
#else
    pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
        solver, shot.sparse_shot.hits);
#endif
    if (!solver.flooder.negative_weight_detection_events.empty())
        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
            solver, solver.flooder.negative_weight_detection_events);
    solver.extract_paths_from_match_edges(
        solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
    // XOR negative weight observables
    for (auto& obs : solver.flooder.negative_weight_observables)
        *(shot.res.obs_crossed.data() + obs) ^= 1;
    // Add negative weight sum to blossom solution weight
    shot.res.weight += solver.flooder.negative_weight_sum;
}

inline void pm::DecodingUnit::extract_obs_mask(pm::Mwpm& solver, pm::ShotContainer& shot, std::vector<uint64_t>& hitsref, size_t tid, std::ostream& t_out) {
#ifdef USE_SHMEM
    pm::MatchingResult bit_packed_res;
    bit_packed_res +=
        pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
            solver, hitsref);

    bit_packed_res +=
        pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
            solver, hitsref);
#else
    pm::MatchingResult bit_packed_res =
        pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
            solver, shot.sparse_shot.hits);
#endif
    if (!solver.flooder.negative_weight_detection_events.empty())
        bit_packed_res +=
            pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                solver, solver.flooder.negative_weight_detection_events);
    // XOR in negative weight observable mask
    bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
    // Translate observable mask into bit vector
    pm::fill_bit_vector_from_obs_mask(
        bit_packed_res.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
    // Add negative weight sum to blossom solution weight
    shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
}

// inline void pm::DecodingUnit::isolate_solution_before_sending(pm::MatchingResult& bit_packed_res, pm::Mwpm& solver, pm::ShotContainer& shot, size_t shot_container_id, std::pair<int, int> p_range, std::pair<pm::GraphFillRegion*, pm::GraphFillRegion*> region_range, std::pair<pm::DetectorNode*, pm::DetectorNode*> node_range, std::ofstream& t_out) {
//     std::vector<std::vector<uint64_t>*> hits;
//     for (size_t p = p_range.first; p < p_range.second - 1; ++p) {
//         if (DEBUG) t_out << "  p" << p << std::flush;
//         hits.emplace_back(&shot.partition_hits[p]);
//         if (DEBUG) t_out << "  vb" << p << std::flush;
//         hits.emplace_back(&shot.virtual_boundary_hits[p]);
//     }
//     if (DEBUG) t_out << "  p" << p_range.second - 1 << std::flush;
//     hits.emplace_back(&shot.partition_hits[p_range.second - 1]);
//     std::vector<pm::GraphFillRegion*> blossom_roots_checked;
//     for (std::vector<uint64_t>* hitsref : hits) {
//         for (uint64_t i : *hitsref) {
//             auto& node_state = solver.flooder.graph.nodes[i].state(shot_container_id);
//             if (node_state.region_that_arrived) {
//                 pm::GraphFillRegion*& blossom_root = node_state.region_that_arrived_top;
//                 if (std::find(blossom_roots_checked.begin(), blossom_roots_checked.end(), blossom_root) == blossom_roots_checked.end()) {
//                     if (DEBUG) t_out << "    checking blossom_root: " << blossom_root << std::endl << std::flush;
//                     if (!check_pointers_for_self_and_all_descendents(blossom_root, region_range, node_range)) {
//                         if (DEBUG) t_out << "      SHATTERING blossom_root\n" << std::flush;
//                         bit_packed_res +=
//                             solver.shatter_blossom_and_extract_matches(blossom_root);
//                     } else  {
//                         blossom_roots_checked.push_back(blossom_root);
//                     }
//                 }
//             }

//         }
//     }

// //         bit_packed_res +=
// //             pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
// //                 solver, hitsref);
// //         if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl;
// //         shot.i_solved_p[p] = false;
// //         if (p < graph.num_virtual_boundaries) {
// //             if (DEBUG) t_out << "  vb" << p;
// //             hitsref = shot.virtual_boundary_hits[p];
// //             bit_packed_res +=
// //                 pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
// //                     solver, hitsref);
// //             if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl;
// //             shot.i_solved_vb[p] = false;
// //         }
// //     }
// //     pm::MatchingResult res;
// //     for (auto& i : detection_events) {
// // #ifdef USE_THREADS
// //         if (mwpm.flooder.graph.nodes[i].state(rotating_buffer_idx).region_that_arrived)
// //             res += mwpm.shatter_blossom_and_extract_matches(mwpm.flooder.graph.nodes[i].state(rotating_buffer_idx).region_that_arrived_top);
// // #else
// //         if (mwpm.flooder.graph.nodes[i].region_that_arrived)
// //             res += mwpm.shatter_blossom_and_extract_matches(mwpm.flooder.graph.nodes[i].region_that_arrived_top);
// // #endif
// //     }
// //     return res;
// }

// Core parallel decoding loop
void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
#ifdef USE_SHMEM
    if (DEBUG) {
        std::string ps = "PE" + std::to_string(pid) + " partitions: ";
        for (auto p : my_partitions)
            ps += std::to_string(p) + " ";
        std::cout << ps << std::endl << std::flush;
    }
    shmem_barrier_all(); // needed to avoid races on symmetric data (signals)
#endif
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        std::ofstream t_out;
        if (DEBUG) {
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
        try {
            while (true) {
                if (DEBUG) {
                    t_out << std::endl
                          << "Starting shot " << shot_id << ", buffer round " << shot_buffer_round << ", buffer_id "
                          << shot_container_id << std::endl;
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
                Task* t = &shot.tasks[tid];
                size_t next_leaf_inc = (num_threads < 2) ? 2 : num_threads;
                size_t next_leaf_id = tid+next_leaf_inc;
                int solver_id = get_solver_id(shot_container_id, t->part, tid);
                if (DEBUG)
                    t_out << "solvers[" << solver_id << "]\n";
                bool stolen = t->try_to_steal_leaf(shot_buffer_round);
                bool i_solved_root = false;
                while (stolen) {  // Got task
                    if (DEBUG) {
                        t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl << std::flush;
                    }
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    // Solve task
                    t->setup();
                    pm::Mwpm& solver = *solvers[solver_id];
                    solver.prepare_for_task(t, shot_id);
#ifdef USE_SHMEM
                    t_out << "  solver bounds: " << solver.flooder.vb_left << "(vb_left) " << solver.flooder.vb_right << " (vb_right)" << std::endl << std::flush;
#endif
                    pm::process_timeline_until_completion(
                        solver,
                        hitsref
#ifdef ENABLE_DRAW_FLAGS
                        ,
                        draw_frames
#endif
                        ,
                        true,
                        tid);
// #ifdef USE_SHMEM
                    if (DEBUG) t_out << "  solved " << (t->is_fusion ? "f" : "p") << t->part << std::endl << std::flush;
                    (t->is_fusion) ? shot.i_solved_vb[t->part] = true : shot.i_solved_p[t->part] = true;
// #endif
                    t->mark_solved();
                    if (t->parent) {
                        bool iamleft = t->child_bit == 1;
                        Task* sibling = (iamleft) ? t->parent->right_child : t->parent->left_child;
                        stolen = t->try_to_steal_parent();
                        t = t->parent;
                        // Try to steal sibling or descendent of sibling
                        if (!stolen && !sibling->is_fusion) {
                            solver_id = get_solver_id(shot_container_id, sibling->part, tid);
                            stolen = sibling->try_to_steal_leaf(shot_buffer_round);
                            t = sibling;
                        }
                        if (!stolen && next_leaf_id < my_partitions.size()) {
                            t = &shot.tasks[next_leaf_id];
                            solver_id = get_solver_id(shot_container_id, t->part, tid);
                            if (DEBUG) t_out << "solvers[" << solver_id << "]\n";
                            stolen = t->try_to_steal_leaf(shot_buffer_round);
                            next_leaf_id += next_leaf_inc;
                        }
                        if (DEBUG) {
                            if (t != nullptr) {
                                t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen
                                      << std::endl << std::flush;
                            }
                        }
                    } else {
                        stolen = false;
                        i_solved_root = true;
                    }
                }
                if (i_solved_root) {
                    auto& solver = *solvers[solver_id];
                    solver.flooder.match_edges.clear();
                    pm::MatchingResult bit_packed_res;
#ifdef USE_SHMEM
//                     // extract solution for partitions in between
//                     // bool even = !(pid % 2);
//                     //   Even PEs map slots (to left, to right), odd PEs (to right, to left)
//                     int vb_left, vb_right;
//                     if (shot.cross_rank_tasks.size() == 2) {
//                         vb_left = shot.cross_rank_tasks[0].vb_right;
//                         vb_right = shot.cross_rank_tasks[1].vb_left;
//                     } else if (shot.cross_rank_tasks[0].iamleft) {
//                         vb_left = -1;
//                         vb_right = shot.cross_rank_tasks[0].vb_left;
//                     } else {
//                         vb_left =  shot.cross_rank_tasks[0].vb_right;
//                         vb_right = graph.num_partitions - 1;
//                     }
// #ifdef ENABLE_DRAW_FLAGS
//                     solver.flooder.vb_left = my_partitions[0]-1;
//                     solver.flooder.vb_right = my_partitions[my_partitions.size()-1];
//                     if (draw_frames) {
//                         draw_frame(solver, pm::MwpmEvent::no_event(), 1000, true, tid);
//                     }
// #endif
//                     t_out << "  extracting from vb_left=" << vb_left << " to vb_right=" << vb_right << std::endl << std::flush;
//                     if (shot.num_observables > sizeof(pm::obs_int) * 8) {
//                         throw std::invalid_argument("More than 64 observables not yet implemented for OpenSHMEM.");
//                         if (vb_left >= 0) {
//                             if (DEBUG) t_out << "  vb" << vb_left;
//                             auto& hitsref = shot.virtual_boundary_hits[vb_left];
//                             pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
//                                 solver, hitsref);
//                             shot.i_solved_vb[vb_left] = false;
//                         }
//                         for (size_t p = vb_left + 1; p <= vb_right; ++p) {
//                             if (DEBUG) t_out << "  p" << p;
//                             auto& hitsref = shot.partition_hits[p];
//                             pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
//                                 solver, hitsref);
//                             shot.i_solved_p[p] = false;
//                             if (p < graph.num_virtual_boundaries) {
//                                 if (DEBUG) t_out << "  vb" << p;
//                                 hitsref = shot.virtual_boundary_hits[p];
//                                 pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
//                                     solver, hitsref);
//                                 shot.i_solved_vb[p] = false;
//                             }
//                         }
//                     } else {
//                         if (vb_left >= 0) {
//                             if (DEBUG) t_out << "  vb" << vb_left;
//                             auto& hitsref = shot.virtual_boundary_hits[vb_left];
//                             bit_packed_res +=
//                                 pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
//                                     solver, hitsref);
//                             if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl;
//                             shot.i_solved_vb[vb_left] = false;
//                         }
//                         for (size_t p = vb_left + 1; p <= vb_right; ++p) {
//                             if (DEBUG) t_out << "  p" << p << std::flush;
//                             auto& hitsref = shot.partition_hits[p];
//                             bit_packed_res +=
//                                 pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
//                                     solver, hitsref);
//                             if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl;
//                             shot.i_solved_p[p] = false;
//                             if (p < graph.num_virtual_boundaries) {
//                                 if (DEBUG) t_out << "  vb" << p;
//                                 hitsref = shot.virtual_boundary_hits[p];
//                                 bit_packed_res +=
//                                     pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
//                                         solver, hitsref);
//                                 if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl;
//                                 shot.i_solved_vb[p] = false;
//                             }
//                         }
//                     }
                    // if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl;
                    if (DEBUG) t_out << "Trying cross-rank fusions" << std::endl << std::flush;
                    for (auto& t : shot.cross_rank_tasks) {
                        if (DEBUG) t_out << "  t.part: " << t.part << std::endl
                                         << "  t.iamleft: " << t.iamleft << std::endl
                                         << "  t.other_pid: " << t.other_pid << std::endl
                                         << "  t.status_shm: " << t.status_shm << std::endl
                                         << "  t.signal_shm: " << t.signal_shm << std::endl
                                         << "  t.fusion_summary_shm: " << t.fusion_summary_shm << std::endl
                                         << "  task_status_ptr: " << task_status_ptr << std::endl
                                         << std::flush;
                        if (t.try_to_steal(pid, t_out)) {
                            if (DEBUG) t_out << "  Stole cross-rank fusion with " << t.other_pid << std::endl << std::flush;
                            t_out << "    entering t.setup();" << std::endl << std::flush;
                            t.setup();
                            t_out << "    exiting t.setup();" << std::endl << std::flush;
                            solver_id = get_solver_id(shot_container_id, t.part);
                            auto& solver = *solvers[solver_id];
                            // Configure solver's flooder bounds
                            //   Hard bounds based on vb_left/vb_right of CrossRankTask
                            solver.prepare_for_task(&t, shot_id);
                            auto& hitsref = shot.virtual_boundary_hits[t.part];
                            // Get remote window
                            get_solution_from_remote_pe(shot_container_id, t, t_out, hitsref);
#ifdef ENABLE_DRAW_FLAGS
                            if (draw_frames) {
                                draw_frame(solver, pm::MwpmEvent::no_event(), 1000, true, tid);
                            }
#endif
                            if (DEBUG) t_out << "  Solving cross-pe fusion " << t.part 
                                                << " bounds " << solver.flooder.vb_left << " " << solver.flooder.vb_right << std::endl << std::flush;
                            pm::process_timeline_until_completion(
                                solver,
                                hitsref
#ifdef ENABLE_DRAW_FLAGS
                                ,
                                draw_frames
#endif
                                ,
                                true,
                                tid);
                            // Mark completion in SHMEM
                            t.mark_solved(pid); 
                            shot.i_solved_vb[t.part] = true;
#ifdef ENABLE_DRAW_FLAGS
                            if (draw_frames) {
                                draw_frame(solver, pm::MwpmEvent::no_event(), 1001, true, tid);
                            }
#endif
                        } else {
                            // I send
                            if (DEBUG) t_out << "  Sending cross-rank fusion data to " << t.other_pid << std::endl;
                            t.setup();
#ifdef ENABLE_DRAW_FLAGS
                            if (draw_frames) {
                                auto& solver = *solvers[get_solver_id(shot_container_id, t.part)];
                                solver.prepare_for_task(&t, shot_id);
                                draw_frame(solver, pm::MwpmEvent::no_event(), 1000, true, tid);
                            }
#endif
                            send_solution_to_remote_pe(shot_container_id, bit_packed_res, t, t_out);
                        }
                    }
#endif
                    if (DEBUG) {
                        t_out << "T" << tid << " extracting solution" << std::endl << std::flush;
                    }
#ifdef ENABLE_DRAW_FLAGS
                    if (draw_frames) {
                        auto& solver = *solvers[solver_id];
                        draw_frame(solver, pm::MwpmEvent::no_event(), 2000, true, tid);
                    }
#endif
                    if (shot.num_observables > sizeof(pm::obs_int) * 8) {
// #ifdef USE_SHMEM
                        size_t i = 0;
                        for (auto& hitsref : shot.partition_hits) {
                            if (shot.i_solved_p[i]) {
                                if (DEBUG) t_out << "  p" << i;
                                pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                    solver, hitsref);
                                shot.i_solved_p[i] = false;
                            }
                            i++;
                        }
                        if (DEBUG) t_out << std::endl;
                        i = 0;
                        for (auto& hitsref : shot.virtual_boundary_hits) {
                            if (shot.i_solved_vb[i]) {
                                if (DEBUG) t_out << "  vb" << i;
                                pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                    solver, hitsref);
                                }
                            shot.i_solved_vb[i] = false;
                            i++;
                        }
                        if (DEBUG) t_out << std::endl;
// #else
//                         pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
//                             solver, shot.sparse_shot.hits);
// #endif
                        if (!solver.flooder.negative_weight_detection_events.empty())
                            shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                solver, solver.flooder.negative_weight_detection_events);
                        solver.extract_paths_from_match_edges(
                            solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
                        // XOR negative weight observables
                        for (auto& obs : solver.flooder.negative_weight_observables)
                            *(shot.res.obs_crossed.data() + obs) ^= 1;
                        // Add negative weight sum to blossom solution weight
                        shot.res.weight += solver.flooder.negative_weight_sum;
                    } else {
// #ifdef USE_SHMEM
                        // pm::MatchingResult& bit_packed_res = shot.obs_mask;
                        size_t i = 0;
                        if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl << std::flush;
                        for (auto& hitsref : shot.partition_hits) {
                            if (shot.i_solved_p[i]) {
                                if (DEBUG) t_out << "  p" << i << std::flush;
                                bit_packed_res +=
                                    pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                        solver, hitsref);
                                if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl << std::flush;
                            }
                            shot.i_solved_p[i] = false;
                            i++;
                        }
                        // if (DEBUG) t_out << std::endl << std::flush;
                        i = 0;
                        for (auto& hitsref : shot.virtual_boundary_hits) {
                            if (shot.i_solved_vb[i]) {
                                if (DEBUG) t_out << "  vb" << i << std::flush;
                                bit_packed_res +=
                                    pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                        solver, hitsref);
                                if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl << std::flush;
                            }
                            shot.i_solved_vb[i] = false;
                            i++;
                        }
                        // if (DEBUG) t_out << std::endl << std::flush;
// #else
//                         pm::MatchingResult bit_packed_res =
//                             pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
//                                 solver, shot.sparse_shot.hits);
// #endif
                        if (!solver.flooder.negative_weight_detection_events.empty()) {
                            if (DEBUG) t_out << "  there are negative detection events" << std::endl << std::flush;
                            bit_packed_res +=
                                shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                    solver, solver.flooder.negative_weight_detection_events);
                        }
                        // XOR in negative weight observable mask
                        bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
                        // Translate observable mask into bit vector
                        pm::fill_bit_vector_from_obs_mask(
                            bit_packed_res.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
                        // Add negative weight sum to blossom solution weight
                        shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
                    }
                    if (DEBUG) t_out << "   obs_mask: " << bit_packed_res.obs_mask << std::endl << std::flush;
#ifdef ENABLE_DRAW_FLAGS
                    if (draw_frames) {
                        auto& solver = *solvers[solver_id];
                        draw_frame(solver, pm::MwpmEvent::no_event(), 2001, true, tid);
                    }
#endif
#ifdef USE_SHMEM
                    if (DEBUG) {
                        for (int p=0; p < graph.num_partitions; ++p) {
                            auto& solver = *solvers[get_solver_id(shot_container_id, p)];
                            bool good = true;
                            for (auto word : solver.flooder.region_arena.shmem_bitmap) {
                                if (word != ~0ULL)
                                    good = false;
                            }
                            if (!good) t_out << "  ERROR: solver for p" << p << " not empty\n" << std::flush; 
                        }
                    }
                    // BARRIER needed to prevent race on extraction/putting mem
                    // Need to make more robust
                    if (DEBUG) t_out << "    entering barrier" << std::endl << std::flush;
                    shmem_barrier_all();
                    if (DEBUG) t_out << "    exiting barrier" << std::endl << std::flush;
#endif
                    shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
                }
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

// Individidual task solver
// void pm::DecodingUnit::solve_task(
//     pm::Mwpm& solver,
//     std::vector<uint64_t>& hits,
//     Task* task,
//     int tid
// #ifdef ENABLE_DRAW_FLAGS
//     ,
//     int draw_frames
// #endif
//     ,
//     int shot_id) {
//     task->setup();
//     solver.prepare_for_task(task, shot_id);
//     // if (DEBUG && tid==0) {
//     //     output_detector_nodes(solver, true);
//     // }
// #ifdef USE_SHMEM
//     if (task->is_cross_pe_fusion) {
//         get_data_from_remote_pe();
//     }
// #endif
//     pm::process_timeline_until_completion(
//         solver,
//         hits
// #ifdef ENABLE_DRAW_FLAGS
//         ,
//         draw_frames
// #endif
//         ,
//         true,
//         tid);
//     task->mark_solved(
//         pid
//     );
//     // if (DEBUG && tid==0) {
//     //     output_solution_state(solver, hits, true);
//     // }
// }
