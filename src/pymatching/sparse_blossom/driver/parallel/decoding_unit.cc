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
#if ENABLE_DRAW_FLAGS
    ,
    bool draw_frames
#endif
    )
    : ensure_search_flooder_included(ensure_search_flooder_included),
      enable_correlations(enable_correlations)
#if ENABLE_DRAW_FLAGS
      ,
      draw_frames(draw_frames)
#endif
{
#ifdef USE_SHMEM
    size_t n_pes = shmem_n_pes();
    pid = shmem_my_pe();
    // --- Create communication contexts ---
    // shmem_ctx_create(0, &summary_ctx);
    // shmem_ctx_create(0, &fallback_ctx);
    // if (!summary_ctx || !fallback_ctx) {
    //     throw std::invalid_argument("Failed to create SHMEM communication contexts.");
    // }
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
    task_status_ptr = static_cast<uint64_t*>(shmem_align(sizeof(uint64_t), (graph.num_partitions-1) * 2 * NUM_BUFFERS_PER_UNIT * sizeof(uint64_t))); // Space for every fusion to have task_status for stealing + put signal
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
    task_fusion_summary_ptr = static_cast<FusionSummary*>(shmem_malloc(task_fusion_summary_size_per_task * graph.num_partitions * NUM_BUFFERS_PER_UNIT));
    if (atomics_ptr == nullptr || task_status_ptr == nullptr || task_fusion_summary_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric atomics buffer.");
    }
    // fusion_summary_ptr = static_cast<FusionSummary*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * 2 * sizeof(FusionSummary)));
    // closest_nb_nelems_per_buffer = graph.node_part_id.size() / graph.num_rounds * SHMEM_INTERSECTION_BUFFER_FACTOR;
    // closest_nb_ptr = static_cast<std::pair<size_t, size_t>*>(
    //     shmem_malloc(closest_nb_nelems_per_buffer * NUM_BUFFERS_PER_UNIT * sizeof(std::pair<size_t, size_t>)));
    // obs_crossed_nelems_per_buffer = graph.graph_ptr->num_observables * sizeof(uint8_t);
    // obs_crossed_nelems_per_buffer = obs_crossed_nelems_per_buffer / sizeof(uint32_t) + (obs_crossed_nelems_per_buffer % sizeof(uint32_t) > 0); // minumum size uint32_t
    // obs_crossed_ptr =
    //     static_cast<uint32_t*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * obs_crossed_nelems_per_buffer * sizeof(uint32_t)));
    // for (size_t i=0; i < NUM_BUFFERS_PER_UNIT*obs_crossed_nelems_per_buffer; ++i) {
    //     obs_crossed_ptr[i] = 0;
    }
#endif
    // --- Create shot buffer ---
    shot_buffer = std::make_shared<pm::ShotBuffer>(
#ifdef USE_SHMEM
        atomics_ptr,
        // (uint8_t*)obs_crossed_ptr,
        // obs_crossed_nelems_per_buffer * sizeof(uint32_t),
#endif
        std::move(reader),
        std::move(writer),
        graph.num_partitions,
        graph.num_virtual_boundaries,
        graph.graph_ptr->num_observables);
    // --- Build tasks ---
    build_tasks_for_round_partitioning(
#ifdef USE_SHMEM
        shmem_n_pes()
#endif
    );
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
    if (max_threads < graph.num_partitions / n_pes) {
        throw std::invalid_argument("The number of threads should be >= the number of partitions / the number of PE's.");
    }
    num_threads = graph.num_partitions / n_pes;  // Current TEST Case for dividing across two PE's
    num_partition_units = 1; // FIX THIS??
    num_solvers_per_buffer = graph.num_partitions;
#else
    num_threads =
        (max_threads > graph.num_partitions) ? graph.num_partitions : max_threads;  // max num_partitions threads
    num_partition_units = graph.num_partitions / num_threads + (graph.num_partitions % num_threads > 0);
    num_solvers_per_buffer = num_threads*num_partition_units;
#endif
    std::cout << "graph.num_partitions = " << graph.num_partitions << "; num_threads: " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
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
#if ENABLE_DRAW_FLAGS
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
    shmem_free(atomics_ptr);
    shmem_free(regions_ptr);
#endif
}

// Builds balanced fusion tree assuming round-based partitioning
void pm::DecodingUnit::build_tasks_for_round_partitioning(
#ifdef USE_SHMEM
    size_t n_pes
#endif
) {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // Build a full fusion tree: less than 2*N tasks. Reserve to keep element addresses stable.
#ifdef USE_SHMEM
    if (n_pes > 2) {
        throw std::invalid_argument("Maximum 2 PEs supported");
    }
    size_t partitions_leftover = graph.num_partitions % n_pes;
    size_t partitions_per_pe = graph.num_partitions / n_pes;
    my_partitions.reserve(partitions_per_pe + (partitions_leftover > pid));
#endif
    for (int shot_container_id=0; shot_container_id < NUM_BUFFERS_PER_UNIT; shot_container_id++) {  // Lazy repeat per shot container -- CAN BE IMPROVED!
        auto& shot_container = shot_buffer->buffer[shot_container_id];
        auto& tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * graph.num_partitions - 1));
        int task_id;
#ifdef USE_SHMEM
        size_t pe = 0;
        size_t pe_partitions_start = 0;
        size_t pe_partitions_end = partitions_per_pe + (partitions_leftover > 0);
#endif
        for (task_id = 0; task_id < graph.num_partitions; ++task_id) {
            tasks.emplace_back(task_id
#ifdef USE_SHMEM
                , pe
#endif
            );
#ifdef USE_SHMEM
            if (pe == pid) {
#endif
            if (shot_container_id == 0) {
                my_partitions.emplace_back(task_id);
            }
#ifdef USE_SHMEM
            }
            if (task_id == pe_partitions_end-1) {
                pe++;
                pe_partitions_start = pe_partitions_end;
                pe_partitions_end += partitions_per_pe + (partitions_leftover - pe > 0);
            }
#endif
        }
        // Save odd trailing task
        int tail_idx = -1;
        if (graph.num_partitions % 2) {
            tail_idx = graph.num_partitions - 1;
        }
        // Add each level of fusions
        int last_step_starts = 0;
        int this_step_starts = task_id;
        int start = 0;
        int step = 2;
        while (start < graph.num_partitions - 1) {
            int counter = 0;
            int i;
            for (i = start; i < graph.num_partitions - 1; i += step) {
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
                tasks.emplace_back(i, left_child, right_child
#ifdef USE_SHMEM
                    , task_status_ptr + (shot_container_id * (graph.num_partitions-1) + i) * 2,
                    task_status_ptr + (shot_container_id * (graph.num_partitions-1) + i) * 2 + 1
#endif
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
    }
    // if (DEBUG) {
    //     std::cout << "DEBUG:" << std::endl;
    //     for (auto& buffer : shot_buffer->buffer) {
    //         std::cout << "Buffer" << std::endl;
    //         for (Task& t : buffer.tasks) {
    //             std::cout << "--part: " << t.part << std::endl
    //                     << "  vb_left: " << t.vb_left << std::endl
    //                     << "  vb_right: " << t.vb_right << std::endl
    //                     << "  is_fusion: " << t.is_fusion << std::endl
    //                     << "  child_bit: " << t.child_bit << std::endl
    //                     << "  left_child: " << t.left_child << std::endl
    //                     << "  right_child: " << t.right_child << std::endl
    //                     << "  parent: " << t.parent << std::endl;
    //         }
    //     }
    // }
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
// void pm::DecodingUnit::solve_cross_process_fusion_and_get_next_shot(
//     size_t shot_container_id, Task* task, size_t tid, size_t num_threads, size_t solver_id, size_t shot_id) {
//     size_t trivial_ctr = 0;
//     // --- setup task and solver ---
//     auto& shot = shot_buffer->buffer[shot_container_id];
//     auto& hitsref = shot.virtual_boundary_hits[task->part];
//     auto& solver = *solvers[solver_id];
//     task->setup();
//     solver.prepare_for_task(task, shot_container_id*NUM_BUFFERS_PER_UNIT);
//     // --- compute lower pid, index ---
//     size_t lower_pid, my_index, other_index;
//     if (pid < other_pid) {
//         lower_pid = pid;
//         my_index = 0;
//         other_index = 1;
//     } else {
//         lower_pid = other_pid;
//         my_index = 1;
//         other_index = 0;
//     }
//     uint64_t* shot_status_base = (uint64_t*) shot.current_buffer_round_shm;
//     if (DEBUG) {
//         std::cout << "Cross PE Fusion: PE" << pid << " shot " << shot_id << std::endl << std::flush;
//     }
//     // --- try to win shot ---
//     uint64_t res = shmem_uint64_atomic_fetch_inc(shot_status_base, lower_pid);
//     if (res == 0) { // I am first, send mem
//         // get nodes (CHANGE TO PUTS!!!)
//         auto& bounds = graph.partition_bounds[other_pid]; // general 2 PE case
//         auto* nodes_base = node_ephemeral_fields_ptr + bounds.first;
//         size_t nodes_nelems = bounds.second - bounds.first + 1;
//         shmem_putmem_signal_nbi(nodes_base, nodes_base, nodes_nelems *
//         sizeof(DetectorNodeEphemeralFields), shot_status_base+1, 1, SHMEM_SIGNAL_ADD, other_pid);
//         std::cout << "PE " << pid << " sent nodes for pe" << other_pid
//             << std::endl;
//         // get regions (Need to add other PE region buffer, CHANGE TO PUTS!!!)
//         auto* regions_base = regions_ptr + shot_container_id*num_threads*regions_nelems_per_solver;
//         size_t regions_nelems = regions_nelems_per_solver * num_threads;
//         shmem_putmem_signal_nbi(regions_base, regions_base, regions_nelems
//             * sizeof(GraphFillRegion), shot_status_base+1, 1, SHMEM_SIGNAL_ADD, other_pid);
//         std::cout << "PE " << pid << " sent regions for pe" << other_pid << std::endl;
//         // Wait for other PE to solve
//         if (pid == lower_pid) {
//             shmem_wait_until(shot_status_base, SHMEM_CMP_EQ, 0);
//         } else {
//             shmem_wait_until(shot_status_base, SHMEM_CMP_EQ, 1);
//             *shot_status_base = 0;
//         }
//         // --- reset atomic synchronizers ---
//         if (DEBUG) {
//             std::cout << "PE" << pid << " freed" << std::endl;
//         }
//     } else {
//         shmem_wait_until(shot_status_base+1, SHMEM_CMP_EQ, 2);
//         if (DEBUG) {
//             std::cout << "PE" << pid << " got nodes and regions" << std::endl;
//         }
//         // --- reset atomic synchronizers ---
//         *(shot_status_base+1) = 0;
//         shmem_uint64_atomic_set(shot_status_base, 0, lower_pid);
//         if (pid == lower_pid) {
//             shmem_uint64_atomic_set(shot_status_base, 1, other_pid);
//         }
//     }
//     // // --- extract solution ---
//     // if (shot.num_observables > sizeof(pm::obs_int) * 8) {
//     //     solver.flooder.match_edges.clear();
//     //     pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(solver, shot.sparse_shot.hits);
//     //     if (!solver.flooder.negative_weight_detection_events.empty())
//     //         shatter_blossoms_for_all_detection_events_and_extract_match_edges(
//     //             solver, solver.flooder.negative_weight_detection_events);
//     //     solver.extract_paths_from_match_edges(
//     //         solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
//     //     // XOR negative weight observables
//     //     for (auto& obs : solver.flooder.negative_weight_observables)
//     //         *(shot.res.obs_crossed.data() + obs) ^= 1;
//     //     // Add negative weight sum to blossom solution weight
//     //     shot.res.weight += solver.flooder.negative_weight_sum;
//     // } else {
//     //     pm::MatchingResult bit_packed_res =
//     //         pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
//     //             solver, shot.sparse_shot.hits);
//     //     if (!solver.flooder.negative_weight_detection_events.empty())
//     //         bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
//     //             solver, solver.flooder.negative_weight_detection_events);
//     //     // XOR in negative weight observable mask
//     //     bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
//     //     // Translate observable mask into bit vector
//     //     pm::fill_bit_vector_from_obs_mask(
//     //         bit_packed_res.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
//     //     // Add negative weight sum to blossom solution weight
//     //     shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
//     // }
// //     uint32_t* obs_crossed_base = obs_crossed_ptr + shot_container_id * obs_crossed_nelems_per_buffer;
// //     if (DEBUG) {
// //         std::cout << "  PE" << pid << " res: " << *obs_crossed_base << std::endl << std::flush;
// //     }
// //     if (pid != lower_pid) {
// //         shmem_ctx_uint32_atomic_xor(summary_ctx, obs_crossed_base, *obs_crossed_base, lower_pid);
// //     }
// // #if ENABLE_DRAW_FRAMES
// //     else if (draw_frames) {
// //         std::filesystem::create_directory("out_parallel/frames/" + std::to_string(shot_id+1));
// //     }
// // #endif
//     // if (DEBUG && pid == lower_pid) {
//     //     std::cout << "  PE" << pid << " res is " << *obs_crossed_base << std::endl << std::flush;
//     // }
//     shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id, pid == lower_pid);
//     // --- reset atomic synchronizers ---
//     // // --- send closest neighbor info ---
//     // std::pair<size_t, size_t>* closest_nb_base = closest_nb_ptr + closest_nb_nelems_per_buffer * shot_container_id;
//     // auto& bounds = graph.vb_bounds[task->part];
//     // for (size_t node = bounds.first; node <= bounds.second; node++) {
//     //     new (closest_nb_base + node) std::pair<size_t, size_t>();
//     // }
//     // shmem_ctx_putmem_nbi(summary_ctx, closest_nb_base, closest_nb_base, closest_nb_nelems_per_buffer, other_pid);
//     // if (shot_status == 1) {  // I perform fusion
//     //     std::cout << "PE " << pid << " to perform fusion " << shot_container_id << "\n" << std::flush;
//     //     // - check fusion_summary
//     //     std::cout << "PE " << pid << " reading FusionSummary("
//     //               << fusion_summary_ptr[other_index].regions_matched_to_vb_size << ", "
//     //               << fusion_summary_ptr[other_index].regions_begin << ")\n"
//     //               << std::flush;
//     //     // - get other PE's solution state -
//     //     // // get nodes (CHANGE TO PUTS!!!)
//     //     // auto& bounds = graph.partition_bounds[other_pid];
//     //     // auto* nodes_base = node_ephemeral_fields_ptr + bounds.first;
//     //     // size_t nodes_nelems = bounds.second - bounds.first + 1;
//     //     // shmem_ctx_getmem_nbi(fallback_ctx, nodes_base, nodes_base, nodes_nelems *
//     //     // sizeof(DetectorNodeEphemeralFields), other_pid); std::cout << "PE " << pid << " got nodes for p" << other_pid
//     //     // << std::endl; get regions (Need to add other PE region buffer, CHANGE TO PUTS!!!) auto* regions_base =
//     //     // regions_ptr + (shot_container_id*num_threads)*regions_nelems_per_solver; size_t regions_nelems =
//     //     // regions_nelems_per_solver * num_threads; shmem_ctx_getmem_nbi(fallback_ctx, regions_base, regions_base,
//     //     // regions_nelems * sizeof(GraphFillRegion), other_pid); std::cout << "PE " << pid << " got regions for p" <<
//     //     // other_pid << std::endl;
//     //     // - fuse -
//     //     if (my_summary.regions_matched_to_vb_size == 0 &&
//     //         other_summary.regions_matched_to_vb_size == 0) {  // trivial case -- no cross-PE changes needed
//     //         // extract solution & write results
//     //     } else {
//     //         shmem_uint64_atomic_set(shot_buffer->shot_container_status + shot_container_id, 3, other_pid);
//     //     }
//     //     // --- write results & get next shot ---
//     //     shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
//     // } else {  // I wait on fusion
//     //     shmem_wait_until(shot_status_base, SHMEM_CMP_GT, 1);
//     //     // Other PE has posted
//     //     std::cout << "PE " << pid << " reading FusionSummary("
//     //               << fusion_summary_ptr[other_index].regions_matched_to_vb_size << ", "
//     //               << fusion_summary_ptr[other_index].regions_begin << ")\n"
//     //               << std::flush;
//     //     if (my_summary.regions_matched_to_vb_size == 0 &&
//     //         other_summary.regions_matched_to_vb_size == 0) {  // trivial case -- no cross-PE changes needed
//     //         // extract solution
//     //     } else {
//     //         shmem_wait_until(shot_status_base, SHMEM_CMP_EQ, 3);
//     //         std::cout << "PE " << pid << " freed " << shot_container_id << "\n" << std::flush;
//     //     }
//     // }
// }

void pm::DecodingUnit::send_solution_to_remote_pe(size_t shot_container_id, Task &t, std::ofstream &t_out) {
    bool iamleft = pid == t.left_pid;
    size_t other_pid = (iamleft) ? t.right_pid : t.left_pid;
    int partition_to_send = (iamleft) ? t.part : t.part+1;
    if (DEBUG) t_out << "    sending partition " << partition_to_send << " data to " << other_pid << std::endl;
    // send FusionSummary
    if (DEBUG) t_out << "    t->regions_to_unmatch.size(): " << t.regions_to_unmatch.size() << std::endl;
    FusionSummary* fusion_summary_base = get_fusion_summary_ptr(shot_container_id, partition_to_send);
    fusion_summary_base->regions_matched_to_vb_size = t.regions_to_unmatch.size();
    fusion_summary_base->regions_ptr_base = regions_ptr;
    // fusion_summary_base->nodes_ptr_base = node_ephemeral_fields_ptr;
    fusion_summary_base->static_nodes_base = graph.graph_ptr->nodes.data();
    for (size_t i=0; i<t.regions_to_unmatch.size(); ++i) {
        fusion_summary_base->regions_matched_to_vb[i] = t.regions_to_unmatch[i];
    }
    //   copy bitmap
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_matched_to_vb + sizeof(GraphFillRegion*) * fusion_summary_base->regions_matched_to_vb_size);
    auto& solver = *solvers[get_solver_id(shot_container_id, partition_to_send)];
    auto& bitmap_src = solver.flooder.region_arena.shmem_bitmap;
    memcpy(bitmap_base, bitmap_src.data(), bitmap_src.size() * sizeof(uint64_t));
    //   construct BlossomChild array
    BlossomChild* child_edges_buff_base = get_child_edges_ptr(shot_container_id, partition_to_send);
    GraphFillRegion* regions_base = get_regions_ptr(shot_container_id, partition_to_send);
    size_t child_edges_counter = 0;
    for (size_t i = 0; i < bitmap_src.size(); ++i) {
        uint64_t word = bitmap_src[i];
         if (word != ~0ULL) {
            for (size_t bit = 0; bit < 64; ++bit) { 
                if (!((word >> bit) & 1ULL)) { // If bit is 0 (taken)
                    size_t index = i * 64 + bit;
                    GraphFillRegion* r = regions_base + index;
                    if (!r->blossom_children.empty()) {
                        for (auto& child_edge : r->blossom_children) {
                            // Store copy of edge but with modified pointers (as offsets)
                            child_edges_buff_base[child_edges_counter] = BlossomChild{index, child_edge};
                            child_edges_counter++;
                        }
                    }
                }
            }
        }
    }
    fusion_summary_base->blossom_children_size = child_edges_counter;
    if (DEBUG) {
        t_out << "    sending FusionSummary\n"
            << "      regions_ptr:" << fusion_summary_base->regions_ptr_base << std::endl
            << "      static_nodes_base: " << fusion_summary_base->static_nodes_base << std::endl
            << "      b_child_size: " << fusion_summary_base->blossom_children_size << std::endl
            << "      r_matched_size: " << fusion_summary_base->regions_matched_to_vb_size << std::endl 
            << "      r_matched: " << fusion_summary_base->regions_matched_to_vb << std::endl
            << "      bitmap: (" << bitmap_base << ")  ";
        for (size_t i=0; i < bitmap_src.size(); ++i) {
            t_out << *(bitmap_base + i) << ", ";
        }
        t_out << std::endl;
    }
    shmem_putmem_signal_nbi(fusion_summary_base, fusion_summary_base, task_fusion_summary_size_per_task, t.signal_shm, 1, SHMEM_SIGNAL_ADD, other_pid);
    if (DEBUG) t_out << "    sent fusion summary to " << other_pid << std::endl;
    // send DetectorNodeEphemeralFields
    auto& bounds = graph.partition_bounds[partition_to_send];
    auto* nodes_base = get_node_fields_ptr(shot_container_id, partition_to_send);
    size_t nodes_nelems = bounds.second - bounds.first + 1;
    if (DEBUG) {
        t_out << "    boundary_partition: " << partition_to_send << "  (" << bounds.first << ", " << bounds.second << ")"
              << "  other " << (iamleft ? graph.partition_bounds[t.part+1].first : graph.partition_bounds[t.part].second) << std::endl
              << "    status_shm: " << t.status_shm << "  signal_shm: " << t.signal_shm << std::endl
              << "    nodes_base: " << nodes_base << " (" << node_ephemeral_fields_ptr << ")  nodes_nelems: " << nodes_nelems << std::endl
              << "    nodes_nelems_per_buffer: " << nodes_nelems_per_buffer << std::endl;
    }
    shmem_putmem_signal_nbi(nodes_base, nodes_base, nodes_nelems *
    sizeof(DetectorNodeEphemeralFields), t.signal_shm, 1, SHMEM_SIGNAL_ADD, other_pid);
    if (DEBUG) t_out << "    sent nodes to pe" << other_pid << std::endl;
    // send GraphFillRegions
    shmem_putmem_signal_nbi(regions_base, regions_base, regions_nelems_per_solver
        * sizeof(GraphFillRegion), t.signal_shm, 1, SHMEM_SIGNAL_ADD, other_pid);
    if (DEBUG) t_out << "    sent regions to pe" << other_pid << std::endl;
    shmem_putmem_signal_nbi(child_edges_buff_base, child_edges_buff_base, child_edges_counter * sizeof(BlossomChild), t.signal_shm, 1, SHMEM_SIGNAL_ADD, other_pid);
    if (DEBUG) t_out << "    sent " << child_edges_counter << " blossom child edges to pe" << other_pid << std::endl;
    // ensure outgoing nodes copied before freeing
    shmem_quiet();
    // free p
    shot_buffer->buffer[shot_container_id].i_solved_p[partition_to_send] = false;
    for (size_t i = 0; i < bitmap_src.size(); ++i) {
        uint64_t word = bitmap_src[i];
         if (word != ~0ULL) {
            for (size_t bit = 0; bit < 64; ++bit) { 
                if (!((word >> bit) & 1ULL)) { // If bit is 0 (taken)
                    size_t index = i * 64 + bit;
                    GraphFillRegion* r = regions_base + index;
                    if (DEBUG) t_out << "  deleting " << r << "  index " << r - regions_base << ", " << r - regions_ptr << std::endl << std::flush;
                    r->cleanup_shell_area();
                    t_out << "  cleaned shell area" << std::endl 
                          << "  r->owner_arena: " << r->owner_arena  << "  " << &solver.flooder.region_arena << std::endl << std::flush;
                    if (r->owner_arena) {
                        r->owner_arena->del(r);
                    } else {
                        throw std::invalid_argument("send_solution_to_remote_pe (rank " + std::to_string(pid) + "): r->owner_arena nullptr");
                    }
                    if (DEBUG) t_out << "  deleted r" << std::endl << std::flush;
                }
            }
        }
    }
    if (DEBUG) t_out << "cleaned up p" << partition_to_send << std::endl << std::flush;
}

bool pm::DecodingUnit::get_solution_from_remote_pe(size_t shot_container_id, Task &t, std::ofstream &t_out, std::vector<uint64_t> &hitsref) {
    bool iamleft = pid == t.left_pid;
    size_t other_pid = (iamleft) ? t.right_pid : t.left_pid;
    int partition_to_get = (iamleft) ? t.part+1 : t.part;
    if (DEBUG) t_out << "    getting partition " << partition_to_get << " data from " << other_pid << std::endl;
    // prep for accesses
    FusionSummary* fusion_summary_base = get_fusion_summary_ptr(shot_container_id, partition_to_get);
    auto& bounds = graph.partition_bounds[partition_to_get];
    auto* nodes_base = get_node_fields_ptr(shot_container_id, partition_to_get);
    size_t nodes_nelems = bounds.second - bounds.first + 1;
    if (DEBUG) {
        t_out << "    boundary_partition: " << partition_to_get << "  (" << bounds.first << ", " << bounds.second << ")"
            << "  other " << (iamleft ? graph.partition_bounds[t.part+1].first : graph.partition_bounds[t.part].second) << std::endl
            << "    status_shm: " << t.status_shm << "  signal_shm: " << t.signal_shm << std::endl
            << "    nodes_base: " << nodes_base << " (" << node_ephemeral_fields_ptr << ")  nodes_nelems: " << nodes_nelems << std::endl
            << "    nodes_nelems_per_buffer: " << nodes_nelems_per_buffer << std::endl;
    }
    auto* regions_base = get_regions_ptr(shot_container_id, partition_to_get);
    auto* child_edges_base = get_child_edges_ptr(shot_container_id, partition_to_get);
    if (DEBUG) {
        t_out << "    regions_base: " << regions_base << std::endl
	      << "    child_edges_base: " << child_edges_base << std::endl;
    }
    // wait for data
    shmem_wait_until(t.signal_shm, SHMEM_CMP_EQ, 4);
    uint64_t* bitmap_base = (uint64_t*)((char*)fusion_summary_base->regions_matched_to_vb + sizeof(GraphFillRegion*)*fusion_summary_base->regions_matched_to_vb_size);
    size_t bitmap_size = regions_nelems_per_solver / 64;
    if (DEBUG) {
        t_out << "    got fusion summary, nodes, & regions from pe" << other_pid <<std::endl
            << "    reading FusionSummary\n"
            << "      regions_ptr: " << fusion_summary_base->regions_ptr_base << std::endl
            << "      static_nodes_base: " << fusion_summary_base->static_nodes_base << std::endl
            << "      b_child_s: " << fusion_summary_base->blossom_children_size << std::endl
            << "      r_matched_size: " << fusion_summary_base->regions_matched_to_vb_size << std::endl
            << "      r_matched: " << fusion_summary_base->regions_matched_to_vb << std::endl
            << "      bitmap_base: " << bitmap_base << std::endl
            << "      bitmap: (" << bitmap_base << ")  ";
    }
    // Copy bitmap
    auto& solver = *solvers[get_solver_id(shot_container_id, partition_to_get)];
    auto& bitmap_src = solver.flooder.region_arena.shmem_bitmap;
    memcpy(bitmap_src.data(), bitmap_base, bitmap_src.size() * sizeof(uint64_t));
    if (DEBUG) {
        for (size_t i=0; i < bitmap_size; ++i) {
            t_out << bitmap_src[i] << ", ";
        }
        t_out << std::endl;
    }
    // Copy regions_matched
    auto& remote_regions_ptr_base = fusion_summary_base->regions_ptr_base;
    // t.regions_to_unmatch.resize(fusion_summary_base->regions_matched_to_vb_size);
    for (size_t i = 0; i < fusion_summary_base->regions_matched_to_vb_size; ++i) {
        t.regions_to_unmatch.push_back((GraphFillRegion*)((char*)fusion_summary_base->regions_matched_to_vb[i] - (char*)remote_regions_ptr_base + (char*)regions_ptr));
    }
    if (DEBUG) t_out << "    t->regions_to_unmatch.size(): " << t.regions_to_unmatch.size() << std::endl;
    // check for trivial case
    auto& remote_static_nodes_base = fusion_summary_base->static_nodes_base;
    auto* local_static_nodes_base = graph.graph_ptr->nodes.data();
    // if (hitsref.size() == 0 && t.regions_matched_to_virtual_boundary.size() == 0 && fusion_summary_base->regions_matched_to_vb_size == 0) {
    //     t_out << "    TRIVIAL CASE!" << std::endl;
    //     return false;
    // } else { // for now, always reconstruct data
        // --- Reconstruct Data ---
        for (size_t word_id = 0; word_id < bitmap_size; ++word_id) {
            uint64_t word = bitmap_src[word_id];
            // Iterate over 0-bits (taken slots)
            if (word != ~0ULL) { // If any bit is 0, it means some slots are taken
                for (size_t bit = 0; bit < 64; ++bit) {
                    if (!((word >> bit) & 1ULL)) { // If bit is 0 (taken)
                        size_t index = word_id * 64 + bit;
                        GraphFillRegion* r = regions_base + index;
                        new (&r->blossom_children) std::vector<RegionEdge>();
                        new (&r->shell_area) std::vector<DetectorNode*>();
                        r->shrink_event_tracker.clear();
                        r->alt_tree_node = nullptr;
                        if (r->blossom_parent)
                            r->blossom_parent = (GraphFillRegion*)((char*)r->blossom_parent - (char*)remote_regions_ptr_base + (char*)regions_ptr);
                        if (r->blossom_parent_top)
                            r->blossom_parent_top = (GraphFillRegion*)((char*)r->blossom_parent_top - (char*)remote_regions_ptr_base + (char*)regions_ptr);
                        if (r->match.region)
                            r->match.region = (GraphFillRegion*)((char*)r->match.region - (char*)remote_regions_ptr_base + (char*)regions_ptr);
                        if (r->match.edge.loc_from)
                            r->match.edge.loc_from = (DetectorNode*)((char*)r->match.edge.loc_from - (char*)remote_static_nodes_base + (char*)local_static_nodes_base);
                        if (r->match.edge.loc_to)
                            r->match.edge.loc_to = (DetectorNode*)((char*)r->match.edge.loc_to - (char*)remote_static_nodes_base + (char*)local_static_nodes_base);
                        r->owner_arena = &solver.flooder.region_arena;
                    }
                }
            }
        }
        if (DEBUG) t_out << "    reset region vectors & rebased points" << std::endl;
        // 2. Reconstruct Blossom Children
        if (DEBUG) t_out << "    reconstructing blossom children (" << fusion_summary_base->blossom_children_size << " edges)" << std::endl;
        for (size_t i=0; i < fusion_summary_base->blossom_children_size; ++i) {
            BlossomChild& child = child_edges_base[i];
            GraphFillRegion* parent = regions_base + child.blossom_parent;

            RegionEdge local_edge = child.region_edge;
            local_edge.region = (GraphFillRegion*)((char*)local_edge.region - (char*)remote_regions_ptr_base + (char*)regions_ptr);
            if (local_edge.edge.loc_from)
                local_edge.edge.loc_from = (DetectorNode*)((char*)local_edge.edge.loc_from - (char*)remote_static_nodes_base + (char*)local_static_nodes_base);
            if (local_edge.edge.loc_to)
                local_edge.edge.loc_to = (DetectorNode*)((char*)local_edge.edge.loc_to - (char*)remote_static_nodes_base + (char*)local_static_nodes_base);
            parent->blossom_children.push_back(local_edge);
        }
        if (DEBUG) t_out << "    reconstructed blossom children" << std::endl;
        // 3. Rebase DetectorNodeEphemeralFields & Reconstruct shell_area
        for (size_t i = 0; i < nodes_nelems; ++i) {
            DetectorNodeEphemeralFields& fields = nodes_base[i];
            // Rebase region_that_arrived
            if (fields.region_that_arrived) {
                fields.region_that_arrived = (GraphFillRegion*)((char*)fields.region_that_arrived - (char*)remote_regions_ptr_base + (char*)regions_ptr);
                // Reconstruct shell_area
                fields.region_that_arrived->shell_area.push_back(&graph.graph_ptr->nodes[bounds.first + i]);
            }
            // Rebase region_that_arrived_top
            if (fields.region_that_arrived_top) {
                fields.region_that_arrived_top = (GraphFillRegion*)((char*)fields.region_that_arrived_top - (char*)remote_regions_ptr_base + (char*)regions_ptr);
            }
            // Rebase reached_from_source
            if (fields.reached_from_source) {
                fields.reached_from_source = (DetectorNode*)((char*)fields.reached_from_source - (char*)remote_static_nodes_base + (char*)local_static_nodes_base);
            }
        }
        if (DEBUG) t_out << "    rebased detector nodes" << std::endl;
        // Claim other partition
        shot_buffer->buffer[shot_container_id].i_solved_p[partition_to_get] = true;
    // }
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
        uint64_t shot_buffer_round =
            shot_buffer->buffer[0].current_buffer_round.load();  // how many times buffer has looped
        uint64_t shot_id = shot_buffer_round * NUM_BUFFERS_PER_UNIT;
        try {
            while (true) {
                if (DEBUG) {
                    t_out << std::endl
                          << "Starting shot " << shot_id << ", buffer round " << shot_buffer_round << ", buffer_id "
                          << shot_container_id << std::endl;
                }
#if ENABLE_DRAW_FLAGS
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
                Task* t = &shot.tasks[my_partitions[tid]];
                size_t next_leaf_id = tid+num_threads;
                int solver_id = get_solver_id(shot_container_id, my_partitions[tid], tid);
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
                    bool need_to_solve = true;
                    if (t->is_cross_pe_fusion) {
                        need_to_solve = get_solution_from_remote_pe(shot_container_id, *t, t_out, hitsref);
                    }
                    if (need_to_solve) {
#endif
                    pm::process_timeline_until_completion(
                        solver,
                        hitsref
#if ENABLE_DRAW_FLAGS
                        ,
                        draw_frames
#endif
                        ,
                        true,
                        tid);
#ifdef USE_SHMEM
                    } // track that I solved this task
                    if (DEBUG) t_out << "  solved " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                    (t->is_fusion) ? shot.i_solved_vb[t->part] = true : shot.i_solved_p[t->part] = true;
#endif
                    t->mark_solved(
#ifdef USE_SHMEM
                        pid
#endif
                    );
//                     solve_task(
//                         ,
//                         hitsref,
//                         t,
//                         tid
// #if ENABLE_DRAW_FLAGS
//                         ,
//                         draw_frames
// #endif
//                         ,
//                         shot_id);
                    if (t->parent) {
                        Task* sibling = (t->child_bit == 1) ? t->parent->right_child : t->parent->left_child;
#ifdef USE_SHMEM
                        if (DEBUG) {
                            t_out << "  trying to steal parent t->parent->is_cross_pe_fusion: " << t->parent->is_cross_pe_fusion << std::endl;
                        }
#endif
                        stolen = t->try_to_steal_parent(
#ifdef USE_SHMEM
                            pid
#endif
                        );
                        t = t->parent;
#ifdef USE_SHMEM
                        if (!stolen && t->is_cross_pe_fusion) {
                            t->setup();
                            if (DEBUG) {
                                t_out << "  sending solution to PE" << ((t->left_pid == pid) ? t->right_pid : t->left_pid) << std::endl;
                            }
                            send_solution_to_remote_pe(shot_container_id, *t, t_out);
                        } else // skip trying to steal sibling
#endif
                        // Try to steal sibling or descendent of sibling
                        if (!stolen && !sibling->is_fusion) {
                            stolen = sibling->try_to_steal_leaf(shot_buffer_round);
                            t = sibling;
                        }
                        if (!stolen && next_leaf_id < my_partitions.size()) {
                            solver_id = get_solver_id(shot_container_id, my_partitions[next_leaf_id], tid);
                            if (DEBUG)
                                t_out << "solvers[" << solver_id << "]\n";
                            t = &shot.tasks[my_partitions[next_leaf_id]];
                            stolen = t->try_to_steal_leaf(shot_buffer_round);
                            next_leaf_id += num_threads;
                        }
                        if (DEBUG) {
                            if (t != nullptr) {
                                t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen
                                      << std::endl;
                            }
                        }
                    } else {
                        stolen = false;
                        i_solved_root = true;
                    }
                }
                if (i_solved_root) {
                    if (DEBUG) {
                        t_out << "T" << tid << " extracting solution" << std::endl;
                    }
                    auto& solver = *solvers[solver_id];
                    if (shot.num_observables > sizeof(pm::obs_int) * 8) {
                        solver.flooder.match_edges.clear();
        #ifdef USE_SHMEM
                        size_t i = 0;
                        for (auto& hitsref : shot.partition_hits) {
                            if (shot.i_solved_p[i]) {
                                if (DEBUG) t_out << "  p" << i;
                                pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                    solver, hitsref);
                            }
                            shot.i_solved_p[i] = false;
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
        #else
                        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                            solver, shot.sparse_shot.hits);
        #endif
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
        #ifdef USE_SHMEM
                        pm::MatchingResult bit_packed_res;
                        size_t i = 0;
                        for (auto& hitsref : shot.partition_hits) {
                            if (shot.i_solved_p[i]) {
                                if (DEBUG) t_out << "  p" << i << std::flush;
                                bit_packed_res +=
                                    pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                        solver, hitsref);
                            }
                            shot.i_solved_p[i] = false;
                            i++;
                        }
                        if (DEBUG) t_out << std::endl << std::flush;
                        i = 0;
                        for (auto& hitsref : shot.virtual_boundary_hits) {
                            if (shot.i_solved_vb[i]) {
                                if (DEBUG) t_out << "  vb" << i << std::flush;
                                bit_packed_res +=
                                    pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                        solver, hitsref);
                            }
                            shot.i_solved_vb[i] = false;
                            i++;
                        }
                        if (DEBUG) t_out << std::endl << std::flush;
        #else
                        pm::MatchingResult bit_packed_res =
                            pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                solver, shot.sparse_shot.hits);
        #endif
                        if (!solver.flooder.negative_weight_detection_events.empty())
                            bit_packed_res +=
                                shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                    solver, solver.flooder.negative_weight_detection_events);
                        // XOR in negative weight observable mask
                        bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
                        // Translate observable mask into bit vector
                        pm::fill_bit_vector_from_obs_mask(
                            bit_packed_res.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
                        // Add negative weight sum to blossom solution weight
                        shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
                    }
#ifdef USE_SHMEM
                    // BARRIER needed to prevent race on extraction/putting mem
                    // Need to make more robust
                    shmem_barrier_all();
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
// #if ENABLE_DRAW_FLAGS
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
// #if ENABLE_DRAW_FLAGS
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
