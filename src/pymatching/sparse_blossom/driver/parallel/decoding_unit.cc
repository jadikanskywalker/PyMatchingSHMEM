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

#define SHMEM_MAX_SYNC_STEPS 5
#define SHMEM_SYNC_SUMMARY 0
#define SHMEM_SYNC_RES 1

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
    pid = shmem_my_pe();
    other_pid = (pid == 0) ? 1 : 0;
    // --- Create communication contexts ---
    shmem_ctx_create(0, &summary_ctx);
    shmem_ctx_create(0, &fallback_ctx);
    if (!summary_ctx || !fallback_ctx) {
        throw std::invalid_argument("Failed to create SHMEM communication contexts.");
    }
#endif
    // --- Create shared matching graph ---
    auto user_graph =
        pm::detector_error_model_to_user_graph(detector_error_model, enable_correlations, num_distinct_weights);
#ifdef USE_SHMEM
    node_ephemeral_fields_ptr = static_cast<DetectorNodeEphemeralFields*>(
        shmem_malloc(NUM_BUFFERS_PER_UNIT * user_graph.nodes.size() * sizeof(DetectorNodeEphemeralFields)));
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
    atomics_ptr = static_cast<uint64_t*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * SHMEM_MAX_SYNC_STEPS * sizeof(uint64_t)));
    if (atomics_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric atomics buffer.");
    }
    fusion_summary_ptr = static_cast<FusionSummary*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * 2 * sizeof(FusionSummary)));
    closest_nb_nelems_per_buffer = graph.node_part_id.size() / graph.num_rounds * SHMEM_INTERSECTION_BUFFER_FACTOR;
    closest_nb_ptr = static_cast<std::pair<size_t, size_t>*>(
        shmem_malloc(closest_nb_nelems_per_buffer * NUM_BUFFERS_PER_UNIT * sizeof(std::pair<size_t, size_t>)));
    obs_crossed_nelems_per_buffer = graph.graph_ptr->num_observables * sizeof(uint8_t);
    obs_crossed_nelems_per_buffer = obs_crossed_nelems_per_buffer / sizeof(uint32_t) + (obs_crossed_nelems_per_buffer % sizeof(uint32_t) > 0); // minumum size uint32_t
    obs_crossed_ptr =
        static_cast<uint32_t*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * obs_crossed_nelems_per_buffer * sizeof(uint32_t)));
    for (size_t i=0; i < NUM_BUFFERS_PER_UNIT*obs_crossed_nelems_per_buffer; ++i) {
        obs_crossed_ptr[i] = 0;
    }
    // regions_matched_to_vb_ptr = static_cast<GraphFillRegionSubstates*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * 2 *
    // SHMEM_MAX_D * sizeof(GraphFillRegionSubstates)));
#endif
    // --- Create shot buffer ---
    shot_buffer = std::make_shared<pm::ShotBuffer>(
#ifdef USE_SHMEM
        atomics_ptr,
        (uint8_t*)obs_crossed_ptr,
        obs_crossed_nelems_per_buffer * sizeof(uint32_t),
#endif
        std::move(reader),
        std::move(writer),
        graph.num_partitions,
        graph.num_virtual_boundaries,
        graph.graph_ptr->num_observables);
    // --- Build tasks ---
    build_tasks_for_round_partitioning();
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
    if (max_threads < graph.num_partitions / 2) {
        throw std::invalid_argument("The number of threads should be >= half the number of partitions.");
    }
    num_threads = graph.num_partitions / 2;  // Current TEST Case for dividing across two PE's
    num_partition_units = 1;
#else
    num_threads =
        (max_threads > graph.num_partitions) ? graph.num_partitions : max_threads;  // max num_partitions threads
    num_partition_units = graph.num_partitions / num_threads + (graph.num_partitions % num_threads > 0);
#endif
    std::cout << "graph.num_partitions = " << graph.num_partitions << "; num_threads: " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
#ifdef USE_SHMEM
    // --- Allocate buffer for GraphFillRegion arenas ---
    if (DEBUG) {
        std::cout << "DEBUG: Allocating Regions\n" << std::flush;
    }
    regions_nelems_per_solver = graph.node_part_id.size() / graph.num_partitions * SHMEM_ARENA_BUFFER_FACTOR;
    if (regions_nelems_per_solver >= 64) {
        regions_nelems_per_solver =
            ((regions_nelems_per_solver * SHMEM_ARENA_BUFFER_FACTOR) / 64) * 64;  // make multiple of 64
    } else {
        regions_nelems_per_solver = 64;
    }
    regions_ptr = static_cast<GraphFillRegion*>(
        shmem_malloc(regions_nelems_per_solver * (num_threads + 1) * NUM_BUFFERS_PER_UNIT * sizeof(GraphFillRegion))); // includes per-buffer buffer for other PE
    if (regions_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric region buffer.");
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
        if (pid < other_pid) {
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
void pm::DecodingUnit::build_tasks_for_round_partitioning() {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // Build a full fusion tree: less than 2*N tasks. Reserve to keep element addresses stable.
    for (auto& shot_container : shot_buffer->buffer) {  // Lazy repeat per shot container -- CAN BE IMPROVED!
        auto& tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * graph.num_partitions - 1));
        // Add tasks for each partitions
        int task_id;
        for (task_id = 0; task_id < graph.num_partitions; ++task_id) {
            tasks.emplace_back(task_id, task_id);
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
                tasks.emplace_back(task_id, i, left_child, right_child);
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
#ifdef USE_SHMEM
        tasks[task_id - 1].is_cross_pe = true;  // for 2 PE simple case
#endif
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
    solvers.reserve(static_cast<size_t>(num_threads * NUM_BUFFERS_PER_UNIT));
    for (int idx = 0; idx < NUM_BUFFERS_PER_UNIT; ++idx) {
        for (int t = 0; t < num_threads * num_partition_units; ++t) {
            // Each solver shares the same MatchingGraph via shared_ptr.
            solvers.emplace_back(
                std::make_shared<pm::Mwpm>(pm::GraphFlooder(
                    graph.graph_ptr,
                    idx
#ifdef USE_SHMEM
                    ,
                    regions_ptr + (idx * (num_threads + 1) + t) * regions_nelems_per_solver,
                    regions_nelems_per_solver
#endif
                    )));
            solvers[t]->flooder.sync_negative_weight_observables_and_detection_events();
        }
    }
}

#ifdef USE_SHMEM
void pm::DecodingUnit::solve_cross_process_fusion_and_get_next_shot(
    size_t shot_container_id, Task* task, size_t tid, size_t num_threads, size_t solver_id, size_t shot_id) {
    size_t trivial_ctr = 0;
    // --- setup task and solver ---
    auto& shot = shot_buffer->buffer[shot_container_id];
    auto& hitsref = shot.virtual_boundary_hits[task->part];
    auto& solver = *solvers[solver_id];
    task->setup();
    solver.prepare_for_task(task, shot_container_id*NUM_BUFFERS_PER_UNIT);
    // --- compute lower pid, index ---
    size_t lower_pid, my_index, other_index;
    if (pid < other_pid) {
        lower_pid = pid;
        my_index = 0;
        other_index = 1;
    } else {
        lower_pid = other_pid;
        my_index = 1;
        other_index = 0;
    }
    uint64_t* shot_status_base = shot_buffer->shot_container_status + shot_container_id*2;
    // uint64_t* shot_counter_base = shot_status_base + 1;
    // shmem_wait_until(shot_counter_base, SHMEM_CMP_EQ, shot_id);
    // shmem_barrier_all();
    // shmem_wait_until(shot_status_base + SHMEM_SYNC_BEGIN, SHMEM_CMP_EQ, shot_id);
    // --- wait until PE ready ---
    if (DEBUG) {
        std::cout << "Cross PE Fusion: PE" << pid << " shot " << shot_id << std::endl << std::flush;
    }
    // --- send fusion summary info ---
    FusionSummary* fusion_summary_base = fusion_summary_ptr + shot_container_id * 2 + my_index;
    new (fusion_summary_base) FusionSummary(
        task->regions_matched_to_virtual_boundary.size(),
        regions_ptr + num_threads * NUM_BUFFERS_PER_UNIT + tid);  // partition_units must == 1
    // std::cout << "  PE" << pid << " putting FusionSummary(" << task->regions_matched_to_virtual_boundary.size() << ", "
    //           << regions_ptr + num_threads * NUM_BUFFERS_PER_UNIT + tid << ")\n"
    //           << std::flush;
    shmem_ctx_putmem_signal(summary_ctx, fusion_summary_base, fusion_summary_base, sizeof(FusionSummary), shot_status_base + SHMEM_SYNC_SUMMARY, shot_id, SHMEM_SIGNAL_SET, other_pid);
    // --- signal ready ---
    shmem_signal_wait_until(shot_status_base + SHMEM_SYNC_SUMMARY, SHMEM_CMP_EQ, shot_id);
    // std::cout << "  PE" << pid << " got summary\n" << std::flush;
    // --- check summaries ---
    auto& my_summary = fusion_summary_ptr[my_index];
    auto& other_summary = fusion_summary_ptr[other_index];
    if (my_summary.regions_matched_to_vb_size == 0 && other_summary.regions_matched_to_vb_size == 0 &&
        hitsref.size() == 0) {
        if (DEBUG && pid == lower_pid) {
            std::cout << "  trivial case " << shot_id << "\n" << std::flush;
        }
        trivial_ctr++;
    }
    // --- extract solution ---
    if (shot.num_observables > sizeof(pm::obs_int) * 8) {
        solver.flooder.match_edges.clear();
        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(solver, shot.sparse_shot.hits);
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
        pm::MatchingResult bit_packed_res =
            pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                solver, shot.sparse_shot.hits);
        if (!solver.flooder.negative_weight_detection_events.empty())
            bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                solver, solver.flooder.negative_weight_detection_events);
        // XOR in negative weight observable mask
        bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
        // Translate observable mask into bit vector
        pm::fill_bit_vector_from_obs_mask(
            bit_packed_res.obs_mask, shot.res.obs_crossed.data(), shot.num_observables);
        // Add negative weight sum to blossom solution weight
        shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
    }
    uint32_t* obs_crossed_base = obs_crossed_ptr + shot_container_id * obs_crossed_nelems_per_buffer;
    if (DEBUG) {
        std::cout << "  PE" << pid << " res: " << *obs_crossed_base << std::endl << std::flush;
    }
    if (pid != lower_pid) {
        shmem_ctx_uint32_atomic_xor(summary_ctx, obs_crossed_base, *obs_crossed_base, lower_pid);
    }
#if ENABLE_DRAW_FRAMES
    else if (draw_frames) {
        std::filesystem::create_directory("out_parallel/frames/" + std::to_string(shot_id+1));
    }
#endif
    shmem_ctx_uint64_atomic_set(summary_ctx, shot_status_base + SHMEM_SYNC_RES, shot_id, other_pid);
    shmem_wait_until(shot_status_base + SHMEM_SYNC_RES, SHMEM_CMP_EQ, shot_id);
    // shmem_barrier_all();

    if (DEBUG && pid == lower_pid) {

        std::cout << "  PE" << pid << " res is " << *obs_crossed_base << std::endl << std::flush;
    }
    shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id, pid == lower_pid);

    // if (pid == lower_pid) { // wait and write
        // shmem_wait_until(shot_status_base, SHMEM_CMP_EQ, PUT_RESULT);
        // std::cout << "  PE" << pid << " got results\n" << std::flush;
        // *shot_counter_base++;
        // shmem_ctx_uint64_atomic_inc(summary_ctx, shot_counter_base, other_pid);
    // } else { // wait to move one
        // shmem_ctx_uint64_atomic_set(summary_ctx, shot_status_base, PUT_RESULT, other_pid);
        // std::cout << "  PE" << pid << " signaled put results\n" << std::flush;
        // shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id, pid == lower_pid);
    // }

    // --- reset atomic synchronizers ---

    // // --- send closest neighbor info ---
    // std::pair<size_t, size_t>* closest_nb_base = closest_nb_ptr + closest_nb_nelems_per_buffer * shot_container_id;
    // auto& bounds = graph.vb_bounds[task->part];
    // for (size_t node = bounds.first; node <= bounds.second; node++) {
    //     new (closest_nb_base + node) std::pair<size_t, size_t>();
    // }
    // shmem_ctx_putmem_nbi(summary_ctx, closest_nb_base, closest_nb_base, closest_nb_nelems_per_buffer, other_pid);

    // if (shot_status == 1) {  // I perform fusion
    //     std::cout << "PE " << pid << " to perform fusion " << shot_container_id << "\n" << std::flush;
    //     // - check fusion_summary
    //     std::cout << "PE " << pid << " reading FusionSummary("
    //               << fusion_summary_ptr[other_index].regions_matched_to_vb_size << ", "
    //               << fusion_summary_ptr[other_index].regions_begin << ")\n"
    //               << std::flush;
    //     // - get other PE's solution state -
    //     // // get nodes (CHANGE TO PUTS!!!)
    //     // auto& bounds = graph.partition_bounds[other_pid];
    //     // auto* nodes_base = node_ephemeral_fields_ptr + bounds.first;
    //     // size_t nodes_nelems = bounds.second - bounds.first + 1;
    //     // shmem_ctx_getmem_nbi(fallback_ctx, nodes_base, nodes_base, nodes_nelems *
    //     // sizeof(DetectorNodeEphemeralFields), other_pid); std::cout << "PE " << pid << " got nodes for p" << other_pid
    //     // << std::endl; get regions (Need to add other PE region buffer, CHANGE TO PUTS!!!) auto* regions_base =
    //     // regions_ptr + (shot_container_id*num_threads)*regions_nelems_per_solver; size_t regions_nelems =
    //     // regions_nelems_per_solver * num_threads; shmem_ctx_getmem_nbi(fallback_ctx, regions_base, regions_base,
    //     // regions_nelems * sizeof(GraphFillRegion), other_pid); std::cout << "PE " << pid << " got regions for p" <<
    //     // other_pid << std::endl;
    //     // - fuse -
    //     if (my_summary.regions_matched_to_vb_size == 0 &&
    //         other_summary.regions_matched_to_vb_size == 0) {  // trivial case -- no cross-PE changes needed
    //         // extract solution & write results
    //     } else {
    //         shmem_uint64_atomic_set(shot_buffer->shot_container_status + shot_container_id, 3, other_pid);
    //     }
    //     // --- write results & get next shot ---
    //     shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
    // } else {  // I wait on fusion
    //     shmem_wait_until(shot_status_base, SHMEM_CMP_GT, 1);
    //     // Other PE has posted
    //     std::cout << "PE " << pid << " reading FusionSummary("
    //               << fusion_summary_ptr[other_index].regions_matched_to_vb_size << ", "
    //               << fusion_summary_ptr[other_index].regions_begin << ")\n"
    //               << std::flush;
    //     if (my_summary.regions_matched_to_vb_size == 0 &&
    //         other_summary.regions_matched_to_vb_size == 0) {  // trivial case -- no cross-PE changes needed
    //         // extract solution
    //     } else {
    //         shmem_wait_until(shot_status_base, SHMEM_CMP_EQ, 3);
    //         std::cout << "PE " << pid << " freed " << shot_container_id << "\n" << std::flush;
    //     }
    // }
}
#endif

// Core parallel decoding loop
void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        std::ofstream t_out;
        if (DEBUG) {
            std::string t_out_name = "out_parallel/";
#ifdef USE_SHMEM
            t_out_name += "p" + std::to_string(pid);
#endif
            t_out_name += "out_parallel/t" + std::to_string(tid) + ".out";
            t_out = (std::ofstream)(t_out_name);
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
                    std::filesystem::create_directory(frames_out);
#endif
                    std::filesystem::create_directory(frames_out + "/t" + std::to_string(tid));
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
                int partition_unit = 0;
                int solver_id = num_threads * (num_partition_units * shot_container_id + partition_unit) + tid;
                if (DEBUG)
                    t_out << "solvers[" << solver_id << "]\n";
                Task* t = &shot.tasks
                               [tid
#ifdef USE_SHMEM
                                + num_threads * pid  // tmp 2 PE case
#endif
                ];
                int next_leaf_id = tid + num_threads;
                bool stolen = t->try_to_steal_leaf(shot_buffer_round);
                while (stolen) {  // Got task
                    if (DEBUG) {
                        t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                    }
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    solve_task(
                        *solvers[solver_id],
                        hitsref,
                        t,
                        tid
#if ENABLE_DRAW_FLAGS
                        ,
                        draw_frames
#endif
                        ,
                        shot_id);
                    if (t->parent) {
#ifdef USE_SHMEM
                        if (t->parent->is_cross_pe) {  // if parent is cross_pe... tmp 2 PE fusion case
                            solve_cross_process_fusion_and_get_next_shot(
                                shot_container_id, t->parent, tid, num_threads, solver_id, shot_id);
                            break;
                        }
#endif
                        Task* sibling = (t->child_bit == 1) ? t->parent->right_child : t->parent->left_child;
                        // Try to steal parent
                        // t = t->try_to_steal_parent_or_descendent(shot_buffer_round);
                        // stolen = (t != nullptr);
                        stolen = t->try_to_steal_parent();
                        t = t->parent;
                        // Try to steal sibling or descendent of sibling
                        if (!stolen && !sibling->is_fusion) {
                            stolen = sibling->try_to_steal_leaf(shot_buffer_round);
                            t = sibling;
                        }
                        if (!stolen && next_leaf_id < graph.num_partitions) {
                            // Try next partition block
                            ++partition_unit;
                            solver_id = num_threads * (num_partition_units * shot_container_id + partition_unit) + tid;
                            if (DEBUG)
                                t_out << "solvers[" << solver_id << "]\n";
                            t = &shot.tasks[next_leaf_id];
                            stolen = t->try_to_steal_leaf(shot_buffer_round);
                            next_leaf_id += num_threads;
                        }
                        if (DEBUG) {
                            if (t != nullptr) {
                                t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen
                                      << std::endl;
                            }
                        }
                    } else {  // I am root, extract solution
                        if (DEBUG) {
                            t_out << "T" << tid << " extracting solution" << std::endl;
                        }
                        auto& solver = *solvers[solver_id];
                        if (shot.num_observables > sizeof(pm::obs_int) * 8) {
                            solver.flooder.match_edges.clear();
                            pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                solver, shot.sparse_shot.hits);
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
                            pm::MatchingResult bit_packed_res =
                                pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                    solver, shot.sparse_shot.hits);
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
                        shot_buffer->write_result_and_get_next_shot(shot_container_id, graph.node_part_id);
                        break;
                    }
                }
                // Move on to next shot buffer
                ++shot_id;
#if NUM_BUFFERS_PER_UNIT > 1
                ++shot_container_id;
                if (shot_container_id >= NUM_BUFFERS_PER_UNIT) {
                    ++shot_buffer_round;
                    shot_container_id = 0;
                }
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
void pm::DecodingUnit::solve_task(
    pm::Mwpm& solver,
    std::vector<uint64_t>& hits,
    Task* task,
    int tid
#if ENABLE_DRAW_FLAGS
    ,
    int draw_frames
#endif
    ,
    int shot_id) {
    task->setup();
    solver.prepare_for_task(task, shot_id);
    // if (DEBUG && tid==0) {
    //     output_detector_nodes(solver, true);
    // }
    pm::process_timeline_until_completion(
        solver,
        hits
#if ENABLE_DRAW_FLAGS
        ,
        draw_frames
#endif
        ,
        true,
        tid);
    task->mark_solved();
    // if (DEBUG && tid==0) {
    //     output_solution_state(solver, hits, true);
    // }
}
