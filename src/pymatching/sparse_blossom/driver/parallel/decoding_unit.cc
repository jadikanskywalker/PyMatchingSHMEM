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

#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"

#ifdef USE_SHMEM
#include <shmem.h>
#endif

void pm::ShotContainer::clear() {
    sparse_shot.clear();
    for (auto& hits : partition_hits) {
        hits.clear();
    }
    for (auto& hits : virtual_boundary_hits) {
        hits.clear();
    }
    res.reset();
}

void pm::DecodingUnit::setup(
#ifdef USE_SHMEM
    void* &regions_ptr,
    uint64_t* atomics_ptr,
#endif
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
    std::unique_ptr<stim::MeasureRecordWriter> writer,
    bool enable_correlations_,
    bool draw_frames_,
    int max_threads,
    stim::DetectorErrorModel dem /* For shot reading */
) {
#ifdef USE_SHMEM
    // Initialize synchronization atomics
    pid = shmem_my_pe();
    other_pid = (pid == 0) ? 1 : 0;
#endif
    // Setup ReaderWriter and shot buffer
    shot_buffer = std::make_shared<pm::ShotBuffer>(
#ifdef USE_SHMEM
        atomics_ptr,
#endif
        std::move(reader), std::move(writer), num_partitions, num_virtual_boundaries, graph_ptr->num_observables);
    build_tasks_for_round_partitioning();
    // Read first NUM_BUFFERS_PER_UNIT shots
    bool shot_read = true;
    int i;
    for (i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot_buffer->buffer[i].sparse_shot);
        if (shot_read) {  // partition shot
            for (auto det : shot_buffer->buffer[i].sparse_shot.hits) {
                int part_id = node_part_id[det];
                if (part_id >= 0) {  // partition
                    shot_buffer->buffer[i].partition_hits[part_id].push_back(det);
                } else {  // virtual_boundary
                    shot_buffer->buffer[i].virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot_buffer->buffer[i].current_shot = i;
        } else {
            break;
        }
    }
    // Setup threads
#ifdef USE_SHMEM
    if (max_threads < num_partitions/2) {
        throw std::invalid_argument("The number of threads should be >= half the number of partitions.");
    }
    int num_threads = num_partitions;
#else
    if (max_threads < num_partitions) {
        throw std::invalid_argument("The number of threads should be >= the number of partitions.");
    }
    int num_threads = num_partitions/2;
#endif
    omp_set_num_threads(num_threads);
#ifdef USE_SHMEM
    // Allocate buffer for GraphFillRegion arenas
    int shmem_buffer_size = graph_ptr->nodes.size() / num_partitions * SHMEM_ARENA_BUFFER_FACTOR;
    if (shmem_buffer_size >= 64) {
        shmem_buffer_size = ((shmem_buffer_size*SHMEM_ARENA_BUFFER_FACTOR)/64)*64; // make multiple of 64
    } else {
        shmem_buffer_size = 64;
    }
    regions_ptr = shmem_malloc(shmem_buffer_size * num_threads * NUM_BUFFERS_PER_UNIT * sizeof(GraphFillRegion));
#endif
    // Setup solvers
    enable_correlations = enable_correlations_;
    draw_frames = draw_frames_;
    build_solvers(
        /*ensure_search_flooder_included=*/enable_correlations,
        /*enable_correlations=*/enable_correlations,
        num_threads
#ifdef USE_SHMEM
        , (GraphFillRegion*) regions_ptr,
        shmem_buffer_size
#endif
    );
    if (draw_frames) {
        auto coords = pm::pick_coords_for_drawing_from_dem(dem, 20);
        for (auto& s : solvers)
            s->coords = coords;
    }
}

// Builds balanced fusion tree assuming ruond-based partitioning
void pm::DecodingUnit::build_tasks_for_round_partitioning() {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // Build a full fusion tree: less than 2*N tasks. Reserve to keep element addresses stable.
    for (auto& shot_container : shot_buffer->buffer) {  // Lazy repeat per shot container -- CAN BE IMPROVED!
        auto& tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * num_partitions - 1));
        // Add tasks for each partitions
        int task_id;
        for (task_id = 0; task_id < num_partitions; ++task_id) {
            tasks.emplace_back(task_id, task_id);
        }
        // Save odd trailing task
        int tail_idx = -1;
        if (num_partitions % 2) {
            tail_idx = num_partitions - 1;
        }
        // Add each level of fusions
        int last_step_starts = 0;
        int this_step_starts = task_id;
        int start = 0;
        int step = 2;
        while (start < num_partitions - 1) {
            int counter = 0;
            int i;
            for (int i = start; i < num_partitions - 1; i += step) {
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
void pm::DecodingUnit::build_solvers(bool ensure_search_flooder_included, bool enable_correlations, int num_threads
#ifdef USE_SHMEM
    , GraphFillRegion* regions_ptr, size_t regions_nelems_per_solver
#endif
) {
    if (num_threads < 1)
        return;
    if (ensure_search_flooder_included || enable_correlations)
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    solvers.clear();
    solvers.reserve(static_cast<size_t>(num_threads * NUM_BUFFERS_PER_UNIT));
    for (int idx = 0; idx < NUM_BUFFERS_PER_UNIT; ++idx) {
        for (int t = 0; t < num_threads; ++t) {
            // Each solver shares the same MatchingGraph via shared_ptr.
            solvers.emplace_back(std::make_shared<pm::Mwpm>(pm::GraphFlooder(graph_ptr, idx
#ifdef USE_SHMEM
                , regions_ptr + idx*regions_nelems_per_solver + t*regions_nelems_per_solver, regions_nelems_per_solver
#endif
            )));
            solvers[t]->flooder.sync_negative_weight_observables_and_detection_events();
        }
    }
}

#ifdef USE_SHMEM
void pm::DecodingUnit::handle_cross_process_fusion_and_get_next_shot(int shot_container_id) {
    // --- acquire thread lock ---
    std::unique_lock<std::mutex> lock(shot_buffer->m);
    shot_buffer->cv.wait(lock, [&] {
        return shot_buffer->next_shot_container_id == shot_container_id;
    });
    auto& shot = shot_buffer->buffer[shot_container_id];
    // --- try to win cross-process fusion ---
    uint64_t shot_status = shot_buffer->shot_container_status[shot_container_id];
    if (shot_status == 0) { // signal other PE
        shmem_atomic_inc(shot_buffer->shot_container_status + shot_container_id, other_pid);
        shmem_quiet();
        if (shot_status == 0) {
            shot_status = shot_buffer->shot_container_status[shot_container_id];
        }
    }
}

#else
// Should only be called by the thread that finishes a shot
//   Assumes round-based partitioning and detection_events given in order of ascending detector node index
void pm::DecodingUnit::write_result_and_get_next_shot(int shot_container_id) {
    // --- acquire ---
    std::unique_lock<std::mutex> lock(shot_buffer->m);
    shot_buffer->cv.wait(lock, [&] {
        return shot_buffer->next_shot_container_id == shot_container_id;
    });
    auto& shot = shot_buffer->buffer[shot_container_id];
    // --- write results ---
    if (DEBUG) {
        std::cout << "T" << omp_get_thread_num() << " writing results for shot buffer " << shot_container_id << std::endl;
    }
    for (size_t k = 0; k < graph_ptr->num_observables; k++) {
        shot_buffer->writer->write_bit(shot.res.obs_crossed[k]);
    }
    shot_buffer->writer->write_end();
    // --- read next shot ---
    shot.clear();
    if (shot_buffer->last_shot_container_id < 0) {
        bool shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot.sparse_shot);
        if (shot_read) {  // partition detection events
            for (auto det : shot.sparse_shot.hits) {
                int part_id = node_part_id[det];
                if (part_id >= 0) {  // partition
                    shot.partition_hits[part_id].push_back(det);
                } else {  // virtual_boundary
                    shot.virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot.current_buffer_round++;
            shot.current_shot += NUM_BUFFERS_PER_UNIT;
        } else {
            if (shot_container_id > 0) {
                shot_buffer->last_shot_container_id = shot_container_id - 1;
            } else {
                shot_buffer->last_shot_container_id = NUM_BUFFERS_PER_UNIT - 1;
            }
        }
    } else if (shot_buffer->last_shot_container_id == shot_container_id) {
        for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
            shot_buffer->buffer[i].current_buffer_round.store(-1);
        }
    }
    // --- increment next_shot_container_id ---
    if (++shot_buffer->next_shot_container_id >= NUM_BUFFERS_PER_UNIT) {
        shot_buffer->next_shot_container_id = 0;
    }
    // --- release ---
    lock.unlock();
    shot_buffer->cv.notify_all();
}
#endif

#ifdef USE_SHMEM
// Core multi-process decoding loop
void pm::DecodingUnit::decode_shots_with_shmem() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
    size_t num_observables = graph_ptr->num_observables;
    shmem_barrier_all();
#pragma omp parallel shared(num_observables)
    {
        const int tid = omp_get_thread_num();
        const int num_threads = omp_get_num_threads();
        std::ofstream t_out;
        if (DEBUG) {
            t_out = (std::ofstream)("out_parallel/t" + std::to_string(tid) + ".out");
            std::cout << "T" << tid << " of " << num_threads << std::endl;
        }
        // Start decoding
        int shot_container_id = 0;
        uint64_t shot_buffer_round = shot_buffer->buffer[0].current_buffer_round.load();
        uint64_t shot_id = shot_buffer_round * NUM_BUFFERS_PER_UNIT;
        try {
            while (true) {
                if (DEBUG) {
                    t_out << std::endl
                          << "Starting shot " << shot_id << ", buffer round " << shot_buffer_round << ", buffer_id "
                          << shot_container_id << std::endl;
                }
                if (draw_frames) {
                    std::filesystem::create_directories(
                        "out_parallel/frames/" + std::to_string(shot_id) + "/p" + std::to_string(pid) + "/t" + std::to_string(tid));
                }
                auto& shot = shot_buffer->buffer[shot_container_id];
                int current_buffer_round = shot.current_buffer_round.load();
                while (current_buffer_round < shot_buffer_round) {  // wait
                    if (current_buffer_round < 0) {
                        break;
                    }
                    _mm_pause();
                    current_buffer_round = shot.current_buffer_round.load();
                }
                if (current_buffer_round < 0) {
                    break;
                }
                pm::Mwpm& solver = *solvers[shot_container_id * num_threads + tid];
                Task* t = &shot.tasks[pid * num_threads + tid];
                bool stolen = t->try_to_steal_leaf(shot_buffer_round);
                while (stolen) {  // Got task
                    if (!t->parent) { // Got root task
                        handle_cross_process_fusion_and_get_next_shot(current_buffer_round);
                        break;
                    }
                    if (DEBUG) {
                        t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                    }
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    solve_task(solver, hitsref, t, tid, draw_frames, shot_id);
                    Task* sibling = (t->child_bit == 1) ? t->parent->right_child : t->parent->left_child;
                    // Try to steal parent
                    t = t->try_to_steal_parent_or_descendent(shot_buffer_round);
                    stolen = (t != nullptr);
                    if (DEBUG) {
                        if (t != nullptr) {
                            t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen
                                  << std::endl;
                        }
                    }
                    //     if (DEBUG) {
                    //         t_out << "T" << tid << " extracting solution" << std::endl;
                    //     }
                    //     if (num_observables > sizeof(pm::obs_int) * 8) {
                    //         solver.flooder.match_edges.clear();
                    //         pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                    //             solver, shot.sparse_shot.hits);
                    //         if (!solver.flooder.negative_weight_detection_events.empty())
                    //             shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                    //                 solver, solver.flooder.negative_weight_detection_events);
                    //         solver.extract_paths_from_match_edges(
                    //             solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
                    //         // XOR negative weight observables
                    //         for (auto& obs : solver.flooder.negative_weight_observables)
                    //             *(shot.res.obs_crossed.data() + obs) ^= 1;
                    //         // Add negative weight sum to blossom solution weight
                    //         shot.res.weight += solver.flooder.negative_weight_sum;
                    //     } else {
                    //         pm::MatchingResult bit_packed_res =
                    //             pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                    //                 solver, shot.sparse_shot.hits);
                    //         if (!solver.flooder.negative_weight_detection_events.empty())
                    //             bit_packed_res +=
                    //                 shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                    //                     solver, solver.flooder.negative_weight_detection_events);
                    //         // XOR in negative weight observable mask
                    //         bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
                    //         // Translate observable mask into bit vector
                    //         pm::fill_bit_vector_from_obs_mask(
                    //             bit_packed_res.obs_mask, shot.res.obs_crossed.data(), num_observables);
                    //         // Add negative weight sum to blossom solution weight
                    //         shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
                    //     }
                }
                // Move on to next shot buffer
                ++shot_id;
                ++shot_container_id;
                if (shot_container_id >= NUM_BUFFERS_PER_UNIT) {
                    ++shot_buffer_round;
                    shot_container_id = 0;
                }
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

#else
// Core multi-threaded multi-active-shot decoding loop
void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
    size_t num_observables = graph_ptr->num_observables;
#pragma omp parallel shared(num_observables)
    {
        const int tid = omp_get_thread_num();
        const int num_threads = omp_get_num_threads();
        std::ofstream t_out;
        if (DEBUG) {
            t_out = (std::ofstream)("out_parallel/t" + std::to_string(tid) + ".out");
            std::cout << "T" << tid << " of " << num_threads << std::endl;
        }
        // Start decoding
        int shot_container_id = 0;
        uint64_t shot_buffer_round = shot_buffer->buffer[0].current_buffer_round.load(); // how many times buffer has looped
        uint64_t shot_id = shot_buffer_round * NUM_BUFFERS_PER_UNIT;
        try {
            while (true) {
                if (DEBUG) {
                    t_out << std::endl
                          << "Starting shot " << shot_id << ", buffer round " << shot_buffer_round << ", buffer_id "
                          << shot_container_id << std::endl;
                }
                if (draw_frames) {
                    std::filesystem::create_directories(
                        "out_parallel/frames/" + std::to_string(shot_id) + "/t" + std::to_string(tid));
                }
                auto& shot = shot_buffer->buffer[shot_container_id];
                int current_buffer_round = shot.current_buffer_round.load();
                while (current_buffer_round < shot_buffer_round) {  // wait
                    if (current_buffer_round < 0) {
                        break;
                    }
                    _mm_pause();
                    current_buffer_round = shot.current_buffer_round.load();
                }
                if (current_buffer_round < 0) {
                    break;
                }
                pm::Mwpm& solver = *solvers[shot_container_id * num_threads + tid];
                Task* t = &shot.tasks[tid];
                bool stolen = t->try_to_steal_leaf(shot_buffer_round);
                while (stolen) {  // Got task
                    if (DEBUG) {
                        t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                    }
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    solve_task(solver, hitsref, t, tid, draw_frames, shot_id);
                    if (t->parent) {
                        Task* sibling = (t->child_bit == 1) ? t->parent->right_child : t->parent->left_child;
                        // Try to steal parent
                        t = t->try_to_steal_parent_or_descendent(shot_buffer_round);
                        stolen = (t != nullptr);
                        // stolen = t->try_to_steal_parent();
                        // t = t->parent;
                        // Try to steal sibling or descendent of sibling
                        // if (!stolen && !sibling->is_fusion) {
                        //     stolen = sibling->try_to_steal_leaf(shot_buffer_round);
                        //     t = sibling;
                        // }
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
                        if (num_observables > sizeof(pm::obs_int) * 8) {
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
                                bit_packed_res.obs_mask, shot.res.obs_crossed.data(), num_observables);
                            // Add negative weight sum to blossom solution weight
                            shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
                        }
                        write_result_and_get_next_shot(shot_container_id);
                        break;
                    }
                }
                // Move on to next shot buffer
                ++shot_id;
                ++shot_container_id;
                if (shot_container_id >= NUM_BUFFERS_PER_UNIT) {
                    ++shot_buffer_round;
                    shot_container_id = 0;
                }
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
#endif

// Individidual task solver
void pm::DecodingUnit::solve_task(
    pm::Mwpm& solver, std::vector<uint64_t>& hits, Task* task, int tid, int draw_frames, int shot_id) {
    task->setup();
    solver.flooder.current_shot = shot_id;
    solver.prepare_for_task(task);
    // if (DEBUG && tid==0) {
    //     output_detector_nodes(solver, true);
    // }
    pm::process_timeline_until_completion(solver, hits, draw_frames, true, tid);
    task->mark_solved();
    // if (DEBUG && tid==0) {
    //     output_solution_state(solver, hits, true);
    // }
}