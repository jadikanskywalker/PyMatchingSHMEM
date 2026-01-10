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
#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"

#include <omp.h>
#include <fstream>
#include <filesystem>

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
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
    std::unique_ptr<stim::MeasureRecordWriter> writer,
    bool enable_correlations_,
    bool draw_frames_,
    int max_threads,
    stim::DetectorErrorModel dem /* For shot reading */
) {
    // Setup ReaderWriter and shot buffer
    shot_buffer = std::make_shared<pm::ShotBuffer>(std::move(reader), std::move(writer), num_partitions, num_virtual_boundaries, graph_ptr->num_observables);
    build_tasks_for_round_partitioning();
    // Read first NUM_ACTIVE_SHOTS_PER_UNIT shots
    bool shot_read = true;
    int i;
    for (i = 0; i < NUM_ACTIVE_SHOTS_PER_UNIT; ++i) {
        shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot_buffer->buffer[i].sparse_shot);
        if (shot_read) { // partition shot
            for (auto det : shot_buffer->buffer[i].sparse_shot.hits) {
                int part_id = node_part_id[det];
                if (part_id >= 0) { // partition
                    shot_buffer->buffer[i].partition_hits[part_id].push_back(det);
                } else { // virtual_boundary
                    shot_buffer->buffer[i].virtual_boundary_hits[-1*part_id-1].push_back(det);
                }
            }
        } else {
            // ++i;
            break;
        }
    }
    // Setup threads
    // num_solver_sets = max_threads / num_partitions;
    if (max_threads < num_partitions) {
        throw std::invalid_argument("The number of threads should be >= the number of partitions.");
    }
    // int num_threads = num_solver_sets * num_partitions;
    int num_threads = num_partitions;
    omp_set_num_threads(num_threads);
    // Setup solvers
    enable_correlations = enable_correlations_;
    draw_frames = draw_frames_;
    build_solvers(
        /*ensure_search_flooder_included=*/enable_correlations,
        /*enable_correlations=*/enable_correlations,
        num_threads
    );
    if (draw_frames) {
        auto coords = pm::pick_coords_for_drawing_from_dem(dem, 20);
        for (auto &s : solvers)
            s->coords = coords;
    }
}

// Builds balanced fusion tree assuming ruond-based partitioning
void pm::DecodingUnit::build_tasks_for_round_partitioning() {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // Build a full fusion tree: at most ~2*N tasks. Reserve to keep element addresses stable.
    for (auto &shot_container : shot_buffer->buffer) { // CAN BE IMPROVED!
        auto &tasks = shot_container.tasks;
        tasks.reserve(static_cast<size_t>(2 * num_partitions - 1));
        // Add tasks for each partitions
        int task_id;
        for (task_id = 0; task_id<num_partitions; ++task_id) {
            tasks.emplace_back(task_id, task_id);
        }
        // Save odd trailing task
        int tail_idx = -1;
        if (num_partitions % 2) {
            tail_idx = num_partitions-1;
        }
        // Add each level of fusions
        int last_step_starts = 0;
        int this_step_starts = task_id;
        int start = 0;
        int step = 2;
        while (start < num_partitions-1) {
            int counter = 0;
            int i;
            for (int i = start; i < num_partitions-1; i += step) {
                Task* left_child = &(tasks)[last_step_starts + 2*counter];
                int right_child_idx;
                if (last_step_starts + 2*counter + 1 < this_step_starts) {
                    right_child_idx = last_step_starts + 2*counter + 1;
                } else {
                    if (tail_idx >= 0) { // fuse last one with tail
                        // std::cout << "Fuse tail_idx: " << tail_idx << " at " << i << std::endl;
                        right_child_idx = tail_idx;
                        tail_idx = -1;
                    } else { // save tail
                        tail_idx = last_step_starts + 2*counter;
                        // std::cout << "Save tail_idx: " << tail_idx << std::endl;
                        break;
                    }
                }
                Task* right_child = &(tasks)[right_child_idx];
                tasks.emplace_back(task_id, i, left_child, right_child);
                ++task_id;
                ++counter;
            }
            if (last_step_starts + 2*counter < this_step_starts) { // save tail
                tail_idx = last_step_starts + 2*counter;
                // std::cout << "Save tail_idx: " << tail_idx << std::endl;
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

// Build per-thread solvers directly from a DEM by constructing a UserGraph once
// and producing independent Mwpm instances via to_mwpm for each thread.
// (Assumes first mwpm has already been built)
void pm::DecodingUnit::build_solvers(
    bool ensure_search_flooder_included,
    bool enable_correlations,
    int num_threads) {
    if (num_threads < 1)
        return;
    if (ensure_search_flooder_included || enable_correlations)
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    solvers.clear();
    solvers.reserve(static_cast<size_t>(num_threads*NUM_ACTIVE_SHOTS_PER_UNIT));
    for (int idx = 0; idx < NUM_ACTIVE_SHOTS_PER_UNIT; ++idx) {
        for (int t = 0; t < num_threads; ++t) {
            // Each solver shares the same MatchingGraph via shared_ptr.
            solvers.emplace_back(std::make_shared<pm::Mwpm>(pm::GraphFlooder(graph_ptr, idx)));
            solvers[t]->flooder.sync_negative_weight_observables_and_detection_events();
        }
    }
}

// Should only be called by thread that finishes a shot
// Assumes round-based partitioning and ascending detection_events
void pm::DecodingUnit::write_result_and_get_next_shot(int shot_buffer_id) {
    // --- acquire ---
    std::unique_lock<std::mutex> lock(shot_buffer->m);
    shot_buffer->cv.wait(lock, [&] {
        return shot_buffer->next_shot_buffer_id == shot_buffer_id;
    });
    auto &shot = shot_buffer->buffer[shot_buffer_id];
    // --- write results ---
    if (DEBUG) {
        std::cout << "T" << omp_get_thread_num() << " writing results for shot buffer " << shot_buffer_id << std::endl;
    }
    for (size_t k = 0; k < graph_ptr->num_observables; k++) {
        shot_buffer->writer->write_bit(shot.res.obs_crossed[k]);
    }
    shot_buffer->writer->write_end();
    // --- read next shot ---
    shot.clear();
    int last_shot_buffer_id = shot_buffer->last_shot_buffer_id.load(std::memory_order_acquire);
    if (last_shot_buffer_id < 0) {
        bool shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot.sparse_shot);
        if (shot_read) { // partition detection events
            for (auto det : shot.sparse_shot.hits) {
                int part_id = node_part_id[det];
                if (part_id >= 0) { // partition
                    shot.partition_hits[part_id].push_back(det);
                } else { // virtual_boundary
                    shot.virtual_boundary_hits[-1*part_id-1].push_back(det);
                }
            }
            ++shot.current_buffer_round;
            // if (DEBUG) {
            //     int i = 0;
            //     std::cout << "T" << omp_get_thread_num() << " read next shot for shot buffer" << shot_buffer_id << std::endl;
            //     for (auto part : shot.partition_hits) {
            //         std::cout << "Part" << i << " hits: ";
            //         for (int det : part) {
            //             std::cout << det << "  ";
            //         }
            //         std::cout << std::endl;
            //         i++;
            //     }
            //     for (auto part : shot.virtual_boundary_hits) {
            //         std::cout << "VB" << i << " hits: ";
            //         for (int det : part) {
            //             std::cout << det << "  ";
            //         }
            //         std::cout << std::endl;
            //         i++;
            //     }
            // }
        } else {
            // shot_buffer->all_shots_read.store(true, std::memory_order_release);
            if (shot_buffer_id > 0) {
                shot_buffer->last_shot_buffer_id = shot_buffer_id-1;
            } else {
                shot_buffer->last_shot_buffer_id = NUM_ACTIVE_SHOTS_PER_UNIT-1;
            }
        }
    } else if (last_shot_buffer_id == shot_buffer_id) {
        done = true;
    }
    // --- increment next_shot_buffer_id ---
    ++shot_buffer->next_shot_buffer_id;
    if (shot_buffer->next_shot_buffer_id >= NUM_ACTIVE_SHOTS_PER_UNIT) {
        shot_buffer->next_shot_buffer_id = 0;
    }
    // --- release ---
    lock.unlock();
    shot_buffer->cv.notify_all();
}

void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
    std::cout << "draw_frames=" << draw_frames << std::endl;
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
        int shot_id = 0;
        int shot_buffer_round = 0; // how many times buffer has looped
        int shot_buffer_id = 0;
        try {
        while (!done) {
            pm::Mwpm& solver = *solvers[shot_buffer_id*num_threads + tid];
            if (DEBUG) {
                t_out << std::endl << "Starting shot " << shot_id
                      << ", buffer round " << shot_buffer_round
                      << ", buffer_id " << shot_buffer_id << std::endl;
            }
            if (draw_frames) {
                std::filesystem::create_directories("out_parallel/frames/" + std::to_string(shot_id) + "/t" + std::to_string(tid));
            }
            auto& shot = shot_buffer->buffer[shot_buffer_id];
            if (shot.current_buffer_round.load(std::memory_order_acquire) < shot_buffer_round) { // wait
                if (DEBUG) {
                    t_out << "Waiting for next shot" << std::endl;
                }
                std::unique_lock<std::mutex> lock(shot_buffer->m);
                if (shot_buffer->last_shot_buffer_id.load(std::memory_order_acquire) >= 0) {
                    lock.unlock();
                    break;
                }
                shot_buffer->cv.wait(lock);
                int last_shot_buffer_id = shot_buffer->last_shot_buffer_id.load(std::memory_order_acquire);
                if (DEBUG) {
                    t_out << "Woke up, last_shot_buffer_id=" << last_shot_buffer_id << std::endl;
                }
                if (last_shot_buffer_id >= 0) {
                    break;
                }
                lock.unlock();
            }
            Task *t = &shot.tasks[tid];
            bool stolen = t->try_to_steal_leaf(shot_buffer_round);
            while (stolen) { // Got task
                if (DEBUG) {
                    t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                }
                auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                solve_task(solver, hitsref, t, tid, draw_frames, shot_id);
                if (t->parent) {
                    Task* sibling = (t->child_bit == 1) ? t->parent->right_child : t->parent->left_child;
                    // Try to steal parent
                    stolen = t->try_to_steal_parent();
                    t = t->parent;
                    // Try to steal sibling or descendent of sibling
                    if (!stolen && !sibling->is_fusion) {
                        stolen = sibling->try_to_steal_leaf(shot_buffer_round);
                        t = sibling;
                    }
                    if (DEBUG) {
                        if (t != nullptr) {
                            t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen << std::endl;
                        }
                    }
                } else { // I am root, extract solution
                    if (DEBUG) {
                        t_out << "T" << tid << " extracting solution" << std::endl;
                    }
                    if (num_observables > sizeof(pm::obs_int) * 8) {
                        solver.flooder.match_edges.clear();
                        pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(solver, shot.sparse_shot.hits);
                        if (!solver.flooder.negative_weight_detection_events.empty())
                            shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                                solver, solver.flooder.negative_weight_detection_events);
                        solver.extract_paths_from_match_edges(solver.flooder.match_edges, shot.res.obs_crossed.data(), shot.res.weight);
                        // XOR negative weight observables
                        for (auto& obs : solver.flooder.negative_weight_observables)
                            *(shot.res.obs_crossed.data() + obs) ^= 1;
                        // Add negative weight sum to blossom solution weight
                        shot.res.weight += solver.flooder.negative_weight_sum;
                    } else {
                        pm::MatchingResult bit_packed_res =
                            pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(solver, shot.sparse_shot.hits);
                        if (!solver.flooder.negative_weight_detection_events.empty())
                            bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                                solver, solver.flooder.negative_weight_detection_events);
                        // XOR in negative weight observable mask
                        bit_packed_res.obs_mask ^= solver.flooder.negative_weight_obs_mask;
                        // Translate observable mask into bit vector
                        pm::fill_bit_vector_from_obs_mask(bit_packed_res.obs_mask, shot.res.obs_crossed.data(), num_observables);
                        // Add negative weight sum to blossom solution weight
                        shot.res.weight = bit_packed_res.weight + solver.flooder.negative_weight_sum;
                    }
                    write_result_and_get_next_shot(shot_buffer_id);
                    break;
                }
            }
            // Move on to next shot buffer
            ++shot_id;
            ++shot_buffer_id;
            if (shot_buffer_id >= NUM_ACTIVE_SHOTS_PER_UNIT) {
                ++shot_buffer_round;
                shot_buffer_id = 0;
            }
        }
        } catch (const std::exception &e) {
            #pragma omp critical
            {
                std::cerr << "ERROR: Shot " << shot_id << " Thread " << tid << " caught exception: " << e.what() << std::endl;
            }
        } catch (...) {
            #pragma omp critical
            {
                std::cerr << "ERROR: Shot " << shot_id << " Thread " << tid << " caught unknown exception." << std::endl;
            }
        }
    }
}

void pm::DecodingUnit::solve_task(pm::Mwpm& solver, std::vector<uint64_t>& hits, Task* task, int tid, int draw_frames, int shot_id) {
    task->setup();
    solver.flooder.current_shot = shot_id;
    solver.prepare_for_task(task);
    // if (DEBUG && tid==0) {
    //     output_detector_nodes(solver, true);
    // }
    // auto& hitsref = (task->is_fusion) ? shot_buffer->buffer[shot_buffer_id].virtual_boundary_hits[task->part] : shot_buffer->buffer[shot_buffer_id].partition_hits[task->part];
    pm::process_timeline_until_completion(solver, hits, draw_frames, true, tid);
    task->mark_solved();
    // if (DEBUG && tid==0) {
    //     output_solution_state(solver, hits, true);
    // }
}