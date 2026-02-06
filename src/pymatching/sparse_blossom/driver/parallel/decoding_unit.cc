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
#include <stdexcept>
#include <omp.h>

#include "pymatching/sparse_blossom/driver/user_graph.h"
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

pm::ShotBuffer::ShotBuffer(
#ifdef USE_SHMEM
    uint64_t* atomics_ptr,
#endif
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
    std::unique_ptr<stim::MeasureRecordWriter> writer,
    int num_partitions,
    int num_virtual_boundaries,
    int num_observables)
    : reader(std::move(reader)), writer(std::move(writer))
{
#ifdef USE_SHMEM
    shot_container_status = atomics_ptr;
    for (int i=0; i<NUM_BUFFERS_PER_UNIT; ++i) {
        shot_container_status[i] = 0;
    }
#endif
    buffer.reserve(static_cast<size_t>(NUM_BUFFERS_PER_UNIT));
    for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        buffer.emplace_back(num_partitions, num_virtual_boundaries, num_observables);
    }
}

pm::SharedMatchingGraph::SharedMatchingGraph() = default;

pm::SharedMatchingGraph::SharedMatchingGraph(
    std::shared_ptr<pm::MatchingGraph> graph_ptr_,
    std::vector<int> node_part_id_,
    int num_partitions_,
    int num_virtual_boundaries_)
: graph_ptr(graph_ptr_),
    node_part_id(node_part_id_),
    num_partitions(num_partitions_),
    num_virtual_boundaries(num_virtual_boundaries_) {
#ifdef USE_SHMEM
    construct_partition_bounds();
#endif
}

#ifdef USE_SHMEM
void pm::SharedMatchingGraph::construct_partition_bounds() {
    if (num_partitions <= 0) {
        throw std::invalid_argument("SharedMatchingGraph requires at least one partition.");
    }
    if (node_part_id.empty()) {
        throw std::invalid_argument("SharedMatchingGraph requires node partition assignments.");
    }
    partition_bounds.resize(num_partitions);
    partition_bounds[0].first = 0;
    size_t partition_bounds_idx = 0;
    int last_part_id = node_part_id[0];
    for (int i = 0; i < node_part_id.size(); ++i) {
        if (node_part_id[i] == last_part_id) {
            continue;
        } else if (node_part_id[i] < 0) { // entering v_b, end of last p
            partition_bounds[partition_bounds_idx].second = i-1;
            last_part_id = node_part_id[i];
        } else { // entering new p
            partition_bounds_idx++;
            partition_bounds[partition_bounds_idx].first = i;
            last_part_id = node_part_id[i];
        }
    }
    partition_bounds[partition_bounds_idx].second = node_part_id.size()-1;
}
#endif

pm::DecodingUnit::DecodingUnit(
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
    std::unique_ptr<stim::MeasureRecordWriter> writer,
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included,
    bool enable_correlations,
    bool draw_frames)
    : ensure_search_flooder_included(ensure_search_flooder_included),
      enable_correlations(enable_correlations),
      draw_frames(draw_frames)
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
        shmem_malloc(user_graph.nodes.size() * NUM_BUFFERS_PER_UNIT * sizeof(DetectorNodeEphemeralFields)));
    if (node_ephemeral_fields_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric detector node buffer.");
    }
#endif
    graph = user_graph.to_shared_matching_graph(
        num_distinct_weights
#ifdef USE_SHMEM
        , node_ephemeral_fields_ptr
#endif
    );
    if (graph.num_partitions <= 0) {
        throw std::invalid_argument("Graph partitioning produced no partitions. Check --rounds_per_partition.");
    }
#ifdef USE_SHMEM
    // --- Allocate sychronization atomics ---    
    atomics_ptr = static_cast<uint64_t*>(shmem_malloc(NUM_BUFFERS_PER_UNIT * sizeof(uint64_t)));
    if (atomics_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric atomics buffer.");
    }
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
    // --- Build tasks ---
    build_tasks_for_round_partitioning();
    // --- Read first NUM_BUFFERS_PER_UNIT shots ---
    bool shot_read = true;
    int i;
    for (i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot_buffer->buffer[i].sparse_shot);
        if (shot_read) {  // partition shot
            for (auto det : shot_buffer->buffer[i].sparse_shot.hits) {
                int part_id = graph.node_part_id[det];
                if (part_id >= 0) {  // partition
                    shot_buffer->buffer[i].partition_hits[part_id].push_back(det);
                } else {  // virtual_boundary
                    shot_buffer->buffer[i].virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot_buffer->buffer[i].current_shot = i;
        } else {
            shot_buffer->last_shot_container_id = i-1;
            if (shot_buffer->last_shot_container_id == NUM_BUFFERS_PER_UNIT) {
                shot_buffer->last_shot_container_id = 0;
            }
            break;
        }
    }
    // --- Set num_threads ---
    int max_threads = omp_get_max_threads();
#ifdef USE_SHMEM
    if (max_threads < graph.num_partitions/2) {
        throw std::invalid_argument("The number of threads should be >= half the number of partitions.");
    }
    num_threads = graph.num_partitions/2; // Current TEST Case for dividing across two PE's
    num_partition_units = 1;
#else
    num_threads = (max_threads > graph.num_partitions) ? graph.num_partitions : max_threads; // max num_partitions threads
    num_partition_units = graph.num_partitions / num_threads + (graph.num_partitions % num_threads > 0);
#endif
    std::cout << "graph.num_partitions = " << graph.num_partitions << "; num_threads: " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
#ifdef USE_SHMEM
    // --- Allocate buffer for GraphFillRegion arenas ---
    regions_nelems_per_solver =
        graph.node_part_id.size() / graph.num_partitions * SHMEM_ARENA_BUFFER_FACTOR;
    if (regions_nelems_per_solver >= 64) {
        regions_nelems_per_solver =
            ((regions_nelems_per_solver * SHMEM_ARENA_BUFFER_FACTOR) / 64) * 64;  // make multiple of 64
    } else {
        regions_nelems_per_solver = 64;
    }
    regions_ptr = static_cast<GraphFillRegion*>(
        shmem_malloc(regions_nelems_per_solver * num_threads * NUM_BUFFERS_PER_UNIT * sizeof(GraphFillRegion)));
    if (regions_ptr == nullptr) {
        throw std::invalid_argument("Failed to allocate symmetric region buffer.");
    }
#endif
    // --- Build solvers ---
    build_solvers();
    if (draw_frames) {
        auto coords = pm::pick_coords_for_drawing_from_dem(detector_error_model, 20);
        for (auto& s : solvers)
            s->coords = coords;
    }
}

pm::DecodingUnit::~DecodingUnit() {
#ifdef USE_SHMEM
    shmem_free(node_ephemeral_fields_ptr);
    shmem_free(atomics_ptr);
    shmem_free(regions_ptr);
#endif
}

// Builds balanced fusion tree assuming ruond-based partitioning
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
            for (int i = start; i < graph.num_partitions - 1; i += step) {
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
void pm::DecodingUnit::build_solvers() {
    if (ensure_search_flooder_included || enable_correlations) {
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    }
    solvers.clear();
    solvers.reserve(static_cast<size_t>(num_threads * NUM_BUFFERS_PER_UNIT));
    for (int idx = 0; idx < NUM_BUFFERS_PER_UNIT; ++idx) {
        for (int t = 0; t < num_threads*num_partition_units; ++t) {
            // Each solver shares the same MatchingGraph via shared_ptr.
            solvers.emplace_back(std::make_shared<pm::Mwpm>(pm::GraphFlooder(graph.graph_ptr, idx
#ifdef USE_SHMEM
                , regions_ptr + (idx*num_threads + t)*regions_nelems_per_solver, regions_nelems_per_solver
#endif
            )));
            solvers[t]->flooder.sync_negative_weight_observables_and_detection_events();
        }
    }
}

#ifdef USE_SHMEM
void pm::DecodingUnit::handle_cross_process_fusion_and_get_next_shot(int shot_container_id, int num_threads) {
    // --- acquire thread lock ---
    std::unique_lock<std::mutex> lock(shot_buffer->m);
    shot_buffer->cv.wait(lock, [&] {
        return shot_buffer->next_shot_container_id == shot_container_id;
    });
    auto& shot = shot_buffer->buffer[shot_container_id];
    // --- prepare fusion summary info ---

    // --- try to win cross-process fusion ---
    int lower_pid = (pid < other_pid) ? pid : other_pid;
    uint64_t* shot_status_ptr = shot_buffer->shot_container_status + shot_container_id;
    uint64_t shot_status = shmem_uint64_atomic_fetch_inc(shot_status_ptr, lower_pid);
    if (shot_status == 1) { // I perform fusion
        std::cout << "PE " << pid << " to perform fusion " << shot_container_id << "\n" << std::flush;
        // - get other PE's solution state -
        // get nodes
        auto& bounds = graph.partition_bounds[other_pid];
        auto* nodes_base = node_ephemeral_fields_ptr + bounds.first;
        size_t nodes_nelems = bounds.second - bounds.first + 1;
        shmem_getmem(nodes_base, nodes_base, nodes_nelems * sizeof(DetectorNodeEphemeralFields), other_pid);
        std::cout << "PE " << pid << " got nodes for p" << other_pid << std::endl;
        // get regions
        auto* regions_base =  regions_ptr + (shot_container_id*num_threads)*regions_nelems_per_solver;
        size_t regions_nelems = regions_nelems_per_solver * num_threads;
        shmem_getmem(regions_base, regions_base, regions_nelems * sizeof(GraphFillRegion), other_pid);
        std::cout << "PE " << pid << " got regions for p" << other_pid << std::endl;
        // - fuse -

        // - write results -
        // if (DEBUG) {
        //     std::cout << "T" << omp_get_thread_num() << " writing results for shot buffer " << shot_container_id << std::endl;
        // }
        // for (size_t k = 0; k < graph.graph_ptr->num_observables; k++) {
        //     shot_buffer->writer->write_bit(shot.res.obs_crossed[k]);
        // }
        // shot_buffer->writer->write_end();
        shmem_uint64_atomic_set(shot_buffer->shot_container_status + shot_container_id, 3, other_pid);
        *shot_status_ptr = 0;
    } else { // I wait on fusion 
        // std::cout << "PE " << pid << " waiting on fusion " << shot_container_id << "\n" << std::flush;
        
        shmem_wait_until(shot_status_ptr, SHMEM_CMP_EQ, 3);
        std::cout << "PE " << pid << " freed " << shot_container_id << "\n" << std::flush;
        *shot_status_ptr = 0;
    }
    // --- read next shot ---
    shot.clear();
    if (shot_buffer->last_shot_container_id < 0) {
        bool shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot.sparse_shot);
        if (shot_read) {  // partition detection events
            for (auto det : shot.sparse_shot.hits) {
                int part_id = graph.node_part_id[det];
                if (part_id >= 0) {  // partition
                    shot.partition_hits[part_id].push_back(det);
                } else {  // virtual_boundary
                    shot.virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot.current_buffer_round++;
            shot.current_shot += NUM_BUFFERS_PER_UNIT;
        } else {
#if NUM_BUFFERS_PER_UNIT > 1
            if (shot_container_id > 0) {
                shot_buffer->last_shot_container_id = shot_container_id - 1;
            } else {
                shot_buffer->last_shot_container_id = NUM_BUFFERS_PER_UNIT - 1;
            }
#else
            shot.current_buffer_round.store(-1);
#endif
        }
    } else if (shot_buffer->last_shot_container_id == shot_container_id) {
        for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
            shot_buffer->buffer[i].current_buffer_round.store(-1);
        }
    }
#if NUM_BUFFERS_PER_UNIT > 1
    // --- increment next_shot_container_id ---
    if (++shot_buffer->next_shot_container_id >= NUM_BUFFERS_PER_UNIT) {
        shot_buffer->next_shot_container_id = 0;
    }
#endif
    // --- release ---
    lock.unlock();
    shot_buffer->cv.notify_all();
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
    for (size_t k = 0; k < graph.graph_ptr->num_observables; k++) {
        shot_buffer->writer->write_bit(shot.res.obs_crossed[k]);
    }
    shot_buffer->writer->write_end();
    // --- read next shot ---
    shot.clear();
    if (shot_buffer->last_shot_container_id < 0) {
        bool shot_read = pm::start_and_read_entire_record_buffered(*shot_buffer->reader, shot.sparse_shot);
        if (shot_read) {  // partition detection events
            for (auto det : shot.sparse_shot.hits) {
                int part_id = graph.node_part_id[det];
                if (part_id >= 0) {  // partition
                    shot.partition_hits[part_id].push_back(det);
                } else {  // virtual_boundary
                    shot.virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot.current_buffer_round++;
            shot.current_shot += NUM_BUFFERS_PER_UNIT;
        } else {
#if NUM_BUFFERS_PER_UNIT > 1
            if (shot_container_id > 0) {
                shot_buffer->last_shot_container_id = shot_container_id - 1;
            } else {
                shot_buffer->last_shot_container_id = NUM_BUFFERS_PER_UNIT - 1;
            }
#else
            shot.current_buffer_round.store(-1);
#endif
        }
    } else if (shot_buffer->last_shot_container_id == shot_container_id) {
        for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
            shot_buffer->buffer[i].current_buffer_round.store(-1);
        }
    }
#if NUM_BUFFERS_PER_UNIT > 1
    // --- increment next_shot_container_id ---
    if (++shot_buffer->next_shot_container_id >= NUM_BUFFERS_PER_UNIT) {
        shot_buffer->next_shot_container_id = 0;
    }
#endif
    // --- release ---
    lock.unlock();
    shot_buffer->cv.notify_all();
}
#endif

// Core multi-process decoding loop
// void decode_shots_with_shmem() {
//     if (enable_correlations) {
//         throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
//     }
//     size_t num_observables = graph.graph_ptr->num_observables;
// #pragma omp parallel shared(num_observables)
//     {
//         const int tid = omp_get_thread_num();
//         std::ofstream t_out;
//         if (DEBUG) {
//             t_out = (std::ofstream)("out_parallel/p" + std::to_string(pid) + "_t" + std::to_string(tid) + ".out");
//             std::cout << "T" << tid << " of " << num_threads << std::endl;
//         }
//         // Start decoding
//         int shot_container_id = 0;
//         uint64_t shot_buffer_round = shot_buffer->buffer[0].current_buffer_round.load();
//         uint64_t shot_id = shot_buffer_round * NUM_BUFFERS_PER_UNIT;
//         try {
//             while (true) {
//                 if (DEBUG) {
//                     t_out << std::endl
//                           << "Starting shot " << shot_id << ", buffer round " << shot_buffer_round << ", buffer_id "
//                           << shot_container_id << std::endl;
//                 }
//                 if (draw_frames) {
//                     std::filesystem::create_directories(
//                         "out_parallel/frames/"  + std::to_string(pid) + "/t" + std::to_string(tid));
//                 }
//                 auto& shot = shot_buffer->buffer[shot_container_id];
//                 int current_buffer_round = shot.current_buffer_round.load();
//                 while (current_buffer_round < shot_buffer_round) {  // wait
//                     if (current_buffer_round < 0) {
//                         break;
//                     }
//                     _mm_pause();
//                     current_buffer_round = shot.current_buffer_round.load();
//                 }
//                 if (current_buffer_round < 0) { // exit condition
//                     break;
//                 }
//                 int partition_unit = 0;
//                 int solver_id = num_threads * (num_partition_units * shot_container_id + partition_unit) + tid;
//                 if (DEBUG) t_out << "solvers[" << solver_id << "]\n";
//                 Task* t = &shot.tasks[tid + num_threads*pid];
//                 // int next_leaf_id = tid + num_threads;
//                 bool stolen = t->try_to_steal_leaf(shot_buffer_round);
//                 while (stolen) {  // Got task
//                     if (DEBUG) {
//                         t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
//                     }
//                     auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
//                     solve_task(*solvers[solver_id], hitsref, t, tid, draw_frames, shot_id);
//                     Task* sibling = (t->child_bit == 1) ? t->parent->right_child : t->parent->left_child;
//                     // Try to steal parent
//                     // t = t->try_to_steal_parent_or_descendent(shot_buffer_round);
//                     // stolen = (t != nullptr);
                    
//                     stolen = t->try_to_steal_parent();
//                     t = t->parent;
//                     // // Try to steal sibling or descendent of sibling
//                     // if (!stolen && !sibling->is_fusion) {
//                     //     stolen = sibling->try_to_steal_leaf(shot_buffer_round);
//                     //     t = sibling;
//                     // }
//                     // if (!stolen && next_leaf_id < graph.num_partitions) {
//                     //     // Try next partition block
//                     //     ++partition_unit;
//                     //     solver_id = num_threads * (num_partition_units * shot_container_id + partition_unit) + tid;
//                     //     if (DEBUG) t_out << "solvers[" << solver_id << "]\n";
//                     //     t = &shot.tasks[next_leaf_id];
//                     //     stolen = t->try_to_steal_leaf(shot_buffer_round);
//                     //     next_leaf_id += num_threads;
//                     // }
//                     if (DEBUG) {
//                         if (t != nullptr) {
//                             t_out << (t->is_fusion ? "  f" : "  p") << t->part << " stolen = " << stolen
//                                     << std::endl;
//                         }
                    
//                     }
//                 }
//                 // Move on to next shot buffer
//                 ++shot_id;
// #if NUM_BUFFERS_PER_UNIT > 1
//                 ++shot_container_id;
//                 if (shot_container_id >= NUM_BUFFERS_PER_UNIT) {
//                     ++shot_buffer_round;
//                     shot_container_id = 0;
//                 }
// #endif
//             }
//         } catch (const std::exception& e) {
// #pragma omp critical
//             {
//                 std::cerr << "ERROR: Shot " << shot_id << " Thread " << tid << " caught exception: " << e.what()
//                           << std::endl;
//             }
//         } catch (...) {
// #pragma omp critical
//             {
//                 std::cerr << "ERROR: Shot " << shot_id << " Thread " << tid << " caught unknown exception."
//                           << std::endl;
//             }
//         }
//     }
// }

// Core multi-threaded multi-active-shot decoding loop
void pm::DecodingUnit::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }
    size_t num_observables = graph.graph_ptr->num_observables;
#pragma omp parallel shared(num_observables)
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
                        "out_parallel/frames/" 
#ifdef USE_SHMEM
                        + std::to_string(shot_id) + "/p"
#endif
                        + std::to_string(shot_id) + "/t" + std::to_string(tid));
                }
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
                if (DEBUG) t_out << "solvers[" << solver_id << "]\n";
                Task* t = &shot.tasks[tid
#ifdef USE_SHMEM
                    + num_threads*pid // tmp 2 PE case
#endif
                ];
                int next_leaf_id = tid + num_threads;
                bool stolen = t->try_to_steal_leaf(shot_buffer_round);
                while (stolen) {  // Got task
                    if (DEBUG) {
                        t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                    }
                    auto& hitsref = (t->is_fusion) ? shot.virtual_boundary_hits[t->part] : shot.partition_hits[t->part];
                    solve_task(*solvers[solver_id], hitsref, t, tid, draw_frames, shot_id);
                    if (t->parent) {
#ifdef USE_SHMEM
                        if (!t->parent->parent) { // if parent is root... tmp 2 PE fusion caser
                            handle_cross_process_fusion_and_get_next_shot(shot_container_id, num_threads);
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
                            if (DEBUG) t_out << "solvers[" << solver_id << "]\n";
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
                    }
#ifndef USE_SHMEM
                    else {  // I am root, extract solution
                        if (DEBUG) {
                            t_out << "T" << tid << " extracting solution" << std::endl;
                        }
                        pm::Mwpm& solver = *solvers[solver_id];
                        // output_solution_state(solver, shot.sparse_shot.hits, true);
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
#endif
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