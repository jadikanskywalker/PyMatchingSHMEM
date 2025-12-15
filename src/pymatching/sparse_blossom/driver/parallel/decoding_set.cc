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
#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_set.h"
#include "stim.h"

#include <filesystem>
#include <fstream>
#include <omp.h>

void DecodingSet::build_solvers(int max_threads, stim::DetectorErrorModel dem) {
    int num_threads = max_threads;
    if (total_partitions < num_threads) {
        num_threads = total_partitions;
    }
    omp_set_num_threads(num_threads); // For one decoding unit at a time...
    // FOR ONE UNIT (no work stealing)
    units[0].build_solvers(
        /*ensure_search_flooder_included=*/enable_correlations,
        /*enable_correlations=*/enable_correlations,
        num_threads
    );
    if (draw_frames) {
        auto coords = pm::pick_coords_for_drawing_from_dem(dem, 20);
        for (auto &s : units[0].solvers)
            s->coords = coords;
    }
}

void DecodingSet::decode_shots() {
    if (enable_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }

    // FOR ONE UNIT PER SET
    auto &unit = units[0];
    auto &mwpm = *unit.solvers[0];
    size_t num_observables = mwpm.flooder.graph.num_observables;

    stim::SparseShot sparse_shot;
    sparse_shot.clear();
    pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);

    int shot_id = 0;
    bool shot_read = pm::start_and_read_entire_record_buffered(*reader, sparse_shot);

    #pragma omp parallel shared(unit, mwpm, num_observables, sparse_shot, res, shot_id, shot_read)
    {
        const int tid = omp_get_thread_num();
        const int num_threads = omp_get_num_threads();
        std::ofstream t_out;
        if (DEBUG) {
            t_out = (std::ofstream)("out_parallel/t" + std::to_string(tid) + ".out");
            std::cout << "T" << tid << " of " << num_threads << std::endl;
        }
        while (shot_read) {
            #pragma omp barrier
            #pragma omp single
            {
                if (DEBUG) {
                    std::cout << std::endl << "Starting shot " << shot_id << std::endl;
                }
                unit.partition_detection_events(sparse_shot.hits);
                if (draw_frames) {
                    std::filesystem::create_directory("out_parallel/frames/" + std::to_string(shot_id));
                }
            }
            try {
                if (DEBUG) {
                    t_out << std::endl << "Starting shot " << shot_id << std::endl;
                }
                if (draw_frames) {
                    std::filesystem::create_directory("out_parallel/frames/" + std::to_string(shot_id) + "/t" + std::to_string(tid));
                }
                std::queue<int> partition_tasks;
                t_out << "  partition_tasks: ";
                for (int i = tid; i < unit.num_partitions; i += num_threads) {
                    t_out << i << "  ";
                    partition_tasks.push(i);
                }
                t_out << std::endl;
                Task *t = &(*unit.tasks)[partition_tasks.front()];
                partition_tasks.pop();
                bool stolen = t->try_to_steal_leaf(shot_id);
                while (true) {
                    if (stolen) {
                        // Got task
                        if (DEBUG) {
                            t_out << "Thread " << tid << " solving " << (t->is_fusion ? "f" : "p") << t->part << std::endl;
                        }
                        unit.solve_task(t, tid, shot_id, draw_frames);
                        t->mark_solved(shot_id);
                        // Try to steal parent
                        if (t->parent) {
                            stolen = t->try_to_steal_parent();
                            t = t->parent;
                            if (DEBUG) {
                                t_out << "  f" << t->part << " stolen = " << stolen << std::endl;
                            }
                        } else {
                            stolen = false;
                        }
                    } else if (!partition_tasks.empty()) {
                        t = &(*unit.tasks)[partition_tasks.front()];
                        partition_tasks.pop();
                        stolen = t->try_to_steal_leaf(shot_id);
                        if (DEBUG) {
                            t_out << "  p" << t->part << " stolen = " << stolen << std::endl;
                        }
                    } else {
                        break;
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
            #pragma omp barrier
            #pragma omp single
            {
                if (num_observables > sizeof(pm::obs_int) * 8) {
                    mwpm.flooder.match_edges.clear();
                    pm::shatter_blossoms_for_all_detection_events_and_extract_match_edges(mwpm, sparse_shot.hits);
                    if (!mwpm.flooder.negative_weight_detection_events.empty())
                        shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                            mwpm, mwpm.flooder.negative_weight_detection_events);
                    mwpm.extract_paths_from_match_edges(mwpm.flooder.match_edges, res.obs_crossed.data(), res.weight);
                        
                    // XOR negative weight observables
                    for (auto& obs : mwpm.flooder.negative_weight_observables)
                        *(res.obs_crossed.data() + obs) ^= 1;
                    // Add negative weight sum to blossom solution weight
                    res.weight += mwpm.flooder.negative_weight_sum;
                } else {
                    pm::MatchingResult bit_packed_res =
                        pm::shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(mwpm, sparse_shot.hits);
                    if (!mwpm.flooder.negative_weight_detection_events.empty())
                        bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                            mwpm, mwpm.flooder.negative_weight_detection_events);
                    // XOR in negative weight observable mask
                    bit_packed_res.obs_mask ^= mwpm.flooder.negative_weight_obs_mask;
                    // Translate observable mask into bit vector
                    pm::fill_bit_vector_from_obs_mask(bit_packed_res.obs_mask, res.obs_crossed.data(), num_observables);
                    // Add negative weight sum to blossom solution weight
                    res.weight = bit_packed_res.weight + mwpm.flooder.negative_weight_sum;
                }
                // Write solution
                for (size_t k = 0; k < num_observables; k++) {
                    writer->write_bit(res.obs_crossed[k]);
                }
                writer->write_end();
                sparse_shot.clear();
                res.reset();
                // Read next shot
                ++shot_id;
                shot_read = pm::start_and_read_entire_record_buffered(*reader, sparse_shot);
            }
        }
    }
}