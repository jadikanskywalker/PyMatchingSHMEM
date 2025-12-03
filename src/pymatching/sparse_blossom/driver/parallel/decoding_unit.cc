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
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_unit.h"

// Build per-thread solvers directly from a DEM by constructing a UserGraph once
// and producing independent Mwpm instances via to_mwpm for each thread.
// (Assumes first mwpm has already been built)
void DecodingUnit::build_solvers(
    bool ensure_search_flooder_included,
    bool enable_correlations,
    int num_threads) {
    if (num_threads < 1)
        return;
    if (ensure_search_flooder_included || enable_correlations)
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    solvers.clear();
    solvers.reserve(static_cast<size_t>(num_threads));
    for (int t = 0; t < num_threads; ++t) {
        // Each solver shares the same MatchingGraph via shared_ptr.
        solvers.emplace_back(std::make_shared<pm::Mwpm>(pm::GraphFlooder(graph_ptr)));
        solvers[t]->flooder.sync_negative_weight_observables_and_detection_events();
    }
}

// Builds balanced fusion tree assuming ruond-based partitioning
void DecodingUnit::build_tasks_for_round_partitioning() {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    int num_partitions = partitions.size();
    // Build a full fusion tree: at most ~2*N tasks. Reserve to keep element addresses stable.
    tasks = std::make_shared<std::vector<Task>>();
    // step_sizes.clear();
    tasks->reserve(static_cast<size_t>(2 * num_partitions - 1));
    // Add tasks for each partitions
    int task_id;
    for (task_id = 0; task_id<num_partitions; ++task_id) {
        tasks->emplace_back(task_id, task_id);
    }
    // Save odd trailing task
    int tail_idx = -1;
    // int tails_vb;
    if (num_partitions % 2) {
        // tails_idx.emplace_back(num_partitions - 1);
        // tails_vb = num_partitions-2;
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
            Task* left_child = &(*tasks)[last_step_starts + 2*counter];
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
            Task* right_child = &(*tasks)[right_child_idx];
            tasks->emplace_back(task_id, i, left_child, right_child);
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
    
    if (DEBUG) {
        std::cout << "DEBUG:" << std::endl;
        for (Task& t : *tasks) {
            std::cout << "--part: " << t.part << std::endl
                      << "  is_fusion: " << t.is_fusion << std::endl
                      << "  left_child: " << t.left_child << std::endl
                      << "  right_child: " << t.right_child << std::endl
                      << "  parent: " << t.parent << std::endl;
        }
    }
}