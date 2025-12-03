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

#ifndef PYMATCHING2_DECODING_UNIT_H
#define PYMATCHING2_DECODING_UNIT_H

#include "../../config_parallel.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"

#include <vector>

namespace pm { class MatchingGraph; class Mwpm; }

// A decoding unit is a connected decoding graph
//   Connected decoding graphs are partitioned to allow parallel solving
//   Decoding task involve solving a partition or fusing partitions along virtual boundaries
struct DecodingUnit {
    const std::vector<std::vector<int>> partitions; 
    const std::vector<std::vector<int>> virtual_boundaries; // start and end index for each boundary
    const std::vector<std::vector<int>> virtual_boundary_partitions; // for each vb, list of partitions it connects

    const std::shared_ptr<pm::MatchingGraph> graph_ptr;
    std::vector<std::shared_ptr<pm::Mwpm>> solvers;

    int num_tasks;
    std::shared_ptr<std::vector<Task>> tasks;
    // std::vector<int> step_sizes;

    DecodingUnit(
        const std::shared_ptr<pm::MatchingGraph> graph_ptr,
        const std::vector<std::vector<int>>& partitions,
        const std::vector<std::vector<int>>& virtual_boundaries,
        const std::vector<std::vector<int>> virtual_boundary_partitions
    ) : graph_ptr(graph_ptr), partitions(partitions), virtual_boundaries(virtual_boundaries),
        virtual_boundary_partitions(virtual_boundary_partitions) {}
    // DecodingUnit(const DecodingUnit&) {
    //     tasks_ptr = tasks_ptr;
    // };
    // DecogingUnit& operator=(const DecodingUnit&) {
    //     tasks_ptr = tasks_ptr;
    // };
    // DecodingUnit(DecodingUnit&& other) noexcept {
    //     status.store(other.status.load());
    // }
    // Task& operator=(Task&& other) noexcept {
    //     status.store(other.status.load());
    //     return *this;
    // }

    void build_solvers(
        bool ensure_search_flooder_included,
        bool enable_correlations,
        int num_threads);

    void build_tasks_for_round_partitioning();
};


#endif // PYMATCHING2_DECODING_UNIT_H