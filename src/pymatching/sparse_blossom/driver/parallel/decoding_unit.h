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

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"

#include <vector>

namespace pm { class MatchingGraph; class Mwpm; }

struct DetectionEventsContainer {
    std::vector<std::vector<uint64_t>> partition_hits;
    std::vector<std::vector<uint64_t>> virtual_boundary_hits;

    // int shot_id{ -1 };

    DetectionEventsContainer(int num_partitions, int num_virtual_boundaries)
     : partition_hits(num_partitions), virtual_boundary_hits(num_virtual_boundaries) {}

    void clear() {
        for (auto& hits : partition_hits) {
            hits.clear();
        }
        for (auto& hits : virtual_boundary_hits) {
            hits.clear();
        }
    }
};

// A decoding unit is a connected decoding graph
//   Connected decoding graphs are partitioned to allow parallel solving
//   Decoding task involve solving a partition or fusing partitions along virtual boundaries
struct DecodingUnit {
    // Graph
    const std::shared_ptr<pm::MatchingGraph> graph_ptr;

    const std::vector<int> node_part_id;
    const int num_partitions;
    const int num_virtual_boundaries;
    // const std::vector<std::vector<int>> partitions; 
    // const std::vector<std::vector<int>> virtual_boundaries; // start and end index for each boundary
    // const std::vector<std::vector<int>> virtual_boundary_partitions; // for each vb, list of partitions it connects

    // Fusion tree
    std::shared_ptr<std::vector<Task>> tasks;
    std::vector<long> virtual_boundary_markers;

    // Solvers
    std::vector<std::shared_ptr<pm::Mwpm>> solvers;

    // Detection Events
    int current_shot;
    // int next_shot;
    DetectionEventsContainer hits;

    DecodingUnit(
        const std::shared_ptr<pm::MatchingGraph> graph_ptr,
        const std::vector<int> node_part_id,
        int num_partitions,
        int num_virtual_boundaries
        // const std::vector<std::vector<int>>& partitions,
        // const std::vector<std::vector<int>>& virtual_boundaries,
        // const std::vector<std::vector<int>> virtual_boundary_partitions
    ) : graph_ptr(graph_ptr), node_part_id(node_part_id), num_partitions(num_partitions), num_virtual_boundaries(num_virtual_boundaries)
        // partitions(partitions), virtual_boundaries(virtual_boundaries), virtual_boundary_partitions(virtual_boundary_partitions)
        , hits(num_partitions, num_virtual_boundaries) {}

    void build_tasks_for_round_partitioning();

    void build_solvers(
        bool ensure_search_flooder_included,
        bool enable_correlations,
        int num_threads);

    void partition_detection_events(const std::vector<uint64_t>& detection_events);

    void solve_task(Task* task, int tid, int shot, int draw_frames);
};


#endif // PYMATCHING2_DECODING_UNIT_H