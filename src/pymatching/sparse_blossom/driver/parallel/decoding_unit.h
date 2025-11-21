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

#include "pymatching/sparse_blossom/driver/parallel/decoding_tree.h"
#include <vector>


struct Partition {
    std::pair<int, int> nodes; // mask range of nodes in partition
};

struct VirtualBoundary {
    std::pair<int, int> nodes; // mask range of nodes in partition
    // std::vector<int> partitions; // partitions this boundary connects
};

// A decoding unit is a connected decoding graph
//   Connected decoding graphs are partitioned to allow parallel solving
//   Decoding task involve solving a partition or fusing partitions along virtual boundaries
struct DecodingUnit {
    const std::vector<int> node_mask; // id mask over all nodes: [ p0 vb0 p1 vb1 p2 ... ]
    const std::vector<Partition> partitions; 
    const std::vector<VirtualBoundary> virtual_boundaries; // start and end index for each boundary

    std::shared_ptr<MatchingGraph> graph_ptr;

    DecodingUnit(
        const std::vector<int>& node_mask,
        const std::vector<Partition>& partitions,
        const std::vector<VirtualBoundary>& virtual_boundaries
    ) : node_mask(node_mask), partitions(partitions), virtual_boundaries(virtual_boundaries) {}
};

#endif // PYMATCHING2_DECODING_UNIT_H