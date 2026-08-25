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

#ifndef PYMATCHING2_SEARCH_DETECTOR_NODE_H
#define PYMATCHING2_SEARCH_DETECTOR_NODE_H

#include "pymatching/sparse_blossom/driver/implied_weights.h"
#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

#include "pymatching/sparse_blossom/config_parallel.h"

#ifdef ENABLE_SHOT_BUFFERS
#include <array>
#endif

namespace pm {

const uint8_t FLIPPED = 1;
const uint8_t WEIGHT_SIGN = 2;

struct SearchDetectorNodeEphemeralFields {
    /// The SearchDetectorNode that this node was reached from in the Dijkstra search
    SearchDetectorNode *reached_from_source = nullptr;

    /// `index_of_predecessor` is the index in `neighbors` of the neighboring detector node that this node was
    /// reached from in the Dijkstra search.
    size_t index_of_predecessor = SIZE_MAX;

    /// The distance from the search source detector node that this node was reached from.
    pm::cumulative_time_int distance_from_source = 0;

    /// Manages the next "look at me" event for the node
    QueuedEventTracker node_event_tracker;

    void reset() {
        reached_from_source = nullptr;
        index_of_predecessor = SIZE_MAX;
        node_event_tracker.clear();
        // distance_from_source intentionally left as-is -- always overwritten before use.
    }
};

class SearchDetectorNode {
   public:
    SearchDetectorNode() = default;

    /// == Ephemeral fields used to track algorithmic state during a search. ==
    /// Never shared cross-PE (unlike DetectorNodeEphemeralFields), so no USE_SHMEM pointer variant is
    /// needed here -- only ENABLE_SHOT_BUFFERS (array of NUM_BUFFERS_PER_UNIT slots) vs. a single slot.
#ifdef ENABLE_SHOT_BUFFERS
    std::array<SearchDetectorNodeEphemeralFields, NUM_BUFFERS_PER_UNIT> ephemeral_fields{};
    inline SearchDetectorNodeEphemeralFields& state(int rotating_buffer_idx) {
        return ephemeral_fields[rotating_buffer_idx];
    }
    inline const SearchDetectorNodeEphemeralFields& state(int rotating_buffer_idx) const {
        return ephemeral_fields[rotating_buffer_idx];
    }
#else
    SearchDetectorNodeEphemeralFields ephemeral_fields;
    inline SearchDetectorNodeEphemeralFields& state(int rotating_buffer_idx) {
        return ephemeral_fields;
    }
    inline const SearchDetectorNodeEphemeralFields& state(int rotating_buffer_idx) const {
        return ephemeral_fields;
    }
#endif // ENABLE_SHOT_BUFFERS

    /// == Permanent fields used to define the structure of the graph. ==
    std::vector<SearchDetectorNode *> neighbors;  /// The node's neighbors.
    std::vector<weight_int> neighbor_weights;     /// Distance crossed by the edge to each neighbor.
    std::vector<std::vector<size_t>>
        neighbor_observable_indices;        /// Indices of observables crossed by the edge to each neighbor.
    std::vector<uint8_t> neighbor_markers;  /// Used to mark edges as "seen" when decoding to an edge list.

    std::vector<std::vector<ImpliedWeight>> neighbor_implied_weights;

    int vb = -1;  // if node is cross partition, virtual boundary it belongs to

    size_t index_of_neighbor(SearchDetectorNode *target) const;

    void reset(int rotating_buffer_idx);
};

}  // namespace pm

#endif  // PYMATCHING2_SEARCH_DETECTOR_NODE_H
