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

#ifndef PYMATCHING_FILL_MATCH_DETECTOR_NODE_H
#define PYMATCHING_FILL_MATCH_DETECTOR_NODE_H

#include <optional>
#include <vector>

#include "pymatching/sparse_blossom/driver/implied_weights.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

#include "pymatching/sparse_blossom/config_parallel.h"

#ifdef ENABLE_SHOT_BUFFERS
#include <array>
#endif

namespace pm {

struct DetectorNodeEphemeralFields {
    GraphFillRegion* region_that_arrived = nullptr;
    GraphFillRegion* region_that_arrived_top = nullptr;
    int32_t wrapped_radius_cached = 0;
    DetectorNode* reached_from_source = nullptr;
    obs_int observables_crossed_from_source = 0;
    cumulative_time_int radius_of_arrival = 0;
    QueuedEventTracker node_event_tracker;

    void reset() {
        region_that_arrived = nullptr;
        region_that_arrived_top = nullptr;
        wrapped_radius_cached = 0;
        reached_from_source = nullptr;
        observables_crossed_from_source = 0;
        radius_of_arrival = 0;
        node_event_tracker.clear();
    }
};

/// A detector node is a location where a detection event might occur.
///
/// It corresponds to a potential symptom that could be seen, and can
/// be used to infer which errors have occurred.
///
/// There is a 1:1 correspondence between detector nodes in a matching
/// graph and the DETECTOR annotations in a Stim circuit.
class DetectorNode {
   public:
    // when using shmem, each neighbors array is stored in shared memory
    DetectorNode() = default;

    /// == Ephemeral fields used to track algorithmic state during matching. ==
#ifdef ENABLE_SHOT_BUFFERS
#ifdef USE_SHMEM
    /// Pointers to shared-memory epmeheral_field structs
    std::array<DetectorNodeEphemeralFields*, NUM_BUFFERS_PER_UNIT> ephemeral_fields{};
    inline DetectorNodeEphemeralFields& state(int rotating_buffer_idx) {
        return *(ephemeral_fields[rotating_buffer_idx]);
    }
    inline const DetectorNodeEphemeralFields& state(int rotating_buffer_idx) const {
        return *(ephemeral_fields[rotating_buffer_idx]);
    }
#else
    /// In-place array of NUM_BUFFERS_PER_UNIT epmeheral_field structs
    std::array<DetectorNodeEphemeralFields, NUM_BUFFERS_PER_UNIT> ephemeral_fields{};
    inline DetectorNodeEphemeralFields& state(int rotating_buffer_idx) {
        return ephemeral_fields[rotating_buffer_idx];
    }
    inline const DetectorNodeEphemeralFields& state(int rotating_buffer_idx) const {
        return ephemeral_fields[rotating_buffer_idx];
    }
#endif // USE_SHMEM
#else
#ifdef USE_SHMEM
    DetectorNodeEphemeralFields* ephemeral_fields=nullptr;
    inline DetectorNodeEphemeralFields& state(int rotating_buffer_idx) {
        return *ephemeral_fields;
    }
    inline const DetectorNodeEphemeralFields& state(int rotating_buffer_idx) const {
        return *ephemeral_fields;
    }
#else
    DetectorNodeEphemeralFields ephemeral_fields;
    inline DetectorNodeEphemeralFields& state(int rotating_buffer_idx) {
        return ephemeral_fields;
    }
    inline const DetectorNodeEphemeralFields& state(int rotating_buffer_idx) const {
        return ephemeral_fields;
    }
#endif
#endif // ENABLE_SHOT_BUFFERS

    /// == Permanent fields used to define the structure of the graph. ==
    std::vector<DetectorNode*> neighbors;       /// The node's neighbors.
    std::vector<weight_int> neighbor_weights;   /// Distance crossed by the edge to each neighbor.
    std::vector<obs_int> neighbor_observables;  /// Observables crossed by the edge to each neighbor.

// ===============
    int vb = -1; // if node is cross partition, virtual boundary it belongs to
    int obs_patch_id = -1;
// ===============

    /// After it reached this node, how much further did the owning search region grow? Also is it currently growing?
    inline VaryingCT local_radius(int rotating_buffer_idx) const {
        const auto& s = state(rotating_buffer_idx);
        if (s.region_that_arrived_top == nullptr) {
            return VaryingCT{0};
        }
        return s.region_that_arrived_top->radius + s.wrapped_radius_cached;
    }

    /// Determines the region that owns this node which is a child of this node's top region.
    GraphFillRegion* heir_region_on_shatter(int rotating_buffer_idx) const;

    /// Check if this node is part the same top-level region as another.
    /// Note that they may have different lower level owners that still merge into the same top level owned.
    inline bool has_same_owner_as(
        const DetectorNode& other,
        int rotating_buffer_idx
    ) const {
        return state(rotating_buffer_idx).region_that_arrived_top ==
               other.state(rotating_buffer_idx).region_that_arrived_top;
    }

    /// Zero all ephemeral fields.
    /// Doesn't free anything or propagate a signal to other objects. Just zeros the fields.
    void reset(int rotating_buffer_idx);

    size_t index_of_neighbor(DetectorNode* neighbor) const;

    /// The 'wrapped radius' is the amount of region growth that has occurred past the point where
    /// this node was hit, up until the last time a blossom was formed around this node. It is the
    /// sum of the radius of (always frozen) not-top-level regions containing this node. It is a
    /// useful intermediate value for quickly computing the local radius of the node, because it
    /// accounts for everything except the (potentially varying) top level region.
    int32_t compute_wrapped_radius(int rotating_buffer_idx) const;

    /// It works out the local radius of the node, but bounds it to not go beyond the given bounding
    /// region. This is useful for understanding the internal structure of the blossoms.
    ///
    /// This is a utility method used to help with drawing the internal state of the algorithm.
    cumulative_time_int compute_local_radius_at_time_bounded_by_region(
        cumulative_time_int time,
        const GraphFillRegion& bounding_region,
        int rotating_buffer_idx
    ) const;

    /// The 'stitch radius' is the place where ownership of an edge stops.
    ///
    /// This is a utility method used to help with drawing the internal state of the algorithm.
    ///
    /// Args:
    ///     time: In case the result is varying, this specifies exactly which time to look at it.
    ///     bounding_region: This says which graph fill region is being "focused on". In particular,
    ///         this may be a frozen region within a blossom in which case the goal is to find the
    ///         exact place along edges where the region froze or where two source nodes were
    ///         bumping up against each other.
    ///     neighbor_index: The edge of interest.
    ///
    /// Returns:
    ///     Empty: This is a fully internal edge. There is no region transition along it.
    ///     An integer: This is the distance, starting from *this node and running along the edge,
    ///         to the transition point. This may be equal to 0, or equal to the full length of the
    ///         edge, or somewhere in between.
    std::optional<float> compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(
        cumulative_time_int time,
        const GraphFillRegion& bounding_region,
        size_t neighbor_index,
        int rotating_buffer_idx
    ) const;

    std::vector<std::vector<ImpliedWeight>> neighbor_implied_weights;
};

}  // namespace pm

#endif  // PYMATCHING_FILL_MATCH_DETECTOR_NODE_H
