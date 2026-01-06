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

#include "pymatching/sparse_blossom/flooder/detector_node.h"

#include <optional>

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

namespace pm {

int32_t DetectorNode::compute_wrapped_radius(
#ifdef USE_THREADS
    int solver_idx
) const {
    const auto& s = state(solver_idx);
    if (s.reached_from_source == nullptr) {
        return 0;
    }
    int32_t total = 0;
    auto r = s.region_that_arrived;
    while (r != s.region_that_arrived_top) {
        total += r->radius.y_intercept();
        r = r->blossom_parent;
    }
    return total - s.radius_of_arrival;
}
#else
) const {
    if (reached_from_source == nullptr) {
        return 0;
    }
    int32_t total = 0;
    auto r = region_that_arrived;
    while (r != region_that_arrived_top) {
        total += r->radius.y_intercept();
        r = r->blossom_parent;
    }
    return total - radius_of_arrival;
}
#endif

void DetectorNode::reset(
#ifdef USE_THREADS
    int solver_idx
) {
    state(solver_idx).reset();
}
#else
) {
    observables_crossed_from_source = 0;
    reached_from_source = nullptr;
    radius_of_arrival = 0;
    region_that_arrived = nullptr;
    region_that_arrived_top = nullptr;
    wrapped_radius_cached = 0;
    node_event_tracker.clear();
}
#endif

size_t DetectorNode::index_of_neighbor(DetectorNode *target) const {
    for (size_t k = 0; k < neighbors.size(); k++) {
        if (neighbors[k] == target) {
            return k;
        }
    }
    throw std::invalid_argument("Failed to find neighbor.");
}

GraphFillRegion *DetectorNode::heir_region_on_shatter(
#ifdef USE_THREADS
    int solver_idx
#endif
) const {
#ifdef USE_THREADS
    GraphFillRegion *r = state(solver_idx).region_that_arrived;
    GraphFillRegion *top = state(solver_idx).region_that_arrived_top;
#else
    GraphFillRegion *r = region_that_arrived;
    GraphFillRegion *top = region_that_arrived_top;
#endif
    while (true) {
        GraphFillRegion *p = r->blossom_parent;
        if (p == top) {
            return r;
        }
        r = p;
    }
}

cumulative_time_int DetectorNode::compute_local_radius_at_time_bounded_by_region(
    cumulative_time_int time,
    const GraphFillRegion &bounding_region
#ifdef USE_THREADS
    , int solver_idx
#endif
) const {
#ifdef USE_THREADS
    const auto& s = state(solver_idx);
    if (s.region_that_arrived == nullptr) {
        return 0;
    }
    if (*s.region_that_arrived > bounding_region) {
        return 0;
    }

    const GraphFillRegion *container = s.region_that_arrived;
    cumulative_time_int container_radius = 0;
    while (true) {
        if (container == nullptr) {
            break;
        }
        if (*container > bounding_region) {
            break;
        }
        container_radius += container->radius.get_distance_at_time(time);
        if (*container == bounding_region) {
            break;
        }
        container = container->blossom_parent;
    }
    return container_radius - s.radius_of_arrival;
#else
    if (region_that_arrived == nullptr) {
        // Nodes not in an region have zero local radius.
        return 0;
    }
    if (*region_that_arrived > bounding_region) {
        // If the target region is a descendant of the region that reached this node, then
        // bounding to that region effectively means this node has not yet been reached. Act
        // like a node not in a region.
        return 0;
    }

    const GraphFillRegion *container = region_that_arrived;
    cumulative_time_int container_radius = 0;
    while (true) {
        if (container == nullptr) {
            // Wasn't inside the bounding region.
            break;
        }
        if (*container > bounding_region) {
            // Was a cousin of some sort of the bounding region, and have just reached the
            // common ancestor.
            break;
        }
        container_radius += container->radius.get_distance_at_time(time);
        if (*container == bounding_region) {
            // Don't go beyond the limits of the bounding region.
            break;
        }
        container = container->blossom_parent;
    }
    return container_radius - radius_of_arrival;
#endif
}

std::optional<float> DetectorNode::compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(
    cumulative_time_int time,
    const GraphFillRegion &bounding_region,
    size_t neighbor_index
#ifdef USE_THREADS
    , int solver_idx
#endif
) const {
    DetectorNode *neighbor = neighbors[neighbor_index];
    cumulative_time_int max_w = neighbor_weights[neighbor_index];
        auto r1 = compute_local_radius_at_time_bounded_by_region(
        time,
        bounding_region
    #ifdef USE_THREADS
        , solver_idx
    #endif
        );
    if (neighbor == nullptr) {
        return (weight_int)std::min(max_w, r1);
    }

        auto r2 = neighbor->compute_local_radius_at_time_bounded_by_region(
        time,
        bounding_region
    #ifdef USE_THREADS
        , solver_idx
    #endif
        );

    // If the nodes at either side of the edge have regions that aren't linked according to the
    // state the mwpm, then the transition must be happening exactly at the local radius.
    if (r1 + r2 < max_w
#ifdef USE_THREADS
        || neighbor->state(solver_idx).region_that_arrived_top != state(solver_idx).region_that_arrived_top
#else
        || neighbor->region_that_arrived_top != region_that_arrived_top
#endif
    ) {
        return (weight_int)r1;
    }
    if (r1 == max_w
#ifdef USE_THREADS
        && *neighbor->state(solver_idx).region_that_arrived > *state(solver_idx).region_that_arrived
#else
        && *neighbor->region_that_arrived > *region_that_arrived
#endif
    ) {
        return max_w;
    }

    // Now we are in the complicated case, where for example we are trying to determine exactly
    // where two child blossom collided along an edge far in the past.

    // There is an additional complication in that *actually* we are focusing on a specific source
    // node, and we want to see the boundaries between the parts of the graph reached by one source
    // node and another source node, even if those nodes are part of the same region.

    // If the edge is between the same two sources, there is no stitch.
#ifdef USE_THREADS
    if (state(solver_idx).reached_from_source == neighbor->state(solver_idx).reached_from_source) {
#else
    if (reached_from_source == neighbor->reached_from_source) {
#endif
        return {};
    }

    // Back up to the time when they collided. Because they must grow in tandem past that point,
    // they each grew half of the extra width.
    float growth_past_collision = r1 + r2 - max_w;
    float growth_past_collision_of_r1 = growth_past_collision / 2;
    return r1 - growth_past_collision_of_r1;
}

}  // namespace pm
