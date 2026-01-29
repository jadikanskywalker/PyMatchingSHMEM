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

#include "pymatching/sparse_blossom/flooder/graph_flooder.h"

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"

#ifdef USE_THREADS
#include <memory>
#include <omp.h>

#include "../config_parallel.h"
#endif

using namespace pm;

#ifdef USE_THREADS
// ===============
GraphFlooder::GraphFlooder()
    : graph_ptr(std::make_shared<MatchingGraph>()),
      graph(*graph_ptr),
      negative_weight_obs_mask(0),
      negative_weight_sum(0) {
}

#ifdef USE_SHMEM
GraphFlooder::GraphFlooder(
    std::shared_ptr<MatchingGraph> graph, int solver_set_idx, GraphFillRegion* shmem_buffer, size_t shmem_buffer_size)
    : graph_ptr(graph),
      graph(*graph_ptr),
      negative_weight_obs_mask(0),
      negative_weight_sum(0),
      rotating_buffer_idx(solver_set_idx),
      region_arena(shmem_buffer, shmem_buffer_size) {
}
#else
GraphFlooder::GraphFlooder(std::shared_ptr<MatchingGraph> graph, int solver_set_idx)
    : graph_ptr(graph),
      graph(*graph_ptr),
      negative_weight_obs_mask(0),
      negative_weight_sum(0),
      rotating_buffer_idx(solver_set_idx) {
}
#endif

GraphFlooder::GraphFlooder(MatchingGraph graph_val)
    : graph_ptr(std::make_shared<MatchingGraph>(std::move(graph_val))),
      graph(*graph_ptr),
      negative_weight_obs_mask(0),
      negative_weight_sum(0) {
}

GraphFlooder::GraphFlooder(GraphFlooder&& flooder) noexcept
    : graph_ptr(std::move(flooder.graph_ptr)),
      graph(*graph_ptr),
      queue(std::move(flooder.queue)),
      region_arena(std::move(flooder.region_arena)),
      match_edges(std::move(flooder.match_edges)),
      negative_weight_detection_events(std::move(flooder.negative_weight_detection_events)),
      negative_weight_observables(std::move(flooder.negative_weight_observables)),
      negative_weight_obs_mask(flooder.negative_weight_obs_mask),
      negative_weight_sum(flooder.negative_weight_sum),
      rotating_buffer_idx(flooder.rotating_buffer_idx) {
}
// ===============
#else
GraphFlooder::GraphFlooder() : negative_weight_obs_mask(0), negative_weight_sum(0) {
}

GraphFlooder::GraphFlooder(MatchingGraph graph)
    : graph(std::move(graph)), negative_weight_obs_mask(0), negative_weight_sum(0) {
}

GraphFlooder::GraphFlooder(GraphFlooder&& flooder) noexcept
    : graph(std::move(flooder.graph)),
      queue(std::move(flooder.queue)),
      region_arena(std::move(flooder.region_arena)),
      match_edges(std::move(flooder.match_edges)),
      negative_weight_detection_events(std::move(flooder.negative_weight_detection_events)),
      negative_weight_observables(std::move(flooder.negative_weight_observables)),
      negative_weight_obs_mask(flooder.negative_weight_obs_mask),
      negative_weight_sum(flooder.negative_weight_sum) {
}
#endif

// #ifdef USE_THREADS
// inline void debug_validate_region_nodes(const GraphFillRegion& r, const std::vector<DetectorNode>& nodes,
// std::set<long> active_parts) {
//     for (DetectorNode *p : r.shell_area) {
//         if (p == nullptr) continue;
//         size_t i = static_cast<size_t>(p - &nodes[0]);
//         if (i >= nodes.size())
//             std::cout << "  ERROR: shell_area nodes out of bounds" << std::endl;
//         if (&nodes[i] != p)
//             std::cout << "  ERROR: could not reference shell_area node" << std::endl;
//         if (!(p->is_active == omp_get_thread_num()))
//             std::cout << "  ERROR: inactive node added to shell_area" << std::endl;
//     }
// }
// #endif

#ifdef USE_THREADS
// ===============
inline bool GraphFlooder::is_active(const DetectorNode* node) const {
    if (node->vb < 0)
        return true;
    return (node->vb > vb_left && node->vb < vb_right);  // DEPENDENT ON ROUND PARTITIONING
}
// ===============
#endif

void GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion& region, DetectorNode& detector_node) {
#ifdef USE_THREADS
    auto& node = detector_node.state(rotating_buffer_idx);
    node.reached_from_source = &detector_node;
    node.radius_of_arrival = 0;
    node.region_that_arrived = &region;
    node.region_that_arrived_top = &region;
    node.wrapped_radius_cached = 0;
#else
    detector_node.reached_from_source = &detector_node;
    detector_node.radius_of_arrival = 0;
    detector_node.region_that_arrived = &region;
    detector_node.region_that_arrived_top = &region;
    detector_node.wrapped_radius_cached = 0;
#endif
    region.shell_area.push_back(&detector_node);
    // #ifdef USE_THREADS
    //     if (DEBUG) {
    //         // std::cout << "    region created at detector node " << &detector_node << std::endl;
    //         debug_validate_region_nodes(region, graph.nodes, active_partitions);
    //     }
    // #endif
    reschedule_events_at_detector_node(detector_node);
}

#ifdef USE_THREADS
std::pair<size_t, cumulative_time_int> GraphFlooder::find_next_event_at_node_not_occupied_by_growing_top_region(
    const DetectorNode& detector_node, VaryingCT rad1) const {
#else
std::pair<size_t, cumulative_time_int> find_next_event_at_node_not_occupied_by_growing_top_region(
    const DetectorNode& detector_node, VaryingCT rad1) {
#endif
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;

    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr)
        start++;

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        auto weight = detector_node.neighbor_weights[i];

        auto neighbor = detector_node.neighbors[i];

#ifdef USE_THREADS
        // ===============
        if (!is_active(neighbor)) {  // skip inactive neighbors
            if (!(neighbor->vb >= 0))
                std::cout << "    NOTE2: neighbor inactive not cross partition (not growing)" << std::endl;
            continue;
        }
// ===============
#endif

#ifdef USE_THREADS
        auto rad2 = neighbor->local_radius(rotating_buffer_idx);
#else
        auto rad2 = neighbor->local_radius();
#endif

        if (rad2.is_growing()) {
            auto collision_time = weight - rad1.y_intercept() - rad2.y_intercept();
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = i;
            }
        }
    }
    return {best_neighbor, best_time};
}

#ifdef USE_THREADS
std::pair<size_t, cumulative_time_int> GraphFlooder::find_next_event_at_node_occupied_by_growing_top_region(
    const DetectorNode& detector_node, const VaryingCT& rad1) const {
#else
std::pair<size_t, cumulative_time_int> find_next_event_at_node_occupied_by_growing_top_region(
    const DetectorNode& detector_node, const VaryingCT& rad1) {
#endif
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;
    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr) {
        // Growing towards boundary
        auto weight = detector_node.neighbor_weights[0];
        auto collision_time = weight - rad1.y_intercept();
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = 0;
        }
        start++;
    }

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        auto weight = detector_node.neighbor_weights[i];

        auto neighbor = detector_node.neighbors[i];

#ifdef USE_THREADS
        // ===============
        // Treat virtual nodes like boundary
        if (!is_active(neighbor)) {
            if (!(neighbor->vb >= 0)) {
                std::cout
                    << "    NOTE3: neighbor inactive not cross partition (growing)" << std::endl
                    << "      node: " << &detector_node << "    vb: " << detector_node.vb
                    << "      vb_left: " << vb_left << "    vb_right: "
                    << vb_right
                    //   << "  shot_marker: " << detector_node.shot_marker << "  flooder.current_shot: " << current_shot
                    << std::endl;
            }
            auto collision_time = weight - rad1.y_intercept();
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = i;
            }
        }
        // Treat active nodes like normal
// ===============
#endif

        if (
#ifdef USE_THREADS
            detector_node.has_same_owner_as(*neighbor, rotating_buffer_idx)
#else
            detector_node.has_same_owner_as(*neighbor)
#endif
        ) {
            continue;
        }
#ifdef USE_THREADS
        auto rad2 = neighbor->local_radius(rotating_buffer_idx);
#else
        auto rad2 = neighbor->local_radius();
#endif
        if (rad2.is_shrinking()) {
            continue;
        }

        auto collision_time = weight - rad1.y_intercept() - rad2.y_intercept();
        if (rad2.is_growing()) {
            collision_time >>= 1;
        }
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = i;
        }
    }
    return {best_neighbor, best_time};
}

std::pair<size_t, cumulative_time_int> GraphFlooder::find_next_event_at_node_returning_neighbor_index_and_time(
    const DetectorNode& detector_node) const {
#ifdef USE_THREADS
    auto rad1 = detector_node.local_radius(rotating_buffer_idx);
#else
    auto rad1 = detector_node.local_radius();
#endif

    if (rad1.is_growing()) {
        return find_next_event_at_node_occupied_by_growing_top_region(detector_node, rad1);
    } else {
        return find_next_event_at_node_not_occupied_by_growing_top_region(detector_node, rad1);
    }
}

void GraphFlooder::reschedule_events_at_detector_node(DetectorNode& detector_node) {
#ifdef USE_THREADS
    // ===============
    if (!is_active(&detector_node)) {
        std::cout << "    NOTE: reschedule called on inactive node" << std::endl;
        return;
    }
// ===============
#endif
    auto x = find_next_event_at_node_returning_neighbor_index_and_time(detector_node);
    if (x.first == SIZE_MAX) {
#ifdef USE_THREADS
        detector_node.state(rotating_buffer_idx).node_event_tracker.set_no_desired_event();
#else
        detector_node.node_event_tracker.set_no_desired_event();
#endif
    } else {
#ifdef USE_THREADS
        detector_node.state(rotating_buffer_idx)
            .node_event_tracker.set_desired_event(
#else
        detector_node.node_event_tracker.set_desired_event(
#endif
                {
                    &detector_node,
                    cyclic_time_int{x.second},
                },
                queue);
    }
}

void GraphFlooder::schedule_tentative_shrink_event(GraphFillRegion& region) {
    cumulative_time_int t;
    if (region.shell_area.empty()) {
        t = region.radius.time_of_x_intercept_for_shrinking();
    } else {
#ifdef USE_THREADS
        t = region.shell_area.back()->local_radius(rotating_buffer_idx).time_of_x_intercept_for_shrinking();
#else
        t = region.shell_area.back()->local_radius().time_of_x_intercept_for_shrinking();
#endif
    }
    region.shrink_event_tracker.set_desired_event(
        {
            &region,
            cyclic_time_int{t},
        },
        queue);
}

void GraphFlooder::do_region_arriving_at_empty_detector_node(
    GraphFillRegion& region, DetectorNode& empty_node, const DetectorNode& from_node, size_t from_to_empty_index) {
#ifdef USE_THREADS
    // ===============
    if (!is_active(&empty_node))
        std::cout << "ERROR: do_region_arriving_at_empty_detector_node -> empty_node is inactive" << std::endl;
// ===============
#endif
#ifdef USE_THREADS
    const auto& from_state = from_node.state(rotating_buffer_idx);
    auto& empty_state = empty_node.state(rotating_buffer_idx);
    empty_state.observables_crossed_from_source =
        (from_state.observables_crossed_from_source ^ from_node.neighbor_observables[from_to_empty_index]);
    empty_state.reached_from_source = from_state.reached_from_source;
    empty_state.radius_of_arrival = region.radius.get_distance_at_time(queue.cur_time);
    empty_state.region_that_arrived = &region;
    empty_state.region_that_arrived_top = region.blossom_parent_top;
    empty_state.wrapped_radius_cached = empty_node.compute_wrapped_radius(rotating_buffer_idx);
#else
    empty_node.observables_crossed_from_source =
        (from_node.observables_crossed_from_source ^ from_node.neighbor_observables[from_to_empty_index]);
    empty_node.reached_from_source = from_node.reached_from_source;
    empty_node.radius_of_arrival = region.radius.get_distance_at_time(queue.cur_time);
    empty_node.region_that_arrived = &region;
    empty_node.region_that_arrived_top = region.blossom_parent_top;
    empty_node.wrapped_radius_cached = empty_node.compute_wrapped_radius();
#endif
    region.shell_area.push_back(&empty_node);
    // #ifdef USE_THREADS
    //     if (DEBUG) debug_validate_region_nodes(region, graph.nodes, active_partitions);
    // #endif
    reschedule_events_at_detector_node(empty_node);
}

MwpmEvent GraphFlooder::do_region_shrinking(GraphFillRegion& region) {
    if (region.shell_area.empty()) {
        return do_blossom_shattering(region);
    } else if (region.shell_area.size() == 1 && region.blossom_children.empty()) {
#ifdef USE_THREADS
        // ===============
        if (!is_active(*region.shell_area.begin()))
            std::cout << "ERROR: do_region_shrinking -> do_generate_implosion on inactive shell_area[0]" << std::endl;
// ===============
#endif
        return do_degenerate_implosion(region);
    } else {
        auto leaving_node = region.shell_area.back();
#ifdef USE_THREADS
        // ===============
        if (!is_active(leaving_node))
            std::cout << "ERROR: do_region_shrinking -> leaving inactive node" << std::endl;
// ===============
#endif
        region.shell_area.pop_back();
#ifdef USE_THREADS
        auto& leave_state = leaving_node->state(rotating_buffer_idx);
        leave_state.region_that_arrived = nullptr;
        leave_state.region_that_arrived_top = nullptr;
        leave_state.wrapped_radius_cached = 0;
        leave_state.reached_from_source = nullptr;
        leave_state.radius_of_arrival = 0;
        leave_state.observables_crossed_from_source = 0;
#else
        leaving_node->region_that_arrived = nullptr;
        leaving_node->region_that_arrived_top = nullptr;
        leaving_node->wrapped_radius_cached = 0;
        leaving_node->reached_from_source = nullptr;
        leaving_node->radius_of_arrival = 0;
        leaving_node->observables_crossed_from_source = 0;
#endif
        reschedule_events_at_detector_node(*leaving_node);
        schedule_tentative_shrink_event(region);
        return MwpmEvent::no_event();
    }
}

MwpmEvent GraphFlooder::do_neighbor_interaction(DetectorNode& src, size_t src_to_dst_index, DetectorNode& dst) {
#ifdef USE_THREADS
    auto& src_state = src.state(rotating_buffer_idx);
    auto& dst_state = dst.state(rotating_buffer_idx);
#endif
    // First check if one region is moving into an empty location
    if (
#ifdef USE_THREADS
        src_state.region_that_arrived && !dst_state.region_that_arrived
#else
        src.region_that_arrived && !dst.region_that_arrived
#endif
    ) {
#ifdef USE_THREADS
        do_region_arriving_at_empty_detector_node(*src_state.region_that_arrived_top, dst, src, src_to_dst_index);
#else
        do_region_arriving_at_empty_detector_node(*src.region_that_arrived_top, dst, src, src_to_dst_index);
#endif
        return MwpmEvent::no_event();
    } else if (
#ifdef USE_THREADS
        dst_state.region_that_arrived && !src_state.region_that_arrived
#else
        dst.region_that_arrived && !src.region_that_arrived
#endif
    ) {
#ifdef USE_THREADS
        do_region_arriving_at_empty_detector_node(
            *dst_state.region_that_arrived_top, src, dst, dst.index_of_neighbor(&src));
#else
        do_region_arriving_at_empty_detector_node(*dst.region_that_arrived_top, src, dst, dst.index_of_neighbor(&src));
#endif
        return MwpmEvent::no_event();
    } else {
        // Two regions colliding
        return RegionHitRegionEventData{
#ifdef USE_THREADS
            src_state.region_that_arrived_top,
            dst_state.region_that_arrived_top,
#else
            src.region_that_arrived_top,
            dst.region_that_arrived_top,
#endif
            CompressedEdge{
#ifdef USE_THREADS
                src_state.reached_from_source,
                dst_state.reached_from_source,
                src_state.observables_crossed_from_source ^ dst_state.observables_crossed_from_source ^
#else
                src.reached_from_source,
                dst.reached_from_source,
                src.observables_crossed_from_source ^ dst.observables_crossed_from_source ^
#endif
                    src.neighbor_observables[src_to_dst_index]},
        };
    }
}

MwpmEvent GraphFlooder::do_region_hit_boundary_interaction(DetectorNode& node) {
    // Drop stale events that fired after shrinking/reset.
    if (
#ifdef USE_THREADS
        node.state(rotating_buffer_idx).reached_from_source == nullptr ||
        node.state(rotating_buffer_idx).region_that_arrived_top == nullptr
#else
        node.reached_from_source == nullptr || node.region_that_arrived_top == nullptr
#endif
    ) {
#ifdef USE_THREADS
        // ===============
        if (DEBUG) {
            std::cout << "  DEBUG: drop stale boundary event at node " << &node
                      << " rfs=" << node.state(rotating_buffer_idx).reached_from_source
                      << " top=" << node.state(rotating_buffer_idx).region_that_arrived_top << std::endl;
        }
// ===============
#endif
        return MwpmEvent::no_event();
    }
    return RegionHitBoundaryEventData{
#ifdef USE_THREADS
        node.state(rotating_buffer_idx).region_that_arrived_top,
        CompressedEdge{
            node.state(rotating_buffer_idx).reached_from_source,
            nullptr,
            node.state(rotating_buffer_idx).observables_crossed_from_source ^ node.neighbor_observables[0]}
#else
        node.region_that_arrived_top,
        CompressedEdge{
            node.reached_from_source, nullptr, node.observables_crossed_from_source ^ node.neighbor_observables[0]}
#endif
    };
}

#ifdef USE_THREADS
// ===============
MwpmEvent GraphFlooder::do_region_hit_virtual_boundary_interaction(DetectorNode& node, size_t virtual_neighbor_index) {
    if (node.state(rotating_buffer_idx).reached_from_source == nullptr ||
        node.state(rotating_buffer_idx).region_that_arrived_top == nullptr) {
        if (DEBUG) {
            std::cout << "  DEBUG: drop stale virtual boundary event at node " << &node
                      << " rfs=" << node.state(rotating_buffer_idx).reached_from_source
                      << " top=" << node.state(rotating_buffer_idx).region_that_arrived_top << std::endl;
        }
        return MwpmEvent::no_event();
    }
    // if (DEBUG)
    //     std::cout << "  DEBUG: region hit virtual boundary" << std::endl
    //               << "    node: " << &node << std::endl
    //               << "    node.reached_from_source: " << node.state(rotating_buffer_idx).reached_from_source <<
    //               std::endl
    //               << "    node.region_that_arrived_top: " << node.state(rotating_buffer_idx).region_that_arrived <<
    //               std::endl
    //               << "    node.neighbors[virtual_neighbor_index]: " << node.neighbors[virtual_neighbor_index] <<
    //               std::endl;
    return RegionHitVirtualBoundaryEventData{
        node.state(rotating_buffer_idx).region_that_arrived_top,
        CompressedEdge{
            node.state(rotating_buffer_idx).reached_from_source,
            node.neighbors[virtual_neighbor_index],
            node.state(rotating_buffer_idx).observables_crossed_from_source ^
                node.neighbor_observables[virtual_neighbor_index]}};
}
// ===============
#endif

MwpmEvent GraphFlooder::do_degenerate_implosion(const GraphFillRegion& region) {
    return RegionHitRegionEventData{
        region.alt_tree_node->parent.alt_tree_node->outer_region,
        region.alt_tree_node->outer_region,
        CompressedEdge{
            region.alt_tree_node->parent.edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.obs_mask ^ region.alt_tree_node->parent.edge.obs_mask}};
}

MwpmEvent GraphFlooder::do_blossom_shattering(GraphFillRegion& region) {
    return BlossomShatterEventData{
        &region,
#ifdef USE_THREADS
        region.alt_tree_node->parent.edge.loc_from->heir_region_on_shatter(region.rotating_buffer_idx),
        region.alt_tree_node->inner_to_outer_edge.loc_from->heir_region_on_shatter(region.rotating_buffer_idx)
#else
        region.alt_tree_node->parent.edge.loc_from->heir_region_on_shatter(),
        region.alt_tree_node->inner_to_outer_edge.loc_from->heir_region_on_shatter()
#endif
    };
}

GraphFillRegion* GraphFlooder::create_blossom(std::vector<RegionEdge>& contained_regions) {
    auto blossom_region = region_arena.alloc_default_constructed();
#ifdef USE_THREADS
    blossom_region->owner_arena = &region_arena;
    blossom_region->rotating_buffer_idx = rotating_buffer_idx;
#endif
    blossom_region->radius = VaryingCT::growing_varying_with_zero_distance_at_time(queue.cur_time);
    blossom_region->blossom_children = std::move(contained_regions);
    for (auto& region_edge : blossom_region->blossom_children) {
        region_edge.region->radius = region_edge.region->radius.then_frozen_at_time(queue.cur_time);
        region_edge.region->wrap_into_blossom(blossom_region);
        region_edge.region->shrink_event_tracker.set_no_desired_event();
    }

    blossom_region->do_op_for_each_node_in_total_area([this](DetectorNode* n) {
        reschedule_events_at_detector_node(*n);
    });

    return blossom_region;
}

bool GraphFlooder::dequeue_decision(FloodCheckEvent ev) {
    switch (ev.tentative_event_type) {
        case FloodCheckEventType::LOOK_AT_NODE: {
            auto& node = *ev.data_look_at_node;
#ifdef USE_THREADS
            return node.state(rotating_buffer_idx).node_event_tracker.dequeue_decision(ev, queue);
#else
            return node.node_event_tracker.dequeue_decision(ev, queue);
#endif
        }
        case FloodCheckEventType::LOOK_AT_SHRINKING_REGION: {
            auto& region = *ev.data_look_at_shrinking_region;
            return region.shrink_event_tracker.dequeue_decision(ev, queue);
        }
        case FloodCheckEventType::NO_FLOOD_CHECK_EVENT:
            return true;
        default:
            throw std::invalid_argument("Unrecognized event type.");
    }
}

FloodCheckEvent GraphFlooder::dequeue_valid() {
    while (true) {
        FloodCheckEvent ev = queue.dequeue();
        if (dequeue_decision(ev)) {
            return ev;
        }
    }
}

void GraphFlooder::set_region_growing(GraphFillRegion& region) {
    region.radius = region.radius.then_growing_at_time(queue.cur_time);

    // No shrinking event can occur while growing.
    region.shrink_event_tracker.set_no_desired_event();

    // Node events can occur while growing, and events in the queue may occur sooner than
    // previously scheduled. Therefore, we must reschedule all the nodes.
    region.do_op_for_each_node_in_total_area([this](DetectorNode* n) {
        reschedule_events_at_detector_node(*n);
    });
}

void GraphFlooder::set_region_frozen(GraphFillRegion& region) {
    bool was_shrinking = region.radius.is_shrinking();
    region.radius = region.radius.then_frozen_at_time(queue.cur_time);

    // No shrinking event can occur while frozen.
    region.shrink_event_tracker.set_no_desired_event();

    // Node events can occur while frozen, from other regions growing into this one.
    // However, those events can only be sooner than the currently scheduled events
    // if the region was previously shrinking (as opposed to growing).
    if (was_shrinking) {
        region.do_op_for_each_node_in_total_area([&](DetectorNode* n) {
            reschedule_events_at_detector_node(*n);
        });
    }
}

void GraphFlooder::set_region_shrinking(GraphFillRegion& region) {
    region.radius = region.radius.then_shrinking_at_time(queue.cur_time);

    // Shrinking events can now occur.
    schedule_tentative_shrink_event(region);

    // No node events can occur while shrinking.
    region.do_op_for_each_node_in_total_area([&](DetectorNode* n) {
#ifdef USE_THREADS
        n->state(rotating_buffer_idx).node_event_tracker.set_no_desired_event();
#else
        n->node_event_tracker.set_no_desired_event();
#endif
    });
}

MwpmEvent GraphFlooder::do_look_at_node_event(DetectorNode& node) {
    auto next = find_next_event_at_node_returning_neighbor_index_and_time(node);
    if (next.second == queue.cur_time) {
        // Need to revisit this node immediately after the mwpm event is handled. There may be an event to handle
        // along another edge, or even along the same edge at a later time (e.g. if this event isn't the *first* event
        // for the neighbor along the current edge when the neighbor is being rescheduled).
#ifdef USE_THREADS
        node.state(rotating_buffer_idx)
            .node_event_tracker.set_desired_event(
#else
        node.node_event_tracker.set_desired_event(
#endif
                {
                    &node,
                    cyclic_time_int{queue.cur_time},
                },
                queue);

        if (node.neighbors[next.first] == nullptr) {
            return do_region_hit_boundary_interaction(node);
        }
#ifdef USE_THREADS
        // ===============
        else if (!is_active(node.neighbors[next.first])) {  // treat virtual nodes like boundary
            if (!(node.neighbors[next.first]->vb >= 0)) {
                std::cout << "ERROR: do_region_hit_virtual_boundary_interaction -> node not cross-partition"
                          << std::endl;
            }
            return do_region_hit_virtual_boundary_interaction(node, next.first);
        }
// ===============
#endif
        auto& neighbor = *node.neighbors[next.first];
        return do_neighbor_interaction(node, next.first, neighbor);
    } else if (next.first != SIZE_MAX) {
        // Need to revisit this node at a later time.
#ifdef USE_THREADS
        node.state(rotating_buffer_idx)
            .node_event_tracker.set_desired_event(
#else
        node.node_event_tracker.set_desired_event(
#endif
                {
                    &node,
                    cyclic_time_int{next.second},
                },
                queue);
    }

    return MwpmEvent::no_event();
}

MwpmEvent GraphFlooder::process_tentative_event_returning_mwpm_event(FloodCheckEvent tentative_event) {
    switch (tentative_event.tentative_event_type) {
        case LOOK_AT_NODE: {
#ifdef USE_THREADS
            // ===============
            if (!is_active(tentative_event.data_look_at_node))
                std::cout << "ERROR: look at node event for inactive node" << std::endl
                          << "  node: " << tentative_event.data_look_at_node << std::endl;
// ===============
#endif
            return do_look_at_node_event(*tentative_event.data_look_at_node);
        }
        case LOOK_AT_SHRINKING_REGION: {
#ifdef USE_THREADS
            // ===============
            if (tentative_event.data_look_at_shrinking_region->shell_area.size() > 0 &&
                !is_active(*tentative_event.data_look_at_shrinking_region->shell_area.begin()))
                std::cout << "ERROR: shrinking region event for inactive source node" << std::endl
                          << "  shell_area[0]: " << *tentative_event.data_look_at_shrinking_region->shell_area.begin()
                          << std::endl;
// ===============
#endif
            return do_region_shrinking(*tentative_event.data_look_at_shrinking_region);
        }
        default:
            throw std::invalid_argument("Unknown tentative event type.");
    }
}

MwpmEvent GraphFlooder::run_until_next_mwpm_notification() {
    while (true) {
        FloodCheckEvent tentative_event = dequeue_valid();
        if (tentative_event.tentative_event_type == NO_FLOOD_CHECK_EVENT) {
            return MwpmEvent::no_event();
        }
        MwpmEvent notification = process_tentative_event_returning_mwpm_event(tentative_event);
        if (notification.event_type != NO_EVENT) {
            return notification;
        }
    }
}

void GraphFlooder::sync_negative_weight_observables_and_detection_events() {
    /// Move set of negative weight detection events into a sorted vector, for faster processing during decoding
    negative_weight_detection_events.clear();
    negative_weight_detection_events.reserve(graph.negative_weight_detection_events_set.size());
    for (auto& det : graph.negative_weight_detection_events_set)
        negative_weight_detection_events.push_back(det);
    /// Move set of negative weight observable indices into a sorted vector, for faster processing during decoding
    negative_weight_observables.clear();
    negative_weight_observables.reserve(graph.negative_weight_observables_set.size());
    for (auto& obs : graph.negative_weight_observables_set)
        negative_weight_observables.push_back(obs);
    negative_weight_obs_mask = 0;
    if (graph.num_observables <= sizeof(pm::obs_int) * 8) {
        // Compute the observable mask for 64 or fewer observables
        for (auto& i : negative_weight_observables)
            negative_weight_obs_mask ^= (pm::obs_int)1 << i;
    }
    negative_weight_sum = graph.negative_weight_sum;
}

#ifdef USE_THREADS
// ===============
// Prepare the flooder to solve a single partition
// void GraphFlooder::prepare_for_solve_partition(int tid, long p) {
//     active_partitions.clear();
//     active_partitions.insert(p);
//     for (DetectorNode& node : graph.nodes)
//         if (node.partition == p && !node.is_cross_partition)
//             node.is_active = tid;
//     current_tid = tid;
// }

// // Prepare the flooder to fuse partitions p1 and p2
// void GraphFlooder::prepare_for_fuse_partitions(int tid, long p_without_virtuals, long p_with_virtuals) {
//     active_partitions.clear();
//     active_partitions.insert(p_without_virtuals);
//     active_partitions.insert(p_with_virtuals);
//     if (DEBUG)
//         std::cout << "  DEBUG: Thread " << tid << " solver preparing to fuse partitions "
//                   << p_without_virtuals << " and " << p_with_virtuals << std::endl;
//     for (DetectorNode& node : graph.nodes)
//         if ((node.partition == p_without_virtuals && !node.is_cross_partition) ||
//             (node.partition == p_with_virtuals))
//             node.is_active = tid;
//     current_tid = tid;
// }
// ===============
#endif
