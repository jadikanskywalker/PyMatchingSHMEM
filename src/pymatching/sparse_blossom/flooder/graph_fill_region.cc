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

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"

using namespace pm;

GraphFillRegionExtractedFields::GraphFillRegionExtractedFields(GraphFillRegion *owning_region)
    : blossom_parent(nullptr), blossom_parent_top(owning_region), radius((0 << 2) + 1) {
#ifdef USE_SHMEM
    owning_region->extracted_fields = this;
#endif
}

void GraphFillRegionExtractedFields::reset(GraphFillRegion *owning_region) {
    blossom_parent = nullptr;
    blossom_parent_top = owning_region;
    radius = pm::VaryingCT((0 << 2) + 1);
    match.clear();
}

GraphFillRegion::GraphFillRegion()
    :
#ifndef USE_SHMEM
      extracted_fields(this),
#endif
      alt_tree_node(nullptr),
      shrink_event_tracker(),
      owner_arena(nullptr),
      rotating_buffer_idx(-1)
{
    // Deliberately not in the initializer list above -- see constructed/destructed's own comment in
    // graph_fill_region.h.
    constructed = 1;
    destructed = 0;
}

void GraphFillRegion::reset() {
    shell_area.clear();
    blossom_children.clear();
    shrink_event_tracker.clear();
    alt_tree_node = nullptr;
    state().reset(this);
    rotating_buffer_idx = -1;
}

#ifndef USE_SHMEM
GraphFillRegion::GraphFillRegion(GraphFillRegion &&other)
    : extracted_fields(this),
      alt_tree_node(std::move(other.alt_tree_node)),
      shrink_event_tracker(std::move(other.shrink_event_tracker)),
      blossom_children(std::move(other.blossom_children)),
      shell_area(std::move(other.shell_area)),
      owner_arena(other.owner_arena),
      rotating_buffer_idx(other.rotating_buffer_idx)
{
    state().blossom_parent = other.state().blossom_parent;
    state().blossom_parent_top = other.state().blossom_parent_top == &other ? this : other.state().blossom_parent_top;
    state().radius = std::move(other.state().radius);
    state().match = std::move(other.state().match);
    constructed = 1;
    destructed = 0;
}
#endif

bool GraphFillRegion::tree_equal(const GraphFillRegion &other) const {
    if (alt_tree_node != other.alt_tree_node || state().radius != other.state().radius ||
        blossom_children.size() != other.blossom_children.size() || shell_area != other.shell_area) {
        return false;
    }
    if (blossom_children.empty())
        return true;
    for (size_t i = 0; i < blossom_children.size(); i++) {
        if (blossom_children[i].edge != other.blossom_children[i].edge)
            return false;
        if (!blossom_children[i].region->tree_equal(*other.blossom_children[i].region))
            return false;
    }
    return true;
}

bool GraphFillRegion::operator==(const GraphFillRegion &rhs) const {
    return tree_equal(rhs);
}

bool GraphFillRegion::operator!=(const GraphFillRegion &rhs) const {
    return !(rhs == *this);
}

void GraphFillRegion::add_match(GraphFillRegion *region, const CompressedEdge &edge) {
    state().match = Match{region, edge};
    region->state().match = Match{this, edge.reversed()};
}

void GraphFillRegion::cleanup_shell_area() {
    for (auto &detector_node : shell_area) {
        detector_node->reset(rotating_buffer_idx);
    }
}

void GraphFillRegion::clear_blossom_parent() {
    state().blossom_parent = nullptr;
    do_op_for_each_descendant_and_self([&](GraphFillRegion *descendant) {
        descendant->state().blossom_parent_top = this;
        for (DetectorNode *n : descendant->shell_area) {
            auto &node_state = n->state(rotating_buffer_idx);
            node_state.region_that_arrived_top = this;
            node_state.wrapped_radius_cached = n->compute_wrapped_radius(rotating_buffer_idx);
        }
    });
}

void GraphFillRegion::clear_blossom_parent_ignoring_wrapped_radius() {
    state().blossom_parent = nullptr;
    do_op_for_each_descendant_and_self([&](GraphFillRegion *descendant) {
        descendant->state().blossom_parent_top = this;
        for (DetectorNode *n : descendant->shell_area) {
            n->state(rotating_buffer_idx).region_that_arrived_top = this;
        }
    });
}

void GraphFillRegion::wrap_into_blossom(GraphFillRegion *new_blossom_parent_and_top) {
    state().blossom_parent = new_blossom_parent_and_top;
    do_op_for_each_descendant_and_self([&](GraphFillRegion *descendant) {
        descendant->state().blossom_parent_top = new_blossom_parent_and_top;
        for (DetectorNode *n : descendant->shell_area) {
            auto &node_state = n->state(rotating_buffer_idx);
            node_state.region_that_arrived_top = new_blossom_parent_and_top;
            node_state.wrapped_radius_cached = n->compute_wrapped_radius(rotating_buffer_idx);
        }
    });
}

bool GraphFillRegion::operator<=(const GraphFillRegion &rhs) const {
    const GraphFillRegion *r = this;
    while (r != nullptr && r != &rhs) {
        r = r->state().blossom_parent;
    }
    return r == &rhs;
}

bool GraphFillRegion::operator<(const GraphFillRegion &rhs) const {
    return this != &rhs && (*this <= rhs);
}

bool GraphFillRegion::operator>(const GraphFillRegion &rhs) const {
    return rhs < *this;
}

bool GraphFillRegion::operator>=(const GraphFillRegion &rhs) const {
    return rhs <= *this;
}
