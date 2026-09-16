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

#ifndef PYMATCHING2_GRAPH_FILL_REGION_H
#define PYMATCHING2_GRAPH_FILL_REGION_H

#include <vector>

#include "pymatching/sparse_blossom/flooder/match.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/region_edge.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

#ifdef USE_SHMEM
#include "pymatching/sparse_blossom/flooder/helpers/shmem_arena.h"
#endif

#include "pymatching/sparse_blossom/arena.h"


namespace pm {

class AltTreeNode;

struct GraphFillRegion {
    /// If this region has merged with others into a blossom, this is that blossom.
    GraphFillRegion* blossom_parent;
    /// The topmost fill region that contains this region. This field must be kept up to date as
    /// the region structure changes.
    GraphFillRegion* blossom_parent_top;
    /// How much this region has grown since its creation (and whether it is currently growing).
    /// For graph fill regions starting from detection events, this is just the actual radius
    /// of the region. For graph fill regions starting from blossom-creation events, this is how
    /// much the blossom has grown since it was created (it's the *extra* radius, not the total
    /// radius starting from the detection events).
    pm::VaryingCT radius;
    /// If the region is matched, as opposed to growing/shrinking, this says what it is matched to.
    pm::Match match;
    /// If this is a top-level region (not a blossom child), this is the alternating tree that
    /// it is part of. Note that it may be a degenerate alternating tree with just a single
    /// graph fill region (this one).
    pm::AltTreeNode* alt_tree_node;
    /// Event tracker for shrink events. Handles ensuring at least, and ideally exactly, one event
    /// to look at the region is in the event queue.
    QueuedEventTracker shrink_event_tracker;
    /// If this region is a blossom, these are its child regions along with the cyclic paths
    /// between the children.
    std::vector<pm::RegionEdge> blossom_children;
    /// The set of nodes directly owned by this region. Note that this vector does not include
    /// the nodes indirectly owned by this region that are owned by the blossom children of this
    /// region (or their children or etc).
    std::vector<pm::DetectorNode*> shell_area;

    void cleanup_shell_area();

#ifdef USE_SHMEM
    SHMEMArena<GraphFillRegion>* owner_arena;
#else
    Arena<GraphFillRegion>* owner_arena;
#endif

    int rotating_buffer_idx = 0;

    // Cumulative construct/destruct counts for whichever memory slot this object currently occupies.
    // GraphFillRegion()'s body sets constructed=1, destructed=0 on a genuine construction (a brand-new
    // slot, or a move-construction -- see below); Arena::del() increments `destructed` directly, and
    // Arena::alloc_default_constructed() increments `constructed` directly when reusing an
    // already-constructed slot (see arena.h -- del() deliberately does NOT reset()/destruct, only
    // alloc_default_constructed() does, right before real reuse, so a stranded region -- del()'d but
    // never recycled, e.g. lost from available by the concurrent-push_back race this exists to catch --
    // keeps every field at its real last-known state for debug printing instead of a blanked one).
    // Deliberately absent from every constructor's initializer list and given no default member
    // initializer, so a placement-new on a recycled slot leaves whatever count was already sitting
    // there untouched until the constructor body's own write -- the counts are tied to the memory
    // slot's whole Arena-driven recycling lifetime, not any one object's. A freshly malloc'd
    // (never-yet-constructed) slot needs both at a well-defined 0, which is why Arena::
    // alloc_unconstructed() zeroes new memory via calloc instead of malloc. A genuinely live object has
    // constructed == destructed + 1; Arena::del()/~Arena() check this immediately before acting, so a
    // repeat/lost construct-destruct pairing is caught at the point of failure instead of inferred
    // later from a lost available.push_back().
    int constructed;
    int destructed;

    GraphFillRegion();
    GraphFillRegion(GraphFillRegion&&);
    GraphFillRegion(const GraphFillRegion&) = delete;

    // Diagnostic alternative to actually destructing+later-reconstructing a recycled Arena slot:
    // clears every vector/pointer field back to the same state the default constructor establishes,
    // WITHOUT ending this object's C++ lifetime and WITHOUT touching constructed/destructed itself --
    // Arena::alloc_default_constructed() calls this (then increments `constructed` itself) only when
    // actually about to reuse an already-constructed slot; Arena::del() never calls this at all (see
    // its own comment in arena.h). A region is therefore never actually destructed mid-run; only
    // ~Arena()'s own final teardown ever calls the real (compiler-generated) destructor.
    void reset();
    bool tree_equal(const pm::GraphFillRegion& other) const;

    void add_match(pm::GraphFillRegion* match, const pm::CompressedEdge& edge);

    template <typename Callable>
    void do_op_for_each_node_in_total_area(const Callable& func);
    template <typename Callable>
    void do_op_for_each_descendant_and_self(const Callable& func);
    void clear_blossom_parent();
    void clear_blossom_parent_ignoring_wrapped_radius();
    void wrap_into_blossom(GraphFillRegion* new_blossom_parent_and_top);

    /// Determines if rhs is an ancestor of, or the same as, lhs.
    /// Caution: The equality used here is reference equality (not tree equality like operator==).
    bool operator<=(const GraphFillRegion& rhs) const;
    /// Determines if rhs is an ancestor of, but not the same as, lhs.
    bool operator<(const GraphFillRegion& rhs) const;
    /// Determines if rhs is a descendent of, or the same as, lhs.
    /// Caution: The equality used here is reference equality (not tree equality like operator==).
    bool operator>=(const GraphFillRegion& rhs) const;
    /// Determines if rhs is a descendent of, but not the same as, lhs.
    bool operator>(const GraphFillRegion& rhs) const;

    bool operator==(const GraphFillRegion& rhs) const;
    bool operator!=(const GraphFillRegion& rhs) const;
};

template <typename Callable>
inline void pm::GraphFillRegion::do_op_for_each_node_in_total_area(const Callable& func) {
    for (size_t i = 0; i < shell_area.size(); i++) {
        func(shell_area[shell_area.size() - i - 1]);
    }
    for (auto& child : blossom_children) {
        child.region->do_op_for_each_node_in_total_area(func);
    }
}

template <typename Callable>
inline void pm::GraphFillRegion::do_op_for_each_descendant_and_self(const Callable& func) {
    func(this);
    for (RegionEdge& child : blossom_children) {
        child.region->do_op_for_each_descendant_and_self(func);
    }
}

#ifdef USE_SHMEM
struct GraphFillRegionSubstates {
    GraphFillRegion* region;
    GraphFillRegion* blossom_parent;
    GraphFillRegion* blossom_parent_top;
    pm::VaryingCT radius; 
    pm::Match match;

    GraphFillRegionSubstates(GraphFillRegion* region_)
      : region(region_),
        blossom_parent(region_->blossom_parent),
        blossom_parent_top(region_->blossom_parent_top),
        radius(region_->radius),
        match(region_->match) {}
};
#endif

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_FILL_REGION_H
