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

#ifndef PYMATCHING2_MWPM_H
#define PYMATCHING2_MWPM_H

#include <fstream>

#include "pymatching/sparse_blossom/flooder/graph_flooder.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"
#include "pymatching/sparse_blossom/search/search_flooder.h"

#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"

namespace pm {

class AltTreeNode;

struct MatchingResult {
    pm::obs_int obs_mask;
    total_weight_int weight;
    MatchingResult();

    bool operator==(const MatchingResult& rhs) const;

    bool operator!=(const MatchingResult& rhs) const;

    MatchingResult(obs_int obs_mask, total_weight_int weight);

    MatchingResult& operator+=(const MatchingResult& rhs);
    MatchingResult operator+(const MatchingResult& rhs) const;
};

struct Mwpm {
    GraphFlooder flooder;
    Arena<AltTreeNode> node_arena;
    SearchFlooder search_flooder;

    std::pair<std::vector<std::pair<float, float>>, std::vector<std::pair<float, float>>> coords;
    TaskBase* task{ nullptr };
    int current_shot{ -1 };

    Mwpm();
    explicit Mwpm(GraphFlooder flooder);
    Mwpm(GraphFlooder flooder, SearchFlooder search_flooder);
    Mwpm(Mwpm&& other) noexcept;
    Mwpm& operator=(Mwpm&& other) noexcept;

    AltTreeNode* make_child(
        AltTreeNode& parent,
        GraphFillRegion* child_inner_region,
        GraphFillRegion* child_outer_region,
        const CompressedEdge& child_inner_to_outer_edge,
        const CompressedEdge& child_compressed_edge);
    void process_event(const pm::MwpmEvent& event);
    void handle_blossom_shattering(const BlossomShatterEventData& event);
    void shatter_descendants_into_matches_and_freeze(AltTreeNode& alt_tree_node);
    void handle_tree_hitting_boundary(const RegionHitBoundaryEventData& event);
    void handle_tree_hitting_virtual_boundary(const RegionHitVirtualBoundaryEventData& event);
    void handle_region_hit_region(const MwpmEvent event);
    void handle_tree_hitting_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_boundary_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_virtual_boundary_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_self(const RegionHitRegionEventData& event, AltTreeNode* common_ancestor);
    void handle_tree_hitting_other_tree(const RegionHitRegionEventData& event);
    // Removes matchings to virtual boundaries, turning matched regions into alternating trees
    void unmatch_virtual_boundaries_between_partitions();
    void prepare_for_task(TaskBase* task, int shot_id);
    // Sets search_flooder's vb/vb_left/vb_right from the given task, bounding
    // extract_paths_from_match_edges' Dijkstra search to the task's own partition (see
    // SearchFlooder::is_active). Separate from prepare_for_task: the task an extraction job is
    // scoped to isn't necessarily the task the solver last grew for, and extraction is a distinct
    // concern from task solving.
    // set_vb_to_part controls whether search_flooder.vb is set to task->vb_marker (default) or left
    // at -1 (never matches a real seam id): a received CRT window is bounded to [vb_left, vb_right]
    // but must NOT re-include the CRT's own seam id, since that seam was already divided -- the
    // local window on this side and the received window on the other side are extracted separately.
    void prepare_for_extraction(TaskBase* task, bool set_vb_to_part = true);
    // Removes region from task->regions_matched_to_virtual_boundary if present (swap-with-back +
    // pop_back -- the same removal handle_tree_hitting_virtual_boundary_match already does when a
    // vb match gets superseded by a region-region match; this is the shared implementation both use).
    // No-op if task is nullptr (the serial, non-parallel mwpm_decoding.cc path, which has no
    // TaskBase/regions_matched_to_virtual_boundary concept at all) or region isn't present (e.g.
    // already removed by handle_tree_hitting_virtual_boundary_match earlier this shot).
    static void remove_from_regions_matched_to_virtual_boundary(TaskBase* task, GraphFillRegion* region);
    // A region shattered here that's matched to a virtual boundary (region->match.edge.loc_to !=
    // nullptr -- the same convention DecodingUnit::send_solution_to_remote_pe already uses) must be
    // pruned from prune_target->regions_matched_to_virtual_boundary at the exact moment it's deleted
    // below: that list is populated during growth (handle_tree_hitting_virtual_boundary) and was
    // otherwise only ever pruned on the "matched to another region instead" path
    // (handle_tree_hitting_virtual_boundary_match) -- this closes the missing removal for the
    // ordinary extraction-time deletion path (the old workaround for that gap,
    // DecodingUnit::prune_stale_regions_matched_to_vb, is removed now that this makes the list never
    // go stale in the first place). prune_target is the caller's own already-correct target
    // (divide_vb's prune_target, send_solution_to_remote_pe's t) -- deliberately NOT Mwpm::task,
    // which is only kept current by prepare_for_task() (called exclusively from the climb loop) and
    // would be stale at every one of these call sites.
    // t_out: optional per-thread debug stream (decode_shots' own already-open one, threaded through
    // rather than opening a separate stream on the same file -- see the call sites in decoding_unit.cc
    // for how it's passed down). When DEBUG and non-null, logs a detected vb-matched region here (see
    // the print inside shatter_blossom_and_extract_matches/_match_edges below).
    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_matches(
        GraphFillRegion* region, MatchingResult& res, TaskBase* prune_target = nullptr,
        std::ofstream* t_out = nullptr);
    MatchingResult shatter_blossom_and_extract_matches(
        GraphFillRegion* region, TaskBase* prune_target = nullptr, std::ofstream* t_out = nullptr);

    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_match_edges(
        GraphFillRegion* region, std::vector<CompressedEdge>& match_edges, TaskBase* prune_target = nullptr,
        std::ofstream* t_out = nullptr);
    void shatter_blossom_and_extract_match_edges(
        GraphFillRegion* region, std::vector<CompressedEdge>& match_edges, TaskBase* prune_target = nullptr,
        std::ofstream* t_out = nullptr);
    void extract_paths_from_match_edges(
        const std::vector<CompressedEdge>& match_edges, uint8_t* obs_begin_ptr, pm::total_weight_int& weight);

    void verify_invariants() const;

    void create_detection_event(DetectorNode* node);
    void reset();
};
}  // namespace pm

#endif  // PYMATCHING2_MWPM_H
