#ifndef PYMATCHING2_SUBSET_PLAN_H
#define PYMATCHING2_SUBSET_PLAN_H

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "pymatching/sparse_blossom/driver/helpers/dem_generator.h"
#include "pymatching/sparse_blossom/driver/parallel/held_set.h"
#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/search/search_graph.h"

namespace pm {

class UserGraph;

/// One cross-observable seam: the two observables it joins and the vb window [vb_left+1, vb_right] of
/// (observable-local) partitions its fusion covers. unit_lo/unit_hi are its extraction-unit span.
struct SeamInfo {
    int oi, oj, vb_left, vb_right;
    int unit_lo{-1}, unit_hi{-1};
};

/// Options for building only the part of the graph a PE needs (see compute_held_set_obs). Observable strategy,
/// SHMEM builds only. `surgeries` must be the same seam spec the graph was generated from.
///
/// limited_memory: the UserGraph is not loaded up front. Instead each PE calls `load_user_graph` inside its turn on
/// its node (PEs of a node take turns, see node_turns.h), builds its partial graphs, and frees the UserGraph before
/// the next PE on the node starts, so a node never holds more than one UserGraph.
///
/// enabled: build only the partitions this PE needs when possible. It quietly falls back to a full build when it
/// cannot help or cannot be done (one PE, not the observable strategy, or a seam spec that does not match the graph).
/// required: the subset build was asked for explicitly, so those last two cases throw instead of falling back.
struct SubsetOptions {
    bool enabled{false};
    bool required{false};
    std::vector<SurgerySpec> surgeries;
    bool poison_unheld{false};
    bool limited_memory{false};
    std::function<UserGraph()> load_user_graph;
};

/// The effective cross-rank fusion window size k for a run (the cap DecodingUnit applies right after graph
/// construction). Throws std::invalid_argument for >2 PEs with fewer than 2 partitions per PE. `note` (if
/// non-null) receives the "NOTE: k ..." message when k is reduced.
int bound_effective_k(
    int k,
    int n_pes,
    size_t num_partitions,
    size_t p_per_obs_patch,
    bool obs_strategy,
    std::string* note = nullptr);

/// Seam infos by walking the built graph's seam-node neighbours (the original discovery). Kept as the
/// oracle for compute_seam_infos_from_spec. Windows are k-padded but NOT yet conflict-resolved.
std::vector<SeamInfo> compute_seam_infos_from_graph(const SharedMatchingGraph& graph, int k);

/// Seam infos derived only from the seam spec + round arithmetic (no graph): same output as
/// compute_seam_infos_from_graph on the graph generated from the same spec.
///   num_rounds: SharedMatchingGraph::num_rounds (distinct rounds of one observable)
///   M: rounds per partition (config_parallel::M)
///   K_p / K_vb: partitions / virtual boundaries per observable (p_per_obs_patch / vb_per_obs_patch)
std::vector<SeamInfo> compute_seam_infos_from_spec(
    const std::vector<SurgerySpec>& surgeries, int num_rounds, int M, int K_p, int K_vb, int k);

/// Global seam window conflict resolution (identical on every PE): narrows overlapping windows of seams
/// that share an observable. Throws std::invalid_argument if a window degenerates.
void resolve_seam_window_conflicts(std::vector<SeamInfo>& seam_infos, int pid);

/// Fills unit_lo/unit_hi (extraction-unit span, L partitions per unit) of every seam.
void assign_extraction_unit_spans(std::vector<SeamInfo>& seam_infos, int L);

/// Fills partition_bounds / vb_bounds ([first,last] node index of each partition / virtual boundary) from a
/// node -> partition (>=0) / virtual boundary (-(vb+1)) assignment. Same result as
/// SharedMatchingGraph::construct_partition_vb_bounds.
void compute_partition_vb_bounds(
    const std::vector<int>& node_part_id,
    size_t num_partitions,
    size_t num_virtual_boundaries,
    std::vector<std::pair<size_t, size_t>>& partition_bounds,
    std::vector<std::pair<size_t, size_t>>& vb_bounds);

/// Partitions / virtual boundaries a PE must have fully built under the OBS strategy (observable-major global
/// ids: observable o owns partitions [o*K_p, (o+1)*K_p) and vbs [o*K_vb, (o+1)*K_vb); seam s is vb K_vb*P+s):
///   - every partition and vb of its own observables [my_obs_start, my_obs_start+my_obs_count);
///   - for every seam with exactly one observable on this PE (a cross-rank task): the seam's own vb, and on the
///     remote observable the window partitions [vb_left+1, vb_right] and the vbs strictly inside the window
///     [vb_left+1, vb_right-1]. The vb just outside each end (vb_left, vb_right) bounds the fusion and is
///     inactive for it, so it only needs to exist as a stub (UserGraph::to_shared_matching_graph builds every
///     unheld neighbour of a held node as a stub with vb stamped);
///   - for a seam with both observables local: just the seam's vb.
/// `seam_infos` must be the k-padded windows BEFORE conflict resolution: a seam's nodes keep edges to every
/// round of the seam's span on both observables, and conflict resolution only narrows the windows the
/// decoding tasks use, not which nodes those edges reach.
HeldSet compute_held_set_obs(
    size_t num_partitions,
    size_t num_virtual_boundaries,
    int num_obs_patches,
    int K_p,
    int K_vb,
    int my_obs_start,
    int my_obs_count,
    const std::vector<SeamInfo>& seam_infos);

/// Result of comparing a partial (held-set) graph against the full graph built from the same UserGraph.
struct PartialGraphCheck {
    std::vector<std::string> errors;  // first max_errors mismatches
    size_t num_errors{0};
    size_t total_nodes{0}, num_held{0}, num_stubs{0};
};

/// Verifies `part` (built with `held`) against `full`: identical global metadata; every held node identical
/// (neighbour identity and order, weights, observable masks, vb, obs_patch_id); every unheld neighbour of a held
/// node is a stub (vb node, vb / obs_patch_id stamped, no neighbours) sitting on the outer end of some window
/// (`unresolved_seam_infos`) or on a seam's vb; every other unheld slot is untouched (all-zero bytes).
PartialGraphCheck verify_partial_against_full(
    const SharedMatchingGraph& full,
    const SharedMatchingGraph& part,
    const HeldSet& held,
    const std::vector<SeamInfo>& unresolved_seam_infos,
    size_t max_errors = 20);

/// Same checks for the SearchGraph pair (`node_part_id` is the MatchingGraph's, shared by both graphs).
PartialGraphCheck verify_partial_search_against_full(
    const SearchGraph& full,
    const SearchGraph& part,
    const std::vector<int>& node_part_id,
    const HeldSet& held,
    int K_vb,
    int num_obs_patches,
    const std::vector<SeamInfo>& unresolved_seam_infos,
    size_t max_errors = 20);

}  // namespace pm

#endif  // PYMATCHING2_SUBSET_PLAN_H
