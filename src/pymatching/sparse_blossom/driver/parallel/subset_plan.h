#ifndef PYMATCHING2_SUBSET_PLAN_H
#define PYMATCHING2_SUBSET_PLAN_H

#include <string>
#include <utility>
#include <vector>

#include "pymatching/sparse_blossom/driver/helpers/dem_generator.h"
#include "pymatching/sparse_blossom/flooder/graph.h"

namespace pm {

/// One cross-observable seam: the two observables it joins and the vb window [vb_left+1, vb_right] of
/// (observable-local) partitions its fusion covers. unit_lo/unit_hi are its extraction-unit span.
struct SeamInfo {
    int oi, oj, vb_left, vb_right;
    int unit_lo{-1}, unit_hi{-1};
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

}  // namespace pm

#endif  // PYMATCHING2_SUBSET_PLAN_H
