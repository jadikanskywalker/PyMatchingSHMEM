#include "pymatching/sparse_blossom/driver/parallel/subset_plan.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <stdexcept>

namespace {

void throw_degenerate_window(int pid, int squeezed, int against, const pm::SeamInfo& s) {
    throw std::invalid_argument(
        "Rank " + std::to_string(pid) + ": global seam window conflict resolution squeezed seam " +
        std::to_string(squeezed) + " to a degenerate window (vb_left=" + std::to_string(s.vb_left) +
        ", vb_right=" + std::to_string(s.vb_right) + ") while resolving against seam " +
        std::to_string(against) +
        " -- 3+ seams sharing overlapping observable sets within too few partitions is not currently "
        "supported.");
}

}  // namespace

int pm::bound_effective_k(
    int k, int n_pes, size_t num_partitions, size_t p_per_obs_patch, bool obs_strategy, std::string* note) {
    size_t p_per_pe = num_partitions / n_pes;
    if (n_pes == 2 && (size_t)k > num_partitions / 2) {
        k = (int)(num_partitions / 2);
        if (note) *note = "NOTE: k set to " + std::to_string(k) + " for 2 PEs";
    } else {
        if (p_per_pe < 2) {
            throw std::invalid_argument("The number of partitions per PE should be >= 2 for more than 2 ranks.");
        } else if (!obs_strategy && (size_t)k > p_per_pe / 2) {
            k = (int)(p_per_pe / 2);
            if (note) *note = "NOTE: k bounded to " + std::to_string(k);
        } else if (obs_strategy && (size_t)k > p_per_obs_patch) {
            k = (int)p_per_obs_patch;
            if (note) *note = "NOTE: k bounded to " + std::to_string(k);
        }
    }
    return k;
}

std::vector<pm::SeamInfo> pm::compute_seam_infos_from_graph(const SharedMatchingGraph& graph, int k) {
    const int K_p = (int)graph.p_per_obs_patch;
    const int K_vb = (int)graph.vb_per_obs_patch;
    const int num_seams = (int)graph.num_virtual_boundaries - K_vb * (int)graph.num_obs_patches;
    std::vector<SeamInfo> seam_infos(num_seams);
    for (int s = 0; s < num_seams; ++s) {
        auto [first, last] = graph.vb_bounds[K_vb * graph.num_obs_patches + s];
        std::map<int, std::pair<int, int>> obs_part_range;  // obs_id -> (min_local_p, max_local_p)
        for (size_t ni = first; ni <= last; ++ni) {
            for (auto* nbr : graph.graph_ptr->nodes[ni].neighbors) {
                if (nbr == nullptr) continue;  // boundary
                int nbr_idx = (int)(nbr - graph.graph_ptr->nodes.data());
                int part_id = graph.node_part_id[nbr_idx];
                int obs_id, low_p, high_p;
                if (part_id >= 0) {  // partition node
                    obs_id = part_id / K_p;
                    low_p = part_id % K_p;
                    high_p = low_p;
                } else {
                    if (K_vb == 0) continue;  // no intra-obs VBs possible
                    int global_vb = -part_id - 1;
                    obs_id = global_vb / K_vb;
                    if (obs_id >= (int)graph.num_obs_patches) continue;  // cross-obs seam node
                    int local_vb = global_vb % K_vb;
                    // vb slot j lies between partitions j and j+1.
                    low_p = local_vb;
                    high_p = local_vb + 1;
                }
                auto [it, inserted] = obs_part_range.try_emplace(obs_id, low_p, high_p);
                if (!inserted) {
                    it->second.first = std::min(it->second.first, low_p);
                    it->second.second = std::max(it->second.second, high_p);
                }
            }
        }
        if (obs_part_range.empty()) {
            throw std::invalid_argument("Seam " + std::to_string(s) + " has no neighbours in any observable.");
        }
        auto it = obs_part_range.begin();
        int oi = it->first;
        int vb_left = std::max(it->second.first - k - 1, -1);
        int vb_right = std::min(it->second.second + k, K_p - 1);
        ++it;
        int oj = (it != obs_part_range.end()) ? it->first : oi;
        seam_infos[s] = {oi, oj, vb_left, vb_right};
    }
    return seam_infos;
}

std::vector<pm::SeamInfo> pm::compute_seam_infos_from_spec(
    const std::vector<SurgerySpec>& surgeries, int num_rounds, int M, int K_p, int K_vb, int k) {
    if (M < 2 || num_rounds < 2) {
        throw std::invalid_argument("compute_seam_infos_from_spec requires M >= 2 and num_rounds >= 2.");
    }
    std::vector<SeamInfo> seam_infos;
    seam_infos.reserve(surgeries.size());
    for (size_t s = 0; s < surgeries.size(); ++s) {
        const SurgerySpec& spec = surgeries[s];
        if (spec.duration < 1) {
            throw std::invalid_argument("Seam " + std::to_string(s) + " has no rounds (duration < 1).");
        }
        int lo = INT32_MAX, hi = INT32_MIN;
        for (int step = 0; step < spec.duration; ++step) {
            // build_seams snaps each target round to the nearest existing round, 0..num_rounds-1.
            int r = std::clamp(spec.start_round + step, 0, num_rounds - 1);
            int b = r / M;
            // A partition round lies in partition b. A virtual-boundary round (last of its M-round block)
            // is tagged as vb b, which lies between partitions b and b+1; a vb that would be the very last
            // round is folded into the final partition (see partition_nodes_by_obs_patch).
            bool is_vb_round = (r % M == M - 1) && (r < num_rounds - 1);
            int low_p = b;
            int high_p = is_vb_round ? b + 1 : b;
            if (is_vb_round && b >= K_vb) {
                throw std::logic_error("compute_seam_infos_from_spec: vb round beyond K_vb.");
            }
            lo = std::min(lo, low_p);
            hi = std::max(hi, high_p);
        }
        // Both observables see the same rounds; the window is sized from the lower-numbered one.
        SeamInfo si;
        si.oi = std::min(spec.obs_a, spec.obs_b);
        si.oj = std::max(spec.obs_a, spec.obs_b);
        si.vb_left = std::max(lo - k - 1, -1);
        si.vb_right = std::min(hi + k, K_p - 1);
        seam_infos.push_back(si);
    }
    return seam_infos;
}

void pm::resolve_seam_window_conflicts(std::vector<SeamInfo>& seam_infos, int pid) {
    const int num_seams = (int)seam_infos.size();
    for (int si_idx = 0; si_idx < num_seams; ++si_idx) {
        auto& si = seam_infos[si_idx];
        for (int sj_idx = si_idx + 1; sj_idx < num_seams; ++sj_idx) {
            auto& sj = seam_infos[sj_idx];
            if (si.oi != sj.oi && si.oi != sj.oj && si.oj != sj.oi && si.oj != sj.oj) continue;
            if (si.vb_left <= sj.vb_left) {
                if (si.vb_right > sj.vb_left) {
                    int mid = (si.vb_right + sj.vb_left + 1) / 2;
                    si.vb_right = mid;
                    sj.vb_left = mid;
                    // A prior adjustment (from a different overlapping pair sharing si or sj) can already
                    // have squeezed the other side of this same seam's window.
                    if (si.vb_left >= si.vb_right) throw_degenerate_window(pid, si_idx, sj_idx, si);
                    if (sj.vb_left >= sj.vb_right) throw_degenerate_window(pid, sj_idx, si_idx, sj);
                }
            } else {
                if (sj.vb_right > si.vb_left) {
                    int mid = (sj.vb_right + si.vb_left + 1) / 2;
                    sj.vb_right = mid;
                    si.vb_left = mid;
                    if (si.vb_left >= si.vb_right) throw_degenerate_window(pid, si_idx, sj_idx, si);
                    if (sj.vb_left >= sj.vb_right) throw_degenerate_window(pid, sj_idx, si_idx, sj);
                }
            }
        }
    }
}

void pm::assign_extraction_unit_spans(std::vector<SeamInfo>& seam_infos, int L) {
    for (auto& si : seam_infos) {
        si.unit_lo = (si.vb_left + 1) / L;
        si.unit_hi = si.vb_right / L;
    }
}

void pm::compute_partition_vb_bounds(
    const std::vector<int>& node_part_id,
    size_t num_partitions,
    size_t num_virtual_boundaries,
    std::vector<std::pair<size_t, size_t>>& partition_bounds,
    std::vector<std::pair<size_t, size_t>>& vb_bounds) {
    if (num_partitions <= 0) {
        throw std::invalid_argument("SharedMatchingGraph requires at least one partition.");
    }
    if (node_part_id.empty()) {
        throw std::invalid_argument("SharedMatchingGraph requires node partition assignments.");
    }
    partition_bounds.resize(num_partitions);
    vb_bounds.resize(num_virtual_boundaries);
    partition_bounds[0].first = 0;
    size_t partition_bounds_idx = 0;
    size_t vb_bounds_idx = 0;
    int last_part_id = node_part_id[0];
    for (int i = 1; i < (int)node_part_id.size(); ++i) {
        if (node_part_id[i] == last_part_id) {
            continue;
        } else if (last_part_id < 0 && node_part_id[i] < 0) {  // vb -> vb
            vb_bounds[vb_bounds_idx++].second = i - 1;
            vb_bounds[vb_bounds_idx].first = i;
            last_part_id = node_part_id[i];
        } else if (node_part_id[i] < 0) {  // p -> vb
            partition_bounds[partition_bounds_idx].second = i - 1;
            vb_bounds[vb_bounds_idx].first = i;
            last_part_id = node_part_id[i];
        } else if (last_part_id < 0) {  // vb -> p
            vb_bounds[vb_bounds_idx++].second = i - 1;
            partition_bounds[++partition_bounds_idx].first = i;
            last_part_id = node_part_id[i];
        } else {  // p -> p (obs boundary)
            partition_bounds[partition_bounds_idx].second = i - 1;
            partition_bounds[++partition_bounds_idx].first = i;
            last_part_id = node_part_id[i];
        }
    }
    // Close the final block (partition or vb if cross-obs nodes are last).
    if (last_part_id < 0) {
        vb_bounds[vb_bounds_idx].second = (int)node_part_id.size() - 1;
    } else {
        partition_bounds[partition_bounds_idx].second = (int)node_part_id.size() - 1;
    }
}

pm::HeldSet pm::compute_held_set_obs(
    size_t num_partitions,
    size_t num_virtual_boundaries,
    int num_obs_patches,
    int K_p,
    int K_vb,
    int my_obs_start,
    int my_obs_count,
    const std::vector<SeamInfo>& seam_infos) {
    HeldSet held;
    held.part_held.assign(num_partitions, false);
    held.vb_held.assign(num_virtual_boundaries, false);
    auto is_local = [&](int o) { return o >= my_obs_start && o < my_obs_start + my_obs_count; };
    for (int o = my_obs_start; o < my_obs_start + my_obs_count; ++o) {
        for (int p = 0; p < K_p; ++p) held.part_held[(size_t)o * K_p + p] = true;
        for (int v = 0; v < K_vb; ++v) held.vb_held[(size_t)o * K_vb + v] = true;
    }
    for (size_t s = 0; s < seam_infos.size(); ++s) {
        const SeamInfo& si = seam_infos[s];
        bool local_i = is_local(si.oi), local_j = is_local(si.oj);
        if (!local_i && !local_j) continue;
        held.vb_held[(size_t)K_vb * num_obs_patches + s] = true;
        for (int remote : {si.oi, si.oj}) {
            if (is_local(remote)) continue;
            for (int p = si.vb_left + 1; p <= si.vb_right; ++p) {
                held.part_held[(size_t)remote * K_p + p] = true;
            }
            // vbs strictly inside the window; the vb just outside each end is inactive for the fusion and is only
            // built as a stub (vb stamped) because it neighbours a window partition.
            for (int v = std::max(si.vb_left + 1, 0); v <= std::min(si.vb_right - 1, K_vb - 1); ++v) {
                held.vb_held[(size_t)remote * K_vb + v] = true;
            }
        }
    }
    return held;
}

namespace {

template <typename Node>
bool is_zero_slot(const Node& n) {
    const unsigned char* b = reinterpret_cast<const unsigned char*>(&n);
    for (size_t i = 0; i < sizeof(Node); ++i)
        if (b[i]) return false;
    return true;
}

// Node-by-node comparison shared by the MatchingGraph and SearchGraph checks. `extra_equal(a, b)` compares the
// node-type-specific per-edge arrays of two held nodes.
template <typename NodeVecT, typename ExtraEqual>
void verify_nodes(
    const NodeVecT& fn,
    const NodeVecT& pn,
    const std::vector<int>& node_part_id,
    const pm::HeldSet& held,
    int K_vb,
    int num_obs,
    const std::vector<pm::SeamInfo>& unresolved,
    pm::PartialGraphCheck& out,
    size_t max_errors,
    ExtraEqual extra_equal) {
    auto fail = [&](const std::string& msg) {
        ++out.num_errors;
        if (out.errors.size() < max_errors) out.errors.push_back(msg);
    };
    if (pn.size() != fn.size()) {
        fail("node count differs");
        return;
    }
    out.total_nodes = fn.size();
    std::vector<bool> is_stub(fn.size(), false);
    for (size_t i = 0; i < fn.size(); ++i) {
        if (!held.node_held(node_part_id[i])) continue;
        for (const auto* nb : fn[i].neighbors) {
            if (nb && !held.node_held(node_part_id[nb - fn.data()])) is_stub[nb - fn.data()] = true;
        }
    }
    auto id = [](size_t i) { return std::to_string(i); };
    for (size_t i = 0; i < fn.size(); ++i) {
        if (!held.node_held(node_part_id[i])) {
            if (!is_stub[i]) {
                if (!is_zero_slot(pn[i])) fail("unheld node " + id(i) + " was touched");
                continue;
            }
            ++out.num_stubs;
            int part_id = node_part_id[i];
            if (part_id >= 0) {
                fail("stub " + id(i) + " is a partition node (held set misses a neighbour)");
                continue;
            }
            int gvb = -part_id - 1;
            if (gvb < K_vb * num_obs) {  // regular vb: must be the vb just outside a window end
                int r = gvb / K_vb, j = gvb % K_vb;
                bool outer = false;
                for (const auto& w : unresolved)
                    if ((w.oi == r || w.oj == r) && (j == w.vb_left || j == w.vb_right)) outer = true;
                if (!outer) fail("regular-vb stub " + id(i) + " is not a window's outer vb");
            }
            if (!pn[i].neighbors.empty()) fail("stub " + id(i) + " has neighbours");
            if (pn[i].vb != fn[i].vb) fail("stub " + id(i) + " vb differs");
            continue;
        }
        ++out.num_held;
        if (pn[i].neighbors.size() != fn[i].neighbors.size()) {
            fail("held node " + id(i) + " neighbour count differs");
            continue;
        }
        if (pn[i].neighbor_weights != fn[i].neighbor_weights) fail("held node " + id(i) + " weights differ");
        if (!extra_equal(fn[i], pn[i])) fail("held node " + id(i) + " per-edge observables/markers differ");
        if (pn[i].vb != fn[i].vb) fail("held node " + id(i) + " vb differs");
        for (size_t j = 0; j < fn[i].neighbors.size(); ++j) {
            const auto* a = fn[i].neighbors[j];
            const auto* b = pn[i].neighbors[j];
            if ((a == nullptr) != (b == nullptr) || (a && (a - fn.data()) != (b - pn.data()))) {
                fail("held node " + id(i) + " neighbour " + std::to_string(j) + " differs");
            }
        }
    }
}

}  // namespace

pm::PartialGraphCheck pm::verify_partial_against_full(
    const SharedMatchingGraph& full,
    const SharedMatchingGraph& part,
    const HeldSet& held,
    const std::vector<SeamInfo>& unresolved,
    size_t max_errors) {
    PartialGraphCheck out;
    auto fail = [&](const std::string& msg) {
        ++out.num_errors;
        if (out.errors.size() < max_errors) out.errors.push_back(msg);
    };
    if (part.node_part_id != full.node_part_id) fail("node_part_id differs");
    if (part.partition_bounds != full.partition_bounds) fail("partition_bounds differ");
    if (part.vb_bounds != full.vb_bounds) fail("vb_bounds differ");
    if (part.graph_ptr->normalising_constant != full.graph_ptr->normalising_constant) fail("normalising_constant differs");
    if (part.graph_ptr->negative_weight_sum != full.graph_ptr->negative_weight_sum) fail("negative_weight_sum differs");
    verify_nodes(
        full.graph_ptr->nodes, part.graph_ptr->nodes, full.node_part_id, held, (int)full.vb_per_obs_patch,
        (int)full.num_obs_patches, unresolved, out, max_errors,
        [](const DetectorNode& a, const DetectorNode& b) { return a.neighbor_observables == b.neighbor_observables; });
    // obs_patch_id is stamped on every built node (held and stub).
    const auto& fn = full.graph_ptr->nodes;
    const auto& pn = part.graph_ptr->nodes;
    for (size_t i = 0; i < fn.size() && out.errors.size() < max_errors; ++i) {
        bool built = held.node_held(full.node_part_id[i]) || !is_zero_slot(pn[i]);
        if (built && pn[i].obs_patch_id != fn[i].obs_patch_id) fail("node " + std::to_string(i) + " obs_patch_id differs");
    }
    return out;
}

pm::PartialGraphCheck pm::verify_partial_search_against_full(
    const SearchGraph& full,
    const SearchGraph& part,
    const std::vector<int>& node_part_id,
    const HeldSet& held,
    int K_vb,
    int num_obs_patches,
    const std::vector<SeamInfo>& unresolved,
    size_t max_errors) {
    PartialGraphCheck out;
    if (part.negative_weight_edges != full.negative_weight_edges) {
        ++out.num_errors;
        out.errors.push_back("negative_weight_edges differ");
    }
    verify_nodes(
        full.nodes, part.nodes, node_part_id, held, K_vb, num_obs_patches, unresolved, out, max_errors,
        [](const SearchDetectorNode& a, const SearchDetectorNode& b) {
            return a.neighbor_observable_indices == b.neighbor_observable_indices &&
                   a.neighbor_markers == b.neighbor_markers;
        });
    return out;
}
