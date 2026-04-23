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

#include "pymatching/sparse_blossom/driver/user_graph.h"

#include "pymatching/rand/rand_gen.h"
#include "pymatching/sparse_blossom/driver/implied_weights.h"

#ifdef USE_THREADS
#include "../config_parallel.h"
#endif

namespace {

double bernoulli_xor(double p1, double p2) {
    return p1 * (1 - p2) + p2 * (1 - p1);
}

}  // namespace


double pm::to_weight_for_correlations(double probability) {
    return std::log((1 - probability) / probability);
}

double pm::merge_weights(double a, double b) {
    auto sgn = std::copysign(1, a) * std::copysign(1, b);
    auto signed_min = sgn * std::min(std::abs(a), std::abs(b));
    return signed_min + std::log(1 + std::exp(-std::abs(a + b))) - std::log(1 + std::exp(-std::abs(a - b)));
}

pm::UserNode::UserNode() : is_boundary(false) {
}

size_t pm::UserNode::index_of_neighbor(size_t node) const {
    auto it = std::find_if(neighbors.begin(), neighbors.end(), [&](const UserNeighbor& neighbor) {
        if (neighbor.pos == 0) {
            return neighbor.edge_it->node1 == node;
        } else if (neighbor.pos == 1) {
            return neighbor.edge_it->node2 == node;
        } else {
            throw std::runtime_error("`neighbor.pos` should be 0 or 1, but got: " + std::to_string(neighbor.pos));
        }
    });
    if (it == neighbors.end())
        return SIZE_MAX;

    return it - neighbors.begin();
}

bool is_valid_probability(double p) {
    return (p >= 0 && p <= 1);
}

void pm::UserGraph::merge_edge_or_boundary_edge(
    size_t node,
    size_t neighbor_index,
    const std::vector<size_t>& parallel_observables,
    double parallel_weight,
    double parallel_error_probability,
    pm::MERGE_STRATEGY merge_strategy) {
    auto& neighbor = nodes[node].neighbors[neighbor_index];
    if (merge_strategy == DISALLOW) {
        throw std::invalid_argument(
            "Edge (" + std::to_string(neighbor.edge_it->node1) + ", " + std::to_string(neighbor.edge_it->node2) +
            ") already exists in the graph. "
            "Parallel edges not permitted with the provided `disallow` `merge_strategy`. Please provide a "
            "different `merge_strategy`.");
    } else if (
        merge_strategy == KEEP_ORIGINAL ||
        (merge_strategy == SMALLEST_WEIGHT && parallel_weight >= neighbor.edge_it->weight)) {
        return;
    } else {
        double new_weight, new_error_probability;
        bool use_new_observables;
        if (merge_strategy == REPLACE || merge_strategy == SMALLEST_WEIGHT) {
            new_weight = parallel_weight;
            new_error_probability = parallel_error_probability;
            use_new_observables = true;
        } else if (merge_strategy == INDEPENDENT) {
            new_weight = pm::merge_weights(parallel_weight, neighbor.edge_it->weight);
            new_error_probability = -1;
            if (is_valid_probability(neighbor.edge_it->error_probability) &&
                is_valid_probability(parallel_error_probability))
                new_error_probability = parallel_error_probability * (1 - neighbor.edge_it->error_probability) +
                                        neighbor.edge_it->error_probability * (1 - parallel_error_probability);
            // We do not need to update the observables. If they do not match up, then the code has distance 2.
            use_new_observables = false;
        } else {
            throw std::invalid_argument("Merge strategy not recognised.");
        }
        // Update the existing edge weight and probability in the adjacency list of `node`
        neighbor.edge_it->weight = new_weight;
        neighbor.edge_it->error_probability = new_error_probability;
        if (use_new_observables)
            neighbor.edge_it->observable_indices = parallel_observables;

        _mwpm_needs_updating = true;
        if (new_error_probability < 0 || new_error_probability > 1)
            _all_edges_have_error_probabilities = false;
    }
}

void pm::UserGraph::add_or_merge_edge(
    size_t node1,
    size_t node2,
    const std::vector<size_t>& observables,
    double weight,
    double error_probability,
    MERGE_STRATEGY merge_strategy) {
    auto max_id = std::max(node1, node2);
    if (max_id + 1 > nodes.size())
        nodes.resize(max_id + 1);

    size_t idx = nodes[node1].index_of_neighbor(node2);

    if (idx == SIZE_MAX) {
        pm::UserEdge edge = {node1, node2, observables, weight, error_probability};
        edges.push_back(edge);
        nodes[node1].neighbors.push_back({std::prev(edges.end()), 1});
        if (node1 != node2)
            nodes[node2].neighbors.push_back({std::prev(edges.end()), 0});

        for (auto& obs : observables) {
            if (obs + 1 > _num_observables)
                _num_observables = obs + 1;
        }
        _mwpm_needs_updating = true;
        if (error_probability < 0 || error_probability > 1)
            _all_edges_have_error_probabilities = false;
    } else {
        merge_edge_or_boundary_edge(node1, idx, observables, weight, error_probability, merge_strategy);
    }
}

void pm::UserGraph::add_or_merge_boundary_edge(
    size_t node,
    const std::vector<size_t>& observables,
    double weight,
    double error_probability,
    MERGE_STRATEGY merge_strategy) {
    if (node + 1 > nodes.size())
        nodes.resize(node + 1);

    size_t idx = nodes[node].index_of_neighbor(SIZE_MAX);

    if (idx == SIZE_MAX) {
        pm::UserEdge edge = {node, SIZE_MAX, observables, weight, error_probability};
        edges.push_back(edge);
        nodes[node].neighbors.push_back({std::prev(edges.end()), 1});

        for (auto& obs : observables) {
            if (obs + 1 > _num_observables)
                _num_observables = obs + 1;
        }
        _mwpm_needs_updating = true;
        if (error_probability < 0 || error_probability > 1)
            _all_edges_have_error_probabilities = false;
    } else {
        merge_edge_or_boundary_edge(node, idx, observables, weight, error_probability, merge_strategy);
    }
}

pm::UserGraph::UserGraph()
    : _num_observables(0), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
}

pm::UserGraph::UserGraph(size_t num_nodes)
    : _num_observables(0), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
    nodes.resize(num_nodes);
}

pm::UserGraph::UserGraph(size_t num_nodes, size_t num_observables)
    : _num_observables(num_observables), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
    nodes.resize(num_nodes);
}

void pm::UserGraph::set_boundary(const std::set<size_t>& boundary) {
    for (auto& n : boundary_nodes)
        nodes[n].is_boundary = false;
    boundary_nodes = boundary;
    for (auto& n : boundary_nodes) {
        if (n >= nodes.size())
            nodes.resize(n + 1);
        nodes[n].is_boundary = true;
    }
    _mwpm_needs_updating = true;
}

std::set<size_t> pm::UserGraph::get_boundary() {
    return boundary_nodes;
}

size_t pm::UserGraph::get_num_observables() {
    return _num_observables;
}

size_t pm::UserGraph::get_num_nodes() {
    return nodes.size();
}

size_t pm::UserGraph::get_num_detectors() {
    return get_num_nodes() - boundary_nodes.size();
}

bool pm::UserGraph::is_boundary_node(size_t node_id) {
    return (node_id == SIZE_MAX) || nodes[node_id].is_boundary;
}

void pm::UserGraph::update_mwpm() {
    _mwpm = to_mwpm(pm::NUM_DISTINCT_WEIGHTS, false);
    _mwpm_needs_updating = false;
}

pm::Mwpm& pm::UserGraph::get_mwpm() {
    if (_mwpm_needs_updating)
        update_mwpm();
    return _mwpm;
}

void pm::UserGraph::add_noise(uint8_t* error_arr, uint8_t* syndrome_arr) const {
    if (!_all_edges_have_error_probabilities)
        return;

    for (auto& e : edges) {
        auto p = e.error_probability;
        if (rand_float(0.0, 1.0) < p) {
            // Flip the observables
            for (auto& obs : e.observable_indices) {
                *(error_arr + obs) ^= 1;
            }
            // Flip the syndrome bits
            *(syndrome_arr + e.node1) ^= 1;
            if (e.node2 != SIZE_MAX)
                *(syndrome_arr + e.node2) ^= 1;
        }
    }

    for (auto& b : boundary_nodes)
        *(syndrome_arr + b) = 0;
}

size_t pm::UserGraph::get_num_edges() {
    return edges.size();
}

bool pm::UserGraph::all_edges_have_error_probabilities() {
    return _all_edges_have_error_probabilities;
}

double pm::UserGraph::max_abs_weight() {
    double max_abs_weight = 0;
    for (auto& e : edges) {
        if (std::abs(e.weight) > max_abs_weight) {
            max_abs_weight = std::abs(e.weight);
        }
    }
    return max_abs_weight;
}

#ifdef USE_THREADS
pm::SharedMatchingGraph pm::UserGraph::to_shared_matching_graph(
    pm::weight_int num_distinct_weights
#ifdef USE_SHMEM
    , DetectorNodeEphemeralFields* node_ephemeral_fields_ptr
#endif
) {
    std::shared_ptr<MatchingGraph> matching_graph_ptr = std::make_shared<pm::MatchingGraph>(nodes.size(), _num_observables);
    pm::MatchingGraph& matching_graph = *matching_graph_ptr;
#ifdef USE_SHMEM
    const int num_nodes = nodes.size();
    for (int i=0; i < num_nodes; ++i) {
        for (int j=0; j < NUM_BUFFERS_PER_UNIT; ++j) {
            matching_graph.nodes[i].ephemeral_fields[j] = node_ephemeral_fields_ptr + i + j*num_nodes;
        }
    }
#endif
    double normalising_constant = to_matching_or_search_graph_helper(
        num_distinct_weights,
        [&](size_t u,
            size_t v,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            matching_graph.add_edge(u, v, weight, observables, implied_weights_for_other_edges);
        },
        [&](size_t u,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            matching_graph.add_boundary_edge(u, weight, observables, implied_weights_for_other_edges);
        });

    matching_graph.normalising_constant = normalising_constant;
    if (boundary_nodes.size() > 0) {
        matching_graph.is_user_graph_boundary_node.clear();
        matching_graph.is_user_graph_boundary_node.resize(nodes.size(), false);
        for (auto& i : boundary_nodes)
            matching_graph.is_user_graph_boundary_node[i] = true;
    }
    matching_graph.convert_implied_weights(normalising_constant);

    size_t num_regular_vb_masks = virtual_boundaries.size();
#ifdef USE_SHMEM
    if (DEBUG) std::cout << "virtual_boundaries.size(): " << virtual_boundaries.size() << "\n" << std::flush;
    if (config_parallel::division_strategy == config_parallel::OBS) {
        num_regular_vb_masks = vb_per_obs_patch;
        if (DEBUG) std::cout << "num_regular_vb_masks: " << num_regular_vb_masks << std::endl << std::flush;
        for (int vb = vb_per_obs_patch; vb < virtual_boundaries.size(); ++vb) {
            size_t vb_marker = vb + vb_per_obs_patch * (num_obs_patches - 1);
            if (DEBUG) std::cout << "  using " << vb_marker << std::endl << std::flush;
            for (int index : virtual_boundaries[vb]) {
                matching_graph.nodes[index].vb = vb_marker;
            }
        }
    }
#endif
    for (int vb=0; vb < num_regular_vb_masks; ++vb) {
        for (int index : virtual_boundaries[vb]) {
            matching_graph.nodes[index].vb = vb;
        }
    }

    return SharedMatchingGraph(matching_graph_ptr, node_part_id, num_partitions, num_virtual_boundaries, num_rounds
#ifdef USE_SHMEM
        , num_obs_patches, p_per_obs_patch, vb_per_obs_patch
#endif
    );
}
#endif

pm::MatchingGraph pm::UserGraph::to_matching_graph(pm::weight_int num_distinct_weights) {
    pm::MatchingGraph matching_graph(nodes.size(), _num_observables);

    double normalising_constant = to_matching_or_search_graph_helper(
        num_distinct_weights,
        [&](size_t u,
            size_t v,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            matching_graph.add_edge(u, v, weight, observables, implied_weights_for_other_edges);
        },
        [&](size_t u,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            matching_graph.add_boundary_edge(u, weight, observables, implied_weights_for_other_edges);
        });

    matching_graph.normalising_constant = normalising_constant;
    if (boundary_nodes.size() > 0) {
        matching_graph.is_user_graph_boundary_node.clear();
        matching_graph.is_user_graph_boundary_node.resize(nodes.size(), false);
        for (auto& i : boundary_nodes)
            matching_graph.is_user_graph_boundary_node[i] = true;
    }

    matching_graph.convert_implied_weights(normalising_constant);

    return matching_graph;
}

pm::SearchGraph pm::UserGraph::to_search_graph(pm::weight_int num_distinct_weights) {
    /// Identical to to_matching_graph but for constructing a pm::SearchGraph
    pm::SearchGraph search_graph(nodes.size());

    double normalizing_constant = to_matching_or_search_graph_helper(
        num_distinct_weights,
        [&](size_t u,
            size_t v,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            search_graph.add_edge(u, v, weight, observables, implied_weights_for_other_edges);
        },
        [&](size_t u,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            search_graph.add_boundary_edge(u, weight, observables, implied_weights_for_other_edges);
        });

    search_graph.convert_implied_weights(normalizing_constant);

// #ifdef USE_THREADS
// // ===============
//     // Propagate partition and virtual metadata
//     search_graph.num_partitions = num_partitions;
//     for (size_t i = 0; i < nodes.size(); ++i) {
//         search_graph.nodes[i].partition = nodes[i].partition;
//         search_graph.nodes[i].is_virtual = nodes[i].is_virtual;
//     }
// // ===============
// #endif

    return search_graph;
}

pm::Mwpm pm::UserGraph::to_mwpm(pm::weight_int num_distinct_weights, bool ensure_search_graph_included) {
    if (_num_observables > sizeof(pm::obs_int) * 8 || ensure_search_graph_included) {
        auto mwpm = pm::Mwpm(
            pm::GraphFlooder(to_matching_graph(num_distinct_weights)),
            pm::SearchFlooder(to_search_graph(num_distinct_weights)));
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
        mwpm.flooder.graph.loaded_from_dem_without_correlations = loaded_from_dem_without_correlations;
        return mwpm;
    } else {
        auto mwpm = pm::Mwpm(pm::GraphFlooder(to_matching_graph(num_distinct_weights)));
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
        mwpm.flooder.graph.loaded_from_dem_without_correlations = loaded_from_dem_without_correlations;
        return mwpm;
    }
}

pm::Mwpm& pm::UserGraph::get_mwpm_with_search_graph() {
    if (!_mwpm_needs_updating && _mwpm.flooder.graph.nodes.size() == _mwpm.search_flooder.graph.nodes.size()) {
        return _mwpm;
    } else {
        _mwpm = to_mwpm(pm::NUM_DISTINCT_WEIGHTS, true);
        _mwpm_needs_updating = false;
        return _mwpm;
    }
}

void pm::UserGraph::handle_dem_instruction(
    double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables) {
    if (detectors.size() == 2) {
        add_or_merge_edge(detectors[0], detectors[1], observables, std::log((1 - p) / p), p, INDEPENDENT);
    } else if (detectors.size() == 1) {
        add_or_merge_boundary_edge(detectors[0], observables, std::log((1 - p) / p), p, INDEPENDENT);
    }
}

void pm::UserGraph::handle_dem_instruction_include_correlations(
    double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables) {
    if (detectors.size() == 2) {
        add_or_merge_edge(detectors[0], detectors[1], observables, pm::to_weight_for_correlations(p), p, INDEPENDENT);
    } else if (detectors.size() == 1) {
        add_or_merge_boundary_edge(detectors[0], observables, pm::to_weight_for_correlations(p), p, INDEPENDENT);
    }
}

void pm::UserGraph::get_nodes_on_shortest_path_from_source(size_t src, size_t dst, std::vector<size_t>& out_nodes) {
    auto& mwpm = get_mwpm_with_search_graph();
    bool src_is_boundary = is_boundary_node(src);
    bool dst_is_boundary = is_boundary_node(dst);
    if (src != SIZE_MAX && src >= nodes.size())
        throw std::invalid_argument("node " + std::to_string(src) + " is not in the graph");
    if (dst != SIZE_MAX && dst >= nodes.size())
        throw std::invalid_argument("node " + std::to_string(dst) + " is not in the graph");
    if (!src_is_boundary) {
        size_t dst_tmp = dst_is_boundary ? SIZE_MAX : dst;
        mwpm.search_flooder.iter_edges_on_shortest_path_from_source(src, dst_tmp, [&](const pm::SearchGraphEdge edge) {
            out_nodes.push_back(edge.detector_node - &mwpm.search_flooder.graph.nodes[0]);
        });
        if (!dst_is_boundary)
            out_nodes.push_back(dst);
    } else if (!dst_is_boundary) {
        std::vector<size_t> temp_out_nodes;
        get_nodes_on_shortest_path_from_source(dst, src, temp_out_nodes);
        for (size_t i = 0; i < temp_out_nodes.size(); i++) {
            out_nodes.push_back(temp_out_nodes[temp_out_nodes.size() - 1 - i]);
        }
    } else {
        throw std::invalid_argument("Both the source and destination vertices provided are boundary nodes");
    }
}

bool pm::UserGraph::has_edge(size_t node1, size_t node2) {
    if (node1 >= nodes.size())
        return false;
    return nodes[node1].index_of_neighbor(node2) != SIZE_MAX;
}

bool pm::UserGraph::has_boundary_edge(size_t node) {
    if (node >= nodes.size())
        return false;
    return nodes[node].index_of_neighbor(SIZE_MAX) != SIZE_MAX;
}

bool pm::UserGraph::get_edge_or_boundary_edge_weight(size_t node1, size_t node2, double& weight_out) {
    if (node1 >= nodes.size()) {
        return false;
    }
    size_t neighbor_idx = nodes[node1].index_of_neighbor(node2);
    if (neighbor_idx == SIZE_MAX) {
        return false;
    }
    weight_out = nodes[node1].neighbors[neighbor_idx].edge_it->weight;
    return true;
}

void pm::UserGraph::set_min_num_observables(size_t num_observables) {
    if (num_observables > _num_observables)
        _num_observables = num_observables;
}

double pm::UserGraph::get_edge_weight_normalising_constant(size_t max_num_distinct_weights) {
    double max_abs_weight = 0;
    bool all_integral_weight = true;
    for (auto& e : edges) {
        if (std::abs(e.weight) > max_abs_weight)
            max_abs_weight = std::abs(e.weight);

        if (round(e.weight) != e.weight) {
            all_integral_weight = false;
        }

        for (auto implied : e.implied_weights_for_other_edges) {
            if (std::abs(implied.implied_weight) > max_abs_weight) {
                max_abs_weight = std::abs(implied.implied_weight);
            }

            if (round(implied.implied_weight) != implied.implied_weight) {
                all_integral_weight = false;
            }

            double current_weight;
            bool has_edge = get_edge_or_boundary_edge_weight(implied.node1, implied.node2, current_weight);
            if (!has_edge) {
                throw std::invalid_argument(
                    "Edge rewrite rule refers to non-existent edge (" + std::to_string(implied.node1) + ", " +
                    std::to_string(implied.node2) + ")");
            }
            bool same_sign = (current_weight * implied.implied_weight) >= 0.;
            if (!same_sign) {
                throw std::invalid_argument(
                    "Edge weight rewrite rules that change the sign of an edge weight are not currently supported.");
            }
        }
    }

    if (max_abs_weight > pm::MAX_USER_EDGE_WEIGHT)
        throw std::invalid_argument(
            "maximum absolute edge weight of " + std::to_string(pm::MAX_USER_EDGE_WEIGHT) + " exceeded.");

    if (all_integral_weight) {
        return 1.0;
    } else {
        pm::weight_int max_half_edge_weight = max_num_distinct_weights - 1;
        return (double)max_half_edge_weight / max_abs_weight;
    }
}

void pm::add_decomposed_error_to_joint_probabilities(
    DecomposedDemError& error,
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>>& joint_probabilites) {
    if (error.components.size() > 1) {
        for (size_t k0 = 0; k0 < error.components.size(); k0++) {
            for (size_t k1 = k0 + 1; k1 < error.components.size(); k1++) {
                auto& c0 = error.components[k0];
                auto& c1 = error.components[k1];
                std::pair<size_t, size_t> e0 = std::minmax(c0.node1, c0.node2);
                std::pair<size_t, size_t> e1 = std::minmax(c1.node1, c1.node2);
                double& p01 = joint_probabilites[e0][e1];
                double& p10 = joint_probabilites[e1][e0];
                p01 = bernoulli_xor(p01, error.probability);
                p10 = bernoulli_xor(p10, error.probability);
            }
        }
    }

    for (auto& e : error.components) {
        double& p = joint_probabilites[std::minmax(e.node1, e.node2)][std::minmax(e.node1, e.node2)];
        p = bernoulli_xor(p, error.probability);
    }
}

pm::UserGraph pm::detector_error_model_to_user_graph(
    const stim::DetectorErrorModel& detector_error_model,
    const bool enable_correlations,
    pm::weight_int num_distinct_weights
    ) {
    pm::UserGraph user_graph(detector_error_model.count_detectors(), detector_error_model.count_observables());
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilites;
    if (enable_correlations) {
        pm::iter_dem_instructions_include_correlations(
            detector_error_model,
            [&](double p, const std::vector<size_t>& detectors, std::vector<size_t>& observables) {
                user_graph.handle_dem_instruction_include_correlations(p, detectors, observables);
            },
            joint_probabilites);

        user_graph.populate_implied_edge_weights(joint_probabilites);
    } else {
        pm::iter_detector_error_model_edges(
            detector_error_model,
            [&](double p, const std::vector<size_t>& detectors, std::vector<size_t>& observables) {
                user_graph.handle_dem_instruction(p, detectors, observables);
            });
        user_graph.loaded_from_dem_without_correlations = true;
    }
#ifdef USE_THREADS
#ifdef USE_SHMEM
    if (config_parallel::division_strategy == config_parallel::OBS) {
        user_graph.partition_nodes_by_obs_patch(detector_error_model);
    } else {
        user_graph.partition_nodes_by_round(detector_error_model);
    }
#else
    user_graph.partition_nodes_by_round(detector_error_model);
#endif
#endif
// reorder_nodes_by_observable() removed: DEM is now generated in observable-major
// order by gen_multi_obs.py so no post-hoc reorder is needed.
    return user_graph;
}

void pm::UserGraph::populate_implied_edge_weights(
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>>& joint_probabilites) {
    for (auto& edge : edges) {
        std::pair<size_t, size_t> current_edge_nodes = std::minmax(edge.node1, edge.node2);
        auto it = joint_probabilites.find(current_edge_nodes);
        if (it != joint_probabilites.end()) {
            const auto& pf = *it;
            std::pair<size_t, size_t> causal_edge = pf.first;
            double marginal_probability = pf.second.at(causal_edge);
            if (marginal_probability == 0)
                continue;

            for (const auto& affected_edge_and_probability : pf.second) {
                std::pair<size_t, size_t> affected_edge = affected_edge_and_probability.first;
                if (affected_edge != causal_edge) {
                    // Since edge weights are computed as std::log((1-p)/p), a probability of more than 0.5 for an
                    // error, would lead to a negatively weighted error. We do not support this (yet), and use a
                    // minimum of 0.5 as an implied probability for an edge to be reweighted.
                    double implied_probability_for_other_edge =
                        std::min(0.5, affected_edge_and_probability.second / marginal_probability);
                    double w = pm::to_weight_for_correlations(implied_probability_for_other_edge);
                    ImpliedWeightUnconverted implied{affected_edge.first, affected_edge.second, w};
                    edge.implied_weights_for_other_edges.push_back(implied);
                }
            }
        }
    }
}

#ifdef USE_THREADS
void pm::UserGraph::partition_nodes_by_round(const stim::DetectorErrorModel& dem) {
    // Query coordinates from stim. Map: det_id -> vector<double> of coords.
    std::set<uint64_t> all_dets;
    size_t num_nodes = nodes.size();
    for (uint64_t k = 0; k < num_nodes; ++k)
        all_dets.emplace_hint(all_dets.end(), k);
    std::map<uint64_t, std::vector<double>> coords_map = dem.get_detector_coordinates(all_dets);
    // M rounds per partition
    int M = config_parallel::M;
    node_part_id.resize(num_nodes);
    virtual_boundaries.clear();
    virtual_boundaries.push_back((std::vector<int>){});
    int p = 0;
    int vb = 1;
    double last_round = -1;
    bool p_or_vb = true; // true = p, false = vb
    int round_counter = 0;
    for (int n=0; n<num_nodes; ++n) {
        auto it = coords_map.find(n);
        if (it == coords_map.end()) {
            // No coordinates available; leave defaults.
            continue;
        }
        const auto& coors = it->second;
        if (coors.empty()) {
            // No coordinate data; skip annotation.
            throw std::invalid_argument("Detector node " + std::to_string(it->first) + " has no coords");
        }
        // Store round
        size_t round_coor = coors.size() - 1;
        nodes[n].round = coors[round_coor];
        // Update current p/vb if necessary
        if (coors[round_coor] > last_round) {
            round_counter++;
            if (round_counter == M) {
                ++p;
                p_or_vb = false; // boundary
            } else if (round_counter > M) {
                round_counter = 1;
                ++vb;
                virtual_boundaries.push_back((std::vector<int>){});
                p_or_vb = true; // partition
            }
            last_round = coors[round_coor];
            num_rounds++;
        }
        // Assign to p or vb
        if (p_or_vb){
            node_part_id[n] = p;
        } else {
            node_part_id[n] = -vb;
            virtual_boundaries[vb-1].push_back(n);
        }
    }
    if (!p_or_vb) { // ended on boundary, reassign to last partition
        --p;
        for (int n : virtual_boundaries[vb-1]) {
            node_part_id[n] = p;
        }
        virtual_boundaries.pop_back();
    }
    num_virtual_boundaries = virtual_boundaries.size();
    num_partitions = p+1;
}
#endif
#ifdef USE_SHMEM
// Assumes nodes are sorted in observable-major order: [obs0 by round][obs1 by round][seam nodes].
// Coords format per detector: [x, [y,] round, obs_id]  where obs_id < 0 encodes seam nodes.
//
// Design:
//   - obs0 creates virtual_boundaries slots; obs1+ ADDS nodes into the same slots (shared).
//   - node_part_id encoding: partition → global_p_offset+local_p; vb → -(global_vb_offset+slot_index+1).
//     obs0 vb0 → -1, obs1 vb0 → -(K_vb+1), etc. (globally unique).
//     node->vb still uses the shared slot index (0..K_vb-1) for GraphFlooder::is_active() gating.
//   - Cross-obs vbs get node->vb slot above K_vb-1. node_part_id increments globally normally
//   - num_virtual_boundaries = K_vb * obs_count + cross_obs_vb_count.
//
// Pending-vb pattern: when in vb, hold vb-round nodes in pending_vb; commit when next partition
// starts. If obs patch ends while pending, assign nodes to the last partition instead.
void pm::UserGraph::partition_nodes_by_obs_patch(const stim::DetectorErrorModel& dem) {
    std::set<uint64_t> all_dets;
    size_t num_nodes = nodes.size();
    for (uint64_t k = 0; k < num_nodes; ++k)
        all_dets.emplace_hint(all_dets.end(), k);
    std::map<uint64_t, std::vector<double>> coords_map = dem.get_detector_coordinates(all_dets);

    int M = config_parallel::M;
    node_part_id.resize(num_nodes);
    virtual_boundaries.clear();

    int local_p        = 0;   // partition index within current observable
    int local_vb_idx   = 0;   // vb slot index (0-based) within current observable
    int global_p_offset = 0;  // cumulative partitions from completed observables
    int global_vb_offset = 0; // cumulative vbs from completed observables

    std::vector<int> pending_vb;
    int pending_vb_slot = -1; // virtual_boundaries index for pending vb nodes

    int K_p  = -1;  // partitions per observable (set at first obs transition)
    int K_vb = -1;  // vb slots per observable

    int obs_count          = 0;
    int cross_obs_vb_count = 0;
    bool in_cross_obs      = false;
    double last_cross_obs_id = std::numeric_limits<double>::max();

    double last_round  = -1.0;
    bool p_or_vb       = true;  // true=partition, false=pending-vb
    int  round_counter = 0;

    if (DEBUG) std::cout << "Entering OBS 0\nStarting p0\n" << std::flush;

    for (int n = 0; n < (int)num_nodes; ++n) {
        auto it = coords_map.find((uint64_t)n);
        if (it == coords_map.end()) continue;
        const auto& coors = it->second;
        if (coors.empty())
            throw std::invalid_argument("Detector " + std::to_string(n) + " has no coords");

        double round      = coors.back();
        double obs_id_val = coors[coors.size() - 2];
        nodes[n].round         = round;
        nodes[n].observable_id = (int)floor(obs_id_val);
        nodes[n].x             = coors[0];
        if (coors.size() > 3) nodes[n].y = coors[1];

        if (round > last_round && !in_cross_obs) {
            if (!p_or_vb) { // end of vb
                // Commit pending vb: this is the first partition round after a vb round.
                if (obs_count == 0)
                    virtual_boundaries.push_back({});  // obs0: create new slot
                for (int ni : pending_vb)
                    virtual_boundaries[pending_vb_slot].push_back(ni);
                pending_vb.clear();
                ++local_vb_idx;
                // start next p
                round_counter = 1;
                p_or_vb = true;
                if (DEBUG) std::cout << "Starting p" << (global_p_offset + local_p) << "\n" << std::flush;
            } else { // end of p
                ++round_counter;
                if (round_counter == M) {
                    // Start pending vb
                    ++local_p;
                    pending_vb_slot = local_vb_idx;
                    p_or_vb = false;
                    if (DEBUG) std::cout << "Starting vb" << -(global_vb_offset + pending_vb_slot + 1) << "\n" << std::flush;
                }
            }
            if (obs_count == 0) ++num_rounds;
        } else if (round < last_round && !in_cross_obs) {
            // Round reset = obs transition (next observable or cross-obs seam section).
            // Record K_p and K_vb at the first obs transition.
            if (K_p == -1) {
                K_p  = p_or_vb ? local_p + 1 : local_p;
                K_vb = (int)virtual_boundaries.size();
            }
            // Cleanup: if obs ended on a pending vb, discard it and keep last partition.
            if (!p_or_vb) {
                int last_part = global_p_offset + local_p - 1;
                for (int ni : pending_vb)
                    node_part_id[ni] = last_part;
                pending_vb.clear();
                p_or_vb = true;
            }
            ++obs_count;

            if (obs_id_val < 0.0) {
                // Entering cross-observable (seam) nodes.
                in_cross_obs = true;
            } else {
                // Entering next observable's patch.
                // round_counter starts at 1 (not 0) because the first round of the new
                // observable is consumed by the transition detection. This keeps vbs aligned.
                global_p_offset  += K_p;
                global_vb_offset += K_vb;
                local_p        = 0;
                local_vb_idx   = 0;
                round_counter  = 1;
                if (DEBUG) std::cout << "Entering OBS " << obs_count << "\nStarting p" << global_p_offset << "\n" << std::flush;
            }
        }
        last_round = round;

        if (in_cross_obs && obs_id_val != last_cross_obs_id) { // new cross-observable node
            ++cross_obs_vb_count;
            virtual_boundaries.push_back({});
            last_cross_obs_id = obs_id_val;
            if (DEBUG) std::cout << "Starting cross-seam vb" << -(K_vb * obs_count + cross_obs_vb_count) << "\n" << std::flush;
        }

        // Assign node_part_id for this node.
        if (in_cross_obs) {
            int slot = K_vb + cross_obs_vb_count - 1;   // virtual_boundaries array index
            node_part_id[n] = -(K_vb * obs_count + cross_obs_vb_count);  // global vb id
            virtual_boundaries[slot].push_back(n);
        } else if (!p_or_vb) {
            node_part_id[n] = -(global_vb_offset + pending_vb_slot + 1);
            pending_vb.push_back(n);
        } else {
            node_part_id[n] = global_p_offset + local_p;
        }
    }

    // Final cleanup for the last observable.
    if (!in_cross_obs) {
        if (K_p == -1) {
            // Only one observable — no round-reset transition was detected.
            K_p  = p_or_vb ? local_p + 1 : local_p;
            K_vb = (int)virtual_boundaries.size();
        }
        if (!p_or_vb) {
            int last_part = global_p_offset + local_p - 1;
            for (int ni : pending_vb)
                node_part_id[ni] = last_part;
            pending_vb.clear();
        }
        ++obs_count;
    }

    num_partitions         = (size_t)K_p * obs_count;
    num_virtual_boundaries = (size_t)K_vb * obs_count + cross_obs_vb_count;
    num_obs_patches            = (size_t)obs_count;
    p_per_obs_patch = K_p;
    vb_per_obs_patch = K_vb;
}
// This function was generated by Claude Haiku 4.5
// void pm::UserGraph::reorder_nodes_by_observable() {
//     // Create a permutation that sorts nodes by observable_id (ascending), then by round
//     size_t num_nodes = nodes.size();
//     std::vector<size_t> permutation(num_nodes);
//     for (size_t i = 0; i < num_nodes; ++i) {
//         permutation[i] = i;
//     }
//     // Sort indices by observable_id, then by round
//     std::sort(permutation.begin(), permutation.end(), [this](size_t a, size_t b) {
//         int obs_a = nodes[a].observable_id;
//         int obs_b = nodes[b].observable_id;
//         if (obs_a != obs_b) {
//             if (obs_a < 0 && obs_b < 0)
//                 return obs_a > obs_b;  // less negative (higher) before more negative
//             if (obs_b < 0)
//                 return true;   // b is seam, a (regular) comes first
//             if (obs_a < 0)
//                 return false;  // a is seam, b (regular) comes first
//             return obs_a < obs_b;
//         }
//         return nodes[a].round < nodes[b].round;
//     });
//     // Build inverse permutation (old_id -> new_id)
//     std::vector<size_t> inverse_perm(num_nodes);
//     for (size_t new_id = 0; new_id < num_nodes; ++new_id) {
//         inverse_perm[permutation[new_id]] = new_id;
//     }
//     // Reorder nodes vector
//     std::vector<UserNode> reordered_nodes(num_nodes);
//     for (size_t new_id = 0; new_id < num_nodes; ++new_id) {
//         reordered_nodes[new_id] = nodes[permutation[new_id]];
//     }
//     nodes = std::move(reordered_nodes);
//     // Update all edges to use new node numbering.
//     // SIZE_MAX is the sentinel for a boundary edge endpoint — leave it unchanged.
//     for (auto& edge : edges) {
//         if (edge.node1 != SIZE_MAX) {
//             if (edge.node1 >= num_nodes)
//                 throw std::invalid_argument(
//                     "Edge references node ID out of range: (" + std::to_string(edge.node1) +
//                     ", " + std::to_string(edge.node2) + ") with num_nodes=" + std::to_string(num_nodes));
//             edge.node1 = inverse_perm[edge.node1];
//         }
//         if (edge.node2 != SIZE_MAX) {
//             if (edge.node2 >= num_nodes)
//                 throw std::invalid_argument(
//                     "Edge references node ID out of range: (" + std::to_string(edge.node1) +
//                     ", " + std::to_string(edge.node2) + ") with num_nodes=" + std::to_string(num_nodes));
//             edge.node2 = inverse_perm[edge.node2];
//         }
//     }
//     // Update boundary_nodes set to use new node numbering
//     std::set<size_t> new_boundary_nodes;
//     for (size_t old_id : boundary_nodes) {
//         new_boundary_nodes.insert(inverse_perm[old_id]);
//     }
//     boundary_nodes = new_boundary_nodes;
//     // Update node_part_id vector if it exists
//     if (!node_part_id.empty()) {
//         std::vector<int> reordered_part_id(num_nodes);
//         for (size_t new_id = 0; new_id < num_nodes; ++new_id) {
//             reordered_part_id[new_id] = node_part_id[permutation[new_id]];
//         }
//         node_part_id = std::move(reordered_part_id);
//     }
// }
#endif
// ===============
// std::set<long> pm::annotate_nodes_with_dem_coordinates(const stim::DetectorErrorModel& dem, pm::UserGraph& g) {
//     // Query coordinates from stim. Map: det_id -> vector<double> of coords.
//     std::set<uint64_t> all_dets;
//     size_t num_nodes = g.nodes.size();
//     for (uint64_t k = 0; k < num_nodes; ++k)
//         all_dets.emplace_hint(all_dets.end(), k);
//     std::map<uint64_t, std::vector<double>> coords_map = dem.get_detector_coordinates(all_dets);
//     // Annotate UserNodes
//     std::set<long> rounds;
//     std::set<long> x;
//     std::set<long> y;
//     for (size_t k = 0; k < g.nodes.size(); ++k) {
//         auto it = coords_map.find(k);
//         if (it == coords_map.end()) {
//             // No coordinates available; leave defaults.
//             continue;
//         }
//         const auto& coors = it->second;
//         if (coors.empty()) {
//             // No coordinate data; skip annotation.
//             throw std::invalid_argument("Detector node " + std::to_string(it->first) + " has no coords");
//         }
//         pm::UserNode& node = g.nodes[k];
//         node.has_coords = true;
//         // Store coordinates if available
//         if (coors.size() >= 2) {
//             node.pos_x = coors[0];
//             node.pos_y = coors[1];
//         } else {
//             // Only one coordinate present; leave x/y defaults
//             node.pos_x = coors[0];
//         }
//         // Use the last value as the round index
//         node.round = (long)lround(coors.back());
//         rounds.insert(node.round);
//         x.insert(node.pos_x);
//         y.insert(node.pos_y);
//     }
//     if (DEBUG)
//         std::cout << "X: " << *x.begin() << " to " << *x.rbegin() << std::endl
//                   << "Y: " << *y.begin() << " to " << *y.rbegin() << std::endl
//                   << "Z: " << *rounds.begin() << " to " << *rounds.rbegin() << std::endl;
//     return rounds;
// }

// void pm::partition_nodes_2d_vertical_split(pm::UserGraph& g, std::set<long> rounds) {
//     throw std::invalid_argument("partition_nodes_2d_vertical_split: not yet implemented");
// }

// ===============
