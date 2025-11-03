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

#ifdef ENABLE_FUSION
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
std::shared_ptr<pm::MatchingGraph> pm::UserGraph::to_matching_graph(pm::weight_int num_distinct_weights) {
    std::shared_ptr<MatchingGraph> matching_graph_ptr = std::make_shared<pm::MatchingGraph>(nodes.size(), _num_observables);
    pm::MatchingGraph& matching_graph = *matching_graph_ptr;
#else
pm::MatchingGraph pm::UserGraph::to_matching_graph(pm::weight_int num_distinct_weights) {
    pm::MatchingGraph matching_graph(nodes.size(), _num_observables);
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

#ifdef ENABLE_FUSION
// ===============
    // Propagate partition and virtual metadata
    // matching_graph.partition_of_node.resize(nodes.size(), -1);
    // matching_graph.is_virtual_node.resize(nodes.size(), false);
    matching_graph.num_partitions = num_partitions;
    for (size_t i = 0; i < nodes.size(); ++i) {
        // matching_graph.partition_of_node[i] = static_cast<int32_t>(nodes[i].partition);
        // matching_graph.is_virtual_node[i] = nodes[i].is_virtual;
        matching_graph.nodes[i].partition = nodes[i].partition;
        matching_graph.nodes[i].is_cross_partition = nodes[i].is_cross_partition;
    }
// ===============
#endif

#ifdef USE_THREADS
    return matching_graph_ptr;
#else
    return matching_graph;
#endif
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

// #ifdef ENABLE_FUSION
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
#ifdef ENABLE_FUSION
// ===============
    // Annotate user nodes with DEM-provided coordinates and virtual node infomation for partitioning.
    const auto& rounds = annotate_nodes_with_dem_coordinates(detector_error_model, user_graph);
    partition_nodes_by_round(user_graph, rounds);
// ===============
#endif
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

#ifdef ENABLE_FUSION
// ===============
std::set<long> pm::annotate_nodes_with_dem_coordinates(const stim::DetectorErrorModel& dem, pm::UserGraph& g) {
    // Query coordinates from stim. Map: det_id -> vector<double> of coords.
    std::set<uint64_t> all_dets;
    all_dets.clear();
    for (uint64_t k = 0; k < dem.count_detectors(); ++k)
        all_dets.insert(k);
    std::map<uint64_t, std::vector<double>> coords_map = dem.get_detector_coordinates(all_dets);
    // Annotate UserNodes
    std::set<long> rounds;
    std::set<long> x;
    std::set<long> y;
    for (size_t k = 0; k < g.nodes.size(); ++k) {
        auto it = coords_map.find(k);
        if (it == coords_map.end()) {
            // No coordinates available; leave defaults.
            continue;
        }
        const auto& coors = it->second;
        if (coors.empty()) {
            // No coordinate data; skip annotation.
            throw std::invalid_argument("Detector node " + std::to_string(it->first) + " has no coords");
        }
        pm::UserNode& node = g.nodes[k];
        node.has_coords = true;
        // Store coordinates if available
        if (coors.size() >= 2) {
            node.pos_x = coors[0];
            node.pos_y = coors[1];
        } else {
            // Only one coordinate present; leave x/y defaults
            node.pos_x = coors[0];
        }
        // Use the last value as the round index
        node.round = (long)lround(coors.back());
        rounds.insert(node.round);
        x.insert(node.pos_x);
        y.insert(node.pos_y);
    }
    if (DEBUG)
        std::cout << "X: " << *x.begin() << " to " << *x.rbegin() << std::endl
                  << "Y: " << *y.begin() << " to " << *y.rbegin() << std::endl
                  << "Z: " << *rounds.begin() << " to " << *rounds.rbegin() << std::endl;
    return rounds;
}


void pm::partition_nodes_by_round(pm::UserGraph& g, std::set<long> rounds) {
    // M rounds per partition
    int M = config_parallel::M;
    int num_p = rounds.size() / M;
    int num_tail = rounds.size() % M;
    int num_tail_per_p = num_tail / num_p;
    g.num_partitions = (long)(long(rounds.size()) / config_parallel::M); // num partitions
    std::vector<long> round_in_partition;
    long p = 0;
    int i = 0;
    int t = 0;
    for (long round : rounds) {
        round_in_partition.emplace_back(p);
        i++;
        if (i == M && num_tail > 0) {
            num_tail--;
            t++;
        } else if (i >= M && t < num_tail_per_p) {
            num_tail--;
            t++;
        } else if (i >= M) {
            p++;
            i = 0;
            t = 0;
        }
    }
    if (DEBUG)
        std::cout << "  M : " << M << "; num_rounds : " << rounds.size() << "; num_partitions : " << g.num_partitions << std::endl;
    // partition nodes
    for (auto& n : g.nodes) {
        if (n.has_coords)
            n.partition = round_in_partition[n.round];
    }
    // mark virtual nodes
    for (const auto& e : g.edges) {
        if (e.node1 == SIZE_MAX || e.node2 == SIZE_MAX) continue; // ignore boundary edges
        const auto& n1 = g.nodes[e.node1];
        const auto& n2 = g.nodes[e.node2];
        if (!n1.has_coords || !n2.has_coords) continue; // need rounds
        if (n1.partition == n2.partition) continue; // not a cross-partition edge
        else if (n1.partition > n2.partition) { // n1 in higher partition
            g.nodes[e.node1].is_cross_partition = true; // mark higher partition node as virtual
            // g.nodes[e.node2].virtual_neighbors.insert(std::make_tuple(e.node1, n1.round)); // save virtual neighbor in lower round node
        } else { // n2 in higher partition
            g.nodes[e.node2].is_cross_partition = true; // mark higher partition node as virtual
            // g.nodes[e.node1].virtual_neighbors.insert(std::make_tuple(e.node2, n2.round)); // save virtual neighbor in lower round node
        }
    }
}

void pm::partition_nodes_2d_vertical_split(pm::UserGraph& g, std::set<long> rounds) {
    g.num_partitions = 2;
    throw std::invalid_argument("partition_nodes_2d_vertical_split: not yet implemented");
}

// ===============
#endif
