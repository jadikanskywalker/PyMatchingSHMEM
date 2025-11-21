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

#ifdef USE_SHMEM
// ===============
#include <shmem.h>
#ifdef I
#undef I
#endif
#endif

#ifdef ENABLE_FUSION
#include "../config_parallel.h"
#include <fstream>
#include <filesystem>
#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"
#endif

#ifdef USE_THREADS
#include <atomic>
#include <omp.h>
// #include "pymatching/sparse_blossom/driver/decoding_task.h"
// ===============
#endif

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

#include "pymatching/sparse_blossom/driver/user_graph.h"
#include <unordered_set>

pm::ExtendedMatchingResult::ExtendedMatchingResult() : obs_crossed(), weight(0) {
}

bool pm::ExtendedMatchingResult::operator==(const ExtendedMatchingResult& rhs) const {
    return (obs_crossed == rhs.obs_crossed) && (weight == rhs.weight);
}

bool pm::ExtendedMatchingResult::operator!=(const ExtendedMatchingResult& rhs) const {
    return !(rhs == *this);
}

pm::ExtendedMatchingResult::ExtendedMatchingResult(size_t num_observables)
    : obs_crossed(num_observables, 0), weight(0) {
}

pm::ExtendedMatchingResult::ExtendedMatchingResult(std::vector<uint8_t> obs_crossed, total_weight_int weight)
    : obs_crossed(std::move(obs_crossed)), weight(weight) {
}

pm::ExtendedMatchingResult& pm::ExtendedMatchingResult::operator+=(const ExtendedMatchingResult& rhs) {
    assert(obs_crossed.size() == rhs.obs_crossed.size());
    for (size_t i = 0; i < obs_crossed.size(); i++) {
        obs_crossed[i] ^= rhs.obs_crossed[i];
    }
    weight += rhs.weight;
    return *this;
}

pm::ExtendedMatchingResult pm::ExtendedMatchingResult::operator+(const ExtendedMatchingResult& rhs) const {
    ExtendedMatchingResult copy = *this;
    copy += rhs;
    return copy;
}

void pm::fill_bit_vector_from_obs_mask(pm::obs_int obs_mask, uint8_t* obs_begin_ptr, size_t num_observables) {
    auto max_obs = sizeof(pm::obs_int) * 8;
    if (num_observables > max_obs)
        throw std::invalid_argument("Too many observables");
    for (size_t i = 0; i < num_observables; i++)
        *(obs_begin_ptr + i) ^= (obs_mask & ((pm::obs_int)1 << i)) >> i;
}

pm::obs_int pm::bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector) {
    auto num_observables = bit_vector.size();
    auto max_obs = sizeof(pm::obs_int) * 8;
    if (num_observables > max_obs)
        throw std::invalid_argument("Too many observables");
    pm::obs_int obs_mask = 0;
    for (size_t i = 0; i < num_observables; i++)
        obs_mask ^= bit_vector[i] << i;
    return obs_mask;
}

pm::Mwpm pm::detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included,
    bool enable_correlations) {
    auto user_graph =
        pm::detector_error_model_to_user_graph(detector_error_model, enable_correlations, num_distinct_weights);
    return user_graph.to_mwpm(num_distinct_weights, ensure_search_flooder_included);
}

void process_timeline_until_completion(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events
#ifdef ENABLE_FUSION
// ===============
    , bool draw_frames=false
// ===============
#endif
#ifdef USE_THREADS
    , bool parallel=false,
    int tid=-1
#endif
) {
    if (!mwpm.flooder.queue.empty()) {
        throw std::invalid_argument("!mwpm.flooder.queue.empty()");
    }
    mwpm.flooder.queue.cur_time = 0;

#ifdef USE_THREADS
    mwpm.unmatch_virtual_boundaries_between_partitions();
#endif

    if (mwpm.flooder.negative_weight_detection_events.empty()) {
        // if (DEBUG) std::cout << "  DEBUG: no negative detection events" << std::endl;
        // Just add detection events if graph has no negative weights
        for (auto& detection : detection_events) {
            if (detection >= mwpm.flooder.graph.nodes.size())
                throw std::invalid_argument(
                    "The detection event with index " + std::to_string(detection) +
                    " does not correspond to a node in the graph, which only has " +
                    std::to_string(mwpm.flooder.graph.nodes.size()) + " nodes.");
            
            if (detection + 1 > mwpm.flooder.graph.is_user_graph_boundary_node.size() ||
                !mwpm.flooder.graph.is_user_graph_boundary_node[detection])
                mwpm.create_detection_event(&mwpm.flooder.graph.nodes[detection]);
                // if (DEBUG) {
                //     std::cout << "  DEBUG: Thread " << tid << " added detection event at node " << &mwpm.flooder.graph.nodes[detection] << std::endl
                //               << "    is_partition: " << mwpm.flooder.graph.nodes[detection].partition << std::endl
                //               << "    is_cross_partition: " << mwpm.flooder.graph.nodes[detection].is_cross_partition << std::endl
                //               << "    is_active: " << mwpm.flooder.graph.nodes[detection].is_active << std::endl;
                // }
        }

    } else {
        // if (DEBUG) std::cout << "  DEBUG: negative detection events" << std::endl;
        // First mark nodes with negative weight detection events
        for (auto& det : mwpm.flooder.negative_weight_detection_events) {
            mwpm.flooder.graph.nodes[det].radius_of_arrival = 1;
        }

        // Now add detection events for unmarked nodes
        for (auto& detection : detection_events) {
            if (detection >= mwpm.flooder.graph.nodes.size())
                throw std::invalid_argument(
                    "Detection event index `" + std::to_string(detection) +
                    "` is larger than any detector node index in the graph.");
            if (!mwpm.flooder.graph.nodes[detection].radius_of_arrival) {
                if (detection + 1 > mwpm.flooder.graph.is_user_graph_boundary_node.size() ||
                    !mwpm.flooder.graph.is_user_graph_boundary_node[detection])
                    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[detection]);
            } else {
                // Unmark node
                mwpm.flooder.graph.nodes[detection].radius_of_arrival = 0;
            }
        }

        for (auto& det : mwpm.flooder.negative_weight_detection_events) {
            if (mwpm.flooder.graph.nodes[det].radius_of_arrival) {
                // Add a detection event if the node is still marked
                mwpm.flooder.graph.nodes[det].radius_of_arrival = 0;
#ifdef USE_THREADS // Only add detection events for active non-virtual nodes
// ===============
                if (mwpm.flooder.is_active(&mwpm.flooder.graph.nodes[det]))
// ===============
#endif
                mwpm.create_detection_event(&mwpm.flooder.graph.nodes[det]);
            }
        }
    }

#ifdef ENABLE_FUSION
// ===============
    int frame = 0;
    if (draw_frames)
        draw_frame(mwpm, pm::MwpmEvent::no_event(), frame, parallel, tid);
    frame++;
// ===============
#endif

    while (true) {
        auto event = mwpm.flooder.run_until_next_mwpm_notification();
#ifdef ENABLE_FUSION
// ===============
        if (draw_frames)
            draw_frame(mwpm, event, frame, parallel, tid);
        frame++;
// ===============
#endif
        if (event.event_type == pm::NO_EVENT)
            break;
        mwpm.process_event(event);
    }

    // If some alternating tree nodes remain, a perfect matching cannot be found
    if (mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size()
    ) {
#ifdef ENABLE_FUSION
// ===============
        if (DEBUG) {
            std::cout << "DEBUG: No perfect matching found" << std::endl;
            const std::unordered_set<pm::AltTreeNode*> freed_nodes(
                mwpm.node_arena.available.begin(), mwpm.node_arena.available.end());
            for (pm::AltTreeNode *alttreenode : mwpm.node_arena.allocated) {
                if (alttreenode == nullptr) continue;
                if (freed_nodes.find(alttreenode) != freed_nodes.end()) continue;
                if (alttreenode->inner_region) {
                    std::cout << "  inner GraphFillRegion: " << alttreenode->inner_region << std::endl;
                    for (auto &detector_node : alttreenode->inner_region->shell_area) {
                        auto index = detector_node - &*mwpm.flooder.graph.nodes.begin();
                        std::cout << "    " << index << " partition=" << mwpm.flooder.graph.nodes[index].partition
                        << " is_active=" << mwpm.flooder.graph.nodes[index].is_active << std::endl;
                    }
                }
                if (alttreenode->outer_region) {
                    std::cout << "  outer GraphFillRegion: " << alttreenode->outer_region << std::endl;
                    for (auto &detector_node : alttreenode->outer_region->shell_area) {
                        auto index = detector_node - &*mwpm.flooder.graph.nodes.begin();
                        std::cout << "    " << index << " partition=" << mwpm.flooder.graph.nodes[index].partition
                        << " is_active=" << mwpm.flooder.graph.nodes[index].is_active << std::endl;
                    }
                }
            }
        }
// ===============
#endif
        mwpm.reset();
        throw std::invalid_argument(
            "No perfect matching could be found. This likely means that the syndrome has odd "
            "parity in the support of a connected component without a boundary.");
    }
}

pm::MatchingResult shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    pm::MatchingResult res;
    for (auto& i : detection_events) {
        if (mwpm.flooder.graph.nodes[i].region_that_arrived)
            res += mwpm.shatter_blossom_and_extract_matches(mwpm.flooder.graph.nodes[i].region_that_arrived_top);
    }
    return res;
}

void shatter_blossoms_for_all_detection_events_and_extract_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    for (auto& i : detection_events) {
        if (mwpm.flooder.graph.nodes[i].region_that_arrived)
            mwpm.shatter_blossom_and_extract_match_edges(
                mwpm.flooder.graph.nodes[i].region_that_arrived_top, mwpm.flooder.match_edges);
    }
}

pm::MatchingResult pm::decode_detection_events_for_up_to_64_observables(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, bool edge_correlations) {
    if (edge_correlations) {
        // Edge correlations might also be called 2-pass matching. This is the slowest and
        // highest-accuracy reweighting rule for correlated decoding in which we first decode to edges,
        // and then reweight the associated edges conditioned on the assumption that an error occurred
        // at that edge.
        std::vector<int64_t> edges;
        decode_detection_events_to_edges(mwpm, detection_events, edges);
        mwpm.flooder.graph.reweight_for_edges(edges);
        mwpm.search_flooder.graph.reweight_for_edges(edges);
    }

    process_timeline_until_completion(mwpm, detection_events);
    auto res = shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(mwpm, detection_events);
    if (!mwpm.flooder.negative_weight_detection_events.empty())
        res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
            mwpm, mwpm.flooder.negative_weight_detection_events);
    res.obs_mask ^= mwpm.flooder.negative_weight_obs_mask;
    res.weight += mwpm.flooder.negative_weight_sum;

    if (edge_correlations) {
        mwpm.flooder.graph.undo_reweights();
        mwpm.search_flooder.graph.undo_reweights();
    }

    return res;
}

void pm::decode_detection_events(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight,
    bool edge_correlations
#ifdef ENABLE_FUSION
// ===============
    , int shot,
    bool draw_frames
// ===============
#endif
    ) {
    if (edge_correlations) {
        // Edge correlations might also be called 2-pass matching. This is the slowest and
        // highest-accuracy reweighting rule for correlated decoding in which we first decode to edges,
        // and then reweight the associated edges conditioned on the assumption that an error occurred
        // at that edge.
        std::vector<int64_t> edges;
        decode_detection_events_to_edges(mwpm, detection_events, edges);
        mwpm.flooder.graph.reweight_for_edges(edges);
        mwpm.search_flooder.graph.reweight_for_edges(edges);
    }

#ifdef USE_SHMEM
// ===============
    int mype = shmem_my_pe();
#else
    int mype = 0;
// ===============
#endif

#ifdef ENABLE_FUSION
// ===============
    if (DEBUG && mype==0)
        output_detection_events(mwpm, detection_events, shot, false);
    if (draw_frames)
        std::filesystem::create_directory("out_serial/frames/" + std::to_string(shot));
    // ===============
#endif

#ifdef USE_THREADS
// ===============
    // for (auto &node : mwpm.flooder.graph.nodes)
    //     node.is_active = 0;
// ===============
#endif

    size_t num_observables = mwpm.flooder.graph.num_observables;
    process_timeline_until_completion(mwpm, detection_events
#ifdef ENABLE_FUSION
// ===============
        , draw_frames
// ===============
#endif
    );

#ifdef ENABLE_FUSION
// ===============
    if (DEBUG && mype==0)
        output_solution_state(mwpm, detection_events, shot, false);
// ===============
#endif

    if (num_observables > sizeof(pm::obs_int) * 8) {
        mwpm.flooder.match_edges.clear();
        shatter_blossoms_for_all_detection_events_and_extract_match_edges(mwpm, detection_events);
        if (!mwpm.flooder.negative_weight_detection_events.empty())
            shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                mwpm, mwpm.flooder.negative_weight_detection_events);
        mwpm.extract_paths_from_match_edges(mwpm.flooder.match_edges, obs_begin_ptr, weight);

        // XOR negative weight observables
        for (auto& obs : mwpm.flooder.negative_weight_observables)
            *(obs_begin_ptr + obs) ^= 1;
        // Add negative weight sum to blossom solution weight
        weight += mwpm.flooder.negative_weight_sum;

    } else {
        pm::MatchingResult bit_packed_res =
            shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(mwpm, detection_events);
        if (!mwpm.flooder.negative_weight_detection_events.empty())
            bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                mwpm, mwpm.flooder.negative_weight_detection_events);
        // XOR in negative weight observable mask
        bit_packed_res.obs_mask ^= mwpm.flooder.negative_weight_obs_mask;
        // Translate observable mask into bit vector
        fill_bit_vector_from_obs_mask(bit_packed_res.obs_mask, obs_begin_ptr, num_observables);
        // Add negative weight sum to blossom solution weight
        weight = bit_packed_res.weight + mwpm.flooder.negative_weight_sum;
    }

    if (edge_correlations) {
        mwpm.flooder.graph.undo_reweights();
        mwpm.search_flooder.graph.undo_reweights();
    }
}

void pm::decode_detection_events_to_match_edges(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    if (mwpm.flooder.negative_weight_sum != 0)
        throw std::invalid_argument(
            "Decoding to matched detection events not supported for graphs containing edges with negative weights.");
    process_timeline_until_completion(mwpm, detection_events);
    mwpm.flooder.match_edges.clear();
    shatter_blossoms_for_all_detection_events_and_extract_match_edges(mwpm, detection_events);
}

void flip_edge(const pm::SearchGraphEdge& edge) {
    edge.detector_node->neighbor_markers[edge.neighbor_index] ^= pm::FLIPPED;
    auto neighbor = edge.detector_node->neighbors[edge.neighbor_index];
    if (neighbor) {
        auto idx_from_neighbor = neighbor->index_of_neighbor(edge.detector_node);
        neighbor->neighbor_markers[idx_from_neighbor] ^= pm::FLIPPED;
    }
}

void pm::decode_detection_events_to_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges) {
    if (mwpm.flooder.graph.nodes.size() != mwpm.search_flooder.graph.nodes.size()) {
        throw std::invalid_argument(
            "Mwpm object does not contain search flooder, which is required to decode to edges.");
    }
    process_timeline_until_completion(mwpm, detection_events);
    mwpm.flooder.match_edges.clear();
    shatter_blossoms_for_all_detection_events_and_extract_match_edges(mwpm, detection_events);
    if (!mwpm.flooder.negative_weight_detection_events.empty())
        shatter_blossoms_for_all_detection_events_and_extract_match_edges(
            mwpm, mwpm.flooder.negative_weight_detection_events);
    // Flip edges with negative weights and add to edges vector.
    for (const auto& neg_node_pair : mwpm.search_flooder.graph.negative_weight_edges) {
        auto node1_ptr = &mwpm.search_flooder.graph.nodes[neg_node_pair.first];
        auto node2_ptr =
            neg_node_pair.second != SIZE_MAX ? &mwpm.search_flooder.graph.nodes[neg_node_pair.second] : nullptr;
        SearchGraphEdge neg_edge = {node1_ptr, node1_ptr->index_of_neighbor(node2_ptr)};
        flip_edge(neg_edge);
        int64_t node1 = neg_edge.detector_node - &mwpm.search_flooder.graph.nodes[0];
        int64_t node2 = node2_ptr ? node2_ptr - &mwpm.search_flooder.graph.nodes[0] : -1;
        edges.push_back(node1);
        edges.push_back(node2);
    }
    // Flip edges along a shortest path between matched detection events and add to edges vector
    for (const auto& match_edge : mwpm.flooder.match_edges) {
        size_t node_from = match_edge.loc_from - &mwpm.flooder.graph.nodes[0];
        size_t node_to = match_edge.loc_to ? match_edge.loc_to - &mwpm.flooder.graph.nodes[0] : SIZE_MAX;
        mwpm.search_flooder.iter_edges_on_shortest_path_from_middle(node_from, node_to, [&](const SearchGraphEdge& e) {
            flip_edge(e);
            int64_t node1 = e.detector_node - &mwpm.search_flooder.graph.nodes[0];
            auto node2_ptr = e.detector_node->neighbors[e.neighbor_index];
            int64_t node2 = node2_ptr ? node2_ptr - &mwpm.search_flooder.graph.nodes[0] : -1;
            edges.push_back(node1);
            edges.push_back(node2);
        });
    }
    // Remove any edges in the edges vector that are no longer flipped. Also unflip edges.
    for (size_t i = 0; i < edges.size() / 2;) {
        int64_t u = edges[2 * i];
        int64_t v = edges[2 * i + 1];
        auto& u_node = mwpm.search_flooder.graph.nodes[u];
        size_t idx;
        if (v == -1) {
            idx = 0;
        } else {
            auto v_ptr = &mwpm.search_flooder.graph.nodes[v];
            idx = u_node.index_of_neighbor(v_ptr);
        }
        if (!(u_node.neighbor_markers[idx] & pm::FLIPPED)) {
            edges[2 * i] = edges[edges.size() - 2];
            edges[2 * i + 1] = edges[edges.size() - 1];
            edges.resize(edges.size() - 2);
        } else {
            flip_edge({&u_node, idx});
            i++;
        }
    }
}

void pm::decode_detection_events_to_edges_with_edge_correlations(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges) {
    decode_detection_events_to_edges(mwpm, detection_events, edges);
    mwpm.flooder.graph.reweight_for_edges(edges);
    mwpm.search_flooder.graph.reweight_for_edges(edges);
    edges.clear();

    decode_detection_events_to_edges(mwpm, detection_events, edges);

    mwpm.flooder.graph.undo_reweights();
    mwpm.search_flooder.graph.undo_reweights();
}

#ifdef ENABLE_FUSION
// ===============
void pm::setup_output_dirs(bool draw_frames=true, bool parallel) {
    std::string out_dir = parallel ? "out_parallel" : "out_serial";
    std::filesystem::create_directory(out_dir);
    if (draw_frames)
            std::filesystem::create_directory(out_dir + "/frames");
}

void pm::output_detector_nodes(pm::Mwpm& mwpm, bool parallel) {
    static int count = 0;
    auto &graph = mwpm.flooder.graph;
    std::string out_dir = parallel ? "out_parallel/" : "out_serial/";
    std::ofstream out(out_dir + "graph_" + std::to_string(count) + ".out");
    count++;
    for (auto &node : graph.nodes) {
        out << "node: " << &(node) << std::endl
            << "  partition : " << node.partition << std::endl;
        if (parallel) 
            out << "  is_active: " << node.is_active << std::endl
                << "  is_cross_partition: " << node.is_cross_partition << std::endl;
        out << "  neighbors : " << std::endl;
        for (int i=0; i<node.neighbors.size(); i++) {
            if (node.neighbors[i] == nullptr)
                continue;
            out << "    neighbor: " << node.neighbors[i] << std::endl
                << "      partition : " << node.neighbors[i]->partition << std::endl;
            if (parallel)
                out << "      is_active: " << node.neighbors[i]->is_active << std::endl
                    << "      is_cross_partition: " << node.neighbors[i]->is_cross_partition << std::endl;
            out << "      weight    : " << node.neighbor_weights[i] << std::endl
                << "      logobs    : " << node.neighbor_observables[i] << std::endl;
            }
    }
    out.close();
}

void pm::output_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, int shot, bool parallel) {
    std::string out_dir = parallel ? "out_parallel/" : "out_serial/";
    std::ofstream out(out_dir + "detection_events_" + std::to_string(shot) + ".out");
    for (uint64_t det : detection_events) {
        out << &mwpm.flooder.graph.nodes[det];
        if (parallel)
            out << " partition=" << mwpm.flooder.graph.nodes[det].partition
                << "  is_cross_partition=" << mwpm.flooder.graph.nodes[det].is_cross_partition;
        out << std::endl;
    }
    out.close();
}

inline bool node_has_detection_event(size_t node_i, const std::vector<uint64_t>& detection_events) {
    auto it = std::find(detection_events.begin(), detection_events.end(), (uint64_t)node_i);
    if (it != detection_events.end())
        return true;
    return false;
}
#endif

#ifdef USE_THREADS
void pm::output_solution_state(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, int shot, bool parallel) {
    std::string out_dir = parallel ? "out_parallel/" : "out_serial/";
    std::string out_name = out_dir + "solution_" + std::to_string(shot);
    if (mwpm.task) {
        out_name += "_" + std::to_string(mwpm.task->p) + (mwpm.task->is_fusion ? std::to_string(mwpm.task->pv) : "");
    }
    std::ofstream out(out_name + ".out");
    // Skip regions that have been freed back to the arena (available). Only walk live objects.
    const std::unordered_set<pm::GraphFillRegion*> freed_regions(
        mwpm.flooder.region_arena.available.begin(), mwpm.flooder.region_arena.available.end());
    auto nodes_begin = &mwpm.flooder.graph.nodes[0];
    auto nodes_end = nodes_begin + mwpm.flooder.graph.nodes.size();
    for (auto region : mwpm.flooder.region_arena.allocated) {
        if (region == nullptr) continue;
        if (freed_regions.find(region) != freed_regions.end()) continue;
        out << "  region: " << region << std::endl
                  << "    radius: " << region->radius << std::endl
                  << "    match.region: " << region->match.region << std::endl
                  << "    match.edge.loc_to: " << region->match.edge.loc_to;
        if (!region->match.edge.loc_to) {
            out << "  -  BOUNDARY  -  NO DETECTION EVENT" << std::endl;
        } else {
            auto ptr = region->match.edge.loc_to;
            if (ptr < nodes_begin || ptr >= nodes_end) {
                out << "  -  INVALID DETECTOR POINTER" << std::endl;
            } else {
                if (parallel && !ptr->is_active) {
                    if (ptr->is_cross_partition)
                        out << "  -  CROSS_PARTITION: partition=" << ptr->partition;
                    else
                        out << "  -  OUT_OF_PARTITION: partition=" << ptr->partition;
                }
                size_t i_to = static_cast<size_t>(ptr - nodes_begin);
                if (&mwpm.flooder.graph.nodes[i_to] == ptr && node_has_detection_event(i_to, detection_events))
                    out << "  -  HAS DETECTION EVENT" << std::endl;
                else
                    out << "  -  NO DETECTION EVENT" << std::endl;
            }
        }

        out << "    match.edge.loc_from: " << region->match.edge.loc_from;
        if (!region->match.edge.loc_from) {
            out << "  -  BOUNDARY  -  NO DETECTION EVENT" << std::endl;
        } else {
            auto ptr = region->match.edge.loc_from;
            if (ptr < nodes_begin || ptr >= nodes_end) {
                out << "  -  INVALID DETECTOR POINTER" << std::endl;
            } else {
                if (parallel && !ptr->is_active) {
                    if (ptr->is_cross_partition)
                        out << "  -  CROSS_PARTITION: partition=" << ptr->partition;
                    else
                        out << "  -  OUT_OF_PARTITION: partition=" << ptr->partition;
                }
                size_t i_from = static_cast<size_t>(ptr - nodes_begin);
                if (&mwpm.flooder.graph.nodes[i_from] == ptr && node_has_detection_event(i_from, detection_events))
                    out << "  -  HAS DETECTION EVENT" << std::endl;
                else
                    out << "  -  NO DETECTION EVENT" << std::endl;
            }
        }
        if (region->shell_area.empty())
            out << "    WARNING: shell_area EMPTY" << std::endl;
        else {
            out << "    shell_area:" << std::endl;
            for (DetectorNode* node : region->shell_area) {
                out << "        " << node;
                if (node == nullptr) continue;
                size_t i = static_cast<size_t>(node - &mwpm.flooder.graph.nodes[0]);
                if (i >= mwpm.flooder.graph.nodes.size())
                    out << "  ERROR: shell_area node out of bounds" << std::endl;
                else if (&mwpm.flooder.graph.nodes[i] != node)
                    out << "  ERROR: could not reference shell_area node" << std::endl;
                else if (node_has_detection_event(i, detection_events))
                    out << "  -  HAS DETECTION EVENT" << std::endl;
                else
                    out << "  -  NO DETECTION EVENT" << std::endl;
            }
        }
        out << "    alt_tree_node: " << region->alt_tree_node << std::endl
                  << "    blossom_parent: " << region->blossom_parent << std::endl
                  << "    blossom_parent_top: " << region->blossom_parent_top << std::endl
                  << "    blossom_children.size(): " << region->blossom_children.size() << std::endl;
    }
    out << "DEBUG: partition 0 detector nodes: " << std::endl;
    for (auto node  = mwpm.flooder.graph.nodes.begin(); node != mwpm.flooder.graph.nodes.end(); ++node) {
        // if (node->partition != 0 && !node->is_virtual) continue;
        out << "  node " << &(*node);
        if (parallel && node->is_active)
            out << "  -  ACTIVE";
        else if (parallel)
            out << "  -  NOT ACTIVE";
        if (parallel && node->is_cross_partition)
            out << "  -  CROSS_PARTITION";
        else if (parallel)
            out << "  -  NOT CROSS_PARTITION";
        if (node_has_detection_event((size_t)(&(*node) - &mwpm.flooder.graph.nodes[0]), detection_events))
            out << "  -  HAS DETECTION EVENT" << std::endl;
        else
            out << "  -  NO DETECTION EVENT" << std::endl;
        if (parallel)
            out << "    partition: " << node->partition << std::endl;
        out << "    region_that_arrived: " << node->region_that_arrived << std::endl
            << "    region_that_arrived_top: " << node->region_that_arrived_top << std::endl;
        if (node->region_that_arrived != node->region_that_arrived_top)
                out << "    NOTE: node.region_that_arrived != region_that_arrived_top" << std::endl;
        out << "    wrapped_radius_cached: " << node->wrapped_radius_cached << std::endl
            << "    reached_from_source: " << node->reached_from_source << std::endl;
        if (node->reached_from_source != 0 && node->reached_from_source != &(*node)) {
            out << "    NOTE: node.reached_from_source != node" << std::endl;
        }
        out << "    observables_crossed_from_source: " << node->observables_crossed_from_source << std::endl
            << "    radius_of_arrival: " << node->radius_of_arrival << std::endl;
    }
    out.close();
}

void pm::draw_frame(pm::Mwpm& mwpm, pm::MwpmEvent ev, int frame_number, bool parallel, int tid) {
    std::string out_dir = parallel ? "out_parallel/" : "out_serial/";
    std::string out_name = out_dir + "frames/" + std::to_string(mwpm.flooder.current_shot) + "/";
    if (tid > -1) {
        out_name += "t" + std::to_string(tid) + "/";
    } 
    out_name +=  "frame"; 
    if (mwpm.task) {
        out_name += "_" + std::to_string(mwpm.task->p) + (mwpm.task->is_fusion ? std::to_string(mwpm.task->pv) : "");
    }
    out_name += "_" + std::to_string(frame_number);
    // Draw decoder state
    std::ofstream out_file(out_name + ".svg");
    if (out_file.fail()) {
        throw std::invalid_argument("Failed to open " + out_name + ".svg to write frame.");
    }
    pm::write_decoder_state_as_svg(mwpm.coords.first, mwpm.coords.second, mwpm, ev, out_file);
    out_file.close();
}

// Define the global per-thread solvers/queues declared in the header.
// std::vector<Task> pm::tasks;
// std::vector<int> pm::partitions_task_id;
std::vector<Task> pm::tasks;
std::vector<long> pm::step_sizes;
std::vector<std::shared_ptr<pm::Mwpm>> pm::solvers;
// WorkStealingDeque pm::task_deque;

// std::vector<std::queue<int>> pm::partition_task_queues;
// std::vector<std::deque<int>> pm::fusion_task_deques;

void pm::init_tasks(int num_partitions) {
    if (DEBUG) {
        std::cout << "DEBUG: initializing tasks" << std::endl;
    }
    // Build a full fusion tree: at most ~2*N tasks. Reserve to keep element addresses stable.
    tasks.clear();
    step_sizes.clear();
    tasks.reserve(static_cast<size_t>(2 * num_partitions - 1));
    // Add tasks for each partitions
    int task_id;
    for (task_id = 0; task_id<num_partitions; ++task_id) {
        tasks.emplace_back(task_id);
        // task_deque.push(new_task);
    }
    step_sizes.emplace_back(num_partitions);

    // Save odd trailing task
    std::vector<int> tails_idx;
    if (num_partitions % 2) {
        tails_idx.emplace_back(num_partitions - 1);
    }
    // Add each step of fusions
    int last_step_starts = 0;
    int num_fusions_this_step = num_partitions / 2 ;
    while (num_fusions_this_step > 0) {
        for (int i = 0; i < num_fusions_this_step; ++i) {
            Task* left_task = &tasks[last_step_starts + 2 * i];
            Task* right_task = &tasks[last_step_starts + 2 * i + 1];
            tasks.emplace_back(left_task, right_task);
            ++task_id;
        }
        if (tails_idx.size() == 2) {
            // Fuse tail pair
            tasks.emplace_back(&tasks[tails_idx[1]], &tasks[tails_idx[0]]);
            ++task_id;
            // task_deque.push(new_task);
            tails_idx.clear();
            num_fusions_this_step++;
        }
        last_step_starts += *step_sizes.rbegin();
        if (num_fusions_this_step % 2) {
            // Add tail
            tails_idx.emplace_back(last_step_starts+num_fusions_this_step-1);
        }
        step_sizes.emplace_back(num_fusions_this_step);
        num_fusions_this_step /= 2;
    }
    // Fuse remaining tail
    if (tails_idx.size() > 1) {
        tasks.emplace_back(&tasks[tails_idx[1]], &tasks[tails_idx[0]]);
        // task_deque.push(new_task);
        step_sizes.emplace_back(1);
    }
    if (DEBUG) {
        std::cout << "DEBUG:" << std::endl;
        for (Task& t : tasks) {
            std::cout << "--p: " << t.p << std::endl
                      << "  pv: " << t.pv << std::endl
                      << "  left_child: " << t.left_child << std::endl
                      << "  right_child: " << t.right_child << std::endl
                      << "  parent: " << t.parent << std::endl
                      << "  is_fusion: " << t.is_fusion << std::endl;
        }
    }
}

// Build per-thread solvers directly from a DEM by constructing a UserGraph once
// and producing independent Mwpm instances via to_mwpm for each thread.
// (Assumes first mwpm has already been built)
void pm::build_thread_solvers(
    pm::Mwpm& mwpm,
    bool ensure_search_flooder_included,
    bool enable_correlations,
    int num_threads) {
    if (num_threads < 1)
        return;
    if (ensure_search_flooder_included || enable_correlations)
        throw std::invalid_argument("Correlations and SearchFlooder are not yet supported with threads");
    pm::solvers.clear();
    pm::solvers.reserve(static_cast<size_t>(num_threads));
    for (int t = 0; t < num_threads; ++t) {
        // Each solver shares the same MatchingGraph via shared_ptr.
        pm::solvers.emplace_back(std::make_shared<pm::Mwpm>(pm::GraphFlooder(mwpm.flooder.graph_ptr)));
        pm::solvers[t]->flooder.sync_negative_weight_observables_and_detection_events();
    }
}

inline std::vector<uint64_t> create_detection_events_mask(pm::Mwpm& solver, const std::vector<uint64_t>& detection_events, Task &t) {
    std::vector<uint64_t> detection_events_mask;
    if (!t.is_fusion) {
        for (auto& det : detection_events) {
            if (t.p == solver.flooder.graph.nodes[det].partition
                && !solver.flooder.graph.nodes[det].is_cross_partition)
                detection_events_mask.push_back(det);
        }
    } else {
        for (auto& det : detection_events) {
            if (t.pv == solver.flooder.graph.nodes[det].partition
                && solver.flooder.graph.nodes[det].is_cross_partition)
                detection_events_mask.push_back(det);
        }     
    }
    return detection_events_mask;
}

inline void solve_task(int tid, int shot, const std::vector<uint64_t>& detection_events, int draw_frames, Task* task) {
    task->setup_regions();
    auto& solver = *pm::solvers[tid];
    solver.prepare_for_task(tid, task);
    auto detection_events_mask = create_detection_events_mask(solver, detection_events, *task);
    if (DEBUG && omp_get_thread_num()==0) {
        output_detector_nodes(solver, true);
    }
    process_timeline_until_completion(solver, detection_events_mask, draw_frames, true, tid);
    if (DEBUG) {
        output_solution_state(solver, detection_events, shot, true);
    }
}

inline void reset_tasks_and_solvers(int tid, int num_threads) {
    // Free solver arena regions
    auto& mwpm = *pm::solvers[tid];
    const std::unordered_set<pm::GraphFillRegion*> freed_regions(
        mwpm.flooder.region_arena.available.begin(), mwpm.flooder.region_arena.available.end());
    for (auto region : mwpm.flooder.region_arena.allocated) {
        if (region == nullptr) continue;
        if (freed_regions.find(region) != freed_regions.end()) continue;
        mwpm.flooder.region_arena.del(region);
    }
    int chunk_size = pm::tasks.size() / num_threads;
    int i;
    for (i = tid*chunk_size; i < tid*chunk_size + chunk_size; ++i) {
        pm::tasks[i].reset();
    }
    if (tid == num_threads-1) {
        for (; i < pm::tasks.size(); ++i) {
            pm::tasks[i].reset();
        }
    }
}
#endif

#ifdef ENABLE_FUSION
void pm::decode_detection_events_in_parallel(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight,
    bool edge_correlations,
    int shot,
    bool draw_frames
) {
    if (edge_correlations) {
        throw std::invalid_argument("Edge correlations are not yet implemented in parallel.");
    }

#ifdef USE_SHMEM
// ===============
    // TODO: Simple Scheduling Algorithm
    int mype = shmem_my_pe();
#else
    int mype = 0;
// ===============
#endif

#ifdef ENABLE_FUSION
// ===============
    if (DEBUG && mype == 0)
        output_detection_events(mwpm, detection_events, shot, true);
    if (draw_frames)
        std::filesystem::create_directory("out_parallel/frames/" + std::to_string(shot));
// ===============
#endif

    size_t num_observables = mwpm.flooder.graph.num_observables;

#if defined(USE_THREADS) && !defined(USE_SHMEM)
// ===============
    // std::cout << "Starting shot " << shot << std::endl;
    int num_partitions = mwpm.flooder.graph.num_partitions;
    mwpm.flooder.graph.reset_active_status_for_all_nodes(); // CAN MAKE THIS MORE EFFICIENT
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        try {
            reset_tasks_and_solvers(tid, num_threads);

            if (draw_frames)
                std::filesystem::create_directory("out_parallel/frames/" + std::to_string(shot) + "/t" + std::to_string(tid));

            #pragma omp barrier

            std::queue<int> partition_tasks;
            std::queue<Task *> fusion_tasks;
            for (int i = tid; i < num_partitions; i += num_threads) {
                partition_tasks.push(i);
            }
            while (true) {
                Task *t = nullptr;
                if (fusion_tasks.empty()) {
                    if (partition_tasks.empty()) {
                        break;
                    } else {
                        t = &tasks[partition_tasks.front()];
                        partition_tasks.pop();
                    }
                } else {
                    t = fusion_tasks.front();
                    fusion_tasks.pop();
                }
                if (t->try_to_steal()) {
                    // Got task
                    // std::cout << "Thread " << tid << " solving " << t->p << (t->is_fusion ? std::to_string(t->pv) : "") << std::endl;
                    solve_task(tid, shot, detection_events, draw_frames, t);
                    t->mark_solved();
                    // Try to steal parent
                    if (t->parent) {
                        fusion_tasks.push(t->parent);
                    }
                }
            }
        } catch (const std::exception &e) {
            #pragma omp critical
            {
                std::cerr << "ERROR: Thread " << tid << " caught exception: " << e.what() << std::endl;
            }
        } catch (...) {
            #pragma omp critical
            {
                std::cerr << "ERROR: Thread " << tid << " caught unknown exception." << std::endl;
            }
        }
    }
#endif

    if (num_observables > sizeof(pm::obs_int) * 8) {
        mwpm.flooder.match_edges.clear();
        shatter_blossoms_for_all_detection_events_and_extract_match_edges(mwpm, detection_events);
        if (!mwpm.flooder.negative_weight_detection_events.empty())
            shatter_blossoms_for_all_detection_events_and_extract_match_edges(
                mwpm, mwpm.flooder.negative_weight_detection_events);
        mwpm.extract_paths_from_match_edges(mwpm.flooder.match_edges, obs_begin_ptr, weight);

        // XOR negative weight observables
        for (auto& obs : mwpm.flooder.negative_weight_observables)
            *(obs_begin_ptr + obs) ^= 1;
        // Add negative weight sum to blossom solution weight
        weight += mwpm.flooder.negative_weight_sum;
    } else {
        pm::MatchingResult bit_packed_res =
            shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(mwpm, detection_events);
        if (!mwpm.flooder.negative_weight_detection_events.empty())
            bit_packed_res += shatter_blossoms_for_all_detection_events_and_extract_obs_mask_and_weight(
                mwpm, mwpm.flooder.negative_weight_detection_events);
        // XOR in negative weight observable mask
        bit_packed_res.obs_mask ^= mwpm.flooder.negative_weight_obs_mask;
        // Translate observable mask into bit vector
        fill_bit_vector_from_obs_mask(bit_packed_res.obs_mask, obs_begin_ptr, num_observables);
        // Add negative weight sum to blossom solution weight
        weight = bit_packed_res.weight + mwpm.flooder.negative_weight_sum;
    }
}
#endif
