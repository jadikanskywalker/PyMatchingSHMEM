#include "pymatching/sparse_blossom/driver/helpers/graph_cache.h"

#include <cstdio>
#include <list>
#include <vector>

namespace {

constexpr uint64_t GRAPH_CACHE_MAGIC = 0x484341474d5024ULL;
constexpr uint32_t GRAPH_CACHE_VERSION = 2;  // bumped: UserNode::vb removed (dead field, never used)

template <typename T>
void write_pod(FILE* f, const T& value) {
    if (fwrite(&value, sizeof(T), 1, f) != 1)
        throw std::runtime_error("graph_cache: failed to write to file");
}

template <typename T>
void read_pod(FILE* f, T& value) {
    if (fread(&value, sizeof(T), 1, f) != 1)
        throw std::runtime_error("graph_cache: failed to read from file (truncated?)");
}

void write_u64_vector(FILE* f, const std::vector<uint64_t>& v) {
    uint64_t n = v.size();
    write_pod(f, n);
    if (n > 0 && fwrite(v.data(), sizeof(uint64_t), n, f) != n)
        throw std::runtime_error("graph_cache: failed to write vector body");
}

std::vector<uint64_t> read_u64_vector(FILE* f) {
    uint64_t n;
    read_pod(f, n);
    std::vector<uint64_t> v(n);
    if (n > 0 && fread(v.data(), sizeof(uint64_t), n, f) != n)
        throw std::runtime_error("graph_cache: failed to read vector body (truncated?)");
    return v;
}

}  // namespace

void pm::write_user_graph_cache(
    pm::UserGraph& user_graph,
    const std::string& path,
    bool enable_correlations,
    int64_t rounds_per_partition,
    int division_strategy
) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f)
        throw std::runtime_error("graph_cache: failed to open " + path + " for writing");

    try {
        // --- Header ---
        write_pod(f, GRAPH_CACHE_MAGIC);
        write_pod(f, GRAPH_CACHE_VERSION);

        uint8_t build_mode = 0;
        // build_mode = 1;
        build_mode = 2;

        write_pod(f, build_mode);

        uint8_t enable_corr_u8 = enable_correlations ? 1 : 0;
        write_pod(f, enable_corr_u8);

        uint8_t loaded_wo_corr_u8 = user_graph.loaded_from_dem_without_correlations ? 1 : 0;
        write_pod(f, loaded_wo_corr_u8);

        uint64_t num_nodes = user_graph.get_num_nodes();
        uint64_t num_edges = user_graph.get_num_edges();
        uint64_t num_observables = user_graph.get_num_observables();
        uint64_t num_boundary_nodes = user_graph.boundary_nodes.size();
        write_pod(f, num_nodes);
        write_pod(f, num_edges);
        write_pod(f, num_observables);
        write_pod(f, num_boundary_nodes);

        write_pod(f, rounds_per_partition);
        uint64_t cache_num_partitions = user_graph.num_partitions;
        uint64_t cache_num_rounds = user_graph.num_rounds;
        uint64_t cache_num_virtual_boundaries = user_graph.num_virtual_boundaries;
        write_pod(f, cache_num_partitions);
        write_pod(f, cache_num_rounds);
        write_pod(f, cache_num_virtual_boundaries);

        int32_t cache_division_strategy = division_strategy;
        write_pod(f, cache_division_strategy);
        uint64_t cache_num_obs_patches = user_graph.num_obs_patches;
        uint64_t cache_p_per_obs_patch = user_graph.p_per_obs_patch;
        uint64_t cache_vb_per_obs_patch = user_graph.vb_per_obs_patch;
        write_pod(f, cache_num_obs_patches);
        write_pod(f, cache_p_per_obs_patch);
        write_pod(f, cache_vb_per_obs_patch);


        // --- Nodes ---
        for (auto& node : user_graph.nodes) {
            uint8_t is_boundary_u8 = node.is_boundary ? 1 : 0;
            write_pod(f, is_boundary_u8);
            uint32_t num_neighbors = (uint32_t)node.neighbors.size();
            write_pod(f, num_neighbors);
            write_pod(f, node.x);
            write_pod(f, node.y);
            write_pod(f, node.round);
            int32_t observable_id = node.observable_id;
            write_pod(f, observable_id);
        }

        // --- Edges (flatten list to stable indices in list order) ---
        for (auto& edge : user_graph.edges) {
            uint64_t node1 = edge.node1;
            uint64_t node2 = edge.node2;
            write_pod(f, node1);
            write_pod(f, node2);
            write_pod(f, edge.weight);
            write_pod(f, edge.error_probability);

            uint32_t num_obs_indices = (uint32_t)edge.observable_indices.size();
            write_pod(f, num_obs_indices);
            for (size_t idx : edge.observable_indices) {
                uint64_t idx64 = idx;
                write_pod(f, idx64);
            }

            uint32_t num_implied = (uint32_t)edge.implied_weights_for_other_edges.size();
            write_pod(f, num_implied);
            for (auto& iw : edge.implied_weights_for_other_edges) {
                uint64_t iw_node1 = iw.node1;
                uint64_t iw_node2 = iw.node2;
                write_pod(f, iw_node1);
                write_pod(f, iw_node2);
                write_pod(f, iw.implied_weight);
            }
        }

        // --- Boundary nodes (std::set iterates in sorted order) ---
        std::vector<uint64_t> boundary_vec(user_graph.boundary_nodes.begin(), user_graph.boundary_nodes.end());
        write_u64_vector(f, boundary_vec);

        // --- node_part_id ---
        write_pod(f, (uint64_t)user_graph.node_part_id.size());
        for (int part_id : user_graph.node_part_id) {
            int32_t part_id_32 = part_id;
            write_pod(f, part_id_32);
        }

        // --- virtual_boundaries (vector<vector<int>>) ---
        uint64_t num_vb_lists = user_graph.virtual_boundaries.size();
        write_pod(f, num_vb_lists);
        for (auto& vb_list : user_graph.virtual_boundaries) {
            uint64_t vb_list_len = vb_list.size();
            write_pod(f, vb_list_len);
            for (int node_idx : vb_list) {
                int32_t node_idx_32 = node_idx;
                write_pod(f, node_idx_32);
            }
        }

        fclose(f);
    } catch (...) {
        fclose(f);
        throw;
    }
}

pm::UserGraph pm::read_user_graph_cache(
    const std::string& path,
    bool expected_enable_correlations,
    int64_t expected_rounds_per_partition,
    int expected_division_strategy
) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f)
        throw std::runtime_error("graph_cache: failed to open " + path + " for reading");

    try {
        uint64_t magic;
        read_pod(f, magic);
        if (magic != GRAPH_CACHE_MAGIC) {
            throw pm::GraphCacheMismatchError("graph_cache: bad magic number in " + path);
        }

        uint32_t format_version;
        read_pod(f, format_version);
        if (format_version != GRAPH_CACHE_VERSION) {
            throw pm::GraphCacheMismatchError(
                "graph_cache: format version mismatch (file has " + std::to_string(format_version) +
                ", expected " + std::to_string(GRAPH_CACHE_VERSION) + ")");
        }

        uint8_t build_mode;
        read_pod(f, build_mode);
        uint8_t expected_build_mode = 0;
        // expected_build_mode = 1;
        expected_build_mode = 2;

        if (build_mode != expected_build_mode) {
            throw pm::GraphCacheMismatchError(
                "graph_cache: build mode mismatch (file has " + std::to_string(build_mode) +
                ", this binary is " + std::to_string(expected_build_mode) + ")");
        }

        uint8_t enable_corr_u8;
        read_pod(f, enable_corr_u8);
        if ((bool)enable_corr_u8 != expected_enable_correlations) {
            throw pm::GraphCacheMismatchError("graph_cache: enable_correlations mismatch");
        }

        uint8_t loaded_wo_corr_u8;
        read_pod(f, loaded_wo_corr_u8);

        uint64_t num_nodes, num_edges, num_observables, num_boundary_nodes;
        read_pod(f, num_nodes);
        read_pod(f, num_edges);
        read_pod(f, num_observables);
        read_pod(f, num_boundary_nodes);

        int64_t cache_rounds_per_partition;
        read_pod(f, cache_rounds_per_partition);
        if (cache_rounds_per_partition != expected_rounds_per_partition) {
            throw pm::GraphCacheMismatchError(
                "graph_cache: rounds_per_partition mismatch (file has " +
                std::to_string(cache_rounds_per_partition) + ", expected " +
                std::to_string(expected_rounds_per_partition) + ")");
        }
        uint64_t cache_num_partitions, cache_num_rounds, cache_num_virtual_boundaries;
        read_pod(f, cache_num_partitions);
        read_pod(f, cache_num_rounds);
        read_pod(f, cache_num_virtual_boundaries);

        int32_t cache_division_strategy;
        read_pod(f, cache_division_strategy);
        if (cache_division_strategy != expected_division_strategy) {
            throw pm::GraphCacheMismatchError(
                "graph_cache: task_division_strategy mismatch (file has " +
                std::to_string(cache_division_strategy) + ", expected " +
                std::to_string(expected_division_strategy) + ")");
        }
        uint64_t cache_num_obs_patches, cache_p_per_obs_patch, cache_vb_per_obs_patch;
        read_pod(f, cache_num_obs_patches);
        read_pod(f, cache_p_per_obs_patch);
        read_pod(f, cache_vb_per_obs_patch);

        pm::UserGraph user_graph(num_nodes, num_observables);
        user_graph.loaded_from_dem_without_correlations = (bool)loaded_wo_corr_u8;

        user_graph.num_partitions = cache_num_partitions;
        user_graph.num_rounds = cache_num_rounds;
        user_graph.num_virtual_boundaries = cache_num_virtual_boundaries;
        user_graph.num_obs_patches = cache_num_obs_patches;
        user_graph.p_per_obs_patch = cache_p_per_obs_patch;
        user_graph.vb_per_obs_patch = cache_vb_per_obs_patch;

        // --- Nodes ---
        // node.neighbors is intentionally left empty here; it is rebuilt below from
        // the edge list once all edges are read (UserNeighbor::edge_it can only be
        // set once the edges live in their final std::list, since it holds an iterator).
        std::vector<uint32_t> expected_num_neighbors(num_nodes, 0);
        for (uint64_t i = 0; i < num_nodes; ++i) {
            auto& node = user_graph.nodes[i];
            uint8_t is_boundary_u8;
            read_pod(f, is_boundary_u8);
            node.is_boundary = (bool)is_boundary_u8;
            uint32_t num_neighbors;
            read_pod(f, num_neighbors);
            expected_num_neighbors[i] = num_neighbors;
            read_pod(f, node.x);
            read_pod(f, node.y);
            read_pod(f, node.round);
            int32_t observable_id;
            read_pod(f, observable_id);
            node.observable_id = observable_id;
        }

        // --- Edges ---
        for (uint64_t i = 0; i < num_edges; ++i) {
            pm::UserEdge edge;
            uint64_t node1, node2;
            read_pod(f, node1);
            read_pod(f, node2);
            edge.node1 = node1;
            edge.node2 = node2;
            read_pod(f, edge.weight);
            read_pod(f, edge.error_probability);

            uint32_t num_obs_indices;
            read_pod(f, num_obs_indices);
            edge.observable_indices.resize(num_obs_indices);
            for (uint32_t j = 0; j < num_obs_indices; ++j) {
                uint64_t idx64;
                read_pod(f, idx64);
                edge.observable_indices[j] = (size_t)idx64;
            }

            uint32_t num_implied;
            read_pod(f, num_implied);
            edge.implied_weights_for_other_edges.resize(num_implied);
            for (uint32_t j = 0; j < num_implied; ++j) {
                uint64_t iw_node1, iw_node2;
                read_pod(f, iw_node1);
                read_pod(f, iw_node2);
                double implied_weight;
                read_pod(f, implied_weight);
                edge.implied_weights_for_other_edges[j] =
                    pm::ImpliedWeightUnconverted{(size_t)iw_node1, (size_t)iw_node2, implied_weight};
            }

            user_graph.edges.push_back(std::move(edge));
        }

        // --- Rebuild neighbors from the now-final edge list ---
        {
            std::vector<std::list<pm::UserEdge>::iterator> edge_iters;
            edge_iters.reserve(num_edges);
            for (auto it = user_graph.edges.begin(); it != user_graph.edges.end(); ++it) {
                edge_iters.push_back(it);
            }
            for (uint64_t i = 0; i < num_edges; ++i) {
                auto& it = edge_iters[i];
                user_graph.nodes[it->node1].neighbors.push_back(pm::UserNeighbor{it, 0});
                // node2 == SIZE_MAX is the sentinel for a boundary edge (no second endpoint) —
                // matches the convention used throughout user_graph.cc (e.g. add_noise()).
                if (it->node2 != SIZE_MAX) {
                    user_graph.nodes[it->node2].neighbors.push_back(pm::UserNeighbor{it, 1});
                }
            }
        }

        // --- Boundary nodes ---
        std::vector<uint64_t> boundary_vec = read_u64_vector(f);
        for (uint64_t b : boundary_vec) {
            user_graph.boundary_nodes.insert((size_t)b);
        }

        // --- node_part_id ---
        uint64_t node_part_id_len;
        read_pod(f, node_part_id_len);
        user_graph.node_part_id.resize(node_part_id_len);
        for (uint64_t i = 0; i < node_part_id_len; ++i) {
            int32_t part_id_32;
            read_pod(f, part_id_32);
            user_graph.node_part_id[i] = part_id_32;
        }

        // --- virtual_boundaries ---
        uint64_t num_vb_lists;
        read_pod(f, num_vb_lists);
        user_graph.virtual_boundaries.resize(num_vb_lists);
        for (uint64_t i = 0; i < num_vb_lists; ++i) {
            uint64_t vb_list_len;
            read_pod(f, vb_list_len);
            user_graph.virtual_boundaries[i].resize(vb_list_len);
            for (uint64_t j = 0; j < vb_list_len; ++j) {
                int32_t node_idx_32;
                read_pod(f, node_idx_32);
                user_graph.virtual_boundaries[i][j] = node_idx_32;
            }
        }

        fclose(f);
        return user_graph;
    } catch (...) {
        fclose(f);
        throw;
    }
}
