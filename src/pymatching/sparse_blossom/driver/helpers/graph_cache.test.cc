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

#include "pymatching/sparse_blossom/driver/helpers/graph_cache.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <tuple>
#include <vector>

#include "pymatching/sparse_blossom/driver/user_graph.h"

#ifdef USE_THREADS
#include "pymatching/sparse_blossom/config_parallel.h"
#endif

namespace {

struct RaiiTempNamedFile {
    std::string path;
    RaiiTempNamedFile() {
        char tmp_filename[] = "/tmp/pm_graph_cache_test_XXXXXX";
        int descriptor = mkstemp(tmp_filename);
        if (descriptor == -1) {
            throw std::runtime_error("Failed to create temporary file.");
        }
        close(descriptor);
        path = tmp_filename;
    }
    ~RaiiTempNamedFile() {
        if (!path.empty()) {
            remove(path.data());
        }
    }
};

using FlatEdge = std::tuple<size_t, size_t, double, double, std::vector<size_t>>;

std::vector<FlatEdge> flatten_edges(pm::UserGraph& g) {
    std::vector<FlatEdge> out;
    for (auto& e : g.edges) {
        out.push_back({e.node1, e.node2, e.weight, e.error_probability, e.observable_indices});
    }
    return out;
}

pm::UserGraph make_test_graph() {
    stim::DetectorErrorModel dem(
        "detector(0, 0, 0) D0\n"
        "detector(0, 0, 1) D1\n"
        "detector(1, 0, 0) D2\n"
        "detector(1, 0, 1) D3\n"
        "error(0.1) D0 D1\n"
        "error(0.2) D1 D2 L0\n"
        "error(0.05) D2 D3 L1\n"
        "error(0.15) D3\n");
    return pm::detector_error_model_to_user_graph(dem, /*enable_correlations=*/false, pm::NUM_DISTINCT_WEIGHTS);
}

}  // namespace

TEST(GraphCache, RoundTripPreservesNodesEdgesAndBoundary) {
#ifdef USE_THREADS
    config_parallel::M = 10;
#ifdef USE_SHMEM
    config_parallel::division_strategy = config_parallel::ROUND;
#endif
#endif

    pm::UserGraph original = make_test_graph();
    RaiiTempNamedFile tmp;

    pm::write_user_graph_cache(
        original, tmp.path, /*enable_correlations=*/false
#ifdef USE_THREADS
        , config_parallel::M
#ifdef USE_SHMEM
        , (int)config_parallel::division_strategy
#endif
#endif
    );

    pm::UserGraph loaded = pm::read_user_graph_cache(
        tmp.path, /*expected_enable_correlations=*/false
#ifdef USE_THREADS
        , config_parallel::M
#ifdef USE_SHMEM
        , (int)config_parallel::division_strategy
#endif
#endif
    );

    ASSERT_EQ(original.get_num_nodes(), loaded.get_num_nodes());
    ASSERT_EQ(original.get_num_observables(), loaded.get_num_observables());
    ASSERT_EQ(original.get_num_edges(), loaded.get_num_edges());
    ASSERT_EQ(original.boundary_nodes, loaded.boundary_nodes);
    ASSERT_EQ(original.loaded_from_dem_without_correlations, loaded.loaded_from_dem_without_correlations);

    ASSERT_EQ(flatten_edges(original), flatten_edges(loaded));

    for (size_t i = 0; i < original.get_num_nodes(); i++) {
        ASSERT_EQ(original.nodes[i].is_boundary, loaded.nodes[i].is_boundary);
        ASSERT_EQ(original.nodes[i].neighbors.size(), loaded.nodes[i].neighbors.size());
    }

    // Functional check: neighbor lookups still resolve to the correct edge data
    // post-load (exercises the edge_it iterator-fixup reconstruction).
    for (size_t i = 0; i < loaded.get_num_nodes(); i++) {
        for (auto& neighbor : loaded.nodes[i].neighbors) {
            size_t other = (neighbor.pos == 0) ? neighbor.edge_it->node2 : neighbor.edge_it->node1;
            EXPECT_NE(loaded.index_of_neighbor(other), SIZE_MAX);
        }
    }

#ifdef USE_THREADS
    ASSERT_EQ(original.node_part_id, loaded.node_part_id);
    ASSERT_EQ(original.num_partitions, loaded.num_partitions);
    ASSERT_EQ(original.num_rounds, loaded.num_rounds);
    ASSERT_EQ(original.num_virtual_boundaries, loaded.num_virtual_boundaries);
    ASSERT_EQ(original.virtual_boundaries, loaded.virtual_boundaries);
#endif
}

TEST(GraphCache, MismatchedEnableCorrelationsThrows) {
#ifdef USE_THREADS
    config_parallel::M = 10;
#ifdef USE_SHMEM
    config_parallel::division_strategy = config_parallel::ROUND;
#endif
#endif

    pm::UserGraph original = make_test_graph();
    RaiiTempNamedFile tmp;

    pm::write_user_graph_cache(
        original, tmp.path, /*enable_correlations=*/false
#ifdef USE_THREADS
        , config_parallel::M
#ifdef USE_SHMEM
        , (int)config_parallel::division_strategy
#endif
#endif
    );

    ASSERT_THROW(
        pm::read_user_graph_cache(
            tmp.path, /*expected_enable_correlations=*/true
#ifdef USE_THREADS
            , config_parallel::M
#ifdef USE_SHMEM
            , (int)config_parallel::division_strategy
#endif
#endif
        ),
        pm::GraphCacheMismatchError);
}

#ifdef USE_THREADS
TEST(GraphCache, MismatchedRoundsPerPartitionThrows) {
    config_parallel::M = 10;
#ifdef USE_SHMEM
    config_parallel::division_strategy = config_parallel::ROUND;
#endif

    pm::UserGraph original = make_test_graph();
    RaiiTempNamedFile tmp;

    pm::write_user_graph_cache(
        original, tmp.path, /*enable_correlations=*/false, config_parallel::M
#ifdef USE_SHMEM
        , (int)config_parallel::division_strategy
#endif
    );

    ASSERT_THROW(
        pm::read_user_graph_cache(
            tmp.path, /*expected_enable_correlations=*/false, /*expected_rounds_per_partition=*/11
#ifdef USE_SHMEM
            , (int)config_parallel::division_strategy
#endif
        ),
        pm::GraphCacheMismatchError);
}
#endif
