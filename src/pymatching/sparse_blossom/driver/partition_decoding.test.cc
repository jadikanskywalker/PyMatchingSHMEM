// Copyright 2025 PyMatching Contributors
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

#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "pymatching/sparse_blossom/config_shmem.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "pymatching/sparse_blossom/flooder/graph.h"
#include "stim.h"

namespace {

std::string find_test_data_file(const char* name) {
    std::vector<std::string> directories_to_check = {
        "data/",
        "../data/",
        "../../data/",
    };
    for (const auto& d : directories_to_check) {
        std::string path = d + std::string(name);
        FILE* f = fopen(path.c_str(), "r");
        if (f != nullptr) {
            fclose(f);
            return path;
        }
    }
    throw std::invalid_argument(std::string("Failed to find test data file ") + name);
}

}

// Verifies that decoding using a single active partition that treats the cross-partition
// neighbor as a virtual boundary yields the same observable outcome as the full decode
// that matches across that edge.
TEST(PartitionDecoding, CrossPartitionEdgeMatchesBoundaryLike) {
#ifndef USE_SHMEM
    GTEST_SKIP() << "Partition gating test requires USE_SHMEM.";
#else
    // Make partitions as fine as possible in time to ensure many cross-partition edges.
    config_shmem::M = 1;

    // Load a DEM with detector coordinates (from a Stim circuit) so partitioning happens.
    std::string dem_path = find_test_data_file("surface_code_rotated_memory_x_13_0.01.dem");
    FILE* dem_file = fopen(dem_path.c_str(), "r");
    ASSERT_NE(dem_file, nullptr);
    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    // Build decoder (this will annotate nodes with rounds and partition them).
    pm::weight_int num_distinct_weights = pm::NUM_DISTINCT_WEIGHTS;
    auto mwpm = pm::detector_error_model_to_mwpm(
        dem,
        num_distinct_weights,
        /*ensure_search_flooder_included=*/false,
        /*enable_correlations=*/false);

    // Find a cross-partition edge u->v where v is in the next partition and is marked virtual.
    size_t u_idx = SIZE_MAX;
    size_t v_idx = SIZE_MAX;
    auto& nodes = mwpm.flooder.graph.nodes;
    ASSERT_GT(mwpm.flooder.graph.num_partitions, 0);
    for (size_t i = 0; i < nodes.size() && u_idx == SIZE_MAX; ++i) {
        const auto& u = nodes[i];
        if (u.neighbors.empty()) continue;
        // Skip boundary neighbor at index 0 if present.
        size_t start = (!u.neighbors.empty() && u.neighbors[0] == nullptr) ? 1 : 0;
        for (size_t k = start; k < u.neighbors.size(); ++k) {
            auto* v = u.neighbors[k];
            if (v == nullptr) continue;
            // Look for edge across a partition boundary to a virtual node in the next partition.
            if (v->partition == u.partition + 1 && v->is_virtual) {
                u_idx = i;
                v_idx = v - &nodes[0];
                break;
            }
        }
    }

    if (u_idx == SIZE_MAX || v_idx == SIZE_MAX) {
        GTEST_SKIP() << "No cross-partition edge with expected virtual semantics found in DEM.";
    }

    // Decode on the full graph with both detections; should match across (u,v).
    std::vector<uint64_t> hits_full = {u_idx, v_idx};
    auto res_full = pm::decode_detection_events_for_up_to_64_observables(
        mwpm, hits_full, /*edge_correlations=*/false);

    // Decode with only the lower partition active, and only seed u; v is in the higher partition.
    std::vector<uint64_t> hits_part = {u_idx};
    mwpm.flooder.active_partitions.clear();
    mwpm.flooder.active_partitions.insert(nodes[u_idx].partition);
    auto res_part = pm::decode_detection_events_for_up_to_64_observables(
        mwpm, hits_part, /*edge_correlations=*/false);

    // The observable mask and weight should match: virtual-boundary decode uses the same edge observables and weight.
    ASSERT_EQ(res_full.obs_mask, res_part.obs_mask);
    ASSERT_EQ(res_full.weight, res_part.weight);
#endif
}
