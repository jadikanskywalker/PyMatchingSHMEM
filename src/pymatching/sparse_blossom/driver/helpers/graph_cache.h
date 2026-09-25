#ifndef _PM_GRAPH_CACHE_H
#define _PM_GRAPH_CACHE_H

#include <stdexcept>
#include <string>

#include "pymatching/sparse_blossom/driver/user_graph.h"

namespace pm {

// Thrown when a graph cache file exists but its embedded metadata (format version,
// enable_correlations, rounds_per_partition, task_division_strategy) does not match
// the parameters of the run trying to read it, or the file fails a basic sanity check.
// Callers should catch this and fall back to rebuilding the graph from the DEM.
struct GraphCacheMismatchError : public std::runtime_error {
    explicit GraphCacheMismatchError(const std::string& what) : std::runtime_error(what) {}
};

// The scalar metadata stored at the start of a graph cache file.
struct GraphCacheHeader {
    uint64_t num_nodes{0}, num_edges{0}, num_observables{0}, num_boundary_nodes{0};
    uint64_t num_partitions{0}, num_rounds{0}, num_virtual_boundaries{0};
    uint64_t num_obs_patches{0}, p_per_obs_patch{0}, vb_per_obs_patch{0};
    bool loaded_from_dem_without_correlations{false};
    // UserGraph::get_num_detectors() of the cached graph.
    uint64_t num_detectors() const { return num_nodes - num_boundary_nodes; }
};

// Takes a non-const reference because UserGraph's read-only accessors
// (get_num_observables(), get_num_nodes(), etc.) are not marked const.
void write_user_graph_cache(
    pm::UserGraph& user_graph,
    const std::string& path,
    bool enable_correlations,
    int64_t rounds_per_partition,
    int division_strategy
);

// Throws GraphCacheMismatchError if the cache's stored parameters don't match the
// expected_* arguments, or std::runtime_error on I/O or file-format errors.
pm::UserGraph read_user_graph_cache(
    const std::string& path,
    bool expected_enable_correlations,
    int64_t expected_rounds_per_partition,
    int expected_division_strategy
);

// Reads only the header (same validation as read_user_graph_cache), without loading the graph.
GraphCacheHeader read_user_graph_cache_header(
    const std::string& path,
    bool expected_enable_correlations,
    int64_t expected_rounds_per_partition,
    int expected_division_strategy
);

}  // namespace pm

#endif  // _PM_GRAPH_CACHE_H
