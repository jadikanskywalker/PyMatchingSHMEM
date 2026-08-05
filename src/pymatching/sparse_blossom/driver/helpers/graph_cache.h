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

}  // namespace pm

#endif  // _PM_GRAPH_CACHE_H
