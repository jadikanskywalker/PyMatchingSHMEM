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

#ifndef PYMATCHING2_MWPM_DECODING_H
#define PYMATCHING2_MWPM_DECODING_H

#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "stim.h"
#include <set>

#ifdef USE_THREADS
#include "pymatching/sparse_blossom/driver/work_stealing_deque.h"
#include <memory>
#endif

namespace pm {

struct ExtendedMatchingResult {
    std::vector<uint8_t> obs_crossed;
    total_weight_int weight;
    ExtendedMatchingResult();
    explicit ExtendedMatchingResult(size_t num_observables);

    bool operator==(const ExtendedMatchingResult& rhs) const;

    bool operator!=(const ExtendedMatchingResult& rhs) const;

    ExtendedMatchingResult(std::vector<uint8_t> obs_crossed, total_weight_int weight);

    void reset();

    ExtendedMatchingResult& operator+=(const ExtendedMatchingResult& rhs);
    ExtendedMatchingResult operator+(const ExtendedMatchingResult& rhs) const;
};

inline void ExtendedMatchingResult::reset() {
    std::fill(obs_crossed.begin(), obs_crossed.end(), 0);
    weight = 0;
}

void fill_bit_vector_from_obs_mask(pm::obs_int obs_mask, uint8_t* obs_begin_ptr, size_t num_observables);
obs_int bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector);

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included = false,
    bool enable_correlations = false);

MatchingResult decode_detection_events_for_up_to_64_observables(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, bool edge_correlations);

/// Used to decode detection events for an existing Mwpm object `mwpm', and a vector of
/// detection event indices `detection_events'. The predicted observables are XOR-ed into an
/// existing uint8_t array with at least `mwpm.flooder.graph.num_observables' elements,
/// the pointer to the first element of which is passed as the `obs_begin_ptr' argument.
/// The weight of the MWPM solution is added to the `weight' argument.
void decode_detection_events(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight,
    bool edge_correlations
#ifdef ENABLE_FUSION
    , int shot = 0,
    bool draw_frames = false
#endif
);

/// Decode detection events using a Mwpm object and vector of detection event indices
/// Returns the compressed edges in the matching: the pairs of detection events that are
/// matched to each other via paths.
void decode_detection_events_to_match_edges(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Decode detection events using a Mwpm object and vector of detection event indices.
/// Returns the edges in the matching: these are pairs of *detectors* forming *edges* in the
/// matching solution (rather than pairs of detection *events* matched via *paths* as returned
/// instead by `decode_detection_events_to_match_edges`).
void decode_detection_events_to_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges);

void decode_detection_events_to_edges_with_edge_correlations(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges);
 
#ifdef ENABLE_FUSION
void setup_output_dirs(bool draw_frames, bool parallel=false);
void output_detector_nodes(pm::Mwpm& mwpm, bool parallel=false);
void output_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, int shot, bool parallel=false);
void output_solution_state(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, int shot, std::set<long>parts, bool parallel=false);
void draw_frame(pm::Mwpm& mwpm, pm::MwpmEvent ev, int shot, int frame_number, bool parallel=false, int thread=-1);
#endif

#ifdef USE_THREADS
// ===============
// Initially, a task is created for each partition and tasks are
//   assigned to thread partition%num_threads.
// Fusion collapses tasks leftward. If p3 & p4 are fused, p4 is
//   added to p3's task and p4's task is marked FINISHED.
extern std::vector<Task> tasks;
extern std::vector<int> partitions_task_id; // Convenience store of which task own each partition
// extern std::vector<std::queue<int>> partition_task_queues;
// extern std::vector<std::deque<int>> fusion_task_deques;
// Stable per-partition solver instances. Use shared_ptr to manage lifetime safely across threads.
extern std::vector<std::shared_ptr<Mwpm>> solvers;
// Build and assign one solver per thread from a DEM, storing stable instances and wiring `solvers`.
void build_partition_solvers(
    pm::Mwpm& mwpm,
    bool ensure_search_flooder_included,
    bool enable_correlations);
void reset_tasks(int num_threads, int num_partitions);
// void init_task_queues(int num_threads, int num_partitions);
#endif

#if defined(USE_THREADS) || defined(USE_SHMEM)
void decode_detection_events_in_parallel(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight,
    bool edge_correlations,
    int shot,
    bool draw_frames);
#endif
// ===============
} // namespace pm

#endif  // PYMATCHING2_MWPM_DECODING_H
