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

#ifndef PYMATCHING2_DECODING_UNIT_H
#define PYMATCHING2_DECODING_UNIT_H

#include <condition_variable>
#include <mutex>
#include <vector>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/shot_buffer.h"

#ifdef USE_SHMEM
#include <shmem.h>
#endif

namespace pm {

#ifdef USE_SHMEM
struct FusionSummary {
    size_t regions_matched_to_vb_size;
    GraphFillRegion* regions_begin; // index where solving partition's arena begins
};
#endif

// A decoding unit is responsible for decoding a connected decoding graph in parallel
//   Connected decoding graphs are partitioned for parallel solving
//   A decoding task involves solving a partition or fusing two solved partition along a virtual boundary
//   Decoding units implement multi-threading and OpenSHMEM behavoior
struct DecodingUnit {
    // Graph
    SharedMatchingGraph graph;

    // Shots
    std::shared_ptr<ShotBuffer> shot_buffer;

    // Params
    bool ensure_search_flooder_included;
    bool enable_correlations;
#if ENABLE_DRAW_FLAGS
    bool draw_frames;
#endif

    // Solvers
    int num_threads;
    int num_partition_units;
    std::vector<std::shared_ptr<Mwpm>> solvers;

#ifdef USE_SHMEM
    int pid;
    int other_pid; // for simple two PE impl.

    // Cross-PE Fusion Contexts
    shmem_ctx_t summary_ctx;
    shmem_ctx_t fallback_ctx;

    // Symmetric memory
    DetectorNodeEphemeralFields* node_ephemeral_fields_ptr{ nullptr };
    GraphFillRegion* regions_ptr{ nullptr };
    // One atomic sychnronizer per shot buffer
    // Each fusing PE fetch_inc's, the second PE wins fusion ()
    uint64_t* atomics_ptr{ nullptr };
    FusionSummary* fusion_summary_ptr { nullptr };
    // GraphFillRegionSubstates* regions_matched_to_vb_ptr { nullptr };
    std::pair<size_t, size_t>* closest_nb_ptr{ nullptr };
    uint32_t* obs_crossed_ptr{ nullptr };

    // Info used for symmetric memory accesses
    size_t regions_nelems_per_solver;
    size_t closest_nb_nelems_per_buffer;
    size_t obs_crossed_nelems_per_buffer;
    GraphFillRegion* regions_base_other_pe{ nullptr }; // base for region SHMEMArena on other PE
#endif

    // Initialization Functions
    DecodingUnit(
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        const stim::DetectorErrorModel& detector_error_model,
        weight_int num_distinct_weights,
        bool ensure_search_flooder_included,
        bool enable_correlations
    #if ENABLE_DRAW_FLAGS
        , bool draw_frames
    #endif
        );

    ~DecodingUnit();

    void build_tasks_for_round_partitioning();

    void build_solvers();

    // Decoding Functions
#ifdef USE_SHMEM
    void solve_cross_process_fusion_and_get_next_shot(size_t shot_container_id, Task* t, size_t tid, size_t num_threads, size_t solver_id, size_t shot_id);
#endif
    void write_result_and_get_next_shot(int shot_container_id);

    void decode_shots();

    void solve_task(Mwpm& solver, std::vector<uint64_t>& hits, Task* task, int tid
#if ENABLE_DRAW_FLAGS
        , int draw_frames
#endif
        , int shot_id);
};

}  // namespace pm

#endif  // PYMATCHING2_DECODING_UNIT_H