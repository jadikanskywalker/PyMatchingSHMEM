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

#include "pymatching/sparse_blossom/flooder/blossom_child.h"

namespace pm {

#ifdef USE_SHMEM
struct FusionSummary {
    GraphFillRegion* regions_ptr_base; // index where solving partition's arena begins
    // DetectorNodeEphemeralFields* nodes_ptr_base; // index where solving partition's nodes begin
    DetectorNode* static_nodes_base; // Base address of graph.nodes on the sender PE
    size_t blossom_children_size;
    size_t regions_matched_to_vb_size;
    GraphFillRegion* regions_matched_to_vb[]; 
    // uint64_t shmem_bitmap[]; // implicitly follows
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

    std::vector<size_t> my_partitions;

#ifdef USE_SHMEM
    int pid;
    // int other_pid; // for simple two PE impl.

    // Symmetric memory
    DetectorNodeEphemeralFields* node_ephemeral_fields_ptr{ nullptr };
    GraphFillRegion* regions_ptr{ nullptr };
    BlossomChild* child_edges_ptr{ nullptr };
    // One atomic sychronizer per shot container
    uint64_t* atomics_ptr{ nullptr };
    uint64_t* task_status_ptr{ nullptr };
    FusionSummary* task_fusion_summary_ptr { nullptr };
    uint32_t* obs_crossed_ptr{ nullptr };

    // Info used for symmetric memory accesses
    size_t nodes_nelems_per_buffer;
    size_t regions_nelems_per_solver;
    size_t child_edges_nelems_per_solver;
    size_t task_fusion_summary_size_per_task;
    size_t regions_matched_to_vb_nelems;
    size_t obs_crossed_nelems_per_buffer;
    GraphFillRegion* regions_base_other_pe{ nullptr }; // base for region SHMEMArena on other PE

    inline FusionSummary* get_fusion_summary_ptr(size_t shot_container_id, int partition_id) {
        return reinterpret_cast<FusionSummary*>(
            reinterpret_cast<char*>(task_fusion_summary_ptr)
            + (shot_container_id * graph.num_partitions + partition_id) * task_fusion_summary_size_per_task
        );
    }

    inline DetectorNodeEphemeralFields* get_node_fields_ptr(size_t shot_container_id, int partition_id) {
        auto& bounds = graph.partition_bounds[partition_id];
        return node_ephemeral_fields_ptr + shot_container_id * nodes_nelems_per_buffer + bounds.first;
    }

    inline GraphFillRegion* get_regions_ptr(size_t shot_container_id, int partition_id) {
        return regions_ptr + (shot_container_id * graph.num_partitions + partition_id) * regions_nelems_per_solver;
    }

    inline BlossomChild* get_child_edges_ptr(size_t shot_container_id, int partition_id) {
        return child_edges_ptr + (shot_container_id * graph.num_partitions + partition_id) * child_edges_nelems_per_solver;
    }

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

    void build_tasks_for_round_partitioning(size_t n_pes);

    void build_solvers();

    inline int get_solver_id(int shot_container_id, int partition, int tid=0) 
    {
#ifdef USE_SHMEM
        return graph.num_partitions * shot_container_id + partition;
#else
        int partition_unit = partition / num_threads;
        return num_threads * (num_partition_units * shot_container_id + partition_unit) + tid;
#endif
    }


    // Decoding Functions
#ifdef USE_SHMEM
    void send_solution_to_remote_pe(size_t shot_container_id, Task &task, std::ofstream &t_out);
    bool get_solution_from_remote_pe(size_t shot_container_id, Task &task, std::ofstream &t_out, std::vector<uint64_t>& hitsref); // returns whether solving is necessary
    void fuse_results_across_pes(size_t shot_container_id);
    void solve_cross_process_fusion_and_get_next_shot(size_t shot_container_id, Task* t, size_t tid, size_t num_threads, size_t solver_id, size_t shot_id);
#endif
    void write_result_and_get_next_shot(int shot_container_id);

    void decode_shots();

//     void solve_task(Mwpm& solver, std::vector<uint64_t>& hits, Task* task, int tid
// #if ENABLE_DRAW_FLAGS
//         , int draw_frames
// #endif
//         , int shot_id);
};

}  // namespace pm

#endif  // PYMATCHING2_DECODING_UNIT_H