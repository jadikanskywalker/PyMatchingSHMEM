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
    DetectorNode* static_nodes_base; // Base address of graph.nodes on the sender PE
    size_t blossom_children_size;
    size_t regions_to_unmatch_size;
    GraphFillRegion* regions_to_unmatch[]; 
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
#ifdef ENABLE_DRAW_FLAGS
    bool draw_frames;
#endif

    // Solvers
    size_t num_threads;
    size_t num_partition_units;
    size_t num_solvers_per_buffer;
    std::vector<std::shared_ptr<Mwpm>> solvers;

    std::vector<size_t> my_partitions;

#ifdef USE_SHMEM
    
    int pid;
    int n_pes;

    // Symmetric memory
    DetectorNodeEphemeralFields* node_ephemeral_fields_ptr{ nullptr }; // ephemeral field buffers for DetectorNodes
    GraphFillRegion* regions_ptr{ nullptr }; // buffers for region SHMEMArenas
    BlossomChild* child_edges_ptr{ nullptr }; // blossom parent-child relationship buffer
    uint64_t* atomics_ptr{ nullptr }; // shot counter for each shot container
    uint64_t* task_status_ptr{ nullptr }; // status and signal for each fusion task
    FusionSummary* task_fusion_summary_ptr { nullptr }; // fusion summary for each fusion task
    // Info used for symmetric memory accesses
    size_t nodes_nelems_per_buffer;
    size_t regions_nelems_per_solver;
    size_t child_edges_nelems_per_solver;
    size_t task_fusion_summary_size_per_task;
    size_t regions_matched_to_vb_nelems;

    inline FusionSummary* get_fusion_summary_ptr(size_t shot_container_id, bool iamleft) {
        bool even = !(pid % 2);
        bool take_second_slot = (even && iamleft) || (!even && !iamleft);
        return reinterpret_cast<FusionSummary*>(
            reinterpret_cast<char*>(task_fusion_summary_ptr)
            + (shot_container_id * SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER + take_second_slot) * task_fusion_summary_size_per_task
        );
    }

    inline uint64_t* get_task_status_ptr(size_t shot_container_id, bool iamleft) {
        // There are two uint64s per slot: status and signal
        //   Simple left/right case: NUM_CROSS_RANK_FUSIONS_PER_BUFFER=2 slots per buffer
        //   Even PEs map slots (to left, to right), odd PEs (to right, to left)
        //   This ensures tasks correspond on each PE
        bool even = !(pid % 2);
        bool take_second_slot = (even && iamleft) || (!even && !iamleft);
        return task_status_ptr + 2*(SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER*shot_container_id + take_second_slot);
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

    inline size_t get_cross_rank_fusion_idx(size_t shot_container_id, size_t index) {
        return SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER * shot_container_id + index;
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
#ifdef ENABLE_DRAW_FLAGS
        , bool draw_frames
#endif
        );

    ~DecodingUnit();

    void build_tasks_for_round_partitioning();

    void build_solvers();

    inline int get_solver_id(int shot_container_id, int partition, int tid=0) 
    {
#ifdef USE_SHMEM
        return num_solvers_per_buffer * shot_container_id + partition;
#else
        int partition_unit = partition / num_threads;
        return num_solvers_per_buffer * shot_container_id + num_threads*partition_unit + tid;
#endif
    }


    // Decoding Functions
#ifdef USE_SHMEM
    void send_solution_to_remote_pe(size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &task, std::ofstream &t_out);
    bool get_solution_from_remote_pe(size_t shot_container_id, pm::MatchingResult& res, CrossRankTask &task, std::ofstream &t_out, std::vector<uint64_t>& hitsref); // returns whether solving is necessary
#endif

    void decode_shots();
};

}  // namespace pm

#endif  // PYMATCHING2_DECODING_UNIT_H