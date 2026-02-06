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
#include <queue>
#include <vector>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"

namespace pm {

#ifdef USE_SHMEM
#define NUM_SYNCHRONIZATION_ATOMICS 2+NUM_BUFFERS_PER_UNIT // number of atomic numbers needed for sychnronization during decoding
#endif

struct ShotContainer {
    std::atomic<int> current_buffer_round{0};  // Used for idle thread spin-wait until new shot read
    uint64_t current_shot;

    stim::SparseShot sparse_shot;

    std::vector<std::vector<uint64_t>> partition_hits;
    std::vector<std::vector<uint64_t>> virtual_boundary_hits;

    pm::ExtendedMatchingResult res;

    std::vector<Task> tasks;  // Tasks handle dynamic fusion tree synchonization

    ShotContainer(
        int num_partitions, int num_virtual_boundaries, int num_observables)
        : partition_hits(num_partitions), virtual_boundary_hits(num_virtual_boundaries), res(num_observables) {
    }

    ShotContainer(const ShotContainer&) = delete;
    ShotContainer& operator=(const ShotContainer&) = delete;

    ShotContainer(ShotContainer&& other) noexcept
        : sparse_shot(std::move(other.sparse_shot)),
          partition_hits(std::move(other.partition_hits)),
          virtual_boundary_hits(std::move(other.virtual_boundary_hits)),
          res(std::move(other.res)),
          tasks(std::move(other.tasks)) {}

    ShotContainer& operator=(ShotContainer&& other) noexcept {
        if (this != &other) {
            sparse_shot = std::move(other.sparse_shot);
            partition_hits = std::move(other.partition_hits);
            virtual_boundary_hits = std::move(other.virtual_boundary_hits);
            res = std::move(other.res);
            tasks = std::move(other.tasks);
        }
        return *this;
    }

    void clear();
};

// rotating shot buffer
struct ShotBuffer {
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader;
    std::unique_ptr<stim::MeasureRecordWriter> writer;

    std::vector<ShotContainer> buffer;

#ifdef USE_SHMEM
    // Symmetric memory pointers
    uint64_t* shot_container_status;
#endif
    int next_shot_container_id{ 0 };
    int last_shot_container_id{ -1 };

    // Lock for result writing & shot reading
    std::mutex m;
    std::condition_variable cv;

    ShotBuffer(
#ifdef USE_SHMEM
        uint64_t* atomics_ptr,
#endif
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        int num_partitions,
        int num_virtual_boundaries,
        int num_observables);
};

struct SharedMatchingGraph {
    std::shared_ptr<pm::MatchingGraph> graph_ptr;
    std::vector<int> node_part_id;
    int num_partitions;
    int num_virtual_boundaries;

#ifdef USE_SHMEM
    std::vector<std::pair<size_t, size_t>> partition_bounds;
#endif

    SharedMatchingGraph();

    SharedMatchingGraph(
        std::shared_ptr<pm::MatchingGraph> graph_ptr_,
        std::vector<int> node_part_id_,
        int num_partitions_,
        int num_virtual_boundaries_);

#ifdef USE_SHMEM
    void construct_partition_bounds();
#endif
};

// A decoding unit is responsible for decoding a connected decoding graph in parallel
//   Connected decoding graphs are partitioned for parallel solving
//   A decoding task involves solving a partition or fusing two solved partition along a virtual boundary
//   Decoding units implement multi-threading and OpenSHMEM behavoior
struct DecodingUnit {
    // Graph
    SharedMatchingGraph graph;
    // std::shared_ptr<pm::MatchingGraph> graph_ptr;

    // std::vector<int> node_part_id;
    // int num_partitions;
    // int num_virtual_boundaries;

    // Shots
    std::shared_ptr<ShotBuffer> shot_buffer;

    // Solvers
    bool ensure_search_flooder_included;
    bool enable_correlations;
    bool draw_frames;
    int num_threads;
    int num_partition_units;
    std::vector<std::shared_ptr<pm::Mwpm>> solvers;

#ifdef USE_SHMEM
    int pid;
    int other_pid; // for simple two PE impl.
    // Symmetric memory
    DetectorNodeEphemeralFields* node_ephemeral_fields_ptr = nullptr;
    uint64_t* atomics_ptr = nullptr;
    GraphFillRegion* regions_ptr = nullptr;
    // Cross-PE Fusion Contexts
    shmem_ctx_t summary_ctx;
    shmem_ctx_t fallback_ctx;
    // Info used for symmetric memory accesses
    size_t regions_nelems_per_solver;
    GraphFillRegion* regions_base_other_pe = nullptr; // base for region SHMEMArena on other PE
#endif

    // Initialization Functions
    DecodingUnit(
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        const stim::DetectorErrorModel& detector_error_model,
        pm::weight_int num_distinct_weights,
        bool ensure_search_flooder_included,
        bool enable_correlations,
        bool draw_frames);

    ~DecodingUnit();

    void build_tasks_for_round_partitioning();

    void build_solvers();

    // Decoding Functions
#ifdef USE_SHMEM
    void handle_cross_process_fusion_and_get_next_shot(int shot_container_id, int num_threads);
#else
    void write_result_and_get_next_shot(int shot_container_id);
#endif

    void decode_shots();

    void solve_task(pm::Mwpm& solver, std::vector<uint64_t>& hits, Task* task, int tid, int draw_frames, int shot_id);
};

}  // namespace pm

#endif  // PYMATCHING2_DECODING_UNIT_H