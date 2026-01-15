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

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"

#include <vector>
#include <mutex>
#include <condition_variable>
#include <queue>

namespace pm { 
// class MatchingGraph; class Mwpm; struct ExtendedMatchingResult; }

struct ShotContainer {
    std::atomic<int> current_buffer_round{0};

    stim::SparseShot sparse_shot;

    std::vector<std::vector<uint64_t>> partition_hits;
    std::vector<std::vector<uint64_t>> virtual_boundary_hits;

    pm::ExtendedMatchingResult res;

    std::vector<Task> tasks;

    ShotContainer(int num_partitions, int num_virtual_boundaries, int num_observables)
     : partition_hits(num_partitions), virtual_boundary_hits(num_virtual_boundaries), res(num_observables) {}

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

struct ShotBuffer {
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader;
    std::unique_ptr<stim::MeasureRecordWriter> writer;

    std::vector<ShotContainer> buffer;
    int next_shot_buffer_id = 0;
    int last_shot_buffer_id = -1; // to mark end

    std::mutex m; // For result writing & shot reading
    std::condition_variable cv;

    ShotBuffer(
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        int num_partitions,
        int num_virtual_boundaries,
        int num_observables
    ) : reader(std::move(reader)), writer(std::move(writer)) {
        buffer.reserve(static_cast<size_t>(NUM_ACTIVE_SHOTS_PER_UNIT));
        for (int i = 0; i < NUM_ACTIVE_SHOTS_PER_UNIT; ++i) {
            buffer.emplace_back(num_partitions, num_virtual_boundaries, num_observables);
        }
    }
};

// A decoding unit is a connected decoding graph
//   Connected decoding graphs are partitioned to allow parallel solving
//   Decoding task involve solving a partition or fusing partitions along virtual boundaries
struct DecodingUnit {
    // Graph
    const std::shared_ptr<pm::MatchingGraph> graph_ptr;

    const std::vector<int> node_part_id;
    const int num_partitions;
    const int num_virtual_boundaries;

    // Shots
    std::shared_ptr<ShotBuffer> shot_buffer;

    // Solvers
    bool enable_correlations;
    bool draw_frames;
    std::vector<std::shared_ptr<pm::Mwpm>> solvers;
    // int num_solver_sets; // num_solvers / num_partitions

    // bool all_shots_read = false;
    bool done = false;

    // Initialization Functions
    DecodingUnit(
        const std::shared_ptr<pm::MatchingGraph> graph_ptr,
        const std::vector<int> node_part_id,
        int num_partitions,
        int num_virtual_boundaries
    ) : graph_ptr(graph_ptr), node_part_id(node_part_id), num_partitions(num_partitions), num_virtual_boundaries(num_virtual_boundaries) {}

    void setup(
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        bool enable_correlations,
        bool draw_frames,
        int max_threads,
        stim::DetectorErrorModel dem /* for drawing coors */
    );

    void build_tasks_for_round_partitioning();

    void build_solvers(
        bool ensure_search_flooder_included,
        bool enable_correlations,
        int num_threads);

    // Decoding Functions
    void write_result_and_get_next_shot(int shot_buffer_id);

    void decode_shots();

    void solve_task(pm::Mwpm& solver, std::vector<uint64_t>& hits, Task* task, int tid, int draw_frames, int shot_id);
};

}

#endif // PYMATCHING2_DECODING_UNIT_H