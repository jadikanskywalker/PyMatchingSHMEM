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

#ifndef PYMATCHING2_SHOT_BUFFER_H
#define PYMATCHING2_SHOT_BUFFER_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <vector>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_task.h"
#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"

namespace pm {

#ifdef USE_SHMEM
enum ShotStatus : uint64_t { READY, PUT_SUMMARY, PUT_RESULT, WROTE_RESULT };
#endif

// ShotContainer isolates everything needed to solve
// a single shot in parallel.
struct ShotContainer {
#ifdef USE_SHMEM
    uint64_t* current_buffer_round_shm;
#endif
    std::atomic<int> current_buffer_round{-1};  // Used for idle thread spin-wait until new shot read

    stim::SparseShot sparse_shot;

    std::vector<std::vector<uint64_t>> partition_hits;
    std::vector<std::vector<uint64_t>> virtual_boundary_hits;

    int num_observables;
    pm::MatchingResult obs_mask;
    pm::ExtendedMatchingResult res;

    std::vector<bool> i_solved_p;
    std::vector<bool> i_solved_vb;
    std::vector<Task> tasks;  // Tasks handle dynamic fusion tree synchonization
#ifdef USE_SHMEM
    std::vector<CrossRankTask> cross_rank_tasks;
#endif

    ShotContainer(
#ifdef USE_SHMEM
        uint64_t* current_buffer_round_ptr,
#endif
        int num_partitions, int num_virtual_boundaries, int num_observables);

    ShotContainer(const ShotContainer&) = delete;
    ShotContainer& operator=(const ShotContainer&) = delete;

    ShotContainer(ShotContainer&& other) noexcept;
    ShotContainer& operator=(ShotContainer&& other) noexcept;

    void clear();
};

// Rotating buffer of shot containers with ordered shot
// reading and result writing.
struct ShotBuffer {
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader;
    std::unique_ptr<stim::MeasureRecordWriter> writer;

    std::vector<ShotContainer> buffer;

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

    // Reads next shot into shot_container_id
    void write_result_and_get_next_shot(int shot_container_id, std::vector<int> &node_part_id
#ifdef USE_SHMEM
        , bool write_results=true
#endif
    );

    void read_shot(int shot_container_id, std::vector<int> &node_part_id);
};

}

#endif // PYMATCHING2_SHOT_BUFFER_H