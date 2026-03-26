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

#include "pymatching/sparse_blossom/driver/parallel/shot_buffer.h"

#include <iostream>
#include <utility>

#include <omp.h>

#ifdef USE_SHMEM
#include <shmem.h>
#endif

namespace pm {

ShotContainer::ShotContainer(
#ifdef USE_SHMEM
    uint64_t* current_buffer_round_ptr,
#endif
    int num_partitions, int num_virtual_boundaries, int num_observables_in)
    : partition_hits(num_partitions),
      virtual_boundary_hits(num_virtual_boundaries),
      num_observables(num_observables_in),
#ifdef USE_SHMEM
      current_buffer_round_shm(current_buffer_round_ptr),
#endif
      i_solved_p(num_partitions, false),
      i_solved_vb(num_virtual_boundaries, false),
      res(num_observables_in) {}

ShotContainer::ShotContainer(ShotContainer&& other) noexcept
    : sparse_shot(std::move(other.sparse_shot)),
      partition_hits(std::move(other.partition_hits)),
      virtual_boundary_hits(std::move(other.virtual_boundary_hits)),
      num_observables(other.num_observables),
      res(std::move(other.res)),
      tasks(std::move(other.tasks)) {}

ShotContainer& ShotContainer::operator=(ShotContainer&& other) noexcept {
    if (this != &other) {
        sparse_shot = std::move(other.sparse_shot);
        partition_hits = std::move(other.partition_hits);
        virtual_boundary_hits = std::move(other.virtual_boundary_hits);
        num_observables = other.num_observables;
        res = std::move(other.res);
        tasks = std::move(other.tasks);
    }
    return *this;
}

void ShotContainer::clear() {
    sparse_shot.clear();
    for (auto& hits : partition_hits) {
        hits.clear();
    }
    for (auto& hits : virtual_boundary_hits) {
        hits.clear();
    }
    res.reset();
}

ShotBuffer::ShotBuffer(
#ifdef USE_SHMEM
    uint64_t* atomics_ptr,
#endif
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader_in,
    std::unique_ptr<stim::MeasureRecordWriter> writer_in,
    int num_partitions,
    int num_virtual_boundaries,
    int num_observables)
    : reader(std::move(reader_in)), writer(std::move(writer_in))
{
    buffer.reserve(static_cast<size_t>(NUM_BUFFERS_PER_UNIT));
    for (size_t i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
        buffer.emplace_back(
#ifdef USE_SHMEM
            atomics_ptr + i,
#endif
            num_partitions, num_virtual_boundaries, num_observables);
    }
}

void ShotBuffer::write_result_and_get_next_shot(
    int shot_container_id,
    std::vector<int>& node_part_id
#ifdef USE_SHMEM
    , bool write_results
#endif
) {
    std::unique_lock<std::mutex> lock(m);
    cv.wait(lock, [&] {
        return next_shot_container_id == shot_container_id;
    });
    auto& shot = buffer[shot_container_id];
#ifdef USE_SHMEM
    if (write_results) {
#endif
        if (DEBUG) {
            std::cout << "  T" << omp_get_thread_num() << " writing results for shot buffer " << shot_container_id << std::endl;
        }
        for (size_t k = 0; k < shot.num_observables; k++) {
            writer->write_bit(shot.res.obs_crossed[k]);
        }
        writer->write_end();
#ifdef USE_SHMEM
    }
#endif
    read_shot(shot_container_id, node_part_id);
#if NUM_BUFFERS_PER_UNIT > 1
    if (++next_shot_container_id >= NUM_BUFFERS_PER_UNIT) {
        next_shot_container_id = 0;
    }
#endif
    lock.unlock();
    cv.notify_all();
}

// Reads a shot
// Should only before multi-threading or when a thread hold the lock
void ShotBuffer::read_shot(int shot_container_id, std::vector<int>& node_part_id) {
    auto& shot = buffer[shot_container_id];
    shot.clear();
#if NUM_BUFFERS_PER_UNIT > 1
    if (last_shot_container_id < 0) {
#endif
        bool shot_read = pm::start_and_read_entire_record_buffered(*reader, shot.sparse_shot);
        if (shot_read) {
            for (auto det : shot.sparse_shot.hits) {
                int part_id = node_part_id[det];
                if (part_id >= 0) {
                    shot.partition_hits[part_id].push_back(det);
                } else {
                    shot.virtual_boundary_hits[-1 * part_id - 1].push_back(det);
                }
            }
            shot.current_buffer_round++;
        } else {
#if NUM_BUFFERS_PER_UNIT > 1
            last_shot_container_id = shot_container_id - 1;
            if (last_shot_container_id < 0) {
                last_shot_container_id = NUM_BUFFERS_PER_UNIT - 1;
            }
#else
            shot.current_buffer_round.store(-1);
#endif
        }
#if NUM_BUFFERS_PER_UNIT > 1
    } else if (last_shot_container_id == shot_container_id) {
        for (int i = 0; i < NUM_BUFFERS_PER_UNIT; ++i) {
            buffer[i].current_buffer_round.store(-1);
        }
    }
#endif
}

}  // namespace pm
