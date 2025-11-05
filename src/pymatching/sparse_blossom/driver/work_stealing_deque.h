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

#ifndef PYMATCHING2_WORK_STEALING_QUEUE_H
#define PYMATCHING2_WORK_STEALING_QUEUE_H

#include <atomic>
#include <memory>
#include <vector>


enum Status { UNSOLVED, BUSY = 0, STEALABLE, IN_TRANSFER, FINISHED };

struct Task {
    int solution_owner_tid{ -1 };
    std::vector<long> partitions{};
    std::atomic<Status> status{ Status::UNSOLVED };

    // Default construct with safe initial values.
    Task() = default;
    // Non-copyable due to atomic members.
    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    // Movable by manually transferring atomic state and moving shared_ptrs via load/store.
    Task(Task&& other) noexcept
        : solution_owner_tid(other.solution_owner_tid),
          partitions(std::move(partitions)),
          status(other.status.load(std::memory_order_relaxed)) {
        // Leave other in a benign state
        other.solution_owner_tid = -1;
        // other.partitions.store(nullptr, std::memory_order_release);
        other.status.store(Status::STEALABLE, std::memory_order_relaxed);
    }

    // void publish_partitions(const std::vector<long>& new_parts) {
    //     auto shared_parts = std::make_shared<const std::vector<long>>(new_parts);
    //     partitions.store(shared_parts, std::memory_order_release);
    // }
};


// If solved=False, partitions should only have one partition
//    -> indicates this partition still needs to be solved
// If solved=True, partitions should have a solved set of one
// or more partitions.
//    -> indicates set is ready to be fused with neighbor
// struct Task {
//     int partition_set_id;
//     // maybe adjacency info or owner ID
// };



// class WorkStealingDeque {
// private:
//     std::vector<Task> buffer;
//     std::atomic<size_t> top;
//     std::atomic<size_t> bottom;
//     const size_t capacity;

// public:
//     WorkStealingDeque(size_t cap = 1024)
//         : buffer(cap), top(0), bottom(0), capacity(cap) {}

//     // Only the thread who owns the queue can push
//     bool push(const Task& task) {
//         size_t b = bottom.load(std::memory_order_relaxed);
//         if (b - top.load(std::memory_order_acquire) >= capacity)
//             return false; // queue full
//         buffer[b % capacity] = task;
//         std::atomic_thread_fence(std::memory_order_release);
//         bottom.store(b + 1, std::memory_order_relaxed);
//         return true;
//     }

//     // Only the thread who owns the queue can pop
//     bool pop(Task& task) {
//         size_t b = bottom.load(std::memory_order_relaxed) - 1;
//         bottom.store(b, std::memory_order_relaxed);
//         std::atomic_thread_fence(std::memory_order_seq_cst);
//         size_t t = top.load(std::memory_order_relaxed);
//         if (t <= b) {
//             task = buffer[b % capacity];
//             return true;
//         } else {
//             bottom.store(t, std::memory_order_relaxed);
//             return false;
//         }
//     }

//     // Other threads can steal
//     bool steal(Task& task) {
//         size_t t = top.load(std::memory_order_acquire);
//         std::atomic_thread_fence(std::memory_order_seq_cst);
//         size_t b = bottom.load(std::memory_order_acquire);
//         if (t < b) {
//             task = buffer[t % capacity];
//             if (!top.compare_exchange_strong(t, t + 1,
//                     std::memory_order_seq_cst,
//                     std::memory_order_relaxed))
//                 return false;
//             return true;
//         }
//         return false;
//     }
// };


#endif // PYMATCHING2_WORK_STEALING_QUEUE_H