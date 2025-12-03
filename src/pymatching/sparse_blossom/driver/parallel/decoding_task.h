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

#ifndef PYMATCHING2_DECODING_TASK_H
#define PYMATCHING2_DECODING_TASK_H

#include <atomic>
#include <memory>
#include <vector>

// Forward declare to avoid cyclic include with mwpm_decoding.h
namespace pm { class GraphFillRegion; }

enum Status { BUSY, FREE };

struct Task {
private:
    std::atomic<Status> status{ Status::FREE };
    int shot_marker{ -1 };

public:
    const int task_id;
    const bool is_fusion;
    // id of partition to solve or virtual boundary to fuse
    const int part;

    Task* left_child{ nullptr };
    Task* right_child{ nullptr };
    Task* parent{ nullptr };

    std::vector<pm::GraphFillRegion *> regions_matched_to_virtual_boundary;
    std::vector<pm::GraphFillRegion *> regions_to_unmatch;

    // Default construct with safe initial values.
    // Task() = default;
    // Task(int task_id, long partition) : task_id(task_id), p(partition), pv(partition), leftmost_p(partition), rightmost_p(partition) {};
    Task(int task_id, int partition) : task_id(task_id), part(partition), is_fusion(false) {}
    Task(int task_id, int vb, Task* left_child, Task* right_child) 
     : task_id(task_id), part(vb), left_child(left_child), right_child(right_child), is_fusion(true) {
        left_child->parent = this;
        right_child->parent = this;
     }
    // Non-copyable due to atomic members.
    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept
     : task_id(other.task_id), part(other.part), is_fusion(other.is_fusion) {
        status.store(other.status.load());
    }
    Task& operator=(Task&& other) noexcept {
        status.store(other.status.load());
        return *this;
    }

    /* Helper Methods */
    void setup() {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        if (is_fusion) {
            for (auto& region : left_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == part)
                    regions_to_unmatch.emplace_back(region);
                else
                    regions_matched_to_virtual_boundary.emplace_back(region);
            }
            for (auto& region : right_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == part)
                    regions_to_unmatch.emplace_back(region);
                else
                    regions_matched_to_virtual_boundary.emplace_back(region);
            }
            left_child->reset();
            right_child->reset();
        }
    };

    /* Sychnorization Methods */
    void mark_solved(int shot) {
        shot_marker = shot;
        status.store(FREE, std::memory_order_release);
    }

    bool is_solved(int shot) {
        return shot_marker == shot;
    }

    // For fusion task, checks if child tasks are solved
    bool is_ready(int shot) {
        if (is_fusion) {
            return (left_child->is_solved(shot) && right_child->is_solved(shot));
        } else {
            return shot_marker < shot;
        }
    }

    // Only child tasks should try to steal their parent
    // Child tasks must mark themselves as SOLVED before trying to steal
    bool try_to_steal(int shot) {
        if (!is_ready(shot)) {
            return false;
        }
        Status expected = FREE;
        return status.compare_exchange_strong(expected, BUSY, std::memory_order_acq_rel);
    }

    void reset() {
        // ++shot_marker;
    }
};

// class WorkStealingDeque {
// private:
//     std::vector<Task> buffer;
//     std::atomic<size_t> top;
//     std::atomic<size_t> bottom;
//     const size_t capacity;

// public:
//     WorkStealingDeque(size_t cap = 1024)
//         : buffer(cap), top(0), bottom(0), capacity(cap) {}

//     // // Only the thread who owns the queue can push
//     // bool push(const Task& task) {
//     //     size_t b = bottom.load(std::memory_order_relaxed);
//     //     if (b - top.load(std::memory_order_acquire) >= capacity)
//     //         return false; // queue full
//     //     buffer[b % capacity] = task;
//     //     std::atomic_thread_fence(std::memory_order_release);
//     //     bottom.store(b + 1, std::memory_order_relaxed);
//     //     return true;
//     // }

//     // // Only the thread who owns the queue can pop
//     // bool pop(Task** task) {
//     //     size_t b = bottom.load(std::memory_order_relaxed) - 1;
//     //     bottom.store(b, std::memory_order_relaxed);
//     //     std::atomic_thread_fence(std::memory_order_seq_cst);
//     //     size_t t = top.load(std::memory_order_relaxed);
//     //     if (t <= b) {
//     //         task = buffer[b % capacity];
//     //         return true;
//     //     } else {
//     //         bottom.store(t, std::memory_order_relaxed);
//     //         return false;
//     //     }
//     // }

//     // // Other threads can steal
//     // bool steal(Task** task) {
//     //     size_t t = top.load(std::memory_order_acquire);
//     //     std::atomic_thread_fence(std::memory_order_seq_cst);
//     //     size_t b = bottom.load(std::memory_order_acquire);
//     //     if (t < b) {
//     //         task = buffer[t % capacity];
//     //         if (task.is_ready) {
//     //             if (!top.compare_exchange_strong(t, t + 1,
//     //                     std::memory_order_seq_cst,
//     //                     std::memory_order_relaxed))
//     //                 return false;
//     //             return true;
//     //         }
//     //         return false;
//     //     }
//     //     return false;
//     // }
// };


#endif // PYMATCHING2_DECODING_TASK_H