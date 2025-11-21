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

enum Status { UNSOLVED, IN_PROGRESS, SOLVED };

struct Task {
private:
    std::atomic<Status> status{ Status::UNSOLVED };

public:
    int task_id;
    bool is_fusion{ 0 };

    // long p{ -1 };
    // long pv{ -1 };
    // long leftmost_p{ -1 };
    // long rightmost_p{ -1 };

    // id of partition to solve or virtual boundary to fuse
    int part;

    Task* left_child{ nullptr };
    Task* right_child{ nullptr };
    Task* parent{ nullptr };

    std::vector<pm::GraphFillRegion *> regions_matched_to_virtual_boundary;
    std::vector<pm::GraphFillRegion *> regions_to_unmatch;

    // Default construct with safe initial values.
    Task() = default;
    // Task(int task_id, long partition) : task_id(task_id), p(partition), pv(partition), leftmost_p(partition), rightmost_p(partition) {};
    Task(int task_id, int partition) : task_id(task_id) {
        part = partition;
    }
    // Task(long partition_without_virtuals, long partition_with_virtuals)
    //   : p(partition_without_virtuals), pv(partition_with_virtuals), p_leftmost(p), p_rightmost(pv), is_fusion(true) {};
    // Task(int t_id, Task* left_child, Task* right_child) : task_id(t_id), left_child(left_child), right_child(right_child), is_fusion(1) {
    //     p = left_child->rightmost_p;
    //     pv = right_child->leftmost_p;
    //     leftmost_p = left_child->leftmost_p;
    //     rightmost_p = right_child->rightmost_p;
    //     left_child->parent = this;
    //     right_child->parent = this;
    //     // status.store(Status::UNSOLVED, std::memory_order_relaxed);
    // }
    Task(int t_id, Task* left_child, Task* right_child, int vb) 
     : task_id(task_id), left_child(left_child), right_child(right_child), is_fusion(true) {
        part = vb;
    }
    // Non-copyable due to atomic members.
    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept {
        status.store(other.status.load());
    }
    Task& operator=(Task&& other) noexcept {
        status.store(other.status.load());
        return *this;
    }

    /* Helper Methods */
    // void setup(long partition) {
    //     p = partition;
    //     pv = partition;
    //     is_fusion = 0;
    // }
    // void setup(Task* left, Task* right) {
    //     left_child = left;
    //     right_child = right;
    //     left_child->parent = this;
    //     right_child->parent = this;
    //     p = left_child->get_rightmost_p();
    //     pv = right_child->get_leftmost_p();
    //     is_fusion = 1;
    // }
    void setup_regions() {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        if (is_fusion) {
            for (auto& region : left_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->partition == pv)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
            for (auto& region : right_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->partition == pv)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
        }
    };

    /* Sychnorization Methods */
    void mark_solved() {
        status.store(SOLVED, std::memory_order_release);
    }

    bool is_solved() {
        return status.load(std::memory_order_acquire) == SOLVED;
    }

    // For fusion task, checks if child tasks are solved
    bool is_ready() {
        return (left_child->is_solved() && right_child->is_solved());
    }

    // Only child tasks should try to steal their parent
    // Child tasks must mark themselves as SOLVED before trying to steal
    bool try_to_steal() {
        if (is_fusion && !is_ready())
            return false;
        Status expected = UNSOLVED;
        return status.compare_exchange_strong(expected, IN_PROGRESS, std::memory_order_acq_rel);
    }

    void reset() {
        status.store(UNSOLVED, std::memory_order_release);
    }
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