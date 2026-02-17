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
namespace pm {
class GraphFillRegion;
}

enum Status { BUSY, FREE };

struct Task {
   private:
    // for partition solves,
    //    status = id of last shot for which task was claimed
    // for fusions,
    //    status = 0 (Unclaimed for current shot)
    //             1 (Right child solved and tried to steal me first)
    //             2 (Left child solved and tried to steal me first)
    //             3 (Both children solved and second one who tried to steal me won)
    // Thus, fusions require resetting status to 0 in mark_solved
    std::atomic<int> status{0};

   public:
    const int task_id;
    // id of partition to solve or virtual boundary to fuse
    const int part;
    const bool is_fusion;

#ifdef USE_SHMEM
    bool is_cross_pe{ false };
#endif

    const int vb_left, vb_right;  // Assumes round-based partitioning

    int child_bit;
    Task* left_child{nullptr};
    Task* right_child{nullptr};
    Task* parent{nullptr};

    std::vector<pm::GraphFillRegion*> regions_to_unmatch; // built from fusion children, unmatched at beginning
    std::vector<pm::GraphFillRegion*> regions_matched_to_virtual_boundary; // saved while solving

    Task(int task_id, int partition)
        : task_id(task_id),
          part(partition),
          is_fusion(false),
          vb_left(partition - 1),
          vb_right(partition)  // Assumes round-based partitioning
    {
        status.store(-1, std::memory_order_release);
    }
    Task(int task_id, int vb, Task* left_child, Task* right_child)
        : task_id(task_id),
          part(vb),
          left_child(left_child),
          right_child(right_child),
          is_fusion(true),
          vb_left(left_child->vb_left), // Assumes round-based partitioning
          vb_right(right_child->vb_right)
    {
        left_child->parent = this;
        left_child->child_bit = 1;
        right_child->parent = this;
        right_child->child_bit = 2;
    }

    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept
        : task_id(other.task_id),
          part(other.part),
          is_fusion(other.is_fusion),
          vb_left(other.vb_left),
          vb_right(other.vb_right) {
        status.store(other.status.load());
    }
    Task& operator=(Task&& other) noexcept {
        status.store(other.status.load());
        return *this;
    }

    /* Helper Methods */
    inline void setup() {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        if (is_fusion) {
            for (auto& region : left_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to->vb == part)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
            for (auto& region : right_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to->vb == part)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
        }
    };

    /* Sychnorization Methods */
    inline void mark_solved() {
        if (is_fusion) {
            status.store(0, std::memory_order_release);
        }
    }

    inline bool try_to_steal_leaf(int next) {
        int expected = next - 1;
        return status.compare_exchange_strong(expected, next, std::memory_order_acq_rel);
    }

    inline bool try_to_steal_parent() {
        int old = parent->status.fetch_or(child_bit, std::memory_order_acq_rel);
        if ((old | child_bit) == 3) {
            return old != 3;
        } else {
            return false;
        }
    }

    inline Task* try_to_steal_parent_or_descendent(int next) {
        int old = parent->status.fetch_or(child_bit, std::memory_order_acq_rel);
        if ((old | child_bit) == 3) {
            if (old == 3) {  // got beat
                return nullptr;
            } else {  // got parent
                return parent;
            }
        }
        // I am first to try to steal parent, try sibling
        Task* sibling = (child_bit == 1) ? parent->right_child : parent->left_child;
        return sibling->try_to_steal_descendent(next);
    }

    inline Task* try_to_steal_descendent(int next) {
        Task* t = nullptr;
        if (!is_fusion) {  // I am leaf
            bool stolen = try_to_steal_leaf(next);
            if (stolen) {
                t = this;
            }
        } else {                                               // Try to steal descendent
            if (status.load(std::memory_order_acquire) > 0) {  // Let other thread have it
                return nullptr;
            }
            t = left_child->try_to_steal_descendent(next);
            if (t == nullptr) {
                t = right_child->try_to_steal_descendent(next);
            }
        }
        return t;
    }
};

#endif  // PYMATCHING2_DECODING_TASK_H