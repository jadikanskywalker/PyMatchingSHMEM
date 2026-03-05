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

#ifdef USE_SHMEM
#include <shmem.h>
#endif

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
#ifdef USE_SHMEM
    uint64_t* status_shm{ nullptr };
    uint64_t* signal_shm{ nullptr };
#endif

    // id of partition to solve or virtual boundary to fuse
    const int part;
    const bool is_fusion;
    const int vb_left, vb_right;  // Assumes round-based partitioning

#ifdef USE_SHMEM
    const size_t partition_assigned_pe{ 0 };
    const bool is_cross_pe_fusion{ false };
    size_t left_pid;
    size_t right_pid;
    bool iamleft;
#endif

    uint64_t child_bit;
    Task* left_child{nullptr};
    Task* right_child{nullptr};
    Task* parent{nullptr};

    std::vector<pm::GraphFillRegion*> regions_to_unmatch; // built from fusion children, unmatched at beginning
    std::vector<pm::GraphFillRegion*> regions_matched_to_virtual_boundary; // saved while solving

    Task(int partition
#ifdef USE_SHMEM
        , size_t assigned_pe
#endif
    )
        : part(partition),
          is_fusion(false),
#ifdef USE_SHMEM
          partition_assigned_pe(assigned_pe),
#endif
          vb_left(partition - 1),
          vb_right(partition)  // Assumes round-based partitioning
    {
        status.store(-1, std::memory_order_release);
    }
    Task(int vb, Task* left_child, Task* right_child
#ifdef USE_SHMEM
        , uint64_t* status_ptr,
        uint64_t* signal_ptr
#endif
    ) : 
#ifdef USE_SHMEM
        status_shm(status_ptr),
        signal_shm(signal_ptr),
#endif
        part(vb),
        is_fusion(true),
#ifdef USE_SHMEM
        partition_assigned_pe(left_child->partition_assigned_pe),
        is_cross_pe_fusion(left_child->partition_assigned_pe != right_child->partition_assigned_pe),
        left_pid(left_child->partition_assigned_pe),
        right_pid(right_child->partition_assigned_pe),
#endif
        vb_left(left_child->vb_left), // Assumes round-based partitioning
        vb_right(right_child->vb_right),
        left_child(left_child),
        right_child(right_child)
    {
        left_child->parent = this;
        left_child->child_bit = 1;
        right_child->parent = this;
        right_child->child_bit = 2;
#ifdef USE_SHMEM
        if (is_cross_pe_fusion && status_shm == nullptr) {
            throw std::invalid_argument("DecodingTask: Cross PE fusion task requires a symmetric status atomic");
        }
        *status_shm = 0;
        *signal_shm = 0;
#endif
    }

    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept
        : part(other.part),
          is_fusion(other.is_fusion),
#ifdef USE_SHMEM
          status_shm(other.status_shm),
          signal_shm(other.signal_shm),
          partition_assigned_pe(other.partition_assigned_pe),
          is_cross_pe_fusion(other.is_cross_pe_fusion),
          left_pid(other.left_pid),
          right_pid(other.right_pid),
#endif
          vb_left(other.vb_left),
          vb_right(other.vb_right) {
        status.store(other.status.load());
        left_child = other.left_child;
        right_child = other.right_child;
        parent = other.parent;
        child_bit = other.child_bit;
        regions_to_unmatch = std::move(other.regions_to_unmatch);
        regions_matched_to_virtual_boundary = std::move(other.regions_matched_to_virtual_boundary);
    }
    Task& operator=(Task&& other) noexcept {
#ifdef USE_SHMEM
        status_shm = other.status_shm;
        signal_shm = other.signal_shm;
        // const members assignments are invalid here, but checks are needed.
        // Since we have const members, this operator= is likely ill-formed if invoked.
        // Fortunately we might not be invoking it if we just use emplace_back on reserve.
#endif
        status.store(other.status.load());
        left_child = other.left_child;
        right_child = other.right_child;
        parent = other.parent;
        child_bit = other.child_bit;
        regions_to_unmatch = std::move(other.regions_to_unmatch);
        regions_matched_to_virtual_boundary = std::move(other.regions_matched_to_virtual_boundary);
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
    // Should be called on both PE's for a cross-PE fusion
    inline void mark_solved(
#ifdef USE_SHMEM
        size_t my_pid
#endif
    ) {
        if (is_fusion) {
#ifdef USE_SHMEM
            if (is_cross_pe_fusion) {
                shmem_uint64_atomic_set(status_shm, 0, left_pid);
                shmem_uint64_atomic_set(signal_shm, 0, my_pid);
                // shmem_quiet();
            } else {
                status.store(0, std::memory_order_release);
            }
#else
            status.store(0, std::memory_order_release);
#endif
        }
    }

    inline bool try_to_steal_leaf(int next) {
        int expected = next - 1;
        return status.compare_exchange_strong(expected, next, std::memory_order_acq_rel);
    }

    inline bool try_to_steal_parent(
#ifdef USE_SHMEM
        size_t my_pid
#endif
    ) {
        int old;
#ifdef USE_SHMEM
        if (parent->is_cross_pe_fusion) {
            old = shmem_uint64_atomic_fetch_or(parent->status_shm, child_bit, parent->left_pid);
        } else {
#endif
            old = parent->status.fetch_or(child_bit, std::memory_order_acq_rel);
#ifdef USE_SHMEM
        }
#endif
        if ((old | child_bit) == 3) {
            return old != 3;
        } else {
            return false;
        }
    }

    // inline Task* try_to_steal_parent_or_descendent(int next) {
    //     int old = parent->status.fetch_or(child_bit, std::memory_order_acq_rel);
    //     if ((old | child_bit) == 3) {
    //         if (old == 3) {  // got beat
    //             return nullptr;
    //         } else {  // got parent
    //             return parent;
    //         }
    //     }
    //     // I am first to try to steal parent, try sibling
    //     Task* sibling = (child_bit == 1) ? parent->right_child : parent->left_child;
    //     return sibling->try_to_steal_descendent(next);
    // }

    // inline Task* try_to_steal_descendent(int next) {
    //     Task* t = nullptr;
    //     if (!is_fusion) {  // I am leaf
    //         bool stolen = try_to_steal_leaf(next);
    //         if (stolen) {
    //             t = this;
    //         }
    //     } else {                                               // Try to steal descendent
    //         if (status.load(std::memory_order_acquire) > 0) {  // Let other thread have it
    //             return nullptr;
    //         }
    //         t = left_child->try_to_steal_descendent(next);
    //         if (t == nullptr) {
    //             t = right_child->try_to_steal_descendent(next);
    //         }
    //     }
    //     return t;
    // }
};

#endif  // PYMATCHING2_DECODING_TASK_H