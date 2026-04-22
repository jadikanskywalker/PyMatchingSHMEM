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
#include <iostream>
#include <shmem.h>
#endif

// Forward declare to avoid cyclic include with mwpm_decoding.h
namespace pm {
class GraphFillRegion;
struct FusionSummary;
}

enum Status { BUSY, FREE };

struct TaskBase {
    int part;
// #ifdef USE_SHMEM
    int vb_marker; // not always equal to part
// #endif
    int vb_left;
    int vb_right;
    bool is_fusion;

    bool is_cross_rank_fusion;

    std::vector<pm::GraphFillRegion*> regions_to_unmatch;
    std::vector<pm::GraphFillRegion*> regions_matched_to_virtual_boundary;

    TaskBase(int part, int vb_left, int vb_right, bool is_fusion, bool is_cross_rank_fusion) 
        : part(part), vb_marker(part), vb_left(vb_left), vb_right(vb_right), is_fusion(is_fusion), is_cross_rank_fusion(is_cross_rank_fusion) {
// #ifdef USE_SHMEM
            // vb_marker = part;
// #endif
    }

    virtual ~TaskBase() = default;

    // TaskBase(const TaskBase&) = delete;
    // TaskBase& operator=(const TaskBase&) = delete;

    // TaskBase(TaskBase&& other) noexcept 
    //     : part(other.part), vb_left(other.vb_left), vb_right(other.vb_right), is_fusion(other.is_fusion),
    //       regions_to_unmatch(std::move(other.regions_to_unmatch)),
    //       regions_matched_to_virtual_boundary(std::move(other.regions_matched_to_virtual_boundary))
    // {}

    // TaskBase& operator=(TaskBase&& other) noexcept {
    //     regions_to_unmatch = std::move(other.regions_to_unmatch);
    //     regions_matched_to_virtual_boundary = std::move(other.regions_matched_to_virtual_boundary);
    //     return *this;
    // }
};

struct Task : public TaskBase {
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

    uint64_t child_bit;
    Task* left_child{nullptr};
    Task* right_child{nullptr};
    Task* parent{nullptr};

#ifdef USE_SHMEM
    bool only_child{ false };
    int vb_solver_offset{ 0 }; // For OBS patch i, this is i (Because no vb between patches, we lose one vb index relative to partition index)
    // int seam_vb_slot{ -1 };             // virtual_boundaries slot index for this seam's VB nodes
    int left_obs_patch_id{-1};   // For OBS local seam tasks: left obs patch
    int right_obs_patch_id{-1};  // For OBS local seam tasks: right obs patch
#endif

    Task(int partition)
        : TaskBase(partition, partition - 1, partition, false, false)
    {
        status.store(-1, std::memory_order_release);
    }
    // Partition leaf with explicit local vb bounds (needed for OBS partitioning)
    Task(int part, int vb_l, int vb_r)
        : TaskBase(part, vb_l, vb_r, false, false)
    {
        status.store(-1, std::memory_order_release);
    }
    Task(int vb, Task* left_child, Task* right_child) : 
        TaskBase(vb, left_child->vb_left, right_child->vb_right, true, false),
        left_child(left_child),
        right_child(right_child)
    {
        left_child->parent = this;
        left_child->child_bit = 1;
        right_child->parent = this;
        right_child->child_bit = 2;
    }

    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept
        : TaskBase(std::move(other))
           {
        status.store(other.status.load());
        left_child = other.left_child;
        right_child = other.right_child;
        parent = other.parent;
        // parent_right = other.parent_right;
        child_bit = other.child_bit;
#ifdef USE_SHMEM
        // seam_vb_slot = other.seam_vb_slot;
        left_obs_patch_id = other.left_obs_patch_id;
        right_obs_patch_id = other.right_obs_patch_id;
#endif
    }
    Task& operator=(Task&& other) noexcept {
        TaskBase::operator=(std::move(other));
        status.store(other.status.load());
        left_child = other.left_child;
        right_child = other.right_child;
        parent = other.parent;
        // parent_right = other.parent_right;
        child_bit = other.child_bit;
#ifdef USE_SHMEM
        // seam_vb_slot = other.seam_vb_slot;
        left_obs_patch_id = other.left_obs_patch_id;
        right_obs_patch_id = other.right_obs_patch_id;
#endif
        return *this;
    }

    /* Helper Methods */
    inline void setup() {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
//         int part_cmp = part;
// #ifdef USE_SHMEM
//         if (seam_vb_slot >= 0 ) part_cmp = seam_vb_slot;
// #endif
        if (is_fusion) {
            for (auto& region : left_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
            for (auto& region : right_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker)
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
        if (is_fusion) return false;
        int expected = next - 1;
        return status.compare_exchange_strong(expected, next, std::memory_order_acq_rel);
    }

    inline bool try_to_steal_parent() {
#ifdef USE_SHMEM
        if (only_child) {
            return true;
        }
#endif
        int old = parent->status.fetch_or(child_bit, std::memory_order_acq_rel);
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

#ifdef USE_SHMEM
struct CrossRankTask : public TaskBase {
public:
    Task* child;
    bool iamleft;

    size_t other_pid{ 0 };

    // SHMEM resources
    uint64_t* status_shm{ nullptr };
    uint64_t* signal_shm{ nullptr };
    uint64_t* done_shm{ nullptr };
    pm::FusionSummary* fusion_summary_shm { nullptr };
    shmem_ctx_t context_shm;
    bool owns_context{ true };

    // For OBS patch fusions
    std::pair<size_t, size_t> left_global_offset{  0, 0 };   // obsA (lower obs) {p_offset, vb_offset}
    std::pair<size_t, size_t> right_global_offset{ 0, 0 };   // obsB (higher obs) {p_offset, vb_offset}
    // int seam_vb_slot{ -1 };             // virtual_boundaries slot index for this seam's VB nodes

    CrossRankTask(
        int vb,
        Task* child,
        bool iamleft,
        int vb_left,
        int vb_right,
        size_t other_pid,
        uint64_t* status_ptr,
        uint64_t* signal_ptr,
        uint64_t* done_ptr,
        pm::FusionSummary* fusion_summary_ptr
    ) : TaskBase(vb, vb_left, vb_right, true, true),
        child(child),
        iamleft(iamleft),
        other_pid(other_pid),
        status_shm(status_ptr),
        signal_shm(signal_ptr),
        done_shm(done_ptr),
        fusion_summary_shm(fusion_summary_ptr)
    {
        if (status_shm == nullptr) {
            throw std::invalid_argument("DecodingTask: Cross PE fusion task requires a symmetric status atomic");
        }
        // Create communication context
        if (shmem_ctx_create(SHMEM_CTX_SERIALIZED, &context_shm)) {
            std::cout << "PE" << shmem_my_pe() << " DecodingTask: Cross PE fusion task failed to create context\n" << std::flush;
            context_shm = SHMEM_CTX_DEFAULT;
        }
        // context_shm = SHMEM_CTX_DEFAULT;
        *status_shm = 0;
        *signal_shm = 0;
        *done_shm = 0;
    }

    CrossRankTask(CrossRankTask&& other) noexcept : TaskBase(std::move(other)) {
        if (DEBUG) std::cout << "PE" << shmem_n_pes() << " CrossRankTask move constructor was called" << std::endl << std::flush;
        child = other.child;
        iamleft = other.iamleft;
        other_pid = other.other_pid;
        status_shm = other.status_shm;
        signal_shm = other.signal_shm;
        done_shm = other.done_shm;
        fusion_summary_shm = other.fusion_summary_shm;
        context_shm = other.context_shm;
        owns_context = other.owns_context;
        left_global_offset  = other.left_global_offset;
        right_global_offset = other.right_global_offset;
        // seam_vb_slot        = other.seam_vb_slot;
        
        // Nullify other's ownership so destructor skips it
        other.owns_context = false;
    }

    ~CrossRankTask() {
        if (owns_context && context_shm != SHMEM_CTX_DEFAULT)
            shmem_ctx_destroy(context_shm);
    }


    /* Helper Methods */
    inline void setup() {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        for (auto& region : child->regions_matched_to_virtual_boundary) {
            if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker) // seam vbs are marked with global part
                regions_to_unmatch.push_back(region);
            // else
            //     regions_matched_to_virtual_boundary.push_back(region);
        }
    };


    /* Sychnorization Methods */
    // Should be called on PE who solved fusion
    inline void mark_solved(size_t my_pid) {
        shmem_ctx_uint64_atomic_set(context_shm, status_shm, 0, (iamleft) ? my_pid : other_pid);
        shmem_ctx_uint64_atomic_set(context_shm, signal_shm, 0, my_pid); // reset my signal
    }

    inline void report_done(int shot_buffer_round) {
        shmem_ctx_uint64_atomic_set(context_shm, done_shm, shot_buffer_round+1, other_pid); // notify other PE we are done
    }

    inline void wait_until_done(size_t my_pid, int shot_buffer_round) {
        shmem_wait_until(done_shm, SHMEM_CMP_GE, shot_buffer_round + 1);
    }

    inline bool try_to_steal(size_t my_pid, std::ostream& t_out) {
        int old, news;
        if (iamleft) {
            if (DEBUG) t_out << "  performing fetch_or on " << status_shm << " my_pid=" << my_pid << " with child_bit=" << 1 << std::endl << std::flush;
            old = shmem_ctx_uint64_atomic_fetch_or(context_shm, status_shm, 1, my_pid);
            news = old | 1;
        } else {
            if (DEBUG) t_out << "  performing fetch_or on " << status_shm << " other_pid=" << other_pid << " with child_bit=" << 2 << std::endl << std::flush;
            old = shmem_ctx_uint64_atomic_fetch_or(context_shm, status_shm, 2, other_pid);
            news = old | 2;
        }
        if (DEBUG) t_out << "  done with fetch_or" << std::endl << std::flush;
        if (news == 3) {
            return old != 3;
        } else {
            return false;
        }
    }

};
#endif

#endif  // PYMATCHING2_DECODING_TASK_H