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

#ifdef SCOREP_USER_ENABLE
#include <scorep/SCOREP_User.h>
#endif

#ifdef USE_SHMEM
#include <iostream>
#include <shmem.h>
#endif

// Forward declare to avoid cyclic include with mwpm_decoding.h
namespace pm {
class GraphFillRegion;
struct FusionSummary;
struct Mwpm;
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

#ifdef USE_SHMEM
    // SHMEM-only: whether this node's parent in the chain is a CrossRankTask. Meaningless without
    // cross-rank fusion, so not worth carrying under plain USE_THREADS.
    bool is_cross_rank_fusion;
#endif

    // Parent in the task chain: local Task fusion or CrossRankTask above this node.
    // nullptr means this node is the chain top (a task graph root).
    TaskBase* parent{nullptr};

    // The solver this task uses when it resolves, chosen once at construction time (see
    // decoding_unit.cc's task-building functions) rather than recomputed from part/vb_solver_offset
    // at use time -- eliminates a whole class of "same partition index, computed two different ways,
    // disagreeing" bugs. Raw, non-owning pointer: lifetime is guaranteed by DecodingUnit's own
    // solvers/remote-arena storage outliving every task built from its shot_buffer.
    pm::Mwpm* solver{nullptr};

    std::vector<pm::GraphFillRegion*> regions_to_unmatch;
    std::vector<pm::GraphFillRegion*> regions_matched_to_virtual_boundary;

    TaskBase(int part, int vb_left, int vb_right, bool is_fusion, bool is_cross_rank_fusion, pm::Mwpm* solver)
        : part(part), vb_left(vb_left), vb_right(vb_right), is_fusion(is_fusion)
#ifdef USE_SHMEM
        , is_cross_rank_fusion(is_cross_rank_fusion)
#endif
        , solver(solver)
    {
        (void)is_cross_rank_fusion;
        if (is_fusion)
            vb_marker = part;
        else
            vb_marker = -1;
    }

    virtual ~TaskBase() = default;

    virtual void setup() = 0;
    virtual void mark_solved(size_t my_pid = 0) = 0;
    virtual bool try_to_steal(size_t val) = 0;
    virtual void reset() = 0;

    TaskBase(TaskBase&& other) noexcept
        : part(other.part), vb_marker(other.vb_marker),
          vb_left(other.vb_left), vb_right(other.vb_right),
          is_fusion(other.is_fusion),
#ifdef USE_SHMEM
          is_cross_rank_fusion(other.is_cross_rank_fusion),
#endif
          parent(other.parent),
          solver(other.solver),
          regions_to_unmatch(std::move(other.regions_to_unmatch)),
          regions_matched_to_virtual_boundary(std::move(other.regions_matched_to_virtual_boundary))
    { other.parent = nullptr; }

    TaskBase& operator=(TaskBase&& other) noexcept {
        part = other.part;
        vb_marker = other.vb_marker;
        vb_left = other.vb_left;
        vb_right = other.vb_right;
        is_fusion = other.is_fusion;
#ifdef USE_SHMEM
        is_cross_rank_fusion = other.is_cross_rank_fusion;
#endif
        parent = other.parent;
        other.parent = nullptr;
        solver = other.solver;
        regions_to_unmatch = std::move(other.regions_to_unmatch);
        regions_matched_to_virtual_boundary = std::move(other.regions_matched_to_virtual_boundary);
        return *this;
    }
};

#ifdef USE_SHMEM
struct CrossRankTask;  // forward declare: Task::left_crt_child needs the pointer type, CrossRankTask
                        // is defined later in this file (and only under USE_SHMEM)
#endif

struct Task : public TaskBase {
   private:
    // for partition solves,
    //    status = id of last shot for which task was claimed
    // for fusions,
    //    status = number of children that have arrived so far this shot (0 = unclaimed). Every
    //             child does the identical fetch_add(1); whoever's fetch_add returns N-1 (last of
    //             N arrivals, N=2 for today's binary fusions) is the winner. Generalizes cleanly to
    //             any future N-way fusion with no per-child bit-position bookkeeping, unlike the
    //             fetch_or/child_bit scheme this replaced.
    // Thus, fusions require resetting status to 0 in mark_solved
    alignas(64) std::atomic<int64_t> status{0};

   public:

    size_t child_bit;
    Task* left_child{nullptr};
    Task* right_child{nullptr};

    // Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md). Only meaningful
    // when is_fusion == true. Two roles, mutually exclusive:
    //   is_extraction_unit_connector: this fusion joins two already-independent units/unit-groups
    //     (a chain link when extract_preemptively, or an inter-unit fusion in the balanced-over-units
    //     tree otherwise). Its vb must be divided before either side is touched independently.
    //   is_extraction_unit_root: this is the top of one standalone extraction unit (a raw leaf, or the
    //     root of a unit's own internal balanced subtree) -- where posting/recursion stops.
    // Both purely tags at construction time, consumed by the extraction-queue logic added separately.
    bool is_extraction_unit_connector{false};
    bool is_extraction_unit_root{false};

    // Preemptive OBS extraction (deferred-chain-with-embedded-seams design, see
    // now-its-time-to-jazzy-turing.md): true on an is_extraction_unit_connector chain-link whose
    // newly-closed unit falls inside an active seam's [unit_lo, unit_hi) span -- tells the decode-loop
    // phase to skip the immediate divide+post (and skip the usual vb_left restriction) until the seam
    // itself resolves and walks back to divide+post every deferred unit at once. False (default) for
    // every ROUND-partitioning fusion and every OBS chain-link outside a seam's span, unchanged from
    // today's immediate-divide behavior.
    bool defer_division{false};

    // bool only_child{ false };
#ifdef USE_SHMEM
    // int seam_vb_slot{ -1 };             // virtual_boundaries slot index for this seam's VB nodes
    int left_obs_patch_id{-1};   // For OBS local seam tasks: left obs patch
    int right_obs_patch_id{-1};  // For OBS local seam tasks: right obs patch

    // A seam fusion is the left child of two different downstream fusions (one continuing each
    // observable's own chain past the seam) -- TaskBase::parent only holds one. By construction-order
    // convention (build the right observable's continuation first, copy its auto-set parent here, then
    // build the left observable's continuation, which overwrites parent to the left value), `parent`
    // always ends up meaning "left observable's continuation" and this field the right's. Doubles as
    // the universally-available "is this a seam fusion" tag (right_obs_parent != nullptr) -- left_obs_
    // patch_id/right_obs_patch_id above only exist under ENABLE_DRAW_FLAGS, not general enough for this.
    Task* right_obs_parent{nullptr};

    // The CRT case needs a fusion whose left operand is a CrossRankTask*, but left_child/right_child
    // above are strictly typed Task* (CrossRankTask is a sibling of Task, not a subclass -- both derive
    // from TaskBase only). Only ever populated on the one fusion immediately following a CRT in a
    // single-local-observable deferred chain (never on right_child's side).
    CrossRankTask* left_crt_child{nullptr};
#endif

    Task(int partition, pm::Mwpm* solver)
        : TaskBase(partition, partition - 1, partition, false, false, solver)
    {
        status.store(-1, std::memory_order_release);
    }
    // Partition leaf with explicit local vb bounds (needed for OBS partitioning)
    Task(int part, int vb_l, int vb_r, pm::Mwpm* solver)
        : TaskBase(part, vb_l, vb_r, false, false, solver)
    {
        status.store(-1, std::memory_order_release);
    }
    Task(int vb, Task* left_child, Task* right_child, pm::Mwpm* solver) :
        TaskBase(vb, left_child->vb_left, right_child->vb_right, true, false, solver),
        left_child(left_child),
        right_child(right_child)
    {
        left_child->parent  = this;
        left_child->child_bit  = 1;
        right_child->parent = this;
        right_child->child_bit = 2;
    }

#ifdef USE_SHMEM
    // The one fusion immediately following a CRT in a single-local-observable deferred chain (see
    // left_crt_child above): left operand is a CrossRankTask*, not a Task*. Declared here, defined out
    // of line after CrossRankTask's own definition later in this file -- CrossRankTask is only
    // forward-declared at this point, so its inherited TaskBase members (vb_left, parent) aren't
    // accessible yet.
    Task(int vb, CrossRankTask* left_crt_child, Task* right_child, pm::Mwpm* solver);
#endif

    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept
        : TaskBase(std::move(other))
           {
        status.store(other.status.load());
        left_child = other.left_child;
        right_child = other.right_child;
        // parent is in TaskBase and moved by TaskBase(std::move(other))
        child_bit = other.child_bit;
        is_extraction_unit_connector = other.is_extraction_unit_connector;
        is_extraction_unit_root = other.is_extraction_unit_root;
        defer_division = other.defer_division;
#ifdef USE_SHMEM
        // seam_vb_slot = other.seam_vb_slot;
        left_obs_patch_id = other.left_obs_patch_id;
        right_obs_patch_id = other.right_obs_patch_id;
        right_obs_parent = other.right_obs_parent;
        left_crt_child = other.left_crt_child;
#endif
    }
    Task& operator=(Task&& other) noexcept {
        TaskBase::operator=(std::move(other));
        status.store(other.status.load());
        left_child = other.left_child;
        right_child = other.right_child;
        // parent is in TaskBase and moved by TaskBase::operator=(std::move(other))
        child_bit = other.child_bit;
        is_extraction_unit_connector = other.is_extraction_unit_connector;
        is_extraction_unit_root = other.is_extraction_unit_root;
        defer_division = other.defer_division;
#ifdef USE_SHMEM
        // seam_vb_slot = other.seam_vb_slot;
        left_obs_patch_id = other.left_obs_patch_id;
        right_obs_patch_id = other.right_obs_patch_id;
        right_obs_parent = other.right_obs_parent;
        left_crt_child = other.left_crt_child;
#endif
        return *this;
    }

    /* Helper Methods */
    void setup() override {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        if (is_fusion) {
            for (auto& region : left_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
            if (left_child == right_child)
                return; // prevent double copy
            for (auto& region : right_child->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker)
                    regions_to_unmatch.push_back(region);
                else
                    regions_matched_to_virtual_boundary.push_back(region);
            }
        }
    };

    /* Sychnorization Methods */
    void mark_solved(size_t my_pid = 0) override {
        (void)my_pid;
        if (is_fusion) {
            status.store(0, std::memory_order_release);
        }
    }

    // for partition leaf, val is the shot id
    // for fusion parent, val is unused -- kept only for interface uniformity with the leaf-claim/
    // CrossRankTask overloads of try_to_steal (see decoding_task.h Design §3 note on Task::status
    // above). The fetch_add-based race needs no per-caller value the way the old fetch_or/child_bit
    // scheme did: every child does the identical fetch_add(1), and whoever's fetch_add returns N-1
    // (last of N arrivals) is the winner.
    bool try_to_steal(size_t val) override {
        if (!is_fusion) { // partition
            int64_t expected = static_cast<int64_t>(val) - 1;
            return status.compare_exchange_strong(expected, static_cast<int64_t>(val), std::memory_order_acq_rel);
        } else { // fusion
#ifdef USE_SHMEM
            if (left_child == right_child) { // one child
                return true;
            }
#endif
            (void)val;
            constexpr int64_t N = 2; // ordinary binary fusion; generalizes to N>2 with no other change
            int64_t old = status.fetch_add(1, std::memory_order_acq_rel);
            return old == N - 1;
        }
    }

    // inline bool try_to_steal_leaf(int next) {
    //     return try_to_steal(static_cast<size_t>(next));
    // }

//     inline bool try_to_steal_parent() {
// #ifdef USE_SHMEM
//         if (only_child) {
//             return true;
//         }
// #endif
//         // parent is guaranteed to be a local Task (caller checks !parent->is_cross_rank_fusion)
//         int old = static_cast<Task*>(parent)->status.fetch_or(child_bit, std::memory_order_acq_rel);
//         if ((old | child_bit) == 3) {
//             return old != 3;
//         } else {
//             return false;
//         }
//     }

    void reset() override {
        status.store((is_fusion) ? 0 : -1, std::memory_order_release);
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
    Task* child; // Want to change this?
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
        pm::FusionSummary* fusion_summary_ptr,
        pm::Mwpm* solver
    ) : TaskBase(vb, vb_left, vb_right, true, true, solver),
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
        // Insert this CRT above the current chain top of child's parent chain.
        // This correctly handles multiple CRTs for the same obs group (they chain sequentially).
        TaskBase* top = child;
        while (top->parent != nullptr) top = top->parent;
        top->parent = this;
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
    void setup() override {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        for (auto& region : child->regions_matched_to_virtual_boundary) {
            if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker) // seam vbs are marked with global part
                regions_to_unmatch.push_back(region);
            // else
            //     regions_matched_to_virtual_boundary.push_back(region);
        }
    };


    /* Sychronization Methods */
    // Should be called on PE who solved fusion
    void mark_solved(size_t my_pid) override {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_BEGIN();
#endif
        shmem_ctx_uint64_atomic_set(context_shm, status_shm, 0, (iamleft) ? my_pid : other_pid);
        shmem_ctx_uint64_atomic_set(context_shm, signal_shm, 0, my_pid); // reset my signal
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_END();
#endif
    }

    inline void report_done(int shot_buffer_round) {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_BEGIN();
#endif
        shmem_ctx_uint64_atomic_set(context_shm, done_shm, shot_buffer_round+1, other_pid); // notify other PE we are done
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_END();
#endif
    }

    inline void wait_until_done(size_t my_pid, int shot_buffer_round) {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_BEGIN();
#endif
        shmem_wait_until(done_shm, SHMEM_CMP_GE, shot_buffer_round + 1);
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_END();
#endif
    }

    bool try_to_steal(size_t my_pid
        // , std::ostream* t_out = nullptr
    ) override {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_BEGIN();
#endif
        uint64_t old, news;
        if (iamleft) {
            // if (t_out && DEBUG) *t_out << "  performing fetch_or on " << status_shm << " my_pid=" << my_pid << " with child_bit=" << 1 << std::endl << std::flush;
            old = shmem_ctx_uint64_atomic_fetch_or(context_shm, status_shm, 1, my_pid);
            news = old | 1;
        } else {
            // if (t_out && DEBUG) *t_out << "  performing fetch_or on " << status_shm << " other_pid=" << other_pid << " with child_bit=" << 2 << std::endl << std::flush;
            old = shmem_ctx_uint64_atomic_fetch_or(context_shm, status_shm, 2, other_pid);
            news = old | 2;
        }
        // if (t_out && DEBUG) *t_out << "  done with fetch_or" << std::endl << std::flush;
        if (news == 3) {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_END();
#endif
            return old != 3;
        } else {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_END();
#endif
            return false;
        }
    }

    void reset() override {
        *status_shm = 0;
        *signal_shm = 0;
        *done_shm = 0;
    }

};

// Out-of-line: CrossRankTask must be a complete type for left_crt_child->vb_left/parent below (only
// forward-declared at the point of Task's own in-class declaration, see left_crt_child's comment).
inline Task::Task(int vb, CrossRankTask* left_crt_child, Task* right_child, pm::Mwpm* solver) :
    TaskBase(vb, left_crt_child->vb_left, right_child->vb_right, true, false, solver),
    left_crt_child(left_crt_child),
    right_child(right_child)
{
    left_crt_child->parent = this;
    right_child->parent = this;
    right_child->child_bit = 2;
}
#endif

#endif  // PYMATCHING2_DECODING_TASK_H