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

#include <algorithm>
#include <atomic>
#include <immintrin.h>
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

// Task graph model (see now-its-time-to-jazzy-turing.md Phase 2 (REVISED)): the binary tree
// (Task::parent/left_child/right_child) is *always* a plain, ordinary tree -- seams and cross-rank
// fusions never sit in it. Instead, each Task carries a vector (special_tasks) of LocalSeamTask
// and/or CrossRankTask instances that trigger when that Task's own newly-closed unit is reached.
// Decode logic processes a solved Task's own special_tasks (in order), then continues the ordinary
// climb via Task::parent -- seam/CRT handling is a self-contained detour, not a structural rewiring
// of the tree.
struct TaskBase {
    int part;
    int vb_marker; // not always equal to part
    int vb_left;
    int vb_right;
    bool is_fusion;



    // Which concrete kind this task is. Replaces the old is_cross_rank_fusion bool with a 3-way tag --
    // scoped (enum class) so the enumerator names can mirror the actual class names (TaskType::Task,
    // TaskType::LocalSeamTask, TaskType::CrossRankTask) without colliding with those class names in
    // the enclosing namespace. Not USE_SHMEM-gated: LocalSeamTask needs it in every build, not just
    // SHMEM ones (only the CrossRankTask value ever goes unused outside USE_SHMEM).
    enum class TaskType { Task, LocalSeamTask, CrossRankTask };
    TaskType type;

    // The solver this task uses when it resolves, chosen once at construction time (see
    // decoding_unit.cc's task-building functions) rather than recomputed from part/vb_solver_offset
    // at use time -- eliminates a whole class of "same partition index, computed two different ways,
    // disagreeing" bugs. Raw, non-owning pointer: lifetime is guaranteed by DecodingUnit's own
    // solvers/remote-arena storage outliving every task built from its shot_buffer.
    pm::Mwpm* solver{nullptr};

    std::vector<pm::GraphFillRegion*> regions_to_unmatch;
    std::vector<pm::GraphFillRegion*> regions_matched_to_virtual_boundary;

    TaskBase(int part, int vb_left, int vb_right, bool is_fusion, TaskType type, pm::Mwpm* solver)
        : part(part), vb_left(vb_left), vb_right(vb_right), is_fusion(is_fusion), type(type), solver(solver)
    {
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

    // Convenience accessors replacing the old stored is_cross_rank_fusion field -- method calls now,
    // not field reads, so every existing call site needs the added parens (mechanical, task-building/
    // decode-logic pass concern, not this one).
    inline bool is_task() const { return type == TaskType::Task; }
    inline bool is_cross_rank_fusion() const { return type == TaskType::CrossRankTask; }
    inline bool is_local_seam_fusion() const { return type == TaskType::LocalSeamTask; }

    TaskBase(TaskBase&& other) noexcept
        : part(other.part), vb_marker(other.vb_marker),
          vb_left(other.vb_left), vb_right(other.vb_right),
          is_fusion(other.is_fusion), type(other.type),
          solver(other.solver),
          regions_to_unmatch(std::move(other.regions_to_unmatch)),
          regions_matched_to_virtual_boundary(std::move(other.regions_matched_to_virtual_boundary))
    {}

    TaskBase& operator=(TaskBase&& other) noexcept {
        part = other.part;
        vb_marker = other.vb_marker;
        vb_left = other.vb_left;
        vb_right = other.vb_right;
        is_fusion = other.is_fusion;
        type = other.type;
        solver = other.solver;
        regions_to_unmatch = std::move(other.regions_to_unmatch);
        regions_matched_to_virtual_boundary = std::move(other.regions_matched_to_virtual_boundary);
        return *this;
    }
};

// Intermediate base shared by LocalSeamTask and CrossRankTask -- both attach to a triggering Task's
// own special_tasks vector rather than participating in the binary tree. Task::special_tasks' own
// element type is SpecialTask*, not TaskBase*, precisely because it can only ever point at one of
// these two kinds, never a plain Task.
struct SpecialTask : public TaskBase {
    using TaskBase::TaskBase;
};

// Forward declarations so Task::add_special_task can be declared here and defined below, once
// LocalSeamTask/CrossRankTask are complete types it needs to reach into (.children / .child).
struct LocalSeamTask;
#ifdef USE_SHMEM
struct CrossRankTask;
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

    // The binary tree is always a plain, ordinary tree of Task nodes -- seams/CRTs never appear here
    // (see special_tasks below). Task*, not TaskBase*: undoes the earlier TaskBase* widening that
    // let a CrossRankTask sit directly as a fusion's own child -- that's no longer how CRTs attach.
    Task* parent{nullptr};
    Task* left_child{nullptr};
    Task* right_child{nullptr};

    // This task's own attached seam(s)/CRT(s), if any -- populated via add_special_task, processed in
    // insertion order by decode logic once this Task is solved, before continuing the ordinary climb
    // via `parent`.
    std::vector<SpecialTask*> special_tasks;

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

    // Preemptive OBS extraction: true on an is_extraction_unit_connector chain-link whose newly-closed
    // unit falls inside an active seam's span -- tells the decode-loop phase to skip the immediate
    // divide+post (and skip the usual vb_left restriction) until the seam itself resolves. False
    // (default) for every ROUND-partitioning fusion and every OBS chain-link outside a seam's span,
    // unchanged from today's immediate-divide behavior. Whether this is still needed in its current
    // shape under the new special_tasks attachment model is an open question for the task-building
    // pass (see now-its-time-to-jazzy-turing.md Phase 2 (REVISED)) -- left as-is here.
    bool defer_division{false};

    // Which observable this task belongs to (OBS partitioning only), set once in
    // build_tasks_for_obs_patch_partitioning. -1 for ROUND partitioning.
    int obs_patch_id{-1};

    // Decentralized shot sync (docs/decentralized_shot_sync_design.md §6): a plain, ever-
    // incrementing counter of how many rounds' worth of extraction this partition leaf has
    // completed. No reset is needed between shots (or even between rounds) -- try_to_steal's CAS
    // already guarantees exactly one thread ever owns a given round for a given partition, so
    // there is no reader/writer race on a stale value to protect against. Only Task::reset()
    // (a full --num_repeats restart) puts it back to 0. Meaningful only for partition leaves
    // (!is_fusion), but kept unconditional on Task for simplicity, mirroring how `status` itself
    // is unconditional despite its two branches meaning different things.
    alignas(64) std::atomic<int64_t> extraction_done{0};

#ifdef USE_SHMEM
    // Set to shot_buffer_round whenever task is part of sent window
    // Reset by finalize_sent_crt_window
    // divide_vb/process_extraction_job skips tasks marked sent or 
    // already extracted
    bool sent_in_crt_window{false};
#endif

    Task(int partition, pm::Mwpm* solver)
        : TaskBase(partition, partition - 1, partition, false, TaskType::Task, solver)
    {
        status.store(-1, std::memory_order_release);
    }
    // Partition leaf with explicit local vb bounds (needed for OBS partitioning)
    Task(int part, int vb_l, int vb_r, pm::Mwpm* solver)
        : TaskBase(part, vb_l, vb_r, false, TaskType::Task, solver)
    {
        status.store(-1, std::memory_order_release);
    }
    // Fusion of two ordinary Task operands -- always Task*, never a CrossRankTask (that no longer sits
    // in the tree as a fusion's own child; see special_tasks above).
    Task(int vb, Task* left, Task* right, pm::Mwpm* solver) :
        TaskBase(vb, left->vb_left, right->vb_right, true, TaskType::Task, solver),
        left_child(left),
        right_child(right)
    {
        left->parent = this;
        right->parent = this;
    }

    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept
        : TaskBase(std::move(other))
    {
        status.store(other.status.load());
        extraction_done.store(other.extraction_done.load());
        parent = other.parent;
        left_child = other.left_child;
        right_child = other.right_child;
        special_tasks = std::move(other.special_tasks);
        is_extraction_unit_connector = other.is_extraction_unit_connector;
        is_extraction_unit_root = other.is_extraction_unit_root;
        defer_division = other.defer_division;
        obs_patch_id = other.obs_patch_id;
        sent_in_crt_window = other.sent_in_crt_window;
    }
    Task& operator=(Task&& other) noexcept {
        TaskBase::operator=(std::move(other));
        status.store(other.status.load());
        extraction_done.store(other.extraction_done.load());
        parent = other.parent;
        left_child = other.left_child;
        right_child = other.right_child;
        special_tasks = std::move(other.special_tasks);
        is_extraction_unit_connector = other.is_extraction_unit_connector;
        is_extraction_unit_root = other.is_extraction_unit_root;
        defer_division = other.defer_division;
        obs_patch_id = other.obs_patch_id;
        sent_in_crt_window = other.sent_in_crt_window;
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

    // Attaches a LocalSeamTask/CrossRankTask to this task: records it on this task's own
    // special_tasks (so decode logic finds it when this task is solved) and, symmetrically, records
    // this task on the special task's own child-holding field (LocalSeamTask::children gets this
    // appended; CrossRankTask::child gets set to this). Defined out-of-line, below, once
    // LocalSeamTask/CrossRankTask are complete types it needs to reach into.
    void add_special_task(SpecialTask* st);

    /* Sychronization Methods */
    void mark_solved(size_t my_pid = 0) override {
        (void)my_pid;
        if (is_fusion) {
            status.store(0, std::memory_order_release);
        }
    }

    // for partition leaf, val is the shot id
    // for fusion parent, val is unused -- kept only for interface uniformity with the leaf-claim/
    // LocalSeamTask/CrossRankTask overloads of try_to_steal. The fetch_add-based race needs no
    // per-caller value the way the old fetch_or/child_bit scheme did: every child does the identical
    // fetch_add(1), and whoever's fetch_add returns N-1 (last of N arrivals) is the winner.
    bool try_to_steal(size_t val) override {
        if (!is_fusion) { // partition
            int64_t expected = static_cast<int64_t>(val) - 1;
            return status.compare_exchange_strong(expected, static_cast<int64_t>(val), std::memory_order_acq_rel);
        } else { // fusion
            if (left_child == right_child) { // one child
                return true;
            }
            (void)val;
            constexpr int64_t N = 2; // ordinary binary fusion; generalizes to N>2 with no other change
            int64_t old = status.fetch_add(1, std::memory_order_acq_rel);
            return old == N - 1;
        }
    }

    void reset() override {
        status.store((is_fusion) ? 0 : -1, std::memory_order_release);
        extraction_done.store(0, std::memory_order_release);
    }

    // Called by whichever thread just finished extracting this partition leaf's own previously-
    // closed round. Plain increment, no round argument: the gating this enables (decode(R) can't
    // start until extract(R-1) is marked done; extract(R) can't happen until decode(R) finishes)
    // already forces every partition's own sequence of these calls to happen in strict round
    // order, so a plain counter can never skip or reorder.
    void mark_extraction_done() {
        extraction_done.fetch_add(1, std::memory_order_release);
    }

    // Called immediately after winning try_to_steal(round) for this leaf, before proceeding to
    // decode -- the CAS succeeding only means this thread now owns the partition going forward,
    // not that its buffer is safe to reuse yet. Pure spin-wait, deliberately not
    // std::atomic::wait()/.notify_all(): measured directly on this exact per-partition site (not
    // just the old global current_buffer_round rendezvous in
    // docs/decentralized_shot_sync_design.md §1) -- even with at most one waiter per atomic,
    // futex wait/wake syscall volume at round-granularity call frequency (10,000 rounds/shot) made
    // decode-only wall time 1.1x-4.2x worse (16->256 threads), dominated by `syscall` self-time
    // (20,050 CPU-s at 256 threads, vs 8,767 CPU-s for the spin it replaced).
    //
    // Exponential _mm_pause() backoff, not a single pause per check: PAUSE's latency is wildly
    // architecture-dependent (measured ~18.5s of self-time on Intel Sapphire Rapids at 32 threads,
    // doesn't even register on AMD Zen4 -- Zen4's PAUSE is only a few cycles vs 100+ on many
    // recent Intel parts), so a fixed one-pause-per-check delay is Intel-tuned by accident. Scaling
    // the *count* of pauses per check instead makes the achieved delay self-calibrating to
    // whatever a PAUSE actually costs on the running core, on either architecture, with no syscall
    // involved (so it can't reproduce the wait/notify regression above). Capped, not unbounded:
    // an uncapped backoff would trade check latency for spin cost on a long wait; 1024 keeps the
    // worst-case delay between checks small relative to a round's own decode work.
    //
    // Exact match against `round` itself, no -1 offset: extraction_done starts at 0, so claiming
    // round 0 checks 0 == 0 and never waits.
    void wait_until_previous_extraction_done(int64_t round) {
        int spin_count = 1;
        constexpr int max_spin_count = 1024;
        while (extraction_done.load(std::memory_order_acquire) != round) {
            for (int i = 0; i < spin_count; ++i) {
                _mm_pause();
            }
            if (spin_count < max_spin_count) {
                spin_count *= 2;
            }
        }
    }
};

// N-ary local convergence -- generalizes the old binary local-seam Task to any number of converging
// observables sharing one boundary, so a 3+-way convergence is one LocalSeamTask, not a cascade of
// nested binary seams (which was the actual source of the child_bit-reliability problems worked
// through earlier -- see now-its-time-to-jazzy-turing.md Phase 2 (REVISED) Context). children starts
// empty and is populated by each converging Task's own Task::add_special_task(this) call (task-
// building pass concern, not this file's) rather than being passed in at construction.
struct LocalSeamTask : public SpecialTask {
   private:
    // Arrival counter, same fetch_add race as Task::status's own fusion branch, generalized from a
    // hardcoded N=2 to children.size(): every converging observable's own thread does the identical
    // fetch_add(1); whoever's fetch_add returns (int64_t)children.size()-1 (last arrival) wins.
    alignas(64) std::atomic<int64_t> status{0};

   public:
    // One triggering predecessor per converging observable: the plain Task itself if this is the
    // first SpecialTask attached to it, or the previously-attached SpecialTask on that same Task if
    // not (see Task::add_special_task) -- so setup() below always pulls from whichever list was most
    // recently pruned for that side, not a stale copy from before an earlier attached special task's
    // own divide_vb call. Populated via Task::add_special_task, one push_back per converging side.
    std::vector<TaskBase*> children;

#ifdef ENABLE_DRAW_FLAGS
    // Generalizes the old Task::left_obs_patch_id/right_obs_patch_id -- one entry per children[i].
    std::vector<int> obs_patch_ids;
#endif

    // Busy-spin signal for the children.size()-1 losing threads to wait on, round-tagged rather than a
    // plain bool: a bool left at `true` after the winner finishes would go stale into the *next*
    // shot (nothing else resets it between ordinary shots -- only mark_solved(), once, by the winner),
    // letting a loser in a later round sail through without ever actually waiting. Storing
    // shot_buffer_round instead means a stale value from any earlier round can never accidentally
    // match the current one, mirroring how Task::try_to_steal's own leaf-claim CAS (`expected = val -
    // 1`) avoids the identical class of bug. -1 = never resolved yet. Set by whichever thread wins
    // this seam's own try_to_steal race, as its last action (after that seam's own mark_solved()).
    // Named decode_done (not the old `ready`) to distinguish this purely intra-shot, inter-thread
    // signal from the unrelated per-partition extraction_done added to Task for shot-advancement
    // gating (docs/decentralized_shot_sync_design.md §7) -- different purpose, different state.
    std::atomic<int64_t> decode_done{-1};

    LocalSeamTask(int vb, int vb_left, int vb_right, pm::Mwpm* solver)
        : SpecialTask(vb, vb_left, vb_right, true, TaskType::LocalSeamTask, solver)
    {}

    LocalSeamTask(const LocalSeamTask&) = delete;
    LocalSeamTask& operator=(const LocalSeamTask&) = delete;
    LocalSeamTask(LocalSeamTask&& other) noexcept
        : SpecialTask(std::move(other))
    {
        status.store(other.status.load());
        children = std::move(other.children);
#ifdef ENABLE_DRAW_FLAGS
        obs_patch_ids = std::move(other.obs_patch_ids);
#endif
        decode_done.store(other.decode_done.load());
    }
    LocalSeamTask& operator=(LocalSeamTask&& other) noexcept {
        SpecialTask::operator=(std::move(other));
        status.store(other.status.load());
        children = std::move(other.children);
#ifdef ENABLE_DRAW_FLAGS
        obs_patch_ids = std::move(other.obs_patch_ids);
#endif
        decode_done.store(other.decode_done.load());
        return *this;
    }

    // Combines regions from every child -- same shape as Task::setup(), just looping children instead
    // of two fixed fields (left_child/right_child). Dedups on push: two children can resolve to the
    // same underlying source (e.g. two seams connecting the same observable pair collapse both
    // sides' own "last attached special task" to one shared predecessor), so the same region can
    // otherwise be discovered more than once. children itself must stay un-deduped -- its size is
    // the arrival count N in try_to_steal().
    void setup() override {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        for (TaskBase* c : children) {
            for (auto& region : c->regions_matched_to_virtual_boundary) {
                if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker) {
                    if (std::find(regions_to_unmatch.begin(), regions_to_unmatch.end(), region) == regions_to_unmatch.end())
                        regions_to_unmatch.push_back(region);
                } else {
                    if (std::find(regions_matched_to_virtual_boundary.begin(), regions_matched_to_virtual_boundary.end(), region) == regions_matched_to_virtual_boundary.end())
                        regions_matched_to_virtual_boundary.push_back(region);
                }
            }
        }
    }

    bool try_to_steal(size_t val) override {
        (void)val;
        int64_t N = (int64_t)children.size();
        int64_t old = status.fetch_add(1, std::memory_order_acq_rel);
        return old == N - 1;
    }

    // Should be called after seam vb is fused
    void mark_solved(size_t my_pid = 0) override {
        (void)my_pid;
        status.store(0, std::memory_order_release);
    }

    // Should only be called after seam vb is divided
    void mark_decode_done(size_t shot_buffer_round) {
        decode_done.store(shot_buffer_round, std::memory_order_release);
    }

    void wait_until_decode_done(size_t shot_buffer_round) {
        while (decode_done.load(std::memory_order_acquire) != shot_buffer_round) {}
    }

    void reset() override {
        status.store(0, std::memory_order_release);
        decode_done.store(-1, std::memory_order_release);
    }
};

#ifdef USE_SHMEM
struct CrossRankTask : public SpecialTask {
public:
    // Cross-rank tasks only ever have one local child -- no N-ary generalization needed here the way
    // LocalSeamTask needed one. TaskBase*, not Task*, for the same reason as LocalSeamTask::children
    // above: the predecessor is the plain Task if this is the first SpecialTask attached to it, or the
    // previously-attached SpecialTask if not. Populated via Task::add_special_task(this), not the
    // constructor.
    TaskBase* child{nullptr};
    bool iamleft;

    size_t other_pid{ 0 };

    // p/vb Tasks inside local window for this CRT
    std::vector<Task*> tasks_in_local_window;

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
        bool iamleft,
        int vb_left,
        int vb_right,
        size_t other_pid,
        uint64_t* status_ptr,
        uint64_t* signal_ptr,
        uint64_t* done_ptr,
        pm::FusionSummary* fusion_summary_ptr,
        pm::Mwpm* solver
    ) : SpecialTask(vb, vb_left, vb_right, true, TaskType::CrossRankTask, solver),
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
        *status_shm = 0;
        *signal_shm = 0;
        *done_shm = 0;
        // No more "walk to top of child's own parent chain and attach there" -- CrossRankTasks no
        // longer sit in the tree at all, and no longer take `child` as a constructor argument: the
        // task-building pass sets it via the triggering Task's own add_special_task(this) call
        // instead (mirroring however LocalSeamTask attaches).
    }

    CrossRankTask(CrossRankTask&& other) noexcept : SpecialTask(std::move(other)) {
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

        // Nullify other's ownership so destructor skips it
        other.owns_context = false;
    }

    ~CrossRankTask() {
        if (owns_context && context_shm != SHMEM_CTX_DEFAULT)
            shmem_ctx_destroy(context_shm);
    }


    /* Helper Methods */
    // Symmetric with LocalSeamTask::setup() -- combines regions from its one predecessor (child),
    // splitting into regions_to_unmatch (matched to this CRT's own vb, to be sent/unmatched) vs
    // regions_matched_to_virtual_boundary (everything else, kept live for a later setup() -- the next
    // attached special task on the same Task, or (via the write-back in decode_shots()) that Task's
    // own future parent). The else branch used to be commented out (this list was never read, since
    // nothing propagated it back to child); restored now that the predecessor chain + write-back make
    // it a genuinely consumed list again.
    void setup() override {
        regions_to_unmatch.clear();
        regions_matched_to_virtual_boundary.clear();
        for (auto& region : child->regions_matched_to_virtual_boundary) {
            if (region->match.edge.loc_to && region->match.edge.loc_to->vb == vb_marker) // seam vbs are marked with global part
                regions_to_unmatch.push_back(region);
            else
                regions_matched_to_virtual_boundary.push_back(region);
        }
    };


    /* Sychronization Methods */

    void mark_solved(size_t my_pid) override {
        throw std::logic_error("CrossRankTask::mark_solved: signal_shm is never reset; this should never be called");
    }

    // Renamed from report_done: under Phase 5, this fires only once this PE's own regions/nodes
    // for this round are actually cleaned up (either the winner's post-extraction, or the loser's
    // deferred finalize_sent_crt_window), so "extraction done" is the accurate name now.
    inline void report_extraction_done(int shot_buffer_round) {
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
        // Named (not just FUNC_BEGIN's auto-generated pm::CrossRankTask::wait_until_done) so
        // analyze_trace_wall_clock.py's WAIT_REGIONS can match it by its literal display name --
        // this is the task-status rendezvous wait (blocks until the other PE has posted its
        // done_shm signal for this shot), one of the pure per-PE ENTER-to-LEAVE cross-PE-cost
        // proxies that stay meaningful even across nodes with unsynced clocks.
        SCOREP_USER_REGION_DEFINE(wait_until_done_named);
        SCOREP_USER_REGION_BEGIN(wait_until_done_named, "wait_until_done", SCOREP_USER_REGION_TYPE_COMMON);
#endif
        shmem_wait_until(done_shm, SHMEM_CMP_GE, shot_buffer_round + 1);
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_END(wait_until_done_named);
        SCOREP_USER_FUNC_END();
#endif
    }

    // Kept only to satisfy TaskBase's pure-virtual interface (mirrors mark_solved()'s own precedent
    // just above) -- every real call site goes through a concrete CrossRankTask*, and CRT's own race
    // needs shot_buffer_round too, so decode_shots() always calls the 2-arg overload below instead.
    bool try_to_steal(size_t val) override {
        throw std::logic_error("CrossRankTask::try_to_steal(size_t): call try_to_steal(my_pid, shot_buffer_round) instead");
    }

    // Left's own copy of status_shm is a monotonically increasing counter, incremented by 1 by BOTH
    // sides every shot, NEVER reset (see mark_race_resolved()'s own comment for why the old
    // reset-every-shot OR-based scheme was unsound). For shot k (0-indexed), the counter reaches
    // 2k+1 on whichever side's own fetch_add is the SECOND to land for that shot -- that side wins.
    //
    // Induction argument for why this can never be "lapped" the way the old scheme was: neither side
    // can even ATTEMPT shot k's vote until wait_until_ready_to_race(k) passes, which for left requires
    // status_shm >= 2k, and for right requires its own local marker == 0 (reset only by shot (k-1)'s
    // winner, which is only known -- synchronously, via this fetch_add's own return value -- at the
    // exact instant status_shm reaches 2k). So neither side's vote for shot k can land before shot
    // k-1's pair of votes has fully landed, by construction -- unlike the old bit-reset scheme, there
    // is no window where a value looks "clean" without actually meaning "ready for the next shot".
    bool try_to_steal(size_t my_pid, int shot_buffer_round) {
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_BEGIN();
#endif
        uint64_t old;
        if (iamleft) {
            shmem_wait_until(status_shm, SHMEM_CMP_GE, (uint64_t)(2 * shot_buffer_round));
            old = shmem_ctx_uint64_atomic_fetch_add(context_shm, status_shm, 1, my_pid);
        } else {
            shmem_wait_until(status_shm, SHMEM_CMP_EQ, shot_buffer_round);
            // Local marker on my own copy first ("I attempted this round") -- consumed only by
            // wait_until_ready_to_race()'s own-copy spin gate, decoupled from the real
            // winner-determination fetch_add below (unchanged target: left's copy). Self first, same
            // ordering discipline as mark_race_resolved.
            // shmem_ctx_uint64_atomic_set(context_shm, status_shm, 1, my_pid);
            old = shmem_ctx_uint64_atomic_fetch_add(context_shm, status_shm, 1, other_pid);
        }
        bool won = (old == 2 * (uint64_t)shot_buffer_round + 1);
        if (won) {
            shmem_ctx_uint64_atomic_inc(context_shm, status_shm, (iamleft) ? other_pid : my_pid);
        }
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_FUNC_END();
#endif
        return won;
    }

    // Called by sender to mark local window tasks as sent
    // divide_vb/process_extraction_job will skip these tasks
    // for the provided shot_buffer_round
    void mark_sent_local_window() {
        for (Task *t : tasks_in_local_window)
            t->sent_in_crt_window = true;
    }

    void unmark_sent_local_window() {
        for (Task *t : tasks_in_local_window)
            t->sent_in_crt_window = false;
    }

    void reset() override {
        *status_shm = 0;
        *signal_shm = 0;
        *done_shm = 0;
    }

};
#endif

inline void Task::add_special_task(SpecialTask* st) {
    // The predecessor is this Task itself only if st is the first special task attached to it; if
    // this Task already has one or more special tasks attached, the previously-attached one is the
    // correct predecessor instead -- its own regions_matched_to_virtual_boundary is what's actually
    // up to date (this Task's own list, populated by its setup() before any special task ran, is
    // stale the moment a first special task's own divide_vb call prunes something from it).
    TaskBase* predecessor = special_tasks.empty()
        ? static_cast<TaskBase*>(this)
        : static_cast<TaskBase*>(special_tasks.back());
    special_tasks.push_back(st);
    if (st->is_local_seam_fusion()) {
        static_cast<LocalSeamTask*>(st)->children.push_back(predecessor);
    }
#ifdef USE_SHMEM
    else if (st->is_cross_rank_fusion()) {
        static_cast<CrossRankTask*>(st)->child = predecessor;
    }
#endif
}

#endif  // PYMATCHING2_DECODING_TASK_H
