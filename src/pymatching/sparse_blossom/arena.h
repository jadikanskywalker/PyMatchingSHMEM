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

#ifndef PYMATCHING_FILL_MATCH_ARENA_H
#define PYMATCHING_FILL_MATCH_ARENA_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "pymatching/sparse_blossom/config_parallel.h"

namespace pm {

/// World's simplest bulk memory owner.
///
/// Memory allocated by the arena is free'd when the arena is destructed.
template <typename T>
struct Arena {
    std::vector<T *> allocated;
    std::vector<T *> available;

    // available_reserve_size: pre-reserves available's capacity so a burst of concurrent del() calls
    // (each doing available.push_back()) can't trigger a concurrent internal buffer growth on it --
    // the actual race this whole Arena has been chasing. 0 (the default) reserves nothing, unchanged
    // from before this parameter existed.
    explicit Arena(size_t available_reserve_size = 0) : allocated(), available() {
        if (available_reserve_size) {
            available.reserve(available_reserve_size);
        }
    }
    Arena(const Arena &) = delete;
    Arena(Arena &&other) : allocated(std::move(other.allocated)), available(std::move(other.available)) {
    }

    T *alloc_unconstructed() {
        if (available.empty()) {
            // calloc, not malloc: constructed/destructed (see graph_fill_region.h) are never written
            // by Arena itself, and their constructor only ever increments -- a brand-new slot that's
            // never been constructed on before must start life zeroed, not garbage, for that first
            // increment (and the constructed==destructed+1 check below) to mean anything.
            T *p = (T *)calloc(1, sizeof(T));
            allocated.push_back(p);
            return p;
        }
        T *result = available.back();
        available.pop_back();
        return result;
    }

    T *alloc_default_constructed() {
        T *result = alloc_unconstructed();
        if constexpr (requires (T t) { t.reset(); t.constructed; }) {
            // del() calls reset() instead of ~T() for this type (see del()'s own comment) -- a
            // region's C++ lifetime is never actually ended by Arena, so a recycled slot is already
            // sitting in a valid, reset (default-equivalent) state from its last del(). Placement-
            // newing a fresh object on top of it would be constructing over an object whose lifetime
            // never ended. constructed==0 is the ground truth for "never actually constructed yet"
            // (calloc leaves a brand-new slot's counters at exactly 0; reset() never touches
            // constructed) -- only that slot needs the real constructor to run, exactly once, ever.
            if (result->constructed == 0) {
                new (result) T();
            } else {
                result->reset();
                result->constructed++;
            }
        } else {
            new (result) T();
        }
        return result;
    }

    void del(T *p) {
        if (DEBUG) {
            // Types without their own constructed/destructed counters (e.g. AltTreeNode) have no
            // protection at all against a double-del() beyond this O(n) scan -- a diagnostic-only
            // safety net to catch that before it corrupts the heap (a second ~T() on an already-
            // destructed vector-holding member is exactly a "double free" trigger). Types with the
            // counters get an extra O(1) check on top: a genuinely live object always has
            // constructed == destructed + 1 (exactly one more construction than destruction so far);
            // anything else means either a lost available.push_back() (the concurrent-vector-growth
            // race -- p was already deleted, but never made it into `available`, so the scan above
            // doesn't catch it either) or a genuine double-del().
            bool already_deleted = std::find(available.begin(), available.end(), p) != available.end();
            if constexpr (requires (T t) { t.constructed; t.destructed; }) {
                if (p->constructed != p->destructed + 1) {
                    std::cout << "ERROR: Arena::del() called on p=" << p
                               << " with constructed=" << p->constructed
                               << " destructed=" << p->destructed
                               << " (expected constructed == destructed + 1)"
                               << std::endl << std::flush;
                    already_deleted = true;
                }
            }
            if (already_deleted) {
                std::cout << "ERROR: Arena::del() called twice on same pointer p=" << p << std::endl << std::flush;
                return;
            }
        }
        if constexpr (requires (T t) { t.destructed; }) {
            p->destructed++;
            // p->reset();
        } else {
            p->~T();
        }
        available.push_back(p);
    }

    ~Arena() {
        std::vector<T *> to_free = std::move(allocated);
        std::vector<T *> not_in_use = std::move(available);

        // Destruct the objects that were still in use.
        std::sort(to_free.begin(), to_free.end());
        std::sort(not_in_use.begin(), not_in_use.end());
        size_t kf = 0;
        size_t kn = 0;
        while (kf < to_free.size()) {
            // if (kf > 0 && to_free[kf] == to_free[kf - 1]) {
            //     // allocated (not available) contains this address twice -- every malloc() should be
            //     // recorded here exactly once, ever, for the arena's whole lifetime, so a duplicate
            //     // means allocated itself was corrupted (the same concurrent-push_back-growth race we
            //     // suspect on available, just landing on this vector instead). If the first occurrence
            //     // already ran ~T() this pass, this second one would silently read/destruct
            //     // already-freed memory -- the constructed/destructed check just below wouldn't catch
            //     // this on its own either, since it also only ever looks at this one address's own
            //     // counters, which the first pass's destruct already left in a "looks fine" state.
            //     std::cout << "ERROR: Arena::~Arena() found p=" << to_free[kf]
            //                << " listed twice in allocated -- skipping the duplicate"
            //                << std::endl << std::flush;
            //     kf++;
            //     continue;
            // }
            if (kn < not_in_use.size() && not_in_use[kn] == to_free[kf]) {
                if constexpr (requires (T t) { t.reset(); t.destructed; }) {
                    to_free[kf]->~T();
                }
                kf++;
                kn++;
            } else {
                if constexpr (requires (T t) { t.constructed; t.destructed; }) {
                    // Treated as still in use (absent from available) -- a genuinely live object
                    // always has constructed == destructed + 1 (exactly one more construction than
                    // destruction so far). Anything else means either a lost available.push_back()
                    // (the concurrent-vector-growth race) or a prior double-del() that already ran
                    // ~T() on this address without ever recording it in available -- destructing it
                    // again here would be exactly the double-free those failure modes lead to, so
                    // skip it and print both counts instead.
                    if (to_free[kf]->constructed != to_free[kf]->destructed + 1) {
                        std::cout << "ERROR: Arena::~Arena() found object p=" << to_free[kf]
                                   << " treated as still in use (absent from available) but "
                                   << "constructed=" << to_free[kf]->constructed
                                   << " destructed=" << to_free[kf]->destructed
                                   << " (expected constructed == destructed + 1)"
                                   << std::endl << std::flush;
                        if constexpr (requires (T t) { t.shell_area; }) {
                            // del() never resets/destructs a GraphFillRegion's fields (only
                            // alloc_default_constructed() does, right before real reuse) -- a stranded
                            // region (del()'d, i.e. destructed incremented, but never recycled, e.g.
                            // lost from available by the very race we're chasing) keeps every field
                            // exactly as it was at that last del(), so this is real, inspectable
                            // last-known state, not a blanked reset() one.
                            T* r = to_free[kf];
                            std::cout << "  blossom_parent=" << r->blossom_parent << std::endl << std::flush;
                            std::cout << "  blossom_parent_top=" << r->blossom_parent_top << std::endl << std::flush;
                            std::cout << "  radius=" << r->radius << std::endl << std::flush;
                            std::cout << "  match.region=" << r->match.region << std::endl << std::flush;
                            std::cout << "  match.edge=" << r->match.edge << std::endl << std::flush;
                            std::cout << "  alt_tree_node=" << r->alt_tree_node << std::endl << std::flush;
                            std::cout << "  blossom_children.size()=" << r->blossom_children.size() << std::endl << std::flush;
                            std::cout << "  shell_area.size()=" << r->shell_area.size() << std::endl << std::flush;
                        }
                        kf++;
                        continue;
                    }
                }
                if constexpr (requires (T t) { t.shell_area; }) {
                    // Identity dump for whatever's about to be destructed here: it passed the
                    // constructed/destructed check above, so it's genuinely still live (never
                    // independently del()'d) -- if this destruct still crashes (shell_area's own
                    // backing buffer corrupted), this is the only place we ever learn anything about
                    // *which* region it was, since nothing else touches it between creation and this
                    // final teardown. Each field gets its own flushed line: if a field access itself
                    // is what's corrupted (e.g. shell_area.front() reading a dangling pointer), a
                    // crash there must not swallow every field printed before it -- std::cout is
                    // buffered, and an abort mid-expression never reaches a later std::flush, so
                    // nothing not already flushed survives.
                    T* r = to_free[kf];
                    std::cout << "NOTE: Arena::~Arena() about to destruct still-live region p=" << r
                               << std::endl << std::flush;
                    std::cout << "  shell_area.size()=" << r->shell_area.size() << std::endl << std::flush;
                    if (!r->shell_area.empty()) {
                        std::cout << "  shell_area.front()=" << r->shell_area.front() << std::endl << std::flush;
                        std::cout << "  shell_area.back()=" << r->shell_area.back() << std::endl << std::flush;
                    }
                    std::cout << "  match.region=" << r->match.region << std::endl << std::flush;
                    std::cout << "  match.edge.loc_to=" << r->match.edge.loc_to << std::endl << std::flush;
                    std::cout << "  alt_tree_node=" << r->alt_tree_node << std::endl << std::flush;
                    std::cout << "  blossom_parent=" << r->blossom_parent << std::endl << std::flush;
                    std::cout << "  blossom_parent_top=" << r->blossom_parent_top << std::endl << std::flush;
                    std::cout << "  blossom_children.size()=" << r->blossom_children.size() << std::endl << std::flush;
                    std::cout << "  rotating_buffer_idx=" << r->rotating_buffer_idx << std::endl << std::flush;
                }
                // No manual bookkeeping needed here (unlike the old allocated-flag scheme): ~T()'s
                // own body increments destructed itself, regardless of which call site (del() or
                // here) triggered it.
                to_free[kf]->~T();
                kf++;
            }
        }

        // Free all allocated memory.
        for (T *v : to_free) {
            free(v);
        }
    }
};
}  // namespace pm

#endif
