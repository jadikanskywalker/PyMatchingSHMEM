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
#ifndef PYMATCHING2_SHMEM_ARENA_H
#define PYMATCHING2_SHMEM_ARENA_H

#include <cstring>
#include <iostream>
#include <limits>

#include <omp.h>
#include <shmem.h>

#include "pymatching/sparse_blossom/config_parallel.h"

// The SHMEMArena, for typename T, :
//   - references a shmem array with SHMEM_ARENA_BUFFER_NELEMS slots for T
//   - assigns new elements at the leftmost available index
//   - tracks which slots are available/taken
//   - falls back to standard heap memory in the case of overflow
template <typename T>
struct SHMEMArena {
    // SHMEM array buffer and availability bitmap
    T* shmem_buffer;
    std::vector<uint64_t> shmem_bitmap;  // 1 = free; 0 = taken
    size_t shmem_buffer_size; // must be a multiple of 64

    // Heap fallback
    std::vector<T*> allocated;
    std::vector<T*> available;

    int num_overflows_tracker{0};

    SHMEMArena() : shmem_buffer(nullptr), shmem_buffer_size(0), allocated(), available(), shmem_bitmap() {
    }
    SHMEMArena(T* shmem_buffer_, size_t shmem_buffer_size_)
        : shmem_buffer(shmem_buffer_), shmem_buffer_size(shmem_buffer_size_), allocated(), available() {
        shmem_bitmap.resize(shmem_buffer_size/64, std::numeric_limits<uint64_t>::max());
        // Could set 0's to support SHMEM_ARENA_BUFFER_NELEMS not a multiple of 64
        // shmem_malloc (like malloc, unlike calloc) makes no zeroing guarantee -- a type with
        // constructed/destructed counters (see graph_fill_region.h) needs every slot in this whole
        // buffer zeroed before its first-ever placement-new, for the same reason Arena::
        // alloc_unconstructed() uses calloc instead of malloc for its own heap fallback below. This
        // constructor runs exactly once per (PE, shot_container, partition) slice for the buffer's
        // entire lifetime, so this is the one place to do it.
        std::memset(shmem_buffer, 0, shmem_buffer_size * sizeof(T));
    }
    SHMEMArena(const SHMEMArena&) = delete;
    SHMEMArena(SHMEMArena&& other)
        : shmem_buffer(other.shmem_buffer),
          shmem_bitmap(std::move(other.shmem_bitmap)),
          shmem_buffer_size(other.shmem_buffer_size),
          allocated(std::move(other.allocated)),
          available(std::move(other.available)) {
    }

    T* alloc_unconstructed() {
        T* result = nullptr;
        for (size_t word_id = 0; word_id < shmem_bitmap.size(); ++word_id) {
            uint64_t word = shmem_bitmap[word_id];
            if (word != 0) {
                size_t available_bit = __builtin_ctzll(word);       // find first bit == 1
                shmem_bitmap[word_id] &= ~(1ULL << available_bit);  // set bit to 0 (taken)
                result = shmem_buffer + word_id * 64 + available_bit;
                break;
            }
        }
        if (result == nullptr) {  // fall back to heap memory
            // std::cout << "ERROR: Fallback to heap" << std::endl
            //           << "  Thread: " << omp_get_thread_num() << "    PE: " << shmem_my_pe() << std::endl;
            if (available.empty()) {
                // calloc, not malloc: see arena.h's own alloc_unconstructed() for why -- a type with
                // constructed/destructed counters needs a brand-new slot zeroed before its first-ever
                // placement-new.
                T* p = (T*)calloc(1, sizeof(T));
                allocated.push_back(p);
                available.push_back(p);
            }
            result = available.back();
            available.pop_back();
        }
        return result;
    }

    T* alloc_default_constructed() {
        T* result = alloc_unconstructed();
        if constexpr (requires (T t) { t.reset(); t.constructed; }) {
            // Mirrors arena.h's own alloc_default_constructed() exactly -- see its comment. Every
            // slot in shmem_buffer is zeroed once up front (this arena's own constructor), and the
            // heap-fallback path callocs, so constructed==0 reliably means "never actually
            // constructed yet" regardless of which of the two backing stores this slot came from.
            // A slot that's been used before is only ever reset() here, right before real reuse --
            // del() itself leaves every field untouched (see del()'s own comment).
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

    void del(T* p) {
        if (DEBUG) {
            // A genuinely live object always has constructed == destructed + 1 (exactly one more
            // construction than destruction so far) -- see arena.h's own del() for the full rationale
            // (this mirrors it exactly). Not atomic, unlike the __atomic_exchange_n check this
            // replaced: a genuine concurrent double-del() can still slip both threads past a plain
            // read here, same as arena.h's own version -- this is a diagnostic, not a fix for the
            // underlying race.
            if constexpr (requires (T t) { t.constructed; t.destructed; }) {
                if (p->constructed != p->destructed + 1) {
                    std::cout << "ERROR: T" << omp_get_thread_num() << " called SHMEMArena::del() on p=" << p
                               << " with constructed=" << p->constructed
                               << " destructed=" << p->destructed
                               << " (expected constructed == destructed + 1)"
                               << std::endl << std::flush;
                    return;
                }
            }
        }
        if constexpr (requires (T t) { t.destructed; }) {
            // Deliberately do NOT reset() or destruct here -- see arena.h's own del() for why: leaves
            // every field as its real last-known state for debug printing, in case this slot gets
            // stranded (del()'d but never recycled). Only alloc_default_constructed() ever calls
            // reset(), right before the slot is about to be reused for real.
            p->destructed++;
        } else {
            p->~T();
        }
        if (p >= shmem_buffer && p < shmem_buffer + shmem_buffer_size) {
            size_t idx = p - shmem_buffer;
            shmem_bitmap[idx/64] |= (1ULL << (idx%64));  // set bit to 1 (free)
        } else {
            available.push_back(p);
        }
    }
};

#endif