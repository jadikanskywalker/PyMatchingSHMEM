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
                T* p = (T*)malloc(sizeof(T));
                allocated.push_back(p);
                available.push_back(p);
            }
            result = available.back();
            available.pop_back();
        }
        if constexpr (requires (T t) { t.allocated; }) {
            result->allocated = true;
        }
        return result;
    }

    T* alloc_default_constructed() {
        T* result = alloc_unconstructed();
        new (result) T();
        if constexpr (requires (T t) { t.allocated; }) {
            result->allocated = true;
        }
        return result;
    }

    void del(T* p) {
        // Mark stale before destroying -- see arena.h's own del() for why (write must land while the
        // object is still formally alive, and before anything else could observe it as reused).
        if constexpr (requires (T t) { t.allocated; }) {
            p->allocated = false;
        }
        p->~T();
        if (p >= shmem_buffer && p < shmem_buffer + shmem_buffer_size) {
            size_t idx = p - shmem_buffer;
            shmem_bitmap[idx/64] |= (1ULL << (idx%64));  // set bit to 1 (free)
        } else {
            available.push_back(p);
        }
    }
};

#endif