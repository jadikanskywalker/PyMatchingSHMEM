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

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

#include <omp.h>
#include <shmem.h>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

namespace pm {

// The GraphFillRegionArena:
//   - references a shmem array of GraphFillRegionExtractedFields (extracted_field) and a same-length
//     process-local array of GraphFillRegions (local_regions); slot i of each belongs to the same region
//   - assigns new slots at the leftmost available index, tracked by one shared bitmap
//   - falls back to heap memory in the case of overflow
struct GraphFillRegionArena {
    GraphFillRegionExtractedFields* extracted_field;
    // Raw (unconstructed, uninitialized) memory: slot i is only valid once i < constructed_count.
    GraphFillRegion* local_regions;
    std::vector<uint64_t> shmem_bitmap;  // 1 = free; 0 = taken
    size_t shmem_buffer_size;            // must be a multiple of 64
    // One past the highest slot index ever constructed. Allocation always takes the lowest free bit, so
    // every slot below this has been constructed at least once and every slot at or above it never has.
    size_t constructed_count{0};

    // Heap fallback. heap_extracted_fields[i] belongs to allocated[i].
    std::vector<GraphFillRegion*> allocated;
    std::vector<GraphFillRegion*> available;
    std::vector<GraphFillRegionExtractedFields*> heap_extracted_fields;

    int num_overflows_tracker{0};

    // True only for a remote-partition shadow arena (decoding_unit.cc's remote_arenas). The sender's put
    // overwrites the whole extracted_field slab, taken and free slots alike, and its bitmap replaces this
    // one, so received regions never pass through alloc_default_constructed(): their constructed/destructed
    // counters are not maintained (and not checked), and a free slot's extracted fields can hold the
    // sender's stale bytes, so they are re-initialized on every alloc here.
    bool is_remote_shadow{false};

    GraphFillRegionArena()
        : extracted_field(nullptr), local_regions(nullptr), shmem_bitmap(), shmem_buffer_size(0) {
    }
    GraphFillRegionArena(
        GraphFillRegionExtractedFields* extracted_field_,
        GraphFillRegion* local_regions_,
        size_t shmem_buffer_size_,
        bool is_remote_shadow_ = false)
        : extracted_field(extracted_field_),
          local_regions(local_regions_),
          shmem_buffer_size(shmem_buffer_size_),
          is_remote_shadow(is_remote_shadow_) {
        shmem_bitmap.resize(shmem_buffer_size / 64, std::numeric_limits<uint64_t>::max());
    }
    GraphFillRegionArena(const GraphFillRegionArena&) = delete;
    GraphFillRegionArena(GraphFillRegionArena&& other)
        : extracted_field(other.extracted_field),
          local_regions(other.local_regions),
          shmem_bitmap(std::move(other.shmem_bitmap)),
          shmem_buffer_size(other.shmem_buffer_size),
          constructed_count(other.constructed_count),
          allocated(std::move(other.allocated)),
          available(std::move(other.available)),
          heap_extracted_fields(std::move(other.heap_extracted_fields)),
          is_remote_shadow(other.is_remote_shadow) {
        // Slot regions point at their owning arena, so they must follow it to its new address.
        for (size_t i = 0; i < constructed_count; ++i) {
            local_regions[i].owner_arena = this;
        }
        for (GraphFillRegion* r : allocated) {
            r->owner_arena = this;
        }
    }

    // Constructs every slot up to and including idx that has never been constructed. keep_extracted
    // leaves the extracted fields bytes alone (they were put there by a remote sender).
    void ensure_constructed(size_t idx, bool keep_extracted) {
        while (constructed_count <= idx) {
            GraphFillRegion* r = new (local_regions + constructed_count) GraphFillRegion();
            if (keep_extracted) {
                r->extracted_fields = extracted_field + constructed_count;
            } else {
                new (extracted_field + constructed_count) GraphFillRegionExtractedFields(r);
            }
            r->owner_arena = this;
            ++constructed_count;
        }
    }

    GraphFillRegion* alloc_default_constructed() {
        GraphFillRegion* result = nullptr;
        for (size_t word_id = 0; word_id < shmem_bitmap.size(); ++word_id) {
            uint64_t word = shmem_bitmap[word_id];
            if (word != 0) {
                size_t available_bit = __builtin_ctzll(word);       // find first bit == 1
                shmem_bitmap[word_id] &= ~(1ULL << available_bit);  // set bit to 0 (taken)
                size_t idx = word_id * 64 + available_bit;
                if (idx >= constructed_count) {
                    ensure_constructed(idx, false);
                } else {
                    local_regions[idx].constructed++;  // already reset by its last del()
                    if (is_remote_shadow) {
                        local_regions[idx].state().reset(local_regions + idx);
                    }
                }
                return local_regions + idx;
            }
        }
        // fall back to heap memory
        if (available.empty()) {
            result = (GraphFillRegion*)calloc(1, sizeof(GraphFillRegion));
            allocated.push_back(result);
            new (result) GraphFillRegion();
            heap_extracted_fields.push_back(new GraphFillRegionExtractedFields(result));
            result->owner_arena = this;
        } else {
            result = available.back();
            available.pop_back();
            result->constructed++;  // already reset by its last del()
        }
        return result;
    }

    void del(GraphFillRegion* p) {
        if (DEBUG && !is_remote_shadow) {
            // A genuinely live region always has constructed == destructed + 1.
            if (p->constructed != p->destructed + 1) {
                std::cout << "ERROR: T" << omp_get_thread_num() << " called GraphFillRegionArena::del() on p=" << p
                          << " with constructed=" << p->constructed
                          << " destructed=" << p->destructed
                          << " (expected constructed == destructed + 1)"
                          << std::endl << std::flush;
                return;
            }
        }
        p->destructed++;
        p->reset();
        if (p >= local_regions && p < local_regions + shmem_buffer_size) {
            size_t idx = p - local_regions;
            shmem_bitmap[idx / 64] |= (1ULL << (idx % 64));  // set bit to 1 (free)
        } else {
            available.push_back(p);
        }
    }
};

}  // namespace pm

#endif
