// Copyright 2025 PyMatching Contributors
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

// This code was adapted from GPT-5-Codex output

#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template <size_t W>
bool read_b8_record_fast(
    stim::MeasureRecordReaderFormatB8<W> &reader,
    stim::SparseShot &shot) {
    const size_t bits = reader.bits_per_record();
    if (bits == 0) {
        return false;
    }
    const size_t bytes_per_record = (bits + 7) >> 3;
    thread_local std::vector<uint8_t> scratch;
    if (scratch.size() < bytes_per_record) {
        scratch.resize(bytes_per_record);
    }

    const size_t bytes_read = fread(scratch.data(), 1, bytes_per_record, reader.in);
    if (bytes_read == 0) {
        return false;
    }
    if (bytes_read != bytes_per_record) {
        throw std::invalid_argument(
            "b8 data ended in middle of record at byte position " +
            std::to_string(bytes_read) +
            ".\nExpected bytes per record was " +
            std::to_string(bytes_per_record) +
            " (" + std::to_string(bits) + " bits padded).");
    }

    shot.hits.clear();
    if (shot.obs_mask.num_bits_padded() < reader.num_observables) {
        shot.obs_mask = stim::simd_bits<64>(reader.num_observables);
    }
    shot.obs_mask.clear();

    const size_t num_detectors_total = reader.num_measurements + reader.num_detectors;
    const size_t total_bits = num_detectors_total + reader.num_observables;

    for (size_t byte_index = 0; byte_index < bytes_per_record; ++byte_index) {
        const uint8_t value = scratch[byte_index];
        if (value == 0) {
            continue;
        }
        const size_t base = byte_index << 3;
        for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset) {
            const size_t bit_index = base + bit_offset;
            if (bit_index >= bits) {
                break;
            }
            if ((value & (1U << bit_offset)) == 0) {
                continue;
            }
            if (bit_index < num_detectors_total) {
                shot.hits.push_back(static_cast<uint64_t>(bit_index));
            } else if (bit_index < total_bits) {
                const size_t observable_index = bit_index - num_detectors_total;
                shot.obs_mask[observable_index] ^= true;
            } else {
                throw std::invalid_argument("Hit index from data is too large.");
            }
        }
    }

    return true;
}

}  // namespace

namespace pm {

bool start_and_read_entire_record_buffered(
    stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH> &reader,
    stim::SparseShot &shot) {
    if (auto *b8_reader = dynamic_cast<stim::MeasureRecordReaderFormatB8<stim::MAX_BITWORD_WIDTH> *>(&reader)) {
        return read_b8_record_fast(*b8_reader, shot);
    }
    return reader.start_and_read_entire_record(shot);
}

}  // namespace pm
