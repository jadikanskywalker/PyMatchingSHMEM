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

#ifndef PYMATCHING2_DRIVER_HELPERS_FAST_B8_READER_H
#define PYMATCHING2_DRIVER_HELPERS_FAST_B8_READER_H

#include "stim/io/measure_record_reader.h"

namespace pm {

/// Reads the next record using a buffered fast path when a B8 formatter is available.
/// Falls back to the underlying Stim implementation for other formats.
bool start_and_read_entire_record_buffered(
    stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH> &reader,
    stim::SparseShot &shot);

}  // namespace pm

#endif  // PYMATCHING2_DRIVER_HELPERS_FAST_B8_READER_H
