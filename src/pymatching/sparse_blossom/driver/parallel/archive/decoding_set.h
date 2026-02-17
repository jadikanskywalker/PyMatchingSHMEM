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

#ifndef PYMATCHING2_DECODING_SET_H
#define PYMATCHING2_DECODING_SET_H

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"
#include "pymatching/sparse_blossom/driver/parallel/decoding_unit.h"

#include "stim.h"

#include <vector>
#include <mutex>

namespace pm { typedef int64_t total_weight_int; };

struct ActiveShot {
    int shot_id;
    stim::SparseShot shot;
    // std::atomic<int> ref_counter{ 0 };
};

// A decoding set is a collection of decoding units,
// where units are the disjoint graphs from the DEM
struct DecodingSet {
    std::vector<DecodingUnit> units;

private:
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader;
    std::unique_ptr<stim::MeasureRecordWriter> writer;

    const bool enable_correlations;
    const bool draw_frames;

    std::mutex read_mutex;

    int total_partitions;

    std::vector<ActiveShot> shot_buffer;
    size_t shot_start = 0;

public:    
    DecodingSet(
        std::vector<DecodingUnit> _units,
        std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> reader,
        std::unique_ptr<stim::MeasureRecordWriter> writer,
        bool enable_correlations,
        bool draw_frames
    ) : units(std::move(_units)), reader(std::move(reader)),
        writer(std::move(writer)), enable_correlations(enable_correlations),
        draw_frames(draw_frames) {
        if (units.size() == 0) {
            throw std::invalid_argument("(DecodingSet) No DecodingUnits were given.");
        }
        shot_buffer.reserve((size_t) NUM_ACTIVE_SHOTS_PER_SET);
        total_partitions = 0;
        for (auto& unit : units) {
            unit.build_tasks_for_round_partitioning();
            total_partitions += unit.num_partitions;
        }
    }

    void build_solvers(int max_threads, stim::DetectorErrorModel dem);

    // ActiveShot& get_shot(int shot_id) {
    //     // return (reader->start_and_read_entire_record());
    // }

    void decode_shots();
};

#endif // PYMATCHING2_DECODING_SET_H