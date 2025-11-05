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

#include "pymatching/sparse_blossom/driver/namespaced_main.h"

#include <cstdlib>
#include <chrono>
#include <cstring>
#include <iostream>
#include <vector>

#include "pymatching/sparse_blossom/diagram/animation_main.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "stim.h"

#ifdef USE_SHMEM
// ===============
#include <omp.h>
#include <shmem.h>
#endif

#ifdef ENABLE_FUSION
#include "../config_parallel.h"
#include "../diagram/mwpm_diagram.h"
#endif

#ifdef USE_THREADS
#include <omp.h>
// ===============
#endif

int main_predict(int argc, const char **argv) {
    stim::check_for_unknown_arguments(
        {
            "--in",
            "--in_format",
            "--in_includes_appended_observables",
            "--out",
            "--out_format",
            "--dem",
            "--enable_correlations"
#if defined(ENABLE_FUSION)
// ===============
            , "--rounds_per_partition",
            "--draw_frames",
            "--parallel"
// ===============
#endif
        },
        {},
        "predict",
        argc,
        argv);

    FILE *shots_in = stim::find_open_file_argument("--in", stdin, "rb", argc, argv);
    FILE *predictions_out = stim::find_open_file_argument("--out", stdout, "wb", argc, argv);
    FILE *dem_file = stim::find_open_file_argument("--dem", nullptr, "r", argc, argv);
    stim::FileFormatData shots_in_format =
        stim::find_enum_argument("--in_format", "b8", stim::format_name_to_enum_map(), argc, argv);
    stim::FileFormatData predictions_out_format =
        stim::find_enum_argument("--out_format", "01", stim::format_name_to_enum_map(), argc, argv);
    bool append_obs = stim::find_bool_argument("--in_includes_appended_observables", argc, argv);
    bool enable_correlations = stim::find_bool_argument("--enable_correlations", argc, argv);

#if defined(ENABLE_FUSION)
// ===============
    config_parallel::M = stim::find_int64_argument("--rounds_per_partition", 10, 1, INT64_MAX, argc, argv);
    bool draw_frames = stim::find_bool_argument("--draw_frames", argc, argv);
#if defined(USE_THREADS) || defined(USE_SHMEM)
    bool parallel = stim::find_bool_argument("--parallel", argc, argv);
#else
    bool parallel = false;
#endif
// ===============
#endif

    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    size_t num_obs = dem.count_observables();
    auto reader = stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>::make(
        shots_in, shots_in_format.id, 0, dem.count_detectors(), append_obs * num_obs);
    auto writer = stim::MeasureRecordWriter::make(predictions_out, predictions_out_format.id);
    writer->begin_result_type('L');

    pm::weight_int num_buckets = pm::NUM_DISTINCT_WEIGHTS;

    auto mwpm = pm::detector_error_model_to_mwpm(
        dem,
        num_buckets,
        /*ensure_search_flooder_included=*/enable_correlations,
        /*enable_correlations=*/enable_correlations);

#ifdef USE_THREADS
// ===============
    int num_threads = omp_get_max_threads();
    if (num_threads > mwpm.flooder.graph.num_partitions) {
        omp_set_num_threads(mwpm.flooder.graph.num_partitions);
        num_threads = mwpm.flooder.graph.num_partitions;
    }
    pm::tasks.resize(static_cast<size_t>(mwpm.flooder.graph.num_partitions));
    pm::partitions_task_id.resize(static_cast<size_t>(mwpm.flooder.graph.num_partitions));

    // pm::partition_task_queues.resize(static_cast<size_t>(num_threads));
    // pm::build_thread_solvers(
    //     mwpm,
    //     /*ensure_search_flooder_included=*/enable_correlations,
    //     /*enable_correlations=*/enable_correlations,
    //     num_threads
    // );
    
    // auto coords = pm::pick_coords_for_drawing_from_dem(dem, 20);
    // mwpm.coords = coords;
    // for (auto &s : pm::solvers)
    //     s->coords = coords;

    // pm::init_tasks(num_threads, mwpm.flooder.graph.num_partitions);
    // pm::init_task_queues(num_threads, mwpm.flooder.graph.num_partitions);
// ===============
#endif
#ifdef ENABLE_FUSION
// ===============
    mwpm.coords = pm::pick_coords_for_drawing_from_dem(dem, 20);
// ===============
#endif

    stim::SparseShot sparse_shot;
    sparse_shot.clear();
    pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);

#ifdef ENABLE_FUSION
// ===============
    pm::setup_output_dirs(draw_frames, parallel);
    if (DEBUG) {
        output_detector_nodes(mwpm, true);
    }
// ===============
#endif

#if defined(USE_THREADS) && !defined(USE_SHMEM)
// ===============
    int i = 0;
    while (reader->start_and_read_entire_record(sparse_shot)) {
        if (parallel)
            pm::decode_detection_events_in_parallel(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations, i, draw_frames);
        else
            pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations, i, draw_frames);
        for (size_t k = 0; k < num_obs; k++) {
            if (DEBUG)
                std::cout << "res.obs_crossed[" << k << "]: " << (static_cast<unsigned int>(res.obs_crossed[k])) << std::endl;
            writer->write_bit(res.obs_crossed[k]);
        }
        writer->write_end();
        sparse_shot.clear();
        res.reset();
        i++;
    }
#elif defined(USE_SHMEM)
    int mype = shmem_my_pe();
    if (DEBUG && mype==0) {
        std::cout << "DEBUG: PE " << mype << " dem.count_detectors() : " << dem.count_detectors() << std::endl;
        output_detector_nodes(mwpm, parallel);
    }

    // Hits buffers
    long *hits_size = (long *)shmem_malloc(sizeof(long));
    uint64_t *hits = (uint64_t *)shmem_malloc(dem.count_detectors() * sizeof(uint64_t));
    // Start processing shots
    int i = 0;
    while (true) {
        // PE 0 reads the shot
        if (mype == 0) {
            if (reader->start_and_read_entire_record(sparse_shot)) {
                *hits_size = (long) sparse_shot.hits.size();
                // Copy hits into buffer
                for (size_t i = 0; i < *hits_size; ++i)
                    hits[i] = sparse_shot.hits[i];
            } else {
                *hits_size = -1;
            }
        }
        shmem_barrier_all();
        // Broadcast hits size
        shmem_long_broadcast(SHMEM_TEAM_WORLD, hits_size, hits_size, 1, 0);
        if (*hits_size < 0) {
            shmem_barrier_all();
            shmem_free(hits_size);
            shmem_free(hits);
            // no more shots
            break;
        }
        shmem_barrier_all();
        // Broadcast hits data
        shmem_uint64_broadcast(SHMEM_TEAM_WORLD, hits, hits, *hits_size, 0);
        // Populate sparse_shot.hits from buffer
        if (mype != 0) {
            sparse_shot.hits.assign(hits, hits + *hits_size);
        }
        // Decode shot
        if (parallel)
            pm::decode_detection_events_in_parallel(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations, i, draw_frames);
        else
            pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations, i, draw_frames);
        if (mype == 0) {
            for (size_t k = 0; k < num_obs; k++) {
                writer->write_bit(res.obs_crossed[k]);
            }
            writer->write_end();
        }

        sparse_shot.clear();
        res.reset();
        i++;
    }
// ===============
#else
    while (reader->start_and_read_entire_record(sparse_shot)) {
        pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations);
        for (size_t k = 0; k < num_obs; k++) {
            writer->write_bit(res.obs_crossed[k]);
        }
        writer->write_end();
        sparse_shot.clear();
        res.reset();
    }
#endif
    if (predictions_out != stdout) {
        fclose(predictions_out);
    }
    if (shots_in != stdin) {
        fclose(shots_in);
    }

    return EXIT_SUCCESS;
}

int main_count_mistakes(int argc, const char **argv) {
    stim::check_for_unknown_arguments(
        {
            "--in",
            "--in_format",
            "--in_includes_appended_observables",
            "--obs_in",
            "--obs_in_format",
            "--out",
            "--dem",
            "--time",
            "--enable_correlations",
        },
        {},
        "count_mistakes",
        argc,
        argv);

    FILE *shots_in = stim::find_open_file_argument("--in", stdin, "rb", argc, argv);
    FILE *obs_in = stim::find_open_file_argument("--obs_in", stdin, "rb", argc, argv);
    FILE *stats_out = stim::find_open_file_argument("--out", stdout, "wb", argc, argv);
    FILE *dem_file = stim::find_open_file_argument("--dem", nullptr, "r", argc, argv);
    stim::FileFormatData shots_in_format =
        stim::find_enum_argument("--in_format", "01", stim::format_name_to_enum_map(), argc, argv);
    stim::FileFormatData obs_in_format =
        stim::find_enum_argument("--obs_in_format", "01", stim::format_name_to_enum_map(), argc, argv);
    bool append_obs = stim::find_bool_argument("--in_includes_appended_observables", argc, argv);
    bool enable_correlations = stim::find_bool_argument("--enable_correlations", argc, argv);

    bool time = stim::find_bool_argument("--time", argc, argv);
    if (!append_obs && obs_in == nullptr) {
        throw std::invalid_argument("Must specify --in_includes_appended_observables or --obs_in.");
    }

    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    size_t num_obs = dem.count_observables();
    std::unique_ptr<stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>> obs_reader;
    if (obs_in != stdin) {
        obs_reader = stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>::make(obs_in, obs_in_format.id, 0, 0, num_obs);
    }
    auto reader = stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>::make(
        shots_in, shots_in_format.id, 0, dem.count_detectors(), append_obs * num_obs);

    pm::weight_int num_buckets = pm::NUM_DISTINCT_WEIGHTS;
    auto mwpm = pm::detector_error_model_to_mwpm(
        dem,
        num_buckets,
        /*ensure_search_flooder_included=*/enable_correlations,
        /*enable_correlations=*/enable_correlations);

    stim::SparseShot sparse_shot;
    stim::SparseShot obs_shot;
    size_t num_mistakes = 0;
    size_t num_shots = 0;
    auto start = std::chrono::steady_clock::now();
    while (reader->start_and_read_entire_record(sparse_shot)) {
        if (obs_reader == nullptr) {
            obs_shot.obs_mask = sparse_shot.obs_mask;
        } else {
            if (!obs_reader->start_and_read_entire_record(obs_shot)) {
                throw std::invalid_argument("Obs data ended before shot data ended.");
            }
        }
        auto res = pm::decode_detection_events_for_up_to_64_observables(mwpm, sparse_shot.hits, enable_correlations);
        if (obs_shot.obs_mask_as_u64() != res.obs_mask) {
            num_mistakes++;
        }
        sparse_shot.clear();
        obs_shot.clear();
        num_shots++;
    }
    fprintf(stats_out, "%zu / %zu\n", num_mistakes, num_shots);
    if (stats_out != stdout) {
        fclose(stats_out);
    }
    if (shots_in != stdin) {
        fclose(shots_in);
    }

    auto end = std::chrono::steady_clock::now();
    auto microseconds = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if (time) {
        std::cerr << "Total decoding time: " << (int)microseconds << "us\n";
        std::cerr << "Decoding time per shot: " << (microseconds / num_shots) << "us\n";
    }

    return EXIT_SUCCESS;
}

int pm::main(int argc, const char **argv) {
    const char *command = "";
    if (argc >= 2) {
        command = argv[1];
    }
    try {
        if (strcmp(command, "predict") == 0) {
#ifdef USE_SHMEM
// ===============
            shmem_init();
            if (DEBUG)
                std::cout << "DEBUG: PE " << shmem_my_pe() << " is running main_predict" << std::endl;
// ===============
#endif
            int status = main_predict(argc, argv);

#ifdef USE_SHMEM
// ===============
            shmem_finalize();
// ===============
#endif
            return status;
        }
        if (strcmp(command, "count_mistakes") == 0) {
            return main_count_mistakes(argc, argv);
        }
        if (strcmp(command, "animate") == 0) {
            return pm::main_animation(argc, argv);
        }
    } catch (std::invalid_argument &ex) {
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    std::stringstream ss;
    ss << "Unrecognized command. Available commands are:\n";
    ss << "    pymatching predict --dem file [--in file] [--out file] [--in_format 01|b8|...] [--out_format 01|b8|...] "
          "[--in_includes_appended_observables]\n";
    ss << "    pymatching count_mistakes --dem file [--in file] [--out file] [--in_format 01|b8|...] [--out_format "
          "01|B8|...] [--in_includes_appended_observables] [--obs_in] [--obs_in_format]\n";
    ss << "    pymatching animate "
          "--dets_in <file> "
          "--dets_in_format 01|b8|... "
          "--out_dir <directory> "
          "--dem_in <file> "
          "--held_frames_per_event # "
          "--held_frames_at_start # "
          "--held_frames_at_end # "
          "--max_growth_between_frames # "
          "--max_edge_weight # "
          "--pixels_per_unit_length # "
          "[--dets_in_includes_appended_observables] "
          "[--quiet]";
    throw std::invalid_argument(ss.str());
}
