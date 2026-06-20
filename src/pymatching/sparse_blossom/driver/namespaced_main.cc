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

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "pymatching/sparse_blossom/diagram/animation_main.h"
#include "pymatching/sparse_blossom/driver/helpers/fast_b8_reader.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "stim.h"

#define OUTPUT_DECODING_TIME

#ifdef USE_THREADS
#include <omp.h>
#include <fstream>
#include "../config_parallel.h"
#include "../diagram/mwpm_diagram.h"
#include "pymatching/sparse_blossom/driver/helpers/dem_generator.h"
#endif

#ifdef USE_SHMEM
#include <shmem.h>
#endif

// #include "profiling/profiling_json.h"

int main_predict(int argc, const char** argv) {
    stim::check_for_unknown_arguments(
        {
            "--in",
            "--in_format",
            "--in_includes_appended_observables",
            "--out",
            "--out_format",
            "--dem",
            "--enable_correlations"
#ifdef OUTPUT_DECODING_TIME
            ,
            "--num_repeats"
#endif
#ifdef USE_THREADS
            ,
            "--rounds_per_partition",
            "--use_threads",
            "--obs_coors_included"
#endif
#ifdef ENABLE_DRAW_FLAGS
            ,
            "--draw_frames"
#endif
#ifdef USE_SHMEM
            ,
            "--cross_rank_fusion_window_size",
            "--task_division_strategy",
            "--gen_code",
            "--gen_task",
            "--gen_distance",
            "--gen_rounds",
            "--gen_num_obs",
            "--gen_depolarization",
            "--gen_surgery_preset",
            "--gen_surgery_spec",
            "--gen_surgery_duration",
            "--gen_p_cross",
            "--gen_boundary_depth",
            "--dem_cache_path",
            "--gen_sample_shots",
            "--gen_sample_seed",
            "--gen_det_out",
            "--gen_obs_out"
#endif
        },
        {},
        "predict",
        argc,
        argv);

    FILE* predictions_out = nullptr;
#ifdef USE_SHMEM
    const char* out_fn_base_const = stim::find_argument("--out", argc, argv);
    std::string out_fn_base = out_fn_base_const ? out_fn_base_const : "";
    if (!out_fn_base.empty() && out_fn_base != "stdout") {
        int pid = shmem_my_pe();
        std::string filename = out_fn_base;
        size_t last_dot = filename.find_last_of(".");
        if (last_dot == std::string::npos) {
            filename += "_pe" + std::to_string(pid);
        } else {
            filename.insert(last_dot, "_pe" + std::to_string(pid));
        }
        predictions_out = fopen(filename.c_str(), "wb");
        if (predictions_out == nullptr) {
             throw std::invalid_argument("Failed to open " + filename);
        }
    } else {
        predictions_out = stim::find_open_file_argument("--out", stdout, "wb", argc, argv);
    }
#else
    predictions_out = stim::find_open_file_argument("--out", stdout, "wb", argc, argv);
#endif
    stim::FileFormatData shots_in_format =
        stim::find_enum_argument("--in_format", "b8", stim::format_name_to_enum_map(), argc, argv);
    stim::FileFormatData predictions_out_format =
        stim::find_enum_argument("--out_format", "01", stim::format_name_to_enum_map(), argc, argv);
    bool append_obs = stim::find_bool_argument("--in_includes_appended_observables", argc, argv);
    bool enable_correlations = stim::find_bool_argument("--enable_correlations", argc, argv);
#ifdef OUTPUT_DECODING_TIME
    int num_repeats = (int)stim::find_int64_argument("--num_repeats", 1, 1, INT64_MAX, argc, argv);
#endif
#ifdef USE_THREADS
    // ===============
    config_parallel::M = stim::find_int64_argument("--rounds_per_partition", 10, 1, INT64_MAX, argc, argv);
    config_parallel::obs_coors_included = stim::find_bool_argument("--obs_coors_included", argc, argv);
    bool use_threads = stim::find_bool_argument("--use_threads", argc, argv);
// ===============
#endif
#ifdef ENABLE_DRAW_FLAGS
    bool draw_frames = stim::find_bool_argument("--draw_frames", argc, argv);
#endif
#ifdef USE_SHMEM
    config_parallel::k = stim::find_int64_argument("--cross_rank_fusion_window_size", 1, 0, INT64_MAX, argc, argv);
    {
        const char* ds = stim::find_argument("--task_division_strategy", argc, argv);
        if (ds == nullptr || strcmp(ds, "round") == 0) {
            config_parallel::division_strategy = config_parallel::ROUND;
        } else if (strcmp(ds, "observable") == 0) {
            config_parallel::division_strategy = config_parallel::OBS;
        } else {
            throw std::invalid_argument(
                std::string("Unknown --task_division_strategy: ") + ds +
                ". Use 'round' or 'observable'.");
        }
    }
#endif

    stim::DetectorErrorModel dem;

    // DEM generation path: build in-memory or load from cache
    const char* gen_code_arg = stim::find_argument("--gen_code", argc, argv);
    const char* dem_cache_arg = stim::find_argument("--dem_cache_path", argc, argv);

    if (gen_code_arg != nullptr) {
        std::string dem_cache_path = dem_cache_arg ? dem_cache_arg : "";

        if (!dem_cache_path.empty()) {
            FILE* cache_check = fopen(dem_cache_path.c_str(), "r");
            if (cache_check) {
                fclose(cache_check);
                dem = pm::MultiObsDemGenerator::read_dem_file(dem_cache_path);
                std::cerr << "Loaded cached DEM from " << dem_cache_path << std::endl;
                goto dem_ready;
            }
        }

        {
            std::string gen_code = gen_code_arg;
            const char* gen_task_arg = stim::find_argument("--gen_task", argc, argv);
            std::string gen_task = gen_task_arg ? gen_task_arg : "memory";
            int gen_distance = (int)stim::find_int64_argument("--gen_distance", 5, 1, INT64_MAX, argc, argv);
            int gen_rounds = (int)stim::find_int64_argument("--gen_rounds", 10, 1, INT64_MAX, argc, argv);
            int gen_num_obs = (int)stim::find_int64_argument("--gen_num_obs", 2, 1, INT64_MAX, argc, argv);
            double gen_depol = stim::find_float_argument("--gen_depolarization", 0.01, 0.0, 1.0, argc, argv);
            int gen_surgery_dur = (int)stim::find_int64_argument("--gen_surgery_duration", 8, 1, INT64_MAX, argc, argv);
            double gen_p_cross_val = stim::find_float_argument("--gen_p_cross", -1.0, -1.0, 1.0, argc, argv);
            double gen_p_cross = (gen_p_cross_val < 0) ? std::min(0.4, 0.2 * gen_depol) : gen_p_cross_val;
            int gen_boundary_depth = (int)stim::find_int64_argument("--gen_boundary_depth", 2, 0, INT64_MAX, argc, argv);

            auto base_dem = pm::MultiObsDemGenerator::generate_base_dem(
                gen_code, gen_task, gen_distance, gen_rounds, gen_depol);

            std::vector<pm::SurgerySpec> surgeries;
            const char* preset_arg = stim::find_argument("--gen_surgery_preset", argc, argv);
            const char* spec_arg = stim::find_argument("--gen_surgery_spec", argc, argv);

            if (preset_arg) {
                std::string preset = preset_arg;
                if (preset == "18obs") surgeries = pm::MultiObsDemGenerator::preset_18obs();
                else if (preset == "24obs") surgeries = pm::MultiObsDemGenerator::preset_24obs();
                else if (preset == "36obs") surgeries = pm::MultiObsDemGenerator::preset_36obs();
                else if (preset == "48obs") surgeries = pm::MultiObsDemGenerator::preset_48obs();
                else if (preset == "64obs") surgeries = pm::MultiObsDemGenerator::preset_64obs();
                else throw std::invalid_argument("Unknown --gen_surgery_preset: " + preset);
            } else if (spec_arg) {
                surgeries = pm::MultiObsDemGenerator::parse_spec(spec_arg, gen_surgery_dur);
            }

            pm::MultiObsDemGenerator generator(
                std::move(base_dem), gen_num_obs, std::move(surgeries), gen_p_cross, gen_boundary_depth);
            dem = generator.generate();

            std::cerr << "Generated DEM: " << dem.count_detectors() << " detectors, "
                      << dem.count_observables() << " observables" << std::endl;

            if (!dem_cache_path.empty()) {
                pm::MultiObsDemGenerator::write_dem_file(dem, dem_cache_path);
                std::cerr << "Cached DEM to " << dem_cache_path << std::endl;
            }

        }
        dem_ready:;

        // Handle sampling if requested (runs whether DEM was generated or loaded from cache)
        const char* gen_det_out = stim::find_argument("--gen_det_out", argc, argv);
        const char* gen_obs_out = stim::find_argument("--gen_obs_out", argc, argv);
        if (gen_det_out && gen_obs_out) {
            size_t sample_shots = (size_t)stim::find_int64_argument("--gen_sample_shots", 1, 1, INT64_MAX, argc, argv);
            uint64_t sample_seed = (uint64_t)stim::find_int64_argument("--gen_sample_seed", 42, 0, INT64_MAX, argc, argv);
            pm::MultiObsDemGenerator::sample_dem(dem, sample_shots, gen_det_out, gen_obs_out, sample_seed);
            std::cerr << "Sampled " << sample_shots << " shots" << std::endl;
        }
    } else {
        FILE* dem_file = stim::find_open_file_argument("--dem", nullptr, "r", argc, argv);
        dem = stim::DetectorErrorModel::from_file(dem_file);
        fclose(dem_file);
    }

    FILE* shots_in = stim::find_open_file_argument("--in", stdin, "rb", argc, argv);
    size_t num_obs = dem.count_observables();
    size_t num_detectors = dem.count_detectors();
    auto reader = stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>::make(
        shots_in, shots_in_format.id, 0, num_detectors, append_obs * num_obs);
    auto writer = stim::MeasureRecordWriter::make(predictions_out, predictions_out_format.id);
    writer->begin_result_type('L');

    pm::weight_int num_buckets = pm::NUM_DISTINCT_WEIGHTS;

#ifdef USE_THREADS
    pm::DecodingUnit decoding_unit(
        std::move(reader),
        std::move(writer),
        dem,
        num_buckets,
        enable_correlations,
        enable_correlations
#ifdef ENABLE_DRAW_FLAGS
        , draw_frames
#endif
    );
#ifdef USE_SHMEM
    shmem_barrier_all();
#endif
    if (DEBUG) {
#ifdef ENABLE_DRAW_FLAGS
        pm::setup_output_dirs(draw_frames, use_threads);
#else
        pm::setup_output_dirs(use_threads);
#endif
    }
#else
    auto mwpm = pm::detector_error_model_to_mwpm(
        dem,
        num_buckets,
        /*ensure_search_flooder_included=*/enable_correlations,
        /*enable_correlations=*/enable_correlations);

    stim::SparseShot sparse_shot;
    sparse_shot.clear();
    pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
#endif

#ifdef USE_THREADS
    decoding_unit.decode_shots();
#else
    while (pm::start_and_read_entire_record_buffered(*reader, sparse_shot)) {
        pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations);
        for (size_t k = 0; k < num_obs; k++) {
            writer->write_bit(res.obs_crossed[k]);
        }
        writer->write_end();
        sparse_shot.clear();
        res.reset();
    }
#endif

#ifdef OUTPUT_DECODING_TIME
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    using std::chrono::steady_clock;

    double total_ms = 0.0;
    for (int rep = 0; rep < num_repeats; rep++) {
        fseek(shots_in, 0, SEEK_SET);
        // reader = stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>::make(
        //     shots_in, shots_in_format.id, 0, num_detectors, append_obs * num_obs);
        fflush(predictions_out);
        fseek(predictions_out, 0, SEEK_SET);
        ftruncate(fileno(predictions_out), 0);
#ifdef USE_THREADS
        // decoding_unit.shot_buffer->reader = std::move(reader);
        decoding_unit.reset();
#endif
        auto t1 = steady_clock::now();
#endif
//---------Decode
#ifdef USE_THREADS
        decoding_unit.decode_shots();
#else
        while (pm::start_and_read_entire_record_buffered(*reader, sparse_shot)) {
            pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight, enable_correlations);
            for (size_t k = 0; k < num_obs; k++) {
                writer->write_bit(res.obs_crossed[k]);
            }
            writer->write_end();
            sparse_shot.clear();
            res.reset();
        }
#endif
//---------Decode
#ifdef OUTPUT_DECODING_TIME
        auto t2 = steady_clock::now();
        total_ms += duration<double, std::milli>(t2 - t1).count();
        std::cout << duration<double, std::milli>(t2 - t1).count() << "ms\n" << std::flush;
    }
    std::cout << "Decoding time: " << total_ms / (double)num_repeats << "ms\n";
#endif

    if (predictions_out != stdout) {
        fclose(predictions_out);
    }
    if (shots_in != stdin) {
        fclose(shots_in);
    }

    return EXIT_SUCCESS;
}

int main_count_mistakes(int argc, const char** argv) {
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

    FILE* shots_in = stim::find_open_file_argument("--in", stdin, "rb", argc, argv);
    FILE* obs_in = stim::find_open_file_argument("--obs_in", stdin, "rb", argc, argv);
    FILE* stats_out = stim::find_open_file_argument("--out", stdout, "wb", argc, argv);
    FILE* dem_file = stim::find_open_file_argument("--dem", nullptr, "r", argc, argv);
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

int pm::main(int argc, const char** argv) {
    const char* command = "";
    if (argc >= 2) {
        command = argv[1];
    }
    try {
        if (strcmp(command, "predict") == 0) {
#ifdef USE_SHMEM
            // ===============
            int provided;
            shmem_init_thread(SHMEM_THREAD_SERIALIZED, &provided);
            if (provided < SHMEM_THREAD_SERIALIZED) {
                std::cerr << "Warning: OpenSHMEM failed to init with SHMEM_THREAD_SERIALED." << std::endl;
            }
#endif
            // profiling_json_init();
            int status = main_predict(argc, argv);
            // profiling_json_finalize();
#ifdef USE_SHMEM
            shmem_finalize();
// ===============
#endif
            return status;
        }
        if (strcmp(command, "count_mistakes") == 0) {
            return main_count_mistakes(argc, argv);
        }
        if (strcmp(command, "animate") == 0) {
            #ifdef USE_SHMEM
            if (DEBUG) std::cout << "calling animate\n";
            #endif
            return pm::main_animation(argc, argv);
        }
    } catch (std::invalid_argument& ex) {
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

// #include "profiling/profiling_json.c"
