#include "pymatching/sparse_blossom/driver/helpers/dem_generator.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "stim/gen/gen_rep_code.h"
#include "stim/gen/gen_surface_code.h"
#include "stim/simulators/dem_sampler.h"
#include "stim/util_top/circuit_to_dem.h"

namespace pm {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

MultiObsDemGenerator::MultiObsDemGenerator(
    stim::DetectorErrorModel base_dem,
    int num_observables,
    std::vector<SurgerySpec> surgeries,
    double p_cross,
    int boundary_depth)
    : base_dem_(std::move(base_dem)),
      num_obs_(num_observables),
      surgeries_(std::move(surgeries)),
      p_cross_(p_cross),
      boundary_depth_(boundary_depth),
      max_x_(0.0),
      K_(base_dem_.count_detectors()) {
}

// ---------------------------------------------------------------------------
// Base DEM generation from Stim circuits
// ---------------------------------------------------------------------------

stim::DetectorErrorModel MultiObsDemGenerator::generate_base_dem(
    const std::string& code,
    const std::string& task,
    int distance,
    int rounds,
    double after_clifford_depolarization) {
    stim::CircuitGenParameters params(rounds, distance, task);
    params.after_clifford_depolarization = after_clifford_depolarization;

    stim::GeneratedCircuit gen = (code == "repetition_code")
        ? stim::generate_rep_code_circuit(params)
        : stim::generate_surface_code_circuit(params);

    stim::DetectorErrorModel dem = stim::circuit_to_dem(gen.circuit, {.decompose_errors = true});
    return dem.flattened();
}

// ---------------------------------------------------------------------------
// Index the base DEM to extract detector positions
// ---------------------------------------------------------------------------

void MultiObsDemGenerator::index_base_dem() {
    max_x_ = 0.0;
    round_dets_.clear();
    x_extremes_.clear();
    all_rounds_.clear();

    std::set<uint64_t> all_det_ids;
    for (const auto& inst : base_dem_.instructions) {
        if (inst.type == stim::DemInstructionType::DEM_DETECTOR) {
            for (const auto& t : inst.target_data) {
                if (t.is_relative_detector_id())
                    all_det_ids.insert(t.val());
            }
        }
    }

    auto coords_map = base_dem_.get_detector_coordinates(all_det_ids);

    for (const auto& [det_id, coords] : coords_map) {
        if (coords.size() < 2)
            continue;
        double x = coords[0];
        double rnd = coords.back();
        int obs_id = 0;  // base DEM is single-obs

        if (x > max_x_)
            max_x_ = x;

        auto key = std::make_pair(obs_id, rnd);
        round_dets_[key].push_back({det_id, coords});

        all_rounds_.insert(rnd);

        auto it = x_extremes_.find(obs_id);
        if (it == x_extremes_.end()) {
            x_extremes_[obs_id] = {x, x};
        } else {
            it->second.first = std::max(it->second.first, x);
            it->second.second = std::min(it->second.second, x);
        }
    }

    for (auto& [key, dets] : round_dets_) {
        std::sort(dets.begin(), dets.end(), [](const DetInfo& a, const DetInfo& b) {
            return a.coords[0] < b.coords[0];
        });
    }
}

// ---------------------------------------------------------------------------
// Boundary detector selection
// ---------------------------------------------------------------------------

std::vector<MultiObsDemGenerator::DetInfo> MultiObsDemGenerator::boundary_dets_at_round(
    const std::map<std::pair<int, double>, std::vector<DetInfo>>& round_dets,
    int obs_id,
    double rnd,
    double x_limit,
    bool is_right_side,
    double x_threshold) {
    auto it = round_dets.find({obs_id, rnd});
    if (it == round_dets.end())
        return {};

    std::vector<DetInfo> filtered;
    for (const auto& d : it->second) {
        double x = d.coords[0];
        if (is_right_side) {
            if (x >= x_limit - x_threshold)
                filtered.push_back(d);
        } else {
            if (x <= x_limit + x_threshold)
                filtered.push_back(d);
        }
    }

    if (filtered.empty())
        return filtered;

    size_t ndim = filtered[0].coords.size();
    if (ndim >= 4) {
        // Surface code: sort by y, then x
        if (is_right_side) {
            std::sort(filtered.begin(), filtered.end(), [](const DetInfo& a, const DetInfo& b) {
                if (a.coords[1] != b.coords[1]) return a.coords[1] < b.coords[1];
                return a.coords[0] > b.coords[0];
            });
        } else {
            std::sort(filtered.begin(), filtered.end(), [](const DetInfo& a, const DetInfo& b) {
                if (a.coords[1] != b.coords[1]) return a.coords[1] < b.coords[1];
                return a.coords[0] < b.coords[0];
            });
        }
    } else {
        // Repetition code: sort by x only
        if (is_right_side) {
            std::sort(filtered.begin(), filtered.end(), [](const DetInfo& a, const DetInfo& b) {
                return a.coords[0] > b.coords[0];
            });
        } else {
            std::sort(filtered.begin(), filtered.end(), [](const DetInfo& a, const DetInfo& b) {
                return a.coords[0] < b.coords[0];
            });
        }
    }

    return filtered;
}

// ---------------------------------------------------------------------------
// Boundary error detection
// ---------------------------------------------------------------------------

bool MultiObsDemGenerator::should_drop_error(
    const stim::DemInstruction& inst,
    const std::set<uint64_t>& seam_adj) {
    int det_count = 0;
    uint64_t sole_det = 0;
    for (const auto& t : inst.target_data) {
        if (t.is_separator()) {
            if (det_count == 1 && seam_adj.count(sole_det))
                return true;
            det_count = 0;
        } else if (t.is_relative_detector_id()) {
            sole_det = t.val();
            det_count++;
        }
    }
    if (det_count == 1 && seam_adj.count(sole_det))
        return true;
    return false;
}

// ---------------------------------------------------------------------------
// Append one observable patch to the DEM
// ---------------------------------------------------------------------------

void MultiObsDemGenerator::append_obs_patch(
    stim::DetectorErrorModel& dem,
    int obs_idx,
    const std::set<uint64_t>& seam_adjacent_dets) {
    uint64_t det_offset = (uint64_t)obs_idx * K_;
    double x_shift = obs_idx * (max_x_ + 2);

    for (const auto& inst : base_dem_.instructions) {
        if (inst.type == stim::DemInstructionType::DEM_DETECTOR) {
            uint64_t det_id = 0;
            for (const auto& t : inst.target_data) {
                if (t.is_relative_detector_id())
                    det_id = t.val();
            }

            std::vector<double> coords(inst.arg_data.begin(), inst.arg_data.end());
            coords[0] += x_shift;
            // Inject obs_idx as second-to-last coordinate
            coords.insert(coords.end() - 1, (double)obs_idx);

            dem.append_detector_instruction(
                coords,
                stim::DemTarget::relative_detector_id(det_id + det_offset),
                "");

        } else if (inst.type == stim::DemInstructionType::DEM_ERROR) {
            // Build offset targets, check for boundary drop
            std::vector<stim::DemTarget> new_targets;
            new_targets.reserve(inst.target_data.size());
            for (const auto& t : inst.target_data) {
                if (t.is_relative_detector_id()) {
                    new_targets.push_back(
                        stim::DemTarget::relative_detector_id(t.val() + det_offset));
                } else if (t.is_observable_id()) {
                    new_targets.push_back(
                        stim::DemTarget::observable_id(t.val() + obs_idx));
                } else {
                    new_targets.push_back(t);  // separator
                }
            }

            // Build a temporary instruction to check for drop
            stim::DemInstruction check_inst{
                inst.arg_data,
                stim::SpanRef<const stim::DemTarget>(new_targets.data(), new_targets.data() + new_targets.size()),
                "",
                stim::DemInstructionType::DEM_ERROR};

            if (!seam_adjacent_dets.empty() && should_drop_error(check_inst, seam_adjacent_dets))
                continue;

            dem.append_error_instruction(inst.arg_data[0], new_targets, "");
        }
        // Skip logical_observable, shift_detectors, repeat_block
    }
}

// ---------------------------------------------------------------------------
// Build seam detectors and edges
// ---------------------------------------------------------------------------

void MultiObsDemGenerator::build_seams(
    stim::DetectorErrorModel& seam_dem,
    uint64_t& next_det_id,
    std::map<int, std::set<uint64_t>>& seam_adjacent_by_obs) {
    // We need per-obs-idx versions of round_dets_ and x_extremes_.
    // The base DEM is single-obs (obs_id=0). For each obs_idx, boundary dets are
    // looked up from the base and then offset.

    std::vector<double> sorted_rounds(all_rounds_.begin(), all_rounds_.end());

    for (size_t seam_idx = 0; seam_idx < surgeries_.size(); ++seam_idx) {
        const auto& gate = surgeries_[seam_idx];
        int obs_a = gate.obs_a;
        int obs_b = gate.obs_b;

        // In the base DEM, all detectors are obs_id=0
        // Obs A's max_x (right boundary) and obs B's min_x (left boundary)
        // after x_shift are:
        double x_shift_a = obs_a * (max_x_ + 2);
        double x_shift_b = obs_b * (max_x_ + 2);
        double max_x_a = max_x_ + x_shift_a;
        double min_x_b = 0.0 + x_shift_b;  // base DEM min_x

        // Get base DEM x_extremes for obs 0
        auto base_ext_it = x_extremes_.find(0);
        if (base_ext_it == x_extremes_.end())
            continue;
        double base_max_x = base_ext_it->second.first;
        double base_min_x = base_ext_it->second.second;

        max_x_a = base_max_x + x_shift_a;
        min_x_b = base_min_x + x_shift_b;

        std::map<int, uint64_t> prev_seam_by_site;

        for (int step = 0; step < gate.duration; ++step) {
            double target_rnd = gate.start_round + step;

            // Snap to nearest available round
            double rnd = *std::min_element(sorted_rounds.begin(), sorted_rounds.end(),
                [target_rnd](double a, double b) {
                    return std::abs(a - target_rnd) < std::abs(b - target_rnd);
                });

            // Get boundary dets from base DEM (obs_id=0)
            auto bdets_a = boundary_dets_at_round(
                round_dets_, 0, rnd, base_max_x, true, boundary_depth_);
            auto bdets_b = boundary_dets_at_round(
                round_dets_, 0, rnd, base_min_x, false, boundary_depth_);

            if (bdets_a.empty() || bdets_b.empty()) {
                prev_seam_by_site.clear();
                continue;
            }

            size_t n_pairs = std::min(bdets_a.size(), bdets_b.size());
            std::map<int, uint64_t> new_prev;

            for (size_t site_i = 0; site_i < n_pairs; ++site_i) {
                const auto& da = bdets_a[site_i];
                const auto& db = bdets_b[site_i];

                // Offset detector IDs for the actual obs patches
                uint64_t d_a_id = da.id + (uint64_t)obs_a * K_;
                uint64_t d_b_id = db.id + (uint64_t)obs_b * K_;

                // Seam detector coordinates
                size_t ndim = da.coords.size();
                double x_seam = max_x_a + (min_x_b - max_x_a) * (site_i + 1.0) / (n_pairs + 1.0);
                std::vector<double> seam_coords = {x_seam};
                for (size_t k = 1; k + 1 < ndim; ++k) {
                    seam_coords.push_back(
                        (da.coords[k] + x_shift_a - x_shift_b + db.coords[k]) / 2.0);
                }
                // Actually, we want the midpoint of the *shifted* coords for spatial dims,
                // but the base coords haven't been shifted yet. Let me reconsider.
                // The base coords are unshifted. We need spatial midpoint in the shifted space.
                // For dims 1..ndim-2, the base coords of a and b are the same (same base DEM).
                // So midpoint of y-coords from the base is fine.
                seam_coords.clear();
                seam_coords.push_back(x_seam);
                for (size_t k = 1; k + 1 < ndim; ++k) {
                    seam_coords.push_back((da.coords[k] + db.coords[k]) / 2.0);
                }
                seam_coords.push_back(-double(seam_idx + 1));  // negative seam index
                seam_coords.push_back(rnd);

                seam_dem.append_detector_instruction(
                    seam_coords,
                    stim::DemTarget::relative_detector_id(next_det_id),
                    "");

                // Track seam-adjacent detectors (using offset IDs)
                seam_adjacent_by_obs[obs_a].insert(d_a_id);
                seam_adjacent_by_obs[obs_b].insert(d_b_id);

                // Lateral edges
                std::vector<stim::DemTarget> lat_a = {
                    stim::DemTarget::relative_detector_id(d_a_id),
                    stim::DemTarget::relative_detector_id(next_det_id)};
                seam_dem.append_error_instruction(p_cross_, lat_a, "");

                std::vector<stim::DemTarget> lat_b = {
                    stim::DemTarget::relative_detector_id(next_det_id),
                    stim::DemTarget::relative_detector_id(d_b_id)};
                seam_dem.append_error_instruction(p_cross_, lat_b, "");

                // Temporal edge
                auto prev_it = prev_seam_by_site.find((int)site_i);
                if (prev_it != prev_seam_by_site.end()) {
                    std::vector<stim::DemTarget> temp_edge = {
                        stim::DemTarget::relative_detector_id(prev_it->second),
                        stim::DemTarget::relative_detector_id(next_det_id)};
                    seam_dem.append_error_instruction(p_cross_, temp_edge, "");
                }

                new_prev[(int)site_i] = next_det_id;
                next_det_id++;
            }

            prev_seam_by_site = std::move(new_prev);
        }
    }
}

// ---------------------------------------------------------------------------
// Main generation
// ---------------------------------------------------------------------------

stim::DetectorErrorModel MultiObsDemGenerator::generate() {
    index_base_dem();

    // Build seams first to determine which boundary dets to strip
    stim::DetectorErrorModel seam_dem;
    uint64_t next_det_id = (uint64_t)num_obs_ * K_;
    std::map<int, std::set<uint64_t>> seam_adjacent_by_obs;

    if (!surgeries_.empty()) {
        build_seams(seam_dem, next_det_id, seam_adjacent_by_obs);
    }

    // Build the full DEM: obs patches in order, then seams
    stim::DetectorErrorModel dem;

    for (int obs_idx = 0; obs_idx < num_obs_; ++obs_idx) {
        auto adj_it = seam_adjacent_by_obs.find(obs_idx);
        const std::set<uint64_t>& adj =
            (adj_it != seam_adjacent_by_obs.end()) ? adj_it->second : std::set<uint64_t>{};
        append_obs_patch(dem, obs_idx, adj);
    }

    // Append seam instructions
    for (const auto& inst : seam_dem.instructions) {
        if (inst.type == stim::DemInstructionType::DEM_DETECTOR) {
            uint64_t det_id = 0;
            for (const auto& t : inst.target_data)
                if (t.is_relative_detector_id())
                    det_id = t.val();
            std::vector<double> coords(inst.arg_data.begin(), inst.arg_data.end());
            dem.append_detector_instruction(
                coords,
                stim::DemTarget::relative_detector_id(det_id),
                "");
        } else if (inst.type == stim::DemInstructionType::DEM_ERROR) {
            std::vector<stim::DemTarget> targets(inst.target_data.begin(), inst.target_data.end());
            dem.append_error_instruction(inst.arg_data[0], targets, "");
        }
    }

    return dem;
}

// ---------------------------------------------------------------------------
// DEM file I/O
// ---------------------------------------------------------------------------

void MultiObsDemGenerator::write_dem_file(
    const stim::DetectorErrorModel& dem,
    const std::string& path) {
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Failed to open " + path + " for writing");
    out << dem;
}

stim::DetectorErrorModel MultiObsDemGenerator::read_dem_file(const std::string& path) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f)
        throw std::runtime_error("Failed to open " + path + " for reading");
    auto dem = stim::DetectorErrorModel::from_file(f);
    fclose(f);
    return dem;
}

// ---------------------------------------------------------------------------
// Detection event sampling
// ---------------------------------------------------------------------------

void MultiObsDemGenerator::sample_dem(
    const stim::DetectorErrorModel& dem,
    size_t num_shots,
    const std::string& det_out_path,
    const std::string& obs_out_path,
    uint64_t seed) {
    std::mt19937_64 rng(seed);

    stim::DemSampler<stim::MAX_BITWORD_WIDTH> sampler(dem, std::move(rng), num_shots);

    FILE* det_out = fopen(det_out_path.c_str(), "wb");
    if (!det_out)
        throw std::runtime_error("Failed to open " + det_out_path);

    FILE* obs_out = fopen(obs_out_path.c_str(), "wb");
    if (!obs_out) {
        fclose(det_out);
        throw std::runtime_error("Failed to open " + obs_out_path);
    }

    sampler.sample_write(
        num_shots,
        det_out, stim::SampleFormat::SAMPLE_FORMAT_B8,
        obs_out, stim::SampleFormat::SAMPLE_FORMAT_01,
        nullptr, stim::SampleFormat::SAMPLE_FORMAT_01,  // no error recording
        nullptr, stim::SampleFormat::SAMPLE_FORMAT_01);  // no replay

    fclose(det_out);
    fclose(obs_out);
}

// ---------------------------------------------------------------------------
// Surgery spec parsing
// ---------------------------------------------------------------------------

std::vector<SurgerySpec> MultiObsDemGenerator::parse_spec(
    const std::string& spec, int default_duration) {
    std::vector<SurgerySpec> gates;
    std::istringstream ss(spec);
    std::string entry;

    while (std::getline(ss, entry, ';')) {
        if (entry.empty())
            continue;

        std::istringstream es(entry);
        std::string tok;
        std::vector<int> parts;
        while (std::getline(es, tok, ',')) {
            parts.push_back(std::stoi(tok));
        }

        if (parts.size() < 3 || parts.size() > 4)
            throw std::invalid_argument(
                "Invalid surgery spec entry '" + entry + "'; expected 'obs_a,obs_b,start[,duration]'");

        int obs_a = parts[0], obs_b = parts[1], start = parts[2];
        int dur = (parts.size() == 4) ? parts[3] : default_duration;

        if (obs_a == obs_b)
            throw std::invalid_argument("obs_a == obs_b in surgery spec '" + entry + "'");
        if (dur < 1)
            throw std::invalid_argument("duration must be >= 1 in surgery spec '" + entry + "'");

        gates.push_back({std::min(obs_a, obs_b), std::max(obs_a, obs_b), start, dur});
    }
    return gates;
}

// ---------------------------------------------------------------------------
// Surgery presets
// ---------------------------------------------------------------------------

static std::vector<SurgerySpec> reduction_adder(int carry, std::vector<int> op_a, std::vector<int> op_b, int t_start) {
    std::vector<SurgerySpec> seams;
    int t = t_start;
    // MAJ forward
    for (int i = 0; i < 4; ++i) {
        seams.push_back({std::min(carry, op_a[i]), std::max(carry, op_a[i]), t, 21}); t += 42;
        seams.push_back({std::min(carry, op_b[i]), std::max(carry, op_b[i]), t, 21}); t += 42;
    }
    // UMA backward
    for (int i = 3; i >= 0; --i) {
        seams.push_back({std::min(carry, op_b[i]), std::max(carry, op_b[i]), t, 21}); t += 42;
        seams.push_back({std::min(carry, op_a[i]), std::max(carry, op_a[i]), t, 21}); t += 42;
    }
    return seams;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_36obs() {
    struct IntraEntry { int carry_off, other_off, round, duration; };
    std::vector<IntraEntry> intra_pattern = {
        {8, 0, 0, 21}, {8, 4, 42, 21}, {8, 1, 84, 21}, {8, 5, 126, 21},
        {8, 2, 168, 21}, {8, 6, 210, 21}, {8, 3, 252, 21}, {8, 7, 294, 21},
        {8, 7, 336, 21}, {8, 3, 378, 21}, {8, 6, 420, 21}, {8, 2, 462, 21},
        {8, 5, 504, 21}, {8, 1, 546, 21}, {8, 4, 588, 21}, {8, 0, 630, 21}
    };

    std::vector<SurgerySpec> gates;
    for (int block = 0; block < 4; ++block) {
        int base = 9 * block;
        for (const auto& e : intra_pattern) {
            int a = base + e.carry_off, b = base + e.other_off;
            gates.push_back({std::min(a, b), std::max(a, b), e.round, e.duration});
        }
    }

    int T_L1 = 672;
    auto adder_E = reduction_adder(8, {4, 5, 6, 7}, {13, 14, 15, 16}, T_L1);
    auto adder_F = reduction_adder(26, {22, 23, 24, 25}, {31, 32, 33, 34}, T_L1);
    int T_L2 = 1386;
    auto adder_G = reduction_adder(17, {13, 14, 15, 16}, {31, 32, 33, 34}, T_L2);

    gates.insert(gates.end(), adder_E.begin(), adder_E.end());
    gates.insert(gates.end(), adder_F.begin(), adder_F.end());
    gates.insert(gates.end(), adder_G.begin(), adder_G.end());

    return gates;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_72obs(int M) {
    // Copy 1: obs 0-35, identical round schedule to preset_36obs().
    std::vector<SurgerySpec> gates = preset_36obs();

    // Copy 2: obs 36-71, same schedule, offset by +36. Independent of copy 1
    // throughout (no cross-copy gates in this phase).
    auto copy2 = preset_36obs();
    for (auto& g : copy2) {
        g.obs_a += 36;
        g.obs_b += 36;
    }
    gates.insert(gates.end(), copy2.begin(), copy2.end());

    // 2*M idle gap: rounds [2058, 2058 + 2*M) have no surgeries at all (both
    // copies just continue plain QEC rounds independently).
    int T_break_end = 2058 + 2 * M;

    // Final reduction: a new top-level carry (obs 72) combining copy 1's and
    // copy 2's adder_G operand groups {13,14,15,16}/{49,50,51,52}. Mirrors
    // adder_G itself (which combined blocks 1 and 3 the same way), just one
    // level higher, combining the two 36obs copies instead of two blocks.
    auto adder_top = reduction_adder(72, {13, 14, 15, 16}, {49, 50, 51, 52}, T_break_end);
    gates.insert(gates.end(), adder_top.begin(), adder_top.end());

    return gates;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_36obs_repeated(int repeats, int M) {
    int period = 2058 + 2 * M;  // one repetition's span + 2*M idle gap before the next
    std::vector<SurgerySpec> gates;
    for (int rep = 0; rep < repeats; ++rep) {
        auto one = preset_36obs();
        int round_offset = rep * period;
        for (auto& g : one) {
            g.start_round += round_offset;
        }
        gates.insert(gates.end(), one.begin(), one.end());
    }
    return gates;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_18obs() {
    std::vector<SurgerySpec> gates = {
        {0, 16, 0, 21}, {8, 16, 42, 21}, {1, 16, 84, 21}, {9, 16, 126, 21},
        {2, 16, 168, 21}, {10, 16, 210, 21}, {3, 16, 252, 21}, {11, 16, 294, 21},
        {4, 16, 336, 21}, {12, 16, 378, 21}, {5, 16, 420, 21}, {13, 16, 462, 21},
        {6, 16, 504, 21}, {14, 16, 546, 21}, {7, 16, 588, 21}, {15, 16, 630, 21},
        {15, 16, 672, 21}, {7, 16, 714, 21}, {14, 16, 756, 21}, {6, 16, 798, 21},
        {5, 16, 840, 21}, {5, 16, 882, 21}, {4, 16, 924, 21}, {4, 16, 966, 21},
        {3, 16, 1008, 21}, {3, 16, 1050, 21}, {2, 16, 1092, 21}, {2, 16, 1134, 21},
        {1, 16, 1176, 21}, {1, 16, 1218, 21}, {0, 16, 1260, 21}, {0, 16, 1302, 21}
    };
    return gates;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_24obs(int M, int duration) {
    struct IntraEntry { int part_idx, off_a, off_b; };
    std::vector<IntraEntry> intra_pattern = {
        {0, 0, 3}, {2, 0, 1}, {4, 0, 1}, {6, 0, 2}, {8, 0, 4}, {10, 0, 5}
    };

    std::vector<SurgerySpec> gates;
    for (int block = 0; block < 4; ++block) {
        int base = 6 * block;
        for (const auto& e : intra_pattern) {
            int a = base + e.off_a, b = base + e.off_b;
            gates.push_back({std::min(a, b), std::max(a, b), e.part_idx * M, duration});
        }
    }

    gates.push_back({5, 11, 12 * M, duration});
    gates.push_back({17, 23, 12 * M, duration});
    gates.push_back({11, 23, 15 * M, duration});

    return gates;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_48obs(int M, int duration) {
    struct IntraEntry { int part_idx, off_a, off_b; };
    std::vector<IntraEntry> intra_pattern = {
        {0, 0, 3}, {2, 0, 1}, {4, 0, 1}, {6, 0, 2}, {8, 0, 4}, {10, 0, 5}
    };

    std::vector<SurgerySpec> gates;
    for (int block = 0; block < 8; ++block) {
        int base = 6 * block;
        for (const auto& e : intra_pattern) {
            int a = base + e.off_a, b = base + e.off_b;
            gates.push_back({std::min(a, b), std::max(a, b), e.part_idx * M, duration});
        }
    }

    gates.push_back({5, 11, 12 * M, duration});
    gates.push_back({17, 23, 12 * M, duration});
    gates.push_back({29, 35, 12 * M, duration});
    gates.push_back({41, 47, 12 * M, duration});
    gates.push_back({11, 23, 15 * M, duration});
    gates.push_back({35, 47, 15 * M, duration});
    gates.push_back({23, 47, 18 * M, duration});

    return gates;
}

std::vector<SurgerySpec> MultiObsDemGenerator::preset_64obs(int M, int duration) {
    struct IntraEntry { int part_idx, off_a, off_b; };
    std::vector<IntraEntry> intra_pattern = {
        {0, 0, 3}, {0, 4, 7}, {2, 0, 1}, {2, 4, 5},
        {4, 0, 1}, {4, 4, 5}, {6, 0, 2}, {6, 4, 6}, {8, 0, 4}
    };

    std::vector<SurgerySpec> gates;
    for (int shift : {0, 15, 30, 45}) {
        int round_shift = shift * M;
        for (int block = 0; block < 8; ++block) {
            int base = 8 * block;
            for (const auto& e : intra_pattern) {
                int a = base + e.off_a, b = base + e.off_b;
                gates.push_back({std::min(a, b), std::max(a, b),
                                 e.part_idx * M + round_shift, duration});
            }
        }

        gates.push_back({7, 15, 11 * M + round_shift, duration});
        gates.push_back({23, 31, 11 * M + round_shift, duration});
        gates.push_back({39, 47, 11 * M + round_shift, duration});
        gates.push_back({55, 63, 11 * M + round_shift, duration});
    }

    gates.push_back({7, 23, 60 * M, duration});
    gates.push_back({39, 55, 60 * M, duration});
    gates.push_back({7, 39, 63 * M, duration});

    return gates;
}

}  // namespace pm
