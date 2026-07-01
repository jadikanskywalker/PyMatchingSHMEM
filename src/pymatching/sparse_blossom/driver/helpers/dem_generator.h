#ifndef _PM_DEM_GENERATOR_H
#define _PM_DEM_GENERATOR_H

#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "stim.h"

namespace pm {

struct SurgerySpec {
    int obs_a, obs_b;
    int start_round;
    int duration;
};

class MultiObsDemGenerator {
   public:
    MultiObsDemGenerator(
        stim::DetectorErrorModel base_dem,
        int num_observables,
        std::vector<SurgerySpec> surgeries,
        double p_cross,
        int boundary_depth = 2);

    stim::DetectorErrorModel generate();

    static std::vector<SurgerySpec> preset_18obs();
    static std::vector<SurgerySpec> preset_24obs(int M = 42, int duration = 21);
    static std::vector<SurgerySpec> preset_36obs();
    static std::vector<SurgerySpec> preset_48obs(int M = 32, int duration = 21);
    static std::vector<SurgerySpec> preset_64obs(int M = 22, int duration = 21);
    // Two independent copies of preset_36obs() (obs 0-35 and obs 36-71, same round
    // schedule as preset_36obs, ending at round 2058), followed by a 2*M idle gap
    // with no surgeries, then one final reduction_adder (a new top-level carry at
    // obs 72) combining each copy's adder_G operand group -- mirrors adder_G itself,
    // just one level higher. Total observables: 73.
    static std::vector<SurgerySpec> preset_72obs(int M = 21);

    static std::vector<SurgerySpec> parse_spec(const std::string& spec, int default_duration = 1);

    static void write_dem_file(const stim::DetectorErrorModel& dem, const std::string& path);
    static stim::DetectorErrorModel read_dem_file(const std::string& path);

    static stim::DetectorErrorModel generate_base_dem(
        const std::string& code,
        const std::string& task,
        int distance,
        int rounds,
        double after_clifford_depolarization);

    static void sample_dem(
        const stim::DetectorErrorModel& dem,
        size_t num_shots,
        const std::string& det_out_path,
        const std::string& obs_out_path,
        uint64_t seed = 42);

   private:
    stim::DetectorErrorModel base_dem_;
    int num_obs_;
    std::vector<SurgerySpec> surgeries_;
    double p_cross_;
    int boundary_depth_;
    double max_x_;
    uint64_t K_;

    struct DetInfo {
        uint64_t id;
        std::vector<double> coords;
    };

    // (obs_id, round) -> dets sorted by x
    std::map<std::pair<int, double>, std::vector<DetInfo>> round_dets_;
    std::map<int, std::pair<double, double>> x_extremes_;
    std::set<double> all_rounds_;

    void index_base_dem();

    void append_obs_patch(
        stim::DetectorErrorModel& dem,
        int obs_idx,
        const std::set<uint64_t>& seam_adjacent_dets);

    void build_seams(
        stim::DetectorErrorModel& seam_dem,
        uint64_t& next_det_id,
        std::map<int, std::set<uint64_t>>& seam_adjacent_by_obs);

    static bool should_drop_error(
        const stim::DemInstruction& inst,
        const std::set<uint64_t>& seam_adj);

    static std::vector<DetInfo> boundary_dets_at_round(
        const std::map<std::pair<int, double>, std::vector<DetInfo>>& round_dets,
        int obs_id,
        double rnd,
        double x_limit,
        bool is_right_side,
        double x_threshold);
};

}  // namespace pm

#endif
