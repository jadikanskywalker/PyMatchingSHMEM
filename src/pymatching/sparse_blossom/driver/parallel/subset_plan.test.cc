#include "pymatching/sparse_blossom/driver/parallel/subset_plan.h"

#include <gtest/gtest.h>

#include <fstream>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/helpers/dem_generator.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"

namespace {

// Sets the process-wide partitioning config for a test and restores it afterwards.
struct ScopedObsConfig {
    int old_M = config_parallel::M;
    config_parallel::div_strgy old_strategy = config_parallel::division_strategy;
    explicit ScopedObsConfig(int M) {
        config_parallel::M = M;
        config_parallel::division_strategy = config_parallel::OBS;
    }
    ~ScopedObsConfig() {
        config_parallel::M = old_M;
        config_parallel::division_strategy = old_strategy;
    }
};

// The (code, task) pairs every equivalence test runs over. Distance 3 keeps the graphs tiny.
const std::vector<std::pair<std::string, std::string>> kCodes = {
    {"repetition_code", "memory"},
    {"surface_code", "rotated_memory_x"},
};

struct BuiltGraph {
    pm::SharedMatchingGraph graph;
#ifdef USE_SHMEM
    std::vector<pm::DetectorNodeEphemeralFields> ephemeral;  // backing store for node ephemeral fields
#endif
};

// A multi-observable repetition-code DEM -> UserGraph (OBS partitioning) -> SharedMatchingGraph.
BuiltGraph build_graph(
    const std::string& code, const std::string& task, int gen_rounds, int num_obs,
    const std::vector<pm::SurgerySpec>& specs) {
    auto base = pm::MultiObsDemGenerator::generate_base_dem(code, task, 3, gen_rounds, 0.01);
    pm::MultiObsDemGenerator generator(std::move(base), num_obs, specs, /*p_cross=*/0.02, /*boundary_depth=*/2);
    stim::DetectorErrorModel dem = generator.generate();
    pm::UserGraph user_graph = pm::detector_error_model_to_user_graph(dem, false, pm::NUM_DISTINCT_WEIGHTS);
    BuiltGraph out;
#ifdef USE_SHMEM
    out.ephemeral.resize(user_graph.nodes.size());
    out.graph = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, out.ephemeral.data());
#else
    out.graph = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
#endif
    return out;
}

}  // namespace

TEST(SubsetPlan, BoundEffectiveK) {
    std::string note;
    // 2 PEs: k capped at num_partitions/2.
    EXPECT_EQ(pm::bound_effective_k(9, 2, 10, 5, true, &note), 5);
    EXPECT_EQ(note, "NOTE: k set to 5 for 2 PEs");
    EXPECT_EQ(pm::bound_effective_k(3, 2, 10, 5, true), 3);
    // >2 PEs: OBS caps at p_per_obs_patch, ROUND at p_per_pe/2.
    note.clear();
    EXPECT_EQ(pm::bound_effective_k(9, 4, 40, 6, true, &note), 6);
    EXPECT_EQ(note, "NOTE: k bounded to 6");
    EXPECT_EQ(pm::bound_effective_k(9, 4, 40, 6, false), 5);
    EXPECT_EQ(pm::bound_effective_k(1, 4, 40, 6, false), 1);
    // Fewer than 2 partitions per PE (with more than 2 PEs) is rejected.
    EXPECT_THROW(pm::bound_effective_k(1, 8, 9, 3, true), std::invalid_argument);
}

TEST(SubsetPlan, ConflictResolutionSplitsOverlap) {
    // Seams 0 and 1 share observable 1; their windows [2,7] and [5,9] overlap.
    std::vector<pm::SeamInfo> seams = {{0, 1, 2, 7}, {1, 2, 5, 9}, {3, 4, 0, 3}};
    pm::resolve_seam_window_conflicts(seams, 0);
    int mid = (7 + 5 + 1) / 2;
    EXPECT_EQ(seams[0].vb_left, 2);
    EXPECT_EQ(seams[0].vb_right, mid);
    EXPECT_EQ(seams[1].vb_left, mid);
    EXPECT_EQ(seams[1].vb_right, 9);
    // Unrelated seam untouched.
    EXPECT_EQ(seams[2].vb_left, 0);
    EXPECT_EQ(seams[2].vb_right, 3);
}

TEST(SubsetPlan, ConflictResolutionDegenerateThrows) {
    // Seam 1's window [1,3] is nested inside seam 0's [0,4]; splitting at the midpoint leaves it empty.
    std::vector<pm::SeamInfo> seams = {{0, 1, 0, 4}, {1, 2, 1, 3}};
    EXPECT_THROW(pm::resolve_seam_window_conflicts(seams, 0), std::invalid_argument);
}

TEST(SubsetPlan, ExtractionUnitSpans) {
    std::vector<pm::SeamInfo> seams = {{0, 1, 3, 9}};
    pm::assign_extraction_unit_spans(seams, 4);
    EXPECT_EQ(seams[0].unit_lo, 1);
    EXPECT_EQ(seams[0].unit_hi, 2);
}

TEST(SubsetPlan, PartitionVbBounds) {
    // partitions 0,1 | vb 0 | partition 2 | seam vbs 1,2
    std::vector<int> node_part_id = {0, 0, 1, -1, -1, 2, 2, 2, -2, -3, -3};
    std::vector<std::pair<size_t, size_t>> pb, vbb;
    pm::compute_partition_vb_bounds(node_part_id, 3, 3, pb, vbb);
    ASSERT_EQ(pb.size(), 3u);
    ASSERT_EQ(vbb.size(), 3u);
    EXPECT_EQ(pb[0], (std::pair<size_t, size_t>{0, 1}));
    EXPECT_EQ(pb[1], (std::pair<size_t, size_t>{2, 2}));
    EXPECT_EQ(pb[2], (std::pair<size_t, size_t>{5, 7}));
    EXPECT_EQ(vbb[0], (std::pair<size_t, size_t>{3, 4}));
    EXPECT_EQ(vbb[1], (std::pair<size_t, size_t>{8, 8}));
    EXPECT_EQ(vbb[2], (std::pair<size_t, size_t>{9, 10}));
}

// The spec-derived seam windows must equal the windows found by walking the graph generated from the
// same spec, for every (rounds, M, k) and a range of seam layouts.
TEST(SubsetPlan, SpecSeamInfosMatchGraphWalk) {
    const std::vector<int> gen_rounds_list = {12, 13, 14, 19};
    const std::vector<int> M_list = {3, 4, 5, 7};
    const int num_obs = 4;
    const std::vector<std::vector<pm::SurgerySpec>> layouts = {
        {{0, 1, 2, 3}},                                   // single seam mid-run
        {{0, 1, 2, 3}, {2, 3, 5, 6}},                     // two independent seams
        {{0, 1, 4, 5}, {1, 2, 4, 5}},                     // shared observable (conflict resolution)
        {{0, 1, 4, 5}, {1, 2, 6, 5}, {2, 3, 9, 4}},       // chain
        {{0, 1, 3, 1}},                                   // single-round seam
        {{0, 1, 11, 6}},                                  // runs past the last round (snapped)
        {{0, 1, 0, 2}},                                   // starts at round 0
    };
    for (const auto& [code, task] : kCodes)
    for (int gen_rounds : gen_rounds_list) {
        for (int M : M_list) {
            ScopedObsConfig cfg(M);
            for (size_t li = 0; li < layouts.size(); ++li) {
                const auto& specs = layouts[li];
                SCOPED_TRACE(code + " gen_rounds=" + std::to_string(gen_rounds) + " M=" + std::to_string(M) +
                             " layout=" + std::to_string(li));
                BuiltGraph built = build_graph(code, task, gen_rounds, num_obs, specs);
                const pm::SharedMatchingGraph& g = built.graph;
                ASSERT_EQ((int)g.num_obs_patches, num_obs);
                ASSERT_EQ((int)(g.num_virtual_boundaries - g.vb_per_obs_patch * g.num_obs_patches), (int)specs.size());
                for (int k : {1, 2}) {
                    auto from_graph = pm::compute_seam_infos_from_graph(g, k);
                    auto from_spec = pm::compute_seam_infos_from_spec(
                        specs, (int)g.num_rounds, M, (int)g.p_per_obs_patch, (int)g.vb_per_obs_patch, k);
                    ASSERT_EQ(from_graph.size(), from_spec.size());
                    for (size_t s = 0; s < from_graph.size(); ++s) {
                        SCOPED_TRACE("k=" + std::to_string(k) + " seam=" + std::to_string(s));
                        EXPECT_EQ(from_graph[s].oi, from_spec[s].oi);
                        EXPECT_EQ(from_graph[s].oj, from_spec[s].oj);
                        EXPECT_EQ(from_graph[s].vb_left, from_spec[s].vb_left);
                        EXPECT_EQ(from_graph[s].vb_right, from_spec[s].vb_right);
                    }
                    // Conflict resolution is a pure function of the windows, so it must agree too
                    // (including whether it throws).
                    bool threw_graph = false, threw_spec = false;
                    try { pm::resolve_seam_window_conflicts(from_graph, 0); } catch (const std::invalid_argument&) { threw_graph = true; }
                    try { pm::resolve_seam_window_conflicts(from_spec, 0); } catch (const std::invalid_argument&) { threw_spec = true; }
                    EXPECT_EQ(threw_graph, threw_spec);
                }
            }
        }
    }
}

TEST(SubsetPlan, SpecSeamInfosRejectsBadArguments) {
    std::vector<pm::SurgerySpec> specs = {{0, 1, 2, 3}};
    EXPECT_THROW(pm::compute_seam_infos_from_spec(specs, 13, 1, 5, 4, 1), std::invalid_argument);
    EXPECT_THROW(pm::compute_seam_infos_from_spec(specs, 1, 4, 5, 4, 1), std::invalid_argument);
    std::vector<pm::SurgerySpec> empty_seam = {{0, 1, 2, 0}};
    EXPECT_THROW(pm::compute_seam_infos_from_spec(empty_seam, 13, 4, 4, 3, 1), std::invalid_argument);
}


// A held-set (partial) build must match a full build node for node on held nodes, leave the rest unconstructed,
// and keep every global quantity identical.
TEST(SubsetPlan, PartialMatchingGraphEqualsFullOnHeldNodes) {
    const int num_obs = 4;
    const std::vector<std::vector<pm::SurgerySpec>> layouts = {
        {{0, 1, 2, 3}},
        {{0, 1, 4, 5}, {1, 2, 4, 5}},
        {{0, 1, 4, 5}, {1, 2, 6, 5}, {2, 3, 9, 4}},
        {{1, 2, 0, 2}, {0, 3, 7, 4}},
    };
    for (const auto& [code, task] : kCodes)
    for (int M : {3, 5}) {
        ScopedObsConfig cfg(M);
        for (size_t li = 0; li < layouts.size(); ++li) {
            const auto& specs = layouts[li];
            auto base = pm::MultiObsDemGenerator::generate_base_dem(code, task, 3, 13, 0.01);
            pm::MultiObsDemGenerator generator(std::move(base), num_obs, specs, 0.02, 2);
            stim::DetectorErrorModel dem = generator.generate();
            pm::UserGraph user_graph = pm::detector_error_model_to_user_graph(dem, false, pm::NUM_DISTINCT_WEIGHTS);
#ifdef USE_SHMEM
            std::vector<pm::DetectorNodeEphemeralFields> eph(user_graph.nodes.size());
            pm::SharedMatchingGraph full = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, eph.data());
#else
            pm::SharedMatchingGraph full = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
#endif
            const int K_p = (int)full.p_per_obs_patch, K_vb = (int)full.vb_per_obs_patch;
            for (int n_pes : {2, 3, 4}) {
                for (int k : {1, 2}) {
                    int k_eff = pm::bound_effective_k(k, n_pes, full.num_partitions, K_p, true);
                    auto seams = pm::compute_seam_infos_from_spec(specs, (int)full.num_rounds, M, K_p, K_vb, k_eff);
                    const auto unresolved = seams;  // the held set covers each seam's full neighbour range
                    try { pm::resolve_seam_window_conflicts(seams, 0); } catch (const std::invalid_argument&) { continue; }
                    const int base_n = num_obs / n_pes, rem = num_obs % n_pes;
                    for (int pid = 0; pid < n_pes; ++pid) {
                        SCOPED_TRACE(code + " M=" + std::to_string(M) + " layout=" + std::to_string(li) + " n_pes=" +
                                     std::to_string(n_pes) + " k=" + std::to_string(k_eff) + " pid=" + std::to_string(pid));
                        int my_start = base_n * pid + std::min(pid, rem);
                        int my_count = base_n + (pid < rem ? 1 : 0);
                        pm::HeldSet held = pm::compute_held_set_obs(
                            full.num_partitions, full.num_virtual_boundaries, num_obs, K_p, K_vb, my_start, my_count, unresolved);
#ifdef USE_SHMEM
                        pm::SharedMatchingGraph part =
                            user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, eph.data(), &held);
#else
                        pm::SharedMatchingGraph part = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, &held);
#endif
                        pm::PartialGraphCheck check = pm::verify_partial_against_full(full, part, held, unresolved);
                        for (const auto& e : check.errors) ADD_FAILURE() << e;
                        EXPECT_EQ(check.num_errors, 0u);
                        EXPECT_LE(check.num_held, check.total_nodes);
                    }
                }
            }
        }
    }
}

// On a longer run each PE should hold far less than the whole graph.
TEST(SubsetPlan, PartialGraphIsSmallerThanFull) {
    ScopedObsConfig cfg(5);
    const int num_obs = 4;
    std::vector<pm::SurgerySpec> specs = {{0, 1, 30, 4}, {2, 3, 30, 4}};
    auto base = pm::MultiObsDemGenerator::generate_base_dem("surface_code", "rotated_memory_x", 3, 60, 0.01);
    pm::MultiObsDemGenerator generator(std::move(base), num_obs, specs, 0.02, 2);
    stim::DetectorErrorModel dem = generator.generate();
    pm::UserGraph user_graph = pm::detector_error_model_to_user_graph(dem, false, pm::NUM_DISTINCT_WEIGHTS);
    const size_t total = user_graph.nodes.size();
#ifdef USE_SHMEM
    std::vector<pm::DetectorNodeEphemeralFields> eph(total);
#endif
    const int K_p = (int)user_graph.p_per_obs_patch, K_vb = (int)user_graph.vb_per_obs_patch;
    auto seams = pm::compute_seam_infos_from_spec(specs, (int)user_graph.num_rounds, 5, K_p, K_vb, 1);
    for (int pid = 0; pid < 4; ++pid) {  // 4 PEs, one observable each
        pm::HeldSet held = pm::compute_held_set_obs(
            user_graph.num_partitions, user_graph.num_virtual_boundaries, num_obs, K_p, K_vb, pid, 1, seams);
#ifdef USE_SHMEM
        pm::SharedMatchingGraph part = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, eph.data(), &held);
#else
        pm::SharedMatchingGraph part = user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, &held);
#endif
        size_t num_held = 0;
        for (size_t i = 0; i < total; ++i) num_held += held.node_held(part.node_part_id[i]);
        EXPECT_LT(num_held, total * 6 / 10) << "pid " << pid;  // ~1 obs of 4 plus a window and seam nodes
    }
}

namespace {
size_t resident_bytes() {
    std::ifstream f("/proc/self/statm");
    size_t size_pages = 0, resident_pages = 0;
    f >> size_pages >> resident_pages;
    return resident_pages * 4096;
}
}  // namespace

// A lazy node array reserves address space for every node but only the constructed nodes become resident.
TEST(SubsetPlan, LazyNodeArrayOnlyResidentWhereConstructed) {
    const size_t n = 20'000'000;  // several GB of address space
    const size_t before = resident_bytes();
    pm::MatchingGraph graph(n, 1, pm::MatchingGraph::LazyNodes{});
    for (size_t i = 5'000'000; i < 5'010'000; ++i) graph.construct_node(i);
    graph.nodes[5'000'000].neighbors.push_back(&graph.nodes[5'000'001]);
    const size_t after = resident_bytes();
    EXPECT_EQ(graph.nodes.size(), n);
    EXPECT_EQ(&graph.nodes[12345] - graph.nodes.data(), 12345);  // contiguous, globally indexed
    EXPECT_EQ(graph.nodes[5'000'005].vb, -1);                     // constructed nodes have default fields
    EXPECT_LT(after - before, size_t(200) << 20) << "resident grew by " << (after - before) << " bytes";
}
