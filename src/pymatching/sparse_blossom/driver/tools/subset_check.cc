// Checks the partial (held-set) graph build against the full build on a real graph cache:
//   verify  : spec-derived seam windows vs graph-walk windows, and per-PE partial graph vs full graph
//   measure : resident memory of one partial (or the full) graph build
// Usage:
//   subset_check verify  --graph_cache_path P --rounds_per_partition M --preset 36obs [--pes 2,4,8] [--k 1] [--search]
//   subset_check measure --graph_cache_path P --rounds_per_partition M --preset 36obs --pe N:PID | --pe full
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "pymatching/sparse_blossom/config_parallel.h"
#include "pymatching/sparse_blossom/driver/helpers/dem_generator.h"
#include "pymatching/sparse_blossom/driver/helpers/graph_cache.h"
#include "pymatching/sparse_blossom/driver/parallel/subset_plan.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"

namespace {

const char* arg(int argc, char** argv, const char* name) {
    for (int i = 1; i + 1 < argc; ++i)
        if (!strcmp(argv[i], name)) return argv[i + 1];
    return nullptr;
}

bool has_flag(int argc, char** argv, const char* name) {
    for (int i = 1; i < argc; ++i)
        if (!strcmp(argv[i], name)) return true;
    return false;
}

double status_mb(const char* key) {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind(key, 0) == 0) {
            std::istringstream ss(line.substr(strlen(key)));
            double kb;
            ss >> kb;
            return kb / 1024.0;
        }
    }
    return -1;
}

std::vector<pm::SurgerySpec> preset(const std::string& name) {
    using G = pm::MultiObsDemGenerator;
    if (name == "18obs") return G::preset_18obs();
    if (name == "24obs") return G::preset_24obs();
    if (name == "36obs") return G::preset_36obs();
    if (name == "48obs") return G::preset_48obs();
    if (name == "64obs") return G::preset_64obs();
    if (name == "72obs") return G::preset_72obs();
    if (name == "144obs") return G::preset_144obs();
    if (name == "128obs") return G::preset_128obs();
    if (name == "256obs") return G::preset_256obs();
    if (name == "36obsx4") return G::preset_36obs_repeated(4);
    throw std::invalid_argument("unknown preset " + name);
}

std::vector<int> parse_list(const std::string& s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) out.push_back(std::stoi(tok));
    return out;
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: subset_check verify|measure ...\n";
        return 2;
    }
    std::string mode = argv[1];
    const char* cache = arg(argc, argv, "--graph_cache_path");
    const char* m_arg = arg(argc, argv, "--rounds_per_partition");
    const char* preset_arg = arg(argc, argv, "--preset");
    if (!cache || !m_arg || !preset_arg) {
        std::cerr << "need --graph_cache_path, --rounds_per_partition, --preset\n";
        return 2;
    }
    const int M = std::stoi(m_arg);
    config_parallel::M = M;
    config_parallel::division_strategy = config_parallel::OBS;
    auto specs = preset(preset_arg);

    auto t0 = std::chrono::steady_clock::now();
    const double rss_start = status_mb("VmRSS:");
    pm::UserGraph user_graph = pm::read_user_graph_cache(cache, false, M, (int)config_parallel::division_strategy);
    const double rss_graph = status_mb("VmRSS:");
    const int P = (int)user_graph.num_obs_patches, K_p = (int)user_graph.p_per_obs_patch,
              K_vb = (int)user_graph.vb_per_obs_patch, num_seams = (int)user_graph.num_virtual_boundaries - K_vb * P;
    std::cout << "loaded in " << seconds_since(t0) << "s: nodes=" << user_graph.nodes.size() << " obs=" << P
              << " K_p=" << K_p << " K_vb=" << K_vb << " seams(graph)=" << num_seams
              << " seams(spec)=" << specs.size() << " rounds=" << user_graph.num_rounds
              << "  UserGraph RSS=" << (rss_graph - rss_start) << " MB" << std::endl;
    if ((int)specs.size() != num_seams) {
        std::cout << "FAIL: spec seam count != graph seam count\n";
        return 1;
    }
#ifdef USE_SHMEM
    std::vector<pm::DetectorNodeEphemeralFields> eph(user_graph.nodes.size());
#endif
    auto build = [&](const pm::HeldSet* held) {
#ifdef USE_SHMEM
        return user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, eph.data(), held);
#else
        return user_graph.to_shared_matching_graph(pm::NUM_DISTINCT_WEIGHTS, held);
#endif
    };
    auto held_for = [&](int n_pes, int pid, int k, std::vector<pm::SeamInfo>* unresolved_out) {
        int k_eff = pm::bound_effective_k(k, n_pes, user_graph.num_partitions, K_p, true);
        auto seams = pm::compute_seam_infos_from_spec(specs, (int)user_graph.num_rounds, M, K_p, K_vb, k_eff);
        if (unresolved_out) *unresolved_out = seams;
        const int base = P / n_pes, rem = P % n_pes;
        int my_start = base * pid + std::min(pid, rem), my_count = base + (pid < rem ? 1 : 0);
        return pm::compute_held_set_obs(
            user_graph.num_partitions, user_graph.num_virtual_boundaries, P, K_p, K_vb, my_start, my_count, seams);
    };

    if (mode == "measure") {
        const char* pe = arg(argc, argv, "--pe");
        const int k = arg(argc, argv, "--k") ? std::stoi(arg(argc, argv, "--k")) : 1;
        if (!pe) return 2;
        const double before = status_mb("VmRSS:");
        auto tb = std::chrono::steady_clock::now();
        pm::HeldSet held;
        bool full = std::string(pe) == "full";
        if (!full) {
            int n_pes = std::stoi(pe), pid = std::stoi(strchr(pe, ':') + 1);
            held = held_for(n_pes, pid, k, nullptr);
        }
        pm::SharedMatchingGraph g = build(full ? nullptr : &held);
        const double after = status_mb("VmRSS:");
        size_t built = 0;
        if (!full)
            for (size_t i = 0; i < user_graph.nodes.size(); ++i) built += held.node_held(g.node_part_id[i]);
        std::cout << "measure pe=" << pe << " build " << seconds_since(tb) << "s  graph RSS delta=" << (after - before)
                  << " MB  (UserGraph " << (rss_graph - rss_start) << " MB)  held nodes="
                  << (full ? user_graph.nodes.size() : built) << "/" << user_graph.nodes.size()
                  << "  VmHWM=" << status_mb("VmHWM:") << " MB" << std::endl;
        return 0;
    }

    // ---- verify ----
    int failures = 0;
    std::vector<int> ks = {1, 2, 3};
    const auto tfull = std::chrono::steady_clock::now();
    pm::SharedMatchingGraph full = build(nullptr);
    std::cout << "full graph built in " << seconds_since(tfull) << "s  RSS=" << status_mb("VmRSS:") << " MB" << std::endl;
    for (int k : ks) {
        auto from_graph = pm::compute_seam_infos_from_graph(full, k);
        auto from_spec = pm::compute_seam_infos_from_spec(specs, (int)full.num_rounds, M, K_p, K_vb, k);
        int bad = 0;
        for (int s = 0; s < num_seams; ++s) {
            const auto &a = from_graph[s], &b = from_spec[s];
            if (a.oi != b.oi || a.oj != b.oj || a.vb_left != b.vb_left || a.vb_right != b.vb_right) {
                if (bad++ < 10)
                    std::cout << "  seam " << s << " graph(oi=" << a.oi << ",oj=" << a.oj << ",l=" << a.vb_left
                              << ",r=" << a.vb_right << ") spec(oi=" << b.oi << ",oj=" << b.oj << ",l=" << b.vb_left
                              << ",r=" << b.vb_right << ")\n";
            }
        }
        std::cout << "seam infos k=" << k << ": " << (num_seams - bad) << "/" << num_seams << " match" << std::endl;
        failures += bad;
    }
    const int k = arg(argc, argv, "--k") ? std::stoi(arg(argc, argv, "--k")) : 1;
    const bool with_search = has_flag(argc, argv, "--search");
    std::shared_ptr<pm::SearchGraph> search_full;
    if (with_search) {
        auto ts = std::chrono::steady_clock::now();
        search_full = user_graph.to_shared_search_graph(pm::NUM_DISTINCT_WEIGHTS);
        std::cout << "full search graph built in " << seconds_since(ts) << "s  RSS=" << status_mb("VmRSS:") << " MB"
                  << std::endl;
    }
    std::vector<int> pes = parse_list(arg(argc, argv, "--pes") ? arg(argc, argv, "--pes") : "2,4,8");
    for (int n_pes : pes) {
        for (int pid = 0; pid < n_pes; ++pid) {
            std::vector<pm::SeamInfo> unresolved;
            pm::HeldSet held;
            try {
                held = held_for(n_pes, pid, k, &unresolved);
            } catch (const std::exception& e) {
                std::cout << "n_pes=" << n_pes << " skipped: " << e.what() << std::endl;
                break;
            }
            auto tp = std::chrono::steady_clock::now();
            pm::SharedMatchingGraph part = build(&held);
            auto check = pm::verify_partial_against_full(full, part, held, unresolved);
            std::cout << "n_pes=" << n_pes << " pid=" << pid << ": held " << check.num_held << "/" << check.total_nodes
                      << " (" << (100.0 * check.num_held / check.total_nodes) << "%) stubs=" << check.num_stubs
                      << " errors=" << check.num_errors << " (" << seconds_since(tp) << "s)" << std::endl;
            for (const auto& e : check.errors) std::cout << "    " << e << std::endl;
            failures += (int)check.num_errors;
            if (with_search) {
                auto spart = user_graph.to_shared_search_graph(pm::NUM_DISTINCT_WEIGHTS, &held);
                auto scheck = pm::verify_partial_search_against_full(
                    *search_full, *spart, full.node_part_id, held, K_vb, P, unresolved);
                std::cout << "   search graph: held " << scheck.num_held << " stubs=" << scheck.num_stubs
                          << " errors=" << scheck.num_errors << std::endl;
                for (const auto& e : scheck.errors) std::cout << "    " << e << std::endl;
                failures += (int)scheck.num_errors;
            }
        }
    }
    std::cout << (failures ? "FAILED" : "ALL OK") << " (" << failures << " mismatches)" << std::endl;
    return failures ? 1 : 0;
}
