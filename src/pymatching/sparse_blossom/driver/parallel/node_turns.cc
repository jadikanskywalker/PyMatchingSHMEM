#include "pymatching/sparse_blossom/driver/parallel/node_turns.h"

#ifdef USE_SHMEM
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

pm::NodeTurns pm::make_node_turns() {
    NodeTurns t;
    t.team = PM_NODE_TEAM;
    t.local_rank = shmem_team_my_pe(t.team);
    t.local_size = shmem_team_n_pes(t.team);
    // Only the launcher's own per-node count is meaningful: SLURM_NTASKS_PER_NODE describes the batch allocation, not
    // how oshrun placed the PEs. A mismatch is fatal for the whole job (a PE that throws alone would leave its peers
    // waiting forever).
    const char* v = std::getenv("OMPI_COMM_WORLD_LOCAL_SIZE");
    if (v && std::atoi(v) > 0 && std::atoi(v) != t.local_size) {
        std::cerr << "PE" << shmem_my_pe() << ": node-local team has " << t.local_size
                  << " PEs but OMPI_COMM_WORLD_LOCAL_SIZE=" << v << "; limited_memory_mode needs a correct node team."
                  << std::endl << std::flush;
        shmem_global_exit(1);
    }
    return t;
}

namespace {
// Memory in use on this node (MemTotal - MemAvailable), in GB.
double node_used_gb() {
    std::ifstream f("/proc/meminfo");
    std::string key;
    long total = 0, avail = 0, value;
    std::string unit;
    while (f >> key >> value >> unit) {
        if (key == "MemTotal:") total = value;
        if (key == "MemAvailable:") avail = value;
    }
    return (total - avail) / (1024.0 * 1024.0);
}
}  // namespace

void pm::run_in_turns(const NodeTurns& turns, const std::function<void()>& work) {
    for (int turn = 0; turn < turns.local_size; ++turn) {
        if (turns.local_rank == turn) {
            try {
                const auto t0 = std::chrono::steady_clock::now();
                work();
                const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                std::cout << "PE" << shmem_my_pe() << " turn " << turn << "/" << turns.local_size << " done in " << secs
                          << "s; node memory in use " << node_used_gb() << " GB" << std::endl << std::flush;
            } catch (const std::exception& e) {
                std::cerr << "PE" << shmem_my_pe() << " failed in its turn: " << e.what() << std::endl << std::flush;
                shmem_global_exit(1);
            }
        }
        shmem_team_sync(turns.team);
    }
}
#endif
