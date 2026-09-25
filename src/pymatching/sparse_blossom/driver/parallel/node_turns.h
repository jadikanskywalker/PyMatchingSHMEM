#ifndef PYMATCHING2_NODE_TURNS_H
#define PYMATCHING2_NODE_TURNS_H

#include <functional>

#ifdef USE_SHMEM
#include <shmem.h>
// Sandia OpenSHMEM provides SHMEMX_TEAM_NODE (built from the runtime's node ranks); other implementations fall back
// to the standard SHMEM_TEAM_SHARED.
#if defined(__has_include) && __has_include(<shmemx-def.h>)
#include <shmemx.h>
#define PM_NODE_TEAM SHMEMX_TEAM_NODE
#else
#define PM_NODE_TEAM SHMEM_TEAM_SHARED
#endif
#endif

namespace pm {

#ifdef USE_SHMEM
/// The PEs sharing this PE's node, used to make them take turns on memory-heavy setup work.
struct NodeTurns {
    shmem_team_t team;
    int local_rank;  // this PE's index among the PEs on its node
    int local_size;  // number of PEs on this node
};

/// Builds the node-local team (SHMEMX_TEAM_NODE) and checks its size against the launcher's per-node PE count when
/// the launcher reports one (a team silently degenerating to size 1 would defeat the mode). A mismatch prints an
/// error and terminates all PEs.
NodeTurns make_node_turns();

/// Runs `work` on each PE of the node in turn (local rank 0 first). Contains no collective operation: PEs on other
/// nodes proceed independently, and the caller synchronises everyone afterwards. If `work` throws, prints the error
/// and terminates all PEs (peers waiting in the node barrier would otherwise hang).
void run_in_turns(const NodeTurns& turns, const std::function<void()>& work);
#endif

}  // namespace pm

#endif  // PYMATCHING2_NODE_TURNS_H
