#pragma once

#include <limits>

#define DEBUG 1 // set for full debugging
#if !DEBUG
#define BARE_DEBUG 0 // set for per shot prints when DEBUG == 0
#else
#define BARE_DEBUG 1
#endif

#define ENABLE_DRAW_FLAGS

#define NUM_BUFFERS_PER_UNIT 1

// #define PROFILE_OMP_BARRIERS

#ifdef USE_SHMEM
// per solver region buffer size in SHMEM = num_nodes_per_partition * SHMEM_ARENA_BUFFER_FACTOR rounded up to a multiple of 64
//    2 would ensure no heap overflow in absolute worst case, can reduce for smaller p & _should_ not segfault
#define SHMEM_ARENA_BUFFER_FACTOR 1 

// region_matched_to_vb buffer size in SHMEM = num_nodes_per_round (eg., d for repetition_code, d^2 for surface_code) * SHMEM_INTERSECTION_BUFFER_FACTOR
#define SHMEM_INTERSECTION_BUFFER_FACTOR 2

// #define SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER 2
// #define SHMEM_USE_FLATTEN_FUSIONS true // not flattening is currently unsupported
// #define SHMEM_MAX_D 30 // used for regions_matched_to_vb buffer in shmem
#endif

#ifdef USE_THREADS
namespace config_parallel {
    inline int M = std::numeric_limits<int>::max();
    inline bool obs_coors_included = false;
    // Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md). L is always
    // active and is the number of partitions per extraction job/unit; must be >= 1 (no power-of-2
    // constraint -- the per-unit balanced-subtree/post-hoc chunking logic handles ragged unit sizes).
    inline int L = 1;
    // Selects tree shape / job-posting timing, decoupled from L's value: false (default) builds the
    // ordinary balanced fusion tree and defers all extraction-job posting until the whole root is
    // solved, then recursively chunks it into <=L-leaf jobs; true builds the chain-of-unit-subtrees
    // and posts one job per checkpoint incrementally, as each one resolves during decode.
    inline bool extract_preemptively = false;
#ifdef USE_SHMEM
    enum div_strgy { ROUND, OBS };
    inline int k = 1;
    inline div_strgy division_strategy = ROUND;
#endif
}
#endif