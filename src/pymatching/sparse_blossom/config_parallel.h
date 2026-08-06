#pragma once

#include <limits>

#define DEBUG 1 // set for full debugging
#if !DEBUG
#define BARE_DEBUG 0 // set for per shot prints when DEBUG == 0
#else
#define BARE_DEBUG 1
#endif

#define ENABLE_DRAW_FLAGS

// #define ENABLE_SHOT_BUFFERS
#ifdef ENABLE_SHOT_BUFFERS
#define NUM_BUFFERS_PER_UNIT 2
#endif

// #define PROFILE_OMP_BARRIERS

#ifdef USE_SHMEM
// per solver region buffer size in SHMEM = num_nodes_per_partition * SHMEM_ARENA_BUFFER_FACTOR rounded up to a multiple of 64.
// 2 would ensure no heap overflow in absolute worst case, can reduce for smaller p
#define SHMEM_ARENA_BUFFER_FACTOR 1 
// region_matched_to_vb buffer size in SHMEM = num_nodes_per_round (eg., d for repetition_code, d^2 for surface_code) * SHMEM_INTERSECTION_BUFFER_FACTOR
#define SHMEM_INTERSECTION_BUFFER_FACTOR 2
#endif

namespace config_parallel {
    inline bool use_threads = true;
    
    inline int M = std::numeric_limits<int>::max(); // rounds per partition
    
    enum div_strgy { ROUND, OBS };
    inline int seam_buffer_size = 1;
    inline int k = 1;
    inline div_strgy division_strategy = ROUND;

    // Unit-checkpointed extraction (see plans/this-is-a-broader-purrfect-crystal.md). L is always
    // active and is the number of partitions per extraction job/unit.
    inline int L = 1;
    // Selects tree shape / job-posting timing, decoupled from L's value: false (default) builds the
    // ordinary balanced fusion tree and defers all extraction-job posting until the whole root is
    // solved, then recursively chunks it into <=L-leaf jobs; true builds the chain-of-unit-subtrees
    // and posts one job per checkpoint incrementally, as each one resolves during decode.
    inline bool extract_preemptively = false;
}