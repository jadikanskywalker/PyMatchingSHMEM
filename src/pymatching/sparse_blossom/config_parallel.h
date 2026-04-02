#pragma once

#define DEBUG 0 // 0 = none; 1 = full;
#define BARE_DEBUG 0

#define ENABLE_DRAW_FLAGS

#define NUM_BUFFERS_PER_UNIT 1

// per solver region buffer size in SHMEM = num_nodes_per_partition * SHMEM_ARENA_BUFFER_FACTOR rounded up to a multiple of 64
//    2 would ensure no heap overflow in absolute worst case, can reduce for smaller p & _should_ not segfault
#define SHMEM_ARENA_BUFFER_FACTOR 1 

// region_matched_to_vb buffer size in SHMEM = num_nodes_per_round (eg., d for repetition_code, d^2 for surface_code) * SHMEM_INTERSECTION_BUFFER_FACTOR
#define SHMEM_INTERSECTION_BUFFER_FACTOR 2

#define SHMEM_NUM_CROSS_RANK_FUSIONS_PER_BUFFER 2
// #define SHMEM_USE_FLATTEN_FUSIONS true // not flattening is currently unsupported
// #define SHMEM_MAX_D 30 // used for regions_matched_to_vb buffer in shmem

#ifdef USE_THREADS
namespace config_parallel {
    inline int M = 0;
#ifdef USE_SHMEM
    inline int k = 1;
#endif
}
#endif