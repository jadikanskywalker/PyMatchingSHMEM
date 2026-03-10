#pragma once

#define DEBUG 1

#define ENABLE_DRAW_FLAGS

#define NUM_BUFFERS_PER_UNIT 1

// per solver region buffer size = num_nodes_per_partition * SHMEM_ARENA_BUFFER_FACTOR
//    1 ensures no heap overflow, can reduce for smaller p, currently could segfault
#define SHMEM_ARENA_BUFFER_FACTOR 1 
#define SHMEM_INTERSECTION_BUFFER_FACTOR 2 // multiplied by d, maximum number of rounds that cross-PE fusion can last

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