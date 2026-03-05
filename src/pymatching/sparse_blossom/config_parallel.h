#pragma once

#define DEBUG 0

#define ENABLE_DRAW_FLAGS 0

#define NUM_BUFFERS_PER_UNIT 1

// per solver region buffer size = num_nodes_per_partition * SHMEM_ARENA_BUFFER_FACTOR
//    1 ensures no heap overflow, can reduce for small enough p
#define SHMEM_ARENA_BUFFER_FACTOR 1 
#define SHMEM_INTERSECTION_BUFFER_FACTOR 2 // multiplied by d, maximum number of rounds that cross-PE fusion can last 
// #define SHMEM_MAX_D 30 // used for regions_matched_to_vb buffer in shmem

#ifdef USE_THREADS
namespace config_parallel {
    inline int M = 0;
}
#endif