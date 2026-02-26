#pragma once

#define DEBUG 1

#define ENABLE_DRAW_FLAGS 1

#define NUM_BUFFERS_PER_UNIT 1

#define SHMEM_ARENA_BUFFER_FACTOR 0.5 // buffer_size = num_nodes_per_partition * SHMEM_ARENA_BUFFER_FACTOR
#define SHMEM_INTERSECTION_BUFFER_FACTOR 2 // multiplied by d, maximum number of rounds that cross-PE fusion can last 
#define SHMEM_MAX_D 30 // used to regions_matched_to_vb buffer in shmem

#ifdef USE_THREADS
namespace config_parallel {
    inline int M = 0;
}
#endif