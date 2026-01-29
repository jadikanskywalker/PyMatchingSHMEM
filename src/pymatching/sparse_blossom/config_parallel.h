#pragma once

#define DEBUG 0

#define NUM_BUFFERS_PER_UNIT 2

#define SHMEM_ARENA_BUFFER_FACTOR 0.5 // buffer_size = num_nodes * SHMEM_ARENA_BUFFER_FACTOR

#ifdef USE_THREADS
namespace config_parallel {
    inline int M = 0;
}
#endif