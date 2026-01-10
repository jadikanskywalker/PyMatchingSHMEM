#pragma once

#define DEBUG 0

#define NUM_ACTIVE_SHOTS_PER_UNIT 4

#ifdef USE_THREADS
namespace config_parallel {
    inline int M = 0;
}
#endif