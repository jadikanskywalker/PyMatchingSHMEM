#pragma once

#define DEBUG 1

#ifdef USE_SHMEM
namespace config_shmem {
    inline int M = 1;
}
#endif