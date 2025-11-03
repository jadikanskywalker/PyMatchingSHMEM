#pragma once

#define DEBUG 1

#ifdef ENABLE_FUSION
namespace config_parallel {
    inline int M = 1;
}
#endif