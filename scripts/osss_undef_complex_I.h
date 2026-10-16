// OSSS header workarounds for C++ (forced into every TU by scripts/build_osss.sh).
// 1) shmem headers include glibc <complex.h>, whose `I` macro breaks stim's GateType::I.
// 2) shmem/api.h includes teams.h outside its extern "C" block, giving team calls C++ linkage.
#ifdef __cplusplus
#include <complex.h>
#include <shmem/defs.h>
extern "C" {
#include <shmem/teams.h>
}
#undef I
#endif
