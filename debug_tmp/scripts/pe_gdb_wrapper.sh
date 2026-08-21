#!/bin/bash
# Same per-rank output redirection as pe_output_wrapper.sh, but runs the real binary under gdb
# so a crash (SIGABRT from glibc's double-free detector, SIGSEGV, etc.) dumps a full backtrace
# for every thread before the process dies, instead of just glibc's own one-line message.
RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
exec gdb -q -batch \
    -ex "handle SIGABRT stop print nopass" \
    -ex "run" \
    -ex "thread apply all bt full" \
    -ex "quit" \
    --args "$@" > "log_shmem_pe${RANK}.out" 2>&1
