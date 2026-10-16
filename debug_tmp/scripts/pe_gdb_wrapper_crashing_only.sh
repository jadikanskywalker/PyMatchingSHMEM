#!/bin/bash
# Same idea as pe_gdb_wrapper.sh, but captures ONLY the crashing thread's own full backtrace
# first (small, reliable) before the full "thread apply all" dump, which can truncate long
# before reaching the actually-interesting thread on a 128-thread run.
RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
exec gdb -q -batch \
    -ex "handle SIGABRT stop print nopass" \
    -ex "run" \
    -ex "echo \n=== CRASHING THREAD ONLY ===\n" \
    -ex "bt full" \
    -ex "echo \n=== ALL THREADS (may truncate) ===\n" \
    -ex "thread apply all bt" \
    -ex "quit" \
    --args "$@" > "log_shmem_pe${RANK}.out" 2>&1
