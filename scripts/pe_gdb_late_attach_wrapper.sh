#!/bin/bash
# Runs the real binary natively (no ptrace overhead at all, avoiding the PRTE-liveness-kill seen
# when wrapping the whole process lifetime in gdb) and only attaches gdb after a delay timed to
# land shortly before the crash window -- a brief one-time ptrace-attach pause instead of
# sustained overhead across the whole run. Per-shot progress ("Starting shot N") only goes to the
# per-thread out_parallel/pN/tM.out trace files (64 of them per PE), not this main log, so a fixed
# delay is simpler/more robust than scanning all of them. Original native crash (job 144389) took
# ~17-18 min total; attach at 12 min to be safely ahead of it with margin to spare.
RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
LOG="log_shmem_pe${RANK}.out"
ATTACH_DELAY_SECONDS=720

"$@" > "$LOG" 2>&1 &
PID=$!

sleep "$ATTACH_DELAY_SECONDS"

if kill -0 "$PID" 2>/dev/null; then
    gdb -q -batch \
        -ex "attach $PID" \
        -ex "handle SIGABRT stop print nopass" \
        -ex "continue" \
        -ex "thread apply all bt full" \
        -ex "quit" \
        > "log_shmem_pe${RANK}_gdb.out" 2>&1
fi

wait "$PID"
