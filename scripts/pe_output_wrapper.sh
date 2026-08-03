#!/bin/bash
# Redirects each PE's own stdout+stderr to log_shmem_pe<id>.out, so per-PE output doesn't get
# interleaved into one shared log. Same rank-detection convention as scorep_wrapper.sh.
RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
exec "$@" > "log_shmem_pe${RANK}.out" 2>&1
