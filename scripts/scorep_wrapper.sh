#!/bin/bash
# Sets a per-PE Score-P experiment directory using the OMPI rank so all PEs
# don't race to create the same directory (--mpp=none treats each PE as serial).
RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
export SCOREP_EXPERIMENT_DIRECTORY=scorep_results_pe${RANK}
exec "$@"
