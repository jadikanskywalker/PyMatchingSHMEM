# PyMatchingSHMEM

PyMatchingSHMEM parallelizes PyMatching v2 using OpenSHMEM. PyMatching v2 is an efficient single core
implementation of the Sparse Blossom algorithm created by Oscar Higgot and Craig Kidney. We parallelize
the program using the same divide-and-conquer strategy employed by Lu Wei in fusion blossom. The 
decoding graph is split into partitions that are solved independently and then fused. However, while 
fusion blossom uses a single multi-threaded process, PyMatchingSHMEM uses OpenSHMEM, a SPMD
specification with rich functionality for atomic and one-sided Remote Memory Access (RMA) operations on
shared memory. Utilizing shared memory provides unique opportunities for optimization...

## Running
...

## Partitioning
This implementation uses detector coordinates from the Stim DetectorErrorModel to partition the
circuit. The case considered here is when N repeated measurement rounds are partioned into 3D
partition of M"<"N rounds each. Each set of M rounds is solved independently and then fused with
neighboring sets.
  - IMPORTANT NOTE: It is assumed that M is sufficiently large that each set of M rounds only
    has edges to the next set of M rounds. Partition i should not connect to partition i+2. (This
    is generally the case, but better to state it).
Additional partition strategies could be implemented to optimize for different cases.

# Work Scheduling


## Notes on Locality
[TENTATIVE] Locality is crucial for optimal performance. This implementation checks enviornment 
variables like SLURM_JOB_NUM_NODES and SLURM_NTASKS_PER_NODE and assumes that processes have a
contiguous layout. This allows the implementation to assign work to PEs more intelligently. For 
example, if the runtime is scaled to two nodes with 16 processes per node, processes 0...15 should
likely be on one node and 16...32 on the other. Half of the graph can be solved on PE's 0...15, and half
on PE's 16...32. This keeps fusions on the same node until the final fusion between the two halves,
minimizing inter-node communication.