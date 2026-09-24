#ifndef PYMATCHING2_HELD_SET_H
#define PYMATCHING2_HELD_SET_H

#include <vector>

namespace pm {

/// Which partitions / virtual boundaries of the global graph a PE needs fully built (nodes with complete
/// neighbour lists). node_part_id[n] >= 0 is partition id, < 0 is virtual boundary -(id+1).
struct HeldSet {
    std::vector<bool> part_held;  // size num_partitions
    std::vector<bool> vb_held;    // size num_virtual_boundaries

    /// True if a node with this node_part_id entry is fully built on this PE.
    inline bool node_held(int part_id) const {
        return part_id >= 0 ? part_held[part_id] : vb_held[-part_id - 1];
    }
};

}  // namespace pm

#endif  // PYMATCHING2_HELD_SET_H
