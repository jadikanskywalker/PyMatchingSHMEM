#ifndef PYMATCHING2_HELD_SET_H
#define PYMATCHING2_HELD_SET_H

#include <vector>

namespace pm {

/// Which partitions / virtual boundaries of the global graph a PE needs fully built (nodes with complete
/// neighbour lists). node_part_id[n] >= 0 is partition id, < 0 is virtual boundary -(id+1).
struct HeldSet {
    std::vector<bool> part_held;  // size num_partitions
    std::vector<bool> vb_held;    // size num_virtual_boundaries
    /// Debug: also construct every node that is neither held nor a stub as a "poisoned" node (vb = INT_MIN and a
    /// bogus neighbour pointer), so any traversal of a node outside the held set crashes instead of silently
    /// reading an empty node.
    bool poison{false};

    /// True if a node with this node_part_id entry is fully built on this PE.
    inline bool node_held(int part_id) const {
        return part_id >= 0 ? part_held[part_id] : vb_held[-part_id - 1];
    }
};

}  // namespace pm

#endif  // PYMATCHING2_HELD_SET_H
