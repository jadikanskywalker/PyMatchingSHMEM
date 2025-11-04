#!/usr/bin/env bash
# filepath: scripts/check_virtual_neighbors.sh
set -euo pipefail
INPUT=${1:-graph.out}

gawk '
function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }

{
    line = trim($0)
    if (line == "") next

    if (line ~ /^node:/) {
        curr_node = trim(substr(line, 6))
        curr_neighbor = ""
    } else if (line ~ /^neighbor:/) {
        split(line, parts, /[[:space:]]*:[[:space:]]*/)
        curr_neighbor = trim(parts[2])
    } else if (line ~ /^partition[[:space:]]*:/) {
        split(line, parts, /[[:space:]]*:[[:space:]]*/)
        val = trim(parts[2]) + 0
        if (curr_neighbor == "")
            node_part[curr_node] = val
        else
            neighbor_part[curr_node, curr_neighbor] = val
    } else if (line ~ /^is_virtual[[:space:]]*:/) {
        split(line, parts, /[[:space:]]*:[[:space:]]*/)
        val = trim(parts[2]) + 0
        if (curr_neighbor == "")
            node_virtual[curr_node] = val
        else
            neighbor_virtual[curr_node, curr_neighbor] = val
    }
}

END {
    violations = 0
    for (node in node_part) {
        node_p = node_part[node]
        node_v = node_virtual[node]

        for (key in neighbor_part) {
            split(key, idx, SUBSEP)
            if (idx[1] != node) continue

            nb = idx[2]
            if (nb == "0") continue

            nb_p = neighbor_part[key]

            if (nb_p == node_p - 1 && node_v != 1) {
                printf("Violation: node %s (partition %d, is_virtual=%d) has neighbour %s in partition %d\n",
                       node, node_p, node_v, nb, nb_p)
                violations++
            }

            if (nb_p == node_p + 1) {
                if (nb in node_virtual) {
                    if (node_virtual[nb] != 1) {
                        printf("Violation: neighbour %s (partition %d, is_virtual=%d) of node %s (partition %d) should be virtual\n",
                               nb, nb_p, node_virtual[nb], node, node_p)
                        violations++
                    }
                } else {
                    printf("Warning: neighbour %s referenced by node %s missing from dump\n", nb, node) > "/dev/stderr"
                }
            }
        }
    }

    if (violations == 0) {
        print "OK: no violations detected"
    } else {
        printf("%d violations found\n", violations)
        exit 1
    }
}
' "$INPUT"