import matplotlib.pyplot as plt

from common_v2 import *

THREADS_PER_NODE = 256
TASKS = [1, 2, 4, 8, 16]


def main():
    apply_style_v2()
    pts = Points()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.3))
    for p in PRESETS:
        rows = {r['ntasks']: r['us'] for r in load_obs(p)
                if r['nodes'] == 1 and r['ntasks'] * r['nthreads'] == THREADS_PER_NODE}
        xs = sorted(rows)
        a1.plot(xs, [rows[x] for x in xs], label=f'{p} obs', **mk(PRESET_COLORS[p]))
        base = t_sp(p)[0]
        a2.plot(xs, [base / rows[x] for x in xs], **mk(PRESET_COLORS[p]))
        for x in xs:
            pts.add('1-node, 256 thr/node', preset=p, ntasks=x, nthreads=THREADS_PER_NODE // x,
                    us=rows[x], T_sp_us=base, speedup_vs_thread_only=base / rows[x])
    log_y(a1, [2, 5, 10, 20, 50])
    a1.set_ylabel('Time per round (µs)')
    a1.legend(ncol=2, loc='upper right', columnspacing=1.0, handlelength=1.6)
    a2.axhline(1, color=MUTED, linewidth=0.8)
    a2.set_ylabel('Speedup vs. best thread-only run')
    for ax in (a1, a2):
        ax.set_xscale('log', base=2)
        ax.set_xticks(TASKS)
        ax.set_xticklabels([f'{t}\n({THREADS_PER_NODE // t})' for t in TASKS])
        ax.minorticks_off()
        ax.set_xlim(0.8, 20)
        ax.set_xlabel('Tasks per node\n(threads per task)')
        for x, lab in ((2, '2 sockets'), (8, '8 NUMA domains')):
            ax.axvline(x, color=MUTED, linestyle=':', linewidth=0.9, zorder=0)
            ax.text(x, 1.01, lab, ha='center', va='bottom', color=MUTED,
                    transform=ax.get_xaxis_transform())
    fig.suptitle('Single-node task decomposition at 256 threads/node (d = 21, p = 0.001)',
                 x=0.01, ha='left', fontweight='bold')
    fig.tight_layout()
    save_v2(fig, 'node_decomposition', pts)


if __name__ == '__main__':
    main()
