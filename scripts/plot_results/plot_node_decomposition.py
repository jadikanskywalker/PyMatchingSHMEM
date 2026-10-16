import matplotlib.pyplot as plt

from common import *
from plot_strong_scaling import COLORS, PRESETS

THREADS_PER_NODE = 256


def main():
    apply_style()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for p in PRESETS:
        rows = {r['ntasks']: r['us'] for r in load_obs(p)
                if r['nodes'] == 1 and r['ntasks'] * r['nthreads'] == THREADS_PER_NODE}
        xs = sorted(rows)
        kw = style_marker_kwargs(COLORS[p])
        a1.plot(xs, [rows[x] for x in xs], label=f'{p} obs', **kw)
        a2.plot(xs, [rows[1] / rows[x] for x in xs], **kw)
    a1.set_yscale('log')
    a1.yaxis.set_major_formatter(plt.ScalarFormatter())
    a1.set_ylabel('Decoding time per round (µs)')
    a1.set_title('Time per round')
    a1.legend(fontsize=8.5, ncol=2)
    a2.set_ylabel('Speedup vs. 1 task per node')
    a2.set_title('Speedup from task decomposition')
    for ax in (a1, a2):
        ax.set_xscale('log', base=2)
        ax.set_xticks([1, 2, 4, 8, 16])
        ax.set_xticklabels(['1', '2', '4', '8', '16'])
        ax.minorticks_off()
        ax.set_xlabel(f'Tasks per node ({THREADS_PER_NODE} threads/node total)')
    fig.suptitle('Single-node task decomposition (d = 21, p = 0.001)', x=0.01, ha='left',
                 fontweight='bold')
    fig.tight_layout()
    save(fig, 'node_decomposition')


if __name__ == '__main__':
    main()
