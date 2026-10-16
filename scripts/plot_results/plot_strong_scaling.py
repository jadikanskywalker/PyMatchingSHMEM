import matplotlib.pyplot as plt

from common import *

PRESETS = [36, 64, 72, 128, 144, 256]
COLORS = dict(zip(PRESETS, [BLUE, ORANGE, AQUA, '#eda100', '#e87ba4', '#4a3aa7']))


def ideal_line(ax, nodes, ys, **kw):
    ax.plot(nodes, ys, color=MUTED, linestyle='--', linewidth=1.2, label='Ideal', **kw)


def main():
    apply_style()
    summary = []
    for p in PRESETS:
        best = best_per_node(p)
        nodes = list(best)
        t = [best[n][0] for n in nodes]
        n0, t0 = nodes[0], t[0]
        sp = best_single_process(p)
        c = COLORS[p]
        fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.2))

        a1.plot(nodes, t, label='Best config', **style_marker_kwargs(c))
        ideal_line(a1, nodes, [t0 * n0 / n for n in nodes])
        a1.set_ylim(0, None)
        a1.set_ylabel('Decoding time per round (µs)')
        a1.set_title('Time per round')

        s_sp, s_n0 = [sp / x for x in t], [t0 / x for x in t]
        a2.plot(nodes, s_sp, label='vs. best single process', **style_marker_kwargs(c))
        a2.plot(nodes, s_n0, color=c, marker='s', markerfacecolor=SURFACE, linestyle=':',
                label=f'vs. best {n0}-node run')
        ideal_line(a2, nodes, nodes)
        a2.set_ylabel('Speedup')
        a2.set_title('Speedup')

        a3.plot(nodes, [100 * s / n for s, n in zip(s_sp, nodes)],
                label='vs. best single process', **style_marker_kwargs(c))
        a3.plot(nodes, [100 * s / n * n0 for s, n in zip(s_n0, nodes)], color=c, marker='s',
                markerfacecolor=SURFACE, linestyle=':', label=f'vs. best {n0}-node run')
        a3.axhline(100, color=MUTED, linestyle='--', linewidth=1.2, label='Ideal')
        a3.set_ylim(0, None)
        a3.set_ylabel('Parallel efficiency (%)')
        a3.set_title('Efficiency')

        for ax in (a1, a2, a3):
            ax.set_xticks(nodes)
            ax.set_xticklabels([str(n) for n in nodes])
            ax.set_xlim(nodes[0] - 0.6, nodes[-1] + 0.6)
            ax.set_xlabel('Nodes')
            ax.legend(loc='best', fontsize=8.5)
        fig.suptitle(f'Strong scaling, {p} observables (d = 21, p = 0.001)',
                     x=0.01, ha='left', fontweight='bold')
        fig.tight_layout()
        save(fig, f'strong_scaling_{p}obs')
        for n, a, b in zip(nodes, s_sp, s_n0):
            summary.append((p, n, best[n][1], best[n][2], round(best[n][0], 3), round(a, 2), round(b, 2)))
    print('preset nodes ntasks nthreads us_per_round speedup_vs_single_process speedup_vs_1node')
    for row in summary:
        print(*row)


if __name__ == '__main__':
    main()
