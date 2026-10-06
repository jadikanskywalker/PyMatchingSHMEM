import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from common_v2 import *

FAMILIES = [[36, 72, 144], [64, 128, 256]]


def main():
    apply_style_v2()
    pts = Points()
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 5.6), sharey='row')
    for col, fam in enumerate(FAMILIES):
        top, bot = axes[0][col], axes[1][col]
        for p in fam:
            best = t_nodes(p)
            nodes = list(best)
            top.plot(nodes, [best[n][0] for n in nodes], label=f'{p} obs', **mk(PRESET_COLORS[p]))
            bot.plot(nodes, [100 * best[1][0] / best[n][0] / n for n in nodes],
                     **mk(PRESET_COLORS[p]))
            sp = t_sp(p)
            top.plot([1], [sp[0]], marker='o', markerfacecolor=SURFACE,
                     markeredgecolor=PRESET_COLORS[p], markeredgewidth=1.5, linestyle='none')
            pts.add('T_sp', preset=p, nodes=1, ntasks=sp[1], nthreads=sp[2], us=sp[0])
            for n in nodes:
                pts.add('T_N', preset=p, nodes=n, ntasks=best[n][1], nthreads=best[n][2],
                        us=best[n][0], efficiency_pct=100 * best[1][0] / best[n][0] / n)
        top.axhline(1, **ONE_US_STYLE)
        top.text(0.99, 0.8, ONE_US_LABEL, ha='right', va='top', color=MUTED,
                 transform=top.get_yaxis_transform())
        bot.axhline(100, **IDEAL_STYLE)
        top.set_title(f'{fam[0]}→{fam[1]}→{fam[2]} observables')
        log_y(top, [1, 2, 5, 10, 20, 50])
        top.set_ylim(0.4, 130)
        bot.set_ylim(0, 110)
        h, l = top.get_legend_handles_labels()
        h.append(Line2D([], [], marker='o', markerfacecolor=SURFACE, markeredgecolor=MUTED,
                        markeredgewidth=1.5, linestyle='none'))
        l.append('Thread-only (1 process)')
        top.legend(h, l, loc='upper right', handlelength=1.6, labelspacing=0.3)
        for ax in (top, bot):
            node_axis(ax)
    axes[0][0].set_ylabel('Time per round, $T_N$ (µs)')
    axes[1][0].set_ylabel('Node efficiency, $E_N$ (%)')
    fig.suptitle('Strong scaling relative to the best 1-node configuration (d = 21, p = 0.001)',
                 x=0.01, ha='left', fontweight='bold')
    fig.tight_layout()
    save_v2(fig, 'scaling_by_observables', pts)


if __name__ == '__main__':
    main()
