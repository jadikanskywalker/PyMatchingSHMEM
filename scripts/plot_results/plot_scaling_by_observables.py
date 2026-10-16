import matplotlib.pyplot as plt

from common import *
from plot_strong_scaling import COLORS

FAMILIES = [[36, 72, 144], [64, 128, 256]]


def main():
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), sharex=False)
    for col, fam in enumerate(FAMILIES):
        top, bot = axes[0][col], axes[1][col]
        for p in fam:
            best = best_per_node(p)
            nodes = list(best)
            sp = best_single_process(p)
            kw = style_marker_kwargs(COLORS[p])
            top.plot(nodes, [sp / best[n][0] for n in nodes], label=f'{p} obs', **kw)
            bot.plot(nodes, [100 * best[nodes[0]][0] / best[n][0] / n for n in nodes], **kw)
        top.plot([1, 16], [1, 16], color=MUTED, linestyle='--', linewidth=1.2, label='Ideal')
        bot.axhline(100, color=MUTED, linestyle='--', linewidth=1.2)
        top.set_title(f'{fam[0]}→{fam[1]}→{fam[2]} observables')
        top.set_ylabel('Speedup vs. best single process')
        bot.set_ylabel('Efficiency vs. best 1-node run (%)')
        bot.set_ylim(0, 110)
        top.legend(fontsize=8.5)
        for ax in (top, bot):
            ax.set_xticks([1, 2, 4, 8, 16])
            ax.tick_params(labelbottom=True)
        top.set_xlabel('Nodes')
        bot.set_xlabel('Nodes')
    fig.suptitle('Scaling by observable count (d = 21, p = 0.001)', x=0.01, ha='left',
                 fontweight='bold')
    fig.tight_layout()
    save(fig, 'scaling_by_observables')


if __name__ == '__main__':
    main()
