import matplotlib.pyplot as plt

from common_v2 import *

SERIES = [('weak_scaling_36_72_144obs', [36, 72, 144]), ('weak_scaling_64_128_256obs', [64, 128, 256])]
# nodes per base-preset unit: level 1 -> 1/2/4 nodes, level 2 -> 2/4/8 nodes
LEVELS = [(1, '-', 'o'), (2, '--', 's')]


def main():
    apply_style_v2()
    for name, presets in SERIES:
        pts = Points()
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.4))
        ymax = 0
        for mult, ls, marker in LEVELS:
            nodes = [mult * p // presets[0] for p in presets]
            runs = [t_nodes(p)[n] for p, n in zip(presets, nodes)]
            t = [r[0] for r in runs]
            eff = [100 * t[0] / x for x in t]
            lab = f'{mult} node{"s" * (mult > 1)} per {presets[0]} obs'
            for ax, ys in ((a1, t), (a2, eff)):
                ax.plot(presets, ys, color=MUTED, linestyle=ls, linewidth=1.2, marker=marker,
                        markerfacecolor=MUTED, markeredgecolor=MUTED, label=lab)
                for p, y in zip(presets, ys):
                    ax.plot([p], [y], marker=marker, color=PRESET_COLORS[p],
                            markeredgecolor=SURFACE, markeredgewidth=1.0, linestyle='none')
            for p, n, r, e in zip(presets, nodes, runs, eff):
                a1.annotate(f'{r[1]}×{r[2]}', (p, r[0]), textcoords='offset points',
                            xytext=(0, 7 if mult == 1 else -14), ha='center', color=INK2)
                pts.add(f'{mult}x', preset=p, nodes=n, ntasks=r[1], nthreads=r[2], us=r[0],
                        efficiency_pct=e)
            ymax = max(ymax, *t)
        a1.set_ylim(0, ymax * 1.25)
        a1.set_ylabel('Time per round (µs)')
        a1.set_title('Time per round (labels: tasks×threads)')
        a2.axhline(100, color=MUTED, linestyle='-', linewidth=0.8)
        a2.set_ylim(0, 115)
        a2.set_ylabel('Weak-scaling efficiency (%)')
        a2.set_title('Efficiency, $T(\\mathrm{base})$ / $T(N)$')
        for ax in (a1, a2):
            ax.set_xscale('log', base=2)
            ax.set_xticks(presets)
            ax.set_xticklabels([str(p) for p in presets])
            ax.minorticks_off()
            ax.set_xlim(presets[0] / 1.25, presets[-1] * 1.25)
            ax.set_xlabel('Observables')
        a1.legend(loc='upper left', handlelength=2.2)
        fig.suptitle(f'Weak scaling, {presets[0]}→{presets[1]}→{presets[2]} observables '
                     f'(d = 21, p = 0.001)', x=0.01, ha='left', fontweight='bold')
        fig.tight_layout()
        save_v2(fig, name, pts)


if __name__ == '__main__':
    main()
