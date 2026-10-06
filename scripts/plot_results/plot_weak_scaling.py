import matplotlib.pyplot as plt

from common import *

# "Proportional" series: nodes grow with observables, starting from BASE_NODES at the smallest preset.
BASE_NODES = 1
SERIES = [('weak_scaling_36_72_144obs', [36, 72, 144]), ('weak_scaling_64_128_256obs', [64, 128, 256])]
COLOR = BLUE


def main():
    apply_style()
    for name, presets in SERIES:
        prop_nodes = [BASE_NODES * p // presets[0] for p in presets]
        best = {p: min(best_per_node(p).items(), key=lambda kv: kv[1][0]) for p in presets}
        t_best = [best[p][1][0] for p in presets]
        t_prop = [best_per_node(p)[n][0] for p, n in zip(presets, prop_nodes)]
        eff = [100 * t_prop[0] / x for x in t_prop]
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
        kw = style_marker_kwargs(COLOR)
        a1.plot(presets, t_best, label='Fastest run at any node count', **kw)
        for p, t, (n, _) in zip(presets, t_best, (best[p] for p in presets)):
            a1.annotate(f'{n} nodes', (p, t), textcoords='offset points', xytext=(0, -16),
                        ha='center', color=INK2, fontsize=9)
        a1.plot(presets, t_prop, color=COLOR, marker='s', markerfacecolor=SURFACE, linestyle=':',
                label=f'Constant observables per node ({"/".join(map(str, prop_nodes))} nodes)')
        for p, t, n in zip(presets, t_prop, prop_nodes):
            a1.annotate(f'{n} node' + 's' * (n > 1), (p, t), textcoords='offset points',
                        xytext=(0, 9), ha='center', color=INK2, fontsize=9)
        a1.set_ylim(0, max(t_prop + t_best) * 1.2)
        a1.set_ylabel('Decoding time per round (µs)')
        a1.set_title('Time per round')
        a1.legend(loc='upper left', fontsize=8.5)
        a2.plot(presets, eff, color=COLOR, marker='o', markeredgecolor=SURFACE, markeredgewidth=1.5)
        a2.axhline(100, color=MUTED, linestyle='--', linewidth=1.2)
        a2.text(presets[-1], 100, 'ideal', ha='right', va='bottom', color=INK2, fontsize=9)
        for p, e in zip(presets, eff):
            a2.annotate(f'{e:.0f}%', (p, e), textcoords='offset points', xytext=(0, -16),
                        ha='center', color=INK2, fontsize=9)
        a2.set_ylim(0, 120)
        a2.set_ylabel('Parallel efficiency (%)')
        a2.set_title('Weak-scaling efficiency, T(base) / T(N)')
        for ax in (a1, a2):
            ax.set_xscale('log', base=2)
            ax.set_xticks(presets)
            ax.set_xticklabels([str(p) for p in presets])
            ax.minorticks_off()
            ax.set_xlim(presets[0] / 1.25, presets[-1] * 1.25)
            ax.set_xlabel('Observables')
        fig.suptitle(f'Weak scaling, {presets[0]}→{presets[1]}→{presets[2]} observables '
                     f'(d = 21, p = 0.001)', x=0.01, ha='left', fontweight='bold')
        fig.tight_layout()
        save(fig, name)
        print(name, [(p, best[p][0], round(b, 2), round(x, 2)) for p, b, x in zip(presets, t_best, t_prop)])


if __name__ == '__main__':
    main()
