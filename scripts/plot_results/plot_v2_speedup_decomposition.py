import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch

from common_v2 import *

N_TOP = 16
NODES = [2, 4, 8, 16]
FAMILIES = [('speedup_decomposition_36_72_144obs', [36, 72, 144]),
            ('speedup_decomposition_64_128_256obs', [64, 128, 256])]


def main():
    apply_style_v2()
    alphas = dict(zip(NODES, [0.12, 0.25, 0.4, 0.6]))
    w = 0.15
    offs = {k: (j - 2) * w for j, k in enumerate(['D'] + NODES)}
    for name, presets in FAMILIES:
        pts = Points()
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.6),
                                     gridspec_kw=dict(width_ratios=[2.2, 1.3]))
        for i, p in enumerate(presets):
            best = t_nodes(p)
            sp = t_sp(p)
            d = sp[0] / best[1][0]
            c = PRESET_COLORS[p]
            a1.bar(i + offs['D'], d - 1, w, bottom=1, color=c, edgecolor=c)
            row = dict(preset=p, T_sp_us=sp[0], T_1n_us=best[1][0], D=d)
            for n in NODES:
                s_n = best[1][0] / best[n][0]
                a1.bar(i + offs[n], s_n - 1, w, bottom=1, color=to_rgba(c, alphas[n]),
                       edgecolor=c, hatch='////', linewidth=0.8)
                row[f'S_{n}'] = s_n
            row['product_D_S_16'] = d * row['S_16']
            pts.add('D_and_S', **row)
        log_y(a1, [0.5, 1, 2, 5, 10, 20])
        a1.set_ylim(0.4, 80)
        a1.axhline(1, color=MUTED, linewidth=0.8)
        a1.set_xlim(-0.6, len(presets) - 0.4)
        a1.set_xticks(range(len(presets)))
        a1.set_xticklabels([str(p) for p in presets])
        a1.set_xlabel('Observables')
        a1.set_ylabel('Speedup factor')
        a1.legend(handles=[Patch(facecolor=MUTED, edgecolor=MUTED,
                                 label='$D$: thread-only → best 1-node')]
                  + [Patch(facecolor=to_rgba(MUTED, alphas[n]), edgecolor=MUTED, hatch='////',
                           label=f'$S_{{{n}}}$: 1 → {n} nodes') for n in NODES],
                  loc='upper left', ncol=2, columnspacing=0.8)
        a1.set_title('Where the speedup comes from')

        big = presets[-1]
        best = t_nodes(big)
        n_best = min(best, key=lambda n: best[n][0])  # node count with the lowest T_N
        vals = [t_sp(big)[0], best[1][0], best[n_best][0]]
        pts.add('waterfall', preset=big, nodes_final=n_best, T_sp_us=vals[0], T_1n_us=vals[1],
                T_final_us=vals[2])
        a2.bar(range(3), vals, 0.4, color=PRESET_COLORS[big])
        for i, v in enumerate(vals):
            a2.text(i, v + vals[0] * 0.02, f'{v:.1f}', ha='center', va='bottom')
        for i in range(2):
            x = i + 0.5
            a2.plot([i + 0.2, x], [vals[i]] * 2, color=MUTED, linestyle=':', linewidth=1.0)
            a2.plot([x, i + 0.8], [vals[i + 1]] * 2, color=MUTED, linestyle=':', linewidth=1.0)
            a2.annotate('', xy=(x, vals[i + 1]), xytext=(x, vals[i]),
                        arrowprops=dict(arrowstyle='->', color=MUTED, linewidth=1.0))
            a2.text(x, (vals[i] + vals[i + 1]) / 2, f'÷{vals[i] / vals[i + 1]:.1f}', ha='center',
                    va='center', color=INK2, bbox=dict(facecolor=SURFACE, edgecolor='none', pad=0.8))
        a2.set_xticks(range(3))
        a2.set_xticklabels(['Thread-\nonly', 'Best\n1-node', f'{n_best}\nnodes'])
        a2.set_ylabel('Time per round (µs)')
        a2.set_ylim(0, vals[0] * 1.2)
        a2.set_title(f'{big} observables')
        fig.tight_layout()
        save_v2(fig, name, pts)


if __name__ == '__main__':
    main()
