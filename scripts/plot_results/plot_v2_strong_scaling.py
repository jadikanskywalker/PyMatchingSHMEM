import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from common_v2 import *


def main():
    apply_style_v2()
    for p in PRESETS:
        pts = Points()
        best = t_nodes(p)
        nodes = list(best)
        t = [best[n][0] for n in nodes]
        t1 = t[0]
        sp = t_sp(p)
        c = PRESET_COLORS[p]
        fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(DOUBLE_COL, 3.4))
        a1.plot(nodes, t, **mk(c))
        a1.plot(nodes, [t1 / n for n in nodes], **IDEAL_STYLE)
        a1.plot([1], [sp[0]], marker='o', markerfacecolor=SURFACE, markeredgecolor=c,
                markeredgewidth=1.5, linestyle='none')
        a1.axhline(1, **ONE_US_STYLE)
        log_y(a1, [0.5, 1, 2, 5, 10, 20, 50, 100])
        a1.set_ylim(0.4, 110)
        a1.set_ylabel('Time per round (µs)')
        s = [t1 / x for x in t]
        a2.plot(nodes, s, **mk(c))
        a2.plot(nodes, nodes, **IDEAL_STYLE)
        log_y(a2, [0.5, 1, 2, 4, 8, 16], base=2)
        a2.set_ylim(0.4, 20)
        a2.set_ylabel('Node speedup, $S_N$')
        a3.plot(nodes, [100 * x / n for x, n in zip(s, nodes)], **mk(c))
        a3.axhline(100, **IDEAL_STYLE)
        a3.set_ylim(0, 110)
        a3.set_ylabel('Node efficiency, $E_N$ (%)')
        for ax in (a1, a2, a3):
            node_axis(ax)
        for n, x in zip(nodes, t):
            pts.add('T_N', preset=p, nodes=n, ntasks=best[n][1], nthreads=best[n][2], us=x,
                    speedup=t1 / x, efficiency_pct=100 * t1 / x / n)
        pts.add('T_sp', preset=p, nodes=1, ntasks=sp[1], nthreads=sp[2], us=sp[0])
        handles = [Line2D([], [], **mk(c)), Line2D([], [], **IDEAL_STYLE),
                   Line2D([], [], marker='o', markerfacecolor=SURFACE, markeredgecolor=c,
                          markeredgewidth=1.5, linestyle='none'),
                   Line2D([], [], **ONE_US_STYLE)]
        labels = ['Best config', 'Ideal', 'Thread-only (1 process)', ONE_US_LABEL]
        fig.legend(handles, labels, loc='lower center', ncol=2, columnspacing=1.5,
                   handlelength=1.8, labelspacing=0.3)
        fig.suptitle(f'Strong scaling, {p} observables (d = 21, p = 0.001)', x=0.01, ha='left',
                     fontweight='bold')
        fig.tight_layout(rect=[0, 0.15, 1, 1])
        save_v2(fig, f'strong_scaling_{p}obs', pts)


if __name__ == '__main__':
    main()
