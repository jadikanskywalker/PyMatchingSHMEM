import matplotlib.pyplot as plt

from common_v2 import *

CASES = [('bench_threads_new', 'threads_p001', 'p = 0.001'),
         ('bench_threads_new_p01', 'threads_p01', 'p = 0.01')]
COLORS = {16: BLUE, 32: ORANGE, 64: AQUA}


def main():
    apply_style_v2()
    for dirname, name, ptitle in CASES:
        pts = Points()
        by_m, serial = load_threads(dirname)
        fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.1))
        best = (1e18, None, None)
        for m, d in sorted(by_m.items()):
            xs = sorted(d)
            ax.plot(xs, [d[x] for x in xs], label=f'M = {m}', **mk(COLORS[m]))
            for x in xs:
                pts.add('threads', M=m, nthreads=x, us=d[x])
                best = min(best, (d[x], m, x))
        pts.add('serial', M='', nthreads=1, us=serial)
        ax.axvline(181, ymax=0.5, color=MUTED, linestyle=':', linewidth=0.9, zorder=0)
        ax.text(181 / 1.05, 0.02, '1 socket', ha='right', va='bottom', color=MUTED,
                transform=ax.get_xaxis_transform())
        ax.text(181 * 1.05, 0.02, '2 sockets', ha='left', va='bottom', color=MUTED,
                transform=ax.get_xaxis_transform())
        ax.text(0.98, 0.98, f'Serial: {serial:.3g} µs\nBest: {best[0]:.2g} µs '
                f'(≈{serial / best[0]:.0f}× serial)', transform=ax.transAxes, ha='right',
                va='top', color=INK2)
        ax.set_xscale('log', base=2)
        ax.set_xticks([16, 32, 64, 128, 256])
        ax.set_xticklabels(['16', '32', '64', '128', '256'])
        ax.minorticks_off()
        ax.set_ylim(0, None)
        ax.set_xlabel('Threads')
        ax.set_ylabel('Time per round (µs)')
        ax.set_title(f'Single-process thread scaling,\nd = 21, {ptitle}')
        ax.legend(loc='upper right', bbox_to_anchor=(1, 0.82), handlelength=1.6,
                  labelspacing=0.3)
        fig.tight_layout()
        save_v2(fig, name, pts)


if __name__ == '__main__':
    main()
