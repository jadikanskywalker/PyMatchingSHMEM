import matplotlib.pyplot as plt

from common import *

CASES = [
    ('bench_threads_new', 'threads_p001', 'p = 0.001'),
    ('bench_threads_new_p01', 'threads_p01', 'p = 0.01'),
]
COLORS = {16: BLUE, 32: ORANGE, 64: AQUA}


def main():
    apply_style()
    for dirname, name, ptitle in CASES:
        by_m, serial = load_threads(dirname)
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        for m, d in sorted(by_m.items()):
            xs = sorted(d)
            ax.plot(xs, [d[x] for x in xs], label=f'M = {m}', **style_marker_kwargs(COLORS[m]))
        fmt_log_axis(ax, 'x', [16, 32, 64, 128, 256])
        ax.set_ylim(0, None)
        ax.text(0.99, 0.97, f'Serial: {serial:.3g} µs per round', transform=ax.transAxes,
                ha='right', va='top', color=INK2, fontsize=9)
        ax.set_xlabel('Threads')
        ax.set_ylabel('Decoding time per round (µs)')
        ax.set_title(f'Single-process thread scaling, d = 21, {ptitle}')
        ax.legend(title='Partition size', title_fontsize=9, loc='upper right', bbox_to_anchor=(1, 0.93))
        save(fig, name)


if __name__ == '__main__':
    main()
