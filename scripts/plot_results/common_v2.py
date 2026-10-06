import csv
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter, NullLocator

from common import *  # noqa: F401,F403 (loaders, ink/palette constants, apply_style)

SUFFIX = '_v2'
NODE_TICKS = [1, 2, 4, 8, 16]
PRESETS = [36, 64, 72, 128, 144, 256]
# One fixed hue per observable count, in every figure
PRESET_COLORS = dict(zip(PRESETS, [BLUE, ORANGE, AQUA, '#eda100', '#e87ba4', '#4a3aa7']))
ONE_US_LABEL = '≈ 1 µs superconducting syndrome round'
ONE_US_STYLE = dict(color=MUTED, linestyle=':', linewidth=1.0)
IDEAL_STYLE = dict(color=MUTED, linestyle='-', linewidth=1.2, alpha=0.6)
SINGLE_COL, DOUBLE_COL = 3.5, 7.16


def apply_style_v2():
    apply_style()
    plt.rcParams.update({
        'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9, 'xtick.labelsize': 9,
        'ytick.labelsize': 9, 'legend.fontsize': 9, 'legend.title_fontsize': 9,
        'figure.titlesize': 9, 'lines.linewidth': 1.6, 'lines.markersize': 5.5,
    })


def mk(color, **kw):
    return dict(color=color, marker='o', markeredgecolor=SURFACE, markeredgewidth=1.0, **kw)


class Points:
    """Traceability sidecar: one row per plotted point with its chosen configuration."""

    def __init__(self):
        self.rows = []

    def add(self, series, **kw):
        self.rows.append(dict(series=series, **kw))


def save_v2(fig, name, points):
    name += SUFFIX
    os.makedirs(PLOTS_DIR, exist_ok=True)
    fig.savefig(os.path.join(PLOTS_DIR, f'{name}.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(PLOTS_DIR, f'{name}.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    cols = []
    for r in points.rows:
        cols += [k for k in r if k not in cols]
    with open(os.path.join(PLOTS_DIR, f'{name}_points.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in points.rows:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})
    print(f'wrote plots/{name}.{{png,pdf}} + _points.csv')


def node_axis(ax):
    ax.set_xscale('log', base=2)
    ax.set_xticks(NODE_TICKS)
    ax.set_xticklabels([str(n) for n in NODE_TICKS])
    ax.minorticks_off()
    ax.set_xlim(0.8, 20)
    ax.tick_params(labelbottom=True)
    ax.set_xlabel('Nodes')


def log_y(ax, ticks, base=10):
    ax.set_yscale('log', base=base)
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_formatter(NullFormatter())


def t_sp(p):
    """(us, ntasks, nthreads) of the best ntasks=1 run (thread-only, 1 process)."""
    return min((r['us'], r['ntasks'], r['nthreads']) for r in load_obs(p) if r['ntasks'] == 1)


def t_nodes(p):
    """{nodes: (us, ntasks, nthreads)}: best run at each node count (T_N); T_1n is entry 1."""
    return best_per_node(p)
