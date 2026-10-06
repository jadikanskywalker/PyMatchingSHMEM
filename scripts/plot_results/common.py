import csv
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
PLOTS_DIR = os.path.join(ROOT, 'plots')
TIME_COL = 'Time per Round (us)'

# Rows to drop, keyed (preset, ntasks, nthreads, ntps); results.csv carries no job id.
# 144obs job 167690 spread 2 64-thread tasks per PE instead of 4.
EXCLUDE = {(144, 16, 64, 2)}

# Fixed hue per entity (categorical slots 1-3 of the reference palette)
INK, INK2, MUTED, GRID, SURFACE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#fcfcfb'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'


def apply_style():
    plt.rcParams.update({
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
        'font.family': ['DejaVu Sans'], 'font.size': 10,
        'text.color': INK, 'axes.labelcolor': INK2, 'axes.titlecolor': INK,
        'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.edgecolor': MUTED,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
        'axes.titlesize': 11, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
        'lines.linewidth': 2, 'lines.markersize': 7,
        'legend.frameon': False, 'legend.labelcolor': INK2,
    })


def style_marker_kwargs(color):
    # 2px surface ring keeps overlapping markers separable
    return dict(color=color, marker='o', markeredgecolor=SURFACE, markeredgewidth=1.5)


def save(fig, name):
    os.makedirs(PLOTS_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(PLOTS_DIR, f'{name}.{ext}'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote plots/{name}.{{png,pdf}}')


def _forward_fill_m(rows):
    m = None
    for r in rows:
        if r['M'] != '':
            m = int(r['M'])
        r['M'] = m
    return rows


def load_obs(preset):
    """Rows of a many-obs results.csv: dicts with ntasks, nodes, nthreads, us (per round)."""
    paths = glob.glob(os.path.join(
        ROOT, f'bench_{preset}obs_new_OSSS', f'bench_{preset}obs_d21_p001_*r_100s_L16', 'results.csv'))
    assert len(paths) == 1, (preset, paths)
    out = []
    with open(paths[0]) as f:
        for r in _forward_fill_m(list(csv.DictReader(f))):
            ntasks, nthreads, ntps = int(r['ntasks']), int(r['nthreads']), int(r['ntps'])
            if (preset, ntasks, nthreads, ntps) in EXCLUDE:
                continue
            out.append(dict(ntasks=ntasks, nthreads=nthreads, nodes=int(float(r['nodes'])),
                            us=float(r[TIME_COL])))
    return out


def best_per_node(preset):
    """{nodes: (min us/round, ntasks, nthreads)} across all configs at that node count."""
    best = {}
    for r in load_obs(preset):
        # the serial row (ntasks=1, nthreads=1) is not a SHMEM config
        if r['nthreads'] < 2 and r['ntasks'] == 1:
            continue
        cur = best.get(r['nodes'])
        if cur is None or r['us'] < cur[0]:
            best[r['nodes']] = (r['us'], r['ntasks'], r['nthreads'])
    return dict(sorted(best.items()))


def best_single_process(preset):
    """Min us/round over all ntasks=1 rows (any thread count)."""
    return min(r['us'] for r in load_obs(preset) if r['ntasks'] == 1)


def load_threads(dirname):
    """({M: {nthreads: us}}, serial_us) from a threads results.csv (M blank on later rows)."""
    with open(os.path.join(ROOT, dirname, 'results.csv')) as f:
        rows = list(csv.DictReader(f))
    serial = next(float(r[TIME_COL]) for r in rows if r['nthreads'] == '1')
    by_m, m = {}, None
    for r in rows:
        if int(r['nthreads']) == 1:
            continue
        if r['M'] != '':
            m = int(r['M'])
        by_m.setdefault(m, {})[int(r['nthreads'])] = float(r[TIME_COL])
    return by_m, serial


def fmt_log_axis(ax, axis, ticks):
    ax.set_xscale('log', base=2) if axis == 'x' else ax.set_yscale('log')
    if axis == 'x':
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.minorticks_off()
