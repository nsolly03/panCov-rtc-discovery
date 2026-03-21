#!/usr/bin/env python3
"""
CGCP Phase 2 Step 7b - Supervisor Presentation Figure
3-panel pharmacophore element summary for NSP12-NSP7
Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step07b_supervisor_figure_NSP12-NSP7.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          9,
    'axes.linewidth':     0.75,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.facecolor':     'white',
    'axes.labelsize':     9,
    'axes.labelpad':      6,
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'xtick.major.pad':    4,
    'ytick.major.pad':    4,
    'axes.grid':          False,
    'legend.fontsize':    8,
    'legend.frameon':     False,
    'figure.facecolor':   'white',
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

EL = {
    'E1': {'color': '#D7191C', 'fill': '#FDAE61',
           'name': 'E1 - Aromatic hydrophobic core'},
    'E2': {'color': '#4DAC26', 'fill': '#A6D96A',
           'name': 'E2 - Electrostatic anchor'},
    'E3': {'color': '#2166AC', 'fill': '#92C5DE',
           'name': 'E3 - Distal aromatic patch'},
}

# id, element, pharm_type, x, y, z, cons, composite
FEATURES = [
    ('PHE440', 'E1', 'AROMATIC',    95.961, 79.936, 117.219, 1.000, 1.000),
    ('PRO412', 'E1', 'HYDROPHOBIC', 96.494, 85.414, 124.025, 1.000, 0.895),
    ('PHE442', 'E1', 'AROMATIC',    99.246, 84.670, 118.641, 1.000, 0.774),
    ('GLU431', 'E2', 'NEG_CHARGE',  99.825, 61.004, 116.860, 1.000, 0.748),
    ('TYR420', 'E1', 'AROMATIC',    88.977, 69.828, 123.848, 1.000, 0.616),
    ('PHE415', 'E1', 'AROMATIC',    90.946, 77.963, 125.821, 1.000, 0.527),
    ('PHE843', 'E3', 'AROMATIC',    90.054, 78.878, 121.234, 1.000, 0.488),
    ('GLY413', 'E1', 'HYDROPHOBIC', 96.386, 83.619, 125.012, 1.000, 0.454),
]

CENTROIDS = {
    'E1': {'x': 94.15, 'y': 80.23, 'z': 122.15, 'radius': 11.74},
    'E2': {'x': 99.83, 'y': 61.00, 'z': 116.86, 'radius': 1.0},
    'E3': {'x': 90.05, 'y': 78.88, 'z': 121.23, 'radius': 2.0},
}

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
OUT_DIR = os.path.join(
    BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-07-pharmacophore')
os.makedirs(OUT_DIR, exist_ok=True)


def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax * 0.04, ymax * 1.12)
    if yticks is not None:
        ax.set_yticks(yticks)


def panel_a(ax):
    for el_id, ct in CENTROIDS.items():
        col  = EL[el_id]['color']
        fill = EL[el_id]['fill']
        circ = Circle(
            (ct['x'], ct['z']),
            radius=ct['radius'],
            facecolor=fill, edgecolor=col,
            linewidth=1.2, alpha=0.25, zorder=1,
        )
        ax.add_patch(circ)
        ax.scatter(ct['x'], ct['z'], marker='+',
                   color=col, s=120, linewidths=1.5, zorder=4)
        ax.text(ct['x'], ct['z'] + ct['radius'] + 0.8,
                el_id, ha='center', va='bottom',
                fontsize=9, color=col, fontweight='bold')

    for fid, el, ftype, x, y, z, cons, comp in FEATURES:
        col  = EL[el]['color']
        fill = EL[el]['fill']
        is_anchor = (fid == 'PHE440')
        ax.scatter(x, z,
                   color=fill, edgecolors=col, linewidths=0.9,
                   s=180 if is_anchor else 70,
                   marker='*' if is_anchor else 'o',
                   zorder=5)
        ox, oz = 1.2, 0.8
        if fid in ('TYR420', 'PHE415', 'PHE843'):
            ox = -3.5
        if fid == 'GLU431':
            oz = -1.5
        ax.text(x + ox, z + oz, fid,
                fontsize=7.5, color=col,
                fontweight='bold' if is_anchor else 'normal',
                va='center')

    ax.set_xlabel('X coordinate (A)', fontsize=9)
    ax.set_ylabel('Z coordinate (A)', fontsize=9)
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=EL[el]['fill'],
                       edgecolor=EL[el]['color'],
                       linewidth=0.9, label=EL[el]['name'])
        for el in ['E1', 'E2', 'E3']
    ]
    handles.append(plt.Line2D([0], [0], marker='*', color='w',
                               markerfacecolor='#D7191C',
                               markersize=10, label='PHE440 anchor'))
    ax.legend(handles=handles, loc='lower right', fontsize=7,
              frameon=True, facecolor='white',
              edgecolor='#CCCCCC', framealpha=0.9)
    ax.set_title(
        'A   2D spatial projection (X-Z plane)\n'
        '(circle = element zone; cross = centroid)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)


def panel_b(ax):
    feats_sorted = sorted(FEATURES, key=lambda f: -f[7])
    labels = [f[0] for f in feats_sorted]
    vals   = [f[7] for f in feats_sorted]
    cols   = [EL[f[1]]['color'] for f in feats_sorted]
    fills  = [EL[f[1]]['fill']  for f in feats_sorted]

    x = np.arange(len(feats_sorted))
    w = 0.38

    for xi, val, fc, ec, f in zip(x, vals, fills, cols, feats_sorted):
        is_anchor = f[0] == 'PHE440'
        ax.bar(xi, val, width=w,
               facecolor=fc, edgecolor=ec,
               linewidth=1.2 if is_anchor else 0.9,
               zorder=3, clip_on=False)
        ax.scatter(xi, val, color=ec,
                   s=22 if is_anchor else 16,
                   zorder=5, clip_on=False)
        ax.text(xi, val + 0.018, f'{val:.3f}',
                ha='center', va='bottom',
                fontsize=7, color=ec,
                fontweight='bold', clip_on=False)
        ax.text(xi, -0.09, f[1],
                ha='center', va='top',
                fontsize=6.5, color=ec, clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.5,
                       rotation=90, ha='center', va='top')
    ax.tick_params(axis='x', pad=2, length=3, width=0.75)
    ax.set_xlim(-0.6, len(feats_sorted) - 0.4)
    ax.set_ylabel('CGCP composite score', fontsize=9)
    prism_axes(ax, ymax=1.25, yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])

    handles = [
        mpatches.Patch(facecolor=EL[el]['fill'],
                       edgecolor=EL[el]['color'],
                       linewidth=0.9, label=EL[el]['name'])
        for el in ['E1', 'E2', 'E3']
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=7,
              frameon=True, facecolor='white',
              edgecolor='#CCCCCC', framealpha=0.9)
    ax.set_title(
        'B   Composite score per pharmacophore residue\n'
        '(color = element; badge below = element ID)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)


def panel_c(ax):
    ax.set_facecolor('white')
    n      = len(FEATURES)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    theta  = np.linspace(0, 2 * np.pi, 100)

    for r_val, col in [(1.00, '#D7191C'),
                       (0.75, '#4DAC26'),
                       (0.50, '#2166AC')]:
        ax.plot(theta, [r_val] * 100,
                color=col, linewidth=0.6,
                linestyle='--', alpha=0.4)

    for i, (fid, el, ftype, x, y, z, cons, comp) in \
            enumerate(FEATURES):
        col  = EL[el]['color']
        fill = EL[el]['fill']
        is_anchor = (fid == 'PHE440')
        ax.bar(angles[i], comp,
               width=2 * np.pi / n * 0.72,
               facecolor=fill, edgecolor=col,
               linewidth=0.9, alpha=0.85, zorder=3)
        ax.scatter(angles[i], comp, color=col,
                   s=80 if is_anchor else 30,
                   marker='*' if is_anchor else 'o',
                   zorder=5, edgecolors='white',
                   linewidths=0.4)
        label_r = comp + 0.13
        rot = np.degrees(angles[i])
        if angles[i] > np.pi:
            rot -= 180
        ax.text(angles[i], label_r, fid,
                ha='center', va='center',
                fontsize=7, color=col,
                fontweight='bold' if is_anchor else 'normal',
                rotation=rot)

    for r_val, lbl in [(0.25, '0.25'), (0.5, '0.50'),
                       (0.75, '0.75'), (1.0, '1.00')]:
        ax.text(0, r_val, lbl, ha='center', va='bottom',
                fontsize=6.5, color='#AAAAAA')

    ax.set_ylim(0, 1.40)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['polar'].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=EL[el]['fill'],
                       edgecolor=EL[el]['color'],
                       linewidth=0.9, label=EL[el]['name'])
        for el in ['E1', 'E2', 'E3']
    ]
    ax.legend(handles=handles, loc='lower center',
              bbox_to_anchor=(0.5, -0.14),
              fontsize=7, frameon=True,
              facecolor='white', edgecolor='#CCCCCC',
              ncol=1, framealpha=0.9)
    ax.set_title(
        'C   Pharmacophore feature wheel\n'
        '(radius = composite score; color = element)',
        loc='center', fontsize=8.5, pad=16, linespacing=1.5)


def make_figure():
    fig = plt.figure(figsize=(16.0, 6.5))
    gs  = fig.add_gridspec(
        1, 3,
        wspace=0.42,
        left=0.06, right=0.97,
        top=0.88, bottom=0.28,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2], polar=True)

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    png_path = os.path.join(
        OUT_DIR, 'Fig_Step07_SupervisorPresentation.png')
    fig.savefig(png_path, dpi=300)
    print(f"  PNG 300dpi: {png_path}")

    png150 = os.path.join(
        OUT_DIR, 'Fig_Step07_SupervisorPresentation_150dpi.png')
    fig.savefig(png150, dpi=150)
    print(f"  PNG 150dpi: {png150}")

    pdf_path = os.path.join(
        OUT_DIR, 'Fig_Step07_SupervisorPresentation.pdf')
    fig.savefig(pdf_path, dpi=300, format='pdf')
    print(f"  PDF:        {pdf_path}")

    plt.close()
    print("Done.")


if __name__ == '__main__':
    print('CGCP Step 7b - Supervisor Presentation Figure')
    make_figure()
    print('\nCopy to desktop:')
    print('  cp ~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/'
          'NSP12-NSP7/step-07-pharmacophore/'
          'Fig_Step07_SupervisorPresentation.png '
          '/mnt/c/Users/nseku/Desktop/')
