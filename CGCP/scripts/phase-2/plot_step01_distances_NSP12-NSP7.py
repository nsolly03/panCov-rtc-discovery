#!/usr/bin/env python3
"""
Phase 2 Step 1 — Prism distance chart
PHE440 distances to catalytic triad vs NSP7 interface residues
Run: python3 plot_step01_distances_NSP12-NSP7.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

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
    'axes.titlesize':     9,
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'xtick.major.pad':    5,
    'ytick.major.pad':    4,
    'xtick.labelsize':    8.5,
    'ytick.labelsize':    8.5,
    'axes.grid':          False,
    'legend.fontsize':    8,
    'legend.frameon':     False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.12,
    'pdf.fonttype':       42,
})

P = {
    'triad':   '#2166AC',
    'triad_f': '#92C5DE',
    'nsp7':    '#4DAC26',
    'nsp7_f':  '#A6D96A',
    'anchor':  '#D7191C',
    'gray':    '#636363',
    'text':    '#1A1A1A',
}

# Distances from structure (Å)
groups = {
    'Catalytic triad\n(negative control)': {
        'labels': ['ASP618', 'SER759', 'ASP760'],
        'values': [22.8, 25.1, 24.1],
        'fill':   P['triad_f'],
        'edge':   P['triad'],
    },
    'NSP7 interface\n(positive control)': {
        'labels': ['LYS2', 'MET3', 'SER4', 'ASP5', 'VAL6'],
        'values': [17.1, 14.0, 11.7, 13.6, 13.2],
        'fill':   P['nsp7_f'],
        'edge':   P['nsp7'],
    },
}

fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2),
                          gridspec_kw={'width_ratios': [3, 5]})
fig.subplots_adjust(wspace=0.50, bottom=0.22,
                    top=0.86, left=0.10, right=0.97)

for ax, (group_name, gdata) in zip(axes, groups.items()):
    labels = gdata['labels']
    vals   = gdata['values']
    fill   = gdata['fill']
    edge   = gdata['edge']
    x      = np.arange(len(labels))
    w      = 0.48

    for xi, val in zip(x, vals):
        ax.bar(xi, val, width=w, facecolor=fill,
               edgecolor=edge, linewidth=0.9,
               zorder=3, clip_on=False)
        ax.text(xi, val + 0.4, f'{val:.1f}',
                ha='center', va='bottom',
                fontsize=8, color=edge,
                fontweight='bold', clip_on=False)

    # Interface threshold line
    ax.axhline(15.0, color=P['nsp7'], linewidth=0.85,
               linestyle=(0, (4, 3)), zorder=1)
    ax.text(len(labels) - 0.42, 15.6,
            '< 15 Å  interface',
            fontsize=7, color=P['nsp7'],
            ha='right', va='bottom', clip_on=False)

    # Non-interface threshold line
    ax.axhline(25.0, color=P['triad'], linewidth=0.85,
               linestyle=(0, (4, 3)), zorder=1)
    ax.text(len(labels) - 0.42, 25.6,
            '> 25 Å  non-interface',
            fontsize=7, color=P['triad'],
            ha='right', va='bottom', clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5,
                       rotation=45, ha='right')
    ax.set_ylabel('Distance from PHE440 Cα (Å)', fontsize=9)
    ax.set_xlim(-0.55, len(labels) - 0.45)

    ymax   = 32
    yticks = [0, 5, 10, 15, 20, 25, 30]
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['left'].set_bounds(0, ymax)
    ax.set_ylim(-0.5, ymax * 1.02)
    ax.set_yticks(yticks)

    panel_id = 'A' if 'triad' in group_name.lower() or \
               'negative' in group_name.lower() else 'B'
    ax.set_title(
        f'{panel_id}   {group_name}\n'
        f'(distances from PHE440 anchor, Cα–Cα)',
        loc='left', fontsize=8.5,
        pad=6, linespacing=1.5,
    )

    # Pass/fail badge
    all_pass = all(
        (v > 25.0 if 'triad' in gdata['edge']
         else v < 15.0)
        for v in vals
    )
    badge_c = P['nsp7'] if all_pass else P['anchor']
    badge_f = P['nsp7_f'] if all_pass else '#FDAE61'
    verdict = 'PASS' if all_pass else 'PARTIAL'
    ax.text(0.97, 0.04, verdict,
            transform=ax.transAxes,
            ha='right', va='bottom',
            fontsize=8.5, fontweight='bold',
            color=badge_c,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor=badge_f,
                      edgecolor=badge_c,
                      linewidth=0.7))

# Shared annotation
fig.text(
    0.53, 0.01,
    'Residue distances measured from PHE440 Cα  '
    '(receptor_NSP12-NSP7_3.pdb, chain A)',
    ha='center', fontsize=7.5,
    color=P['gray'], style='italic',
)

OUT = os.path.expanduser(
    '~/projects/rtc-pan-coronavirus/CGCP/'
    '02-deep-dive/NSP12-NSP7/step-01-structure/visuals'
)
os.makedirs(OUT, exist_ok=True)
out_path = os.path.join(OUT, 'Fig_Step01_Distances_Prism.png')
fig.savefig(out_path)
plt.close()
print(f"Saved: {out_path}")
