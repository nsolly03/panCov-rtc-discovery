#!/usr/bin/env python3
"""
CGCP Phase 2 Step 5 - Conservation Overlay: NSP12-NSP8
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from prism_style import (
    apply_prism, prism_axes,
    make_legend, panel_title,
    COLORS
)

apply_prism()

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
S4_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP8/step-04-clusters')
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP8/step-05-conservation')
os.makedirs(OUT_DIR, exist_ok=True)

TIERS = {'identical':1.0,'high':0.8,'moderate':0.6,'variable':0.0}

TIER_COLORS = {
    'identical': '#D7191C',
    'high': '#F46D43',
    'moderate': '#FEE090',
    'variable': '#AAAAAA',
}
TIER_FILLS = {
    'identical': '#FDAE61',
    'high': '#FDD0A2',
    'moderate': '#FFFFBF',
    'variable': '#DDDDDD',
}

def get_tier(cons):
    if cons == 1.0: return 'identical'
    if cons >= 0.8: return 'high'
    if cons >= 0.6: return 'moderate'
    return 'variable'


def load_clusters():
    path = os.path.join(S4_DIR, 'clusters_NSP12-NSP8.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            records.append({
                'residue': row['residue'],
                'chain': row['chain'],
                'position': int(row['position']),
                'aa': row['aa'],
                'primary_feature': row['primary_feature'],
                'conservation': float(row['conservation']),
                'composite': float(row['composite']),
                'cluster': int(row['cluster']),
                'x': float(row['x']),
                'y': float(row['y']),
                'z': float(row['z']),
            })
    return records


def overlay_conservation(records):
    for r in records:
        r['cons_tier'] = get_tier(r['conservation'])

    records.sort(key=lambda x: (
        0 if x['primary_feature']=='anchor_primary' else
        1 if x['primary_feature']=='anchor_secondary' else 2,
        -x['conservation'],
        -x['composite']
    ))
    return records


def make_figure(records):
    fig = plt.figure(figsize=(18.5, 6.5))
    gs  = fig.add_gridspec(1, 1,
                           left=0.06, right=0.97,
                           top=0.88, bottom=0.38)

    ax = fig.add_subplot(gs[0,0])

    # ── Top residues (FIXED) ─────────────────────────
    top = sorted(records, key=lambda x: -x['composite'])[:30]

    x = np.arange(len(top))
    labels = [f"{r['aa'][:3]}{r['position']}" for r in top]

    for xi, r in enumerate(top):
        fc = TIER_FILLS[r['cons_tier']]
        ec = TIER_COLORS[r['cons_tier']]

        is_p = r['primary_feature']=='anchor_primary'
        is_s = r['primary_feature']=='anchor_secondary'

        lw = 2.0 if is_p else (1.5 if is_s else 1.2)

        ax.bar(xi, r['conservation'], width=0.38,
               facecolor=fc, edgecolor=ec,
               linewidth=lw)

    # ── X labels (FIXED CLEAN) ───────────────────────
    reduced_labels = [
        labels[i] if i % 2 == 0 else ''
        for i in range(len(labels))
    ]

    ax.set_xticks(x)
    ax.set_xticklabels(reduced_labels)

    for lbl in ax.get_xticklabels():
        lbl.set_rotation(90)
        lbl.set_fontsize(6)
        lbl.set_ha('center')

    # ── Y axis ───────────────────────────────────────
    ax.set_ylabel('Conservation score',
                  fontsize=11, fontweight='bold')

    prism_axes(ax, ymax=1.1)

    # ── Axis styling ────────────────────────────────
    ax.spines['left'].set_linewidth(1.8)
    ax.spines['bottom'].set_linewidth(1.8)
    ax.tick_params(width=1.6, length=7)

    panel_title(ax, 'A',
                'Conservation (top residues)',
                'Reduced labels for clarity')

    path = os.path.join(OUT_DIR, 'Fig_Step05_FIXED.png')
    fig.savefig(path)
    plt.close()

    print("Figure saved:", path)


if __name__ == '__main__':
    records = load_clusters()
    records = overlay_conservation(records)
    make_figure(records)
