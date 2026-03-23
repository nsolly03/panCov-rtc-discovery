#!/usr/bin/env python3
"""
CGCP Phase 2 Step 4 - DBSCAN Spatial Clustering: NSP12-NSP8
Dual anchor: LEU387 (primary) + LYS332 (secondary electrostatic)
"""

import os
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from sklearn.cluster import DBSCAN
from prism_style import (
    apply_prism, prism_axes,
    set_xticklabels_vertical,
    make_legend, panel_title,
    COLORS, CLUSTER_COLORS, CLUSTER_FILLS
)

apply_prism()

# ── JSON encoder ─────────────────────────────────────────────
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):  return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray):  return obj.tolist()
        if isinstance(obj, np.bool_):    return bool(obj)
        return super().default(obj)

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
S3_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP8/step-03-features')
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP8/step-04-clusters')
os.makedirs(OUT_DIR, exist_ok=True)

EPS         = 8.0
MIN_SAMPLES = 2
NOISE_COL   = COLORS['gray']
NOISE_FILL  = COLORS['gray_fill']

# ── Short decision labels ────────────────────────────────────
def short_label(decision):
    return {
        'PRIMARY PHARMACOPHORE': 'PRIMARY',
        'ELECTROSTATIC ANCHOR': 'ELECTRO',
        'SECONDARY PHARMACOPHORE': 'SECOND',
        'SUPPORTING': 'SUPPORT',
        'DEPRIORITIZE': 'LOW'
    }.get(decision, decision[:10])


# ── Load features ────────────────────────────────────────────
def load_features():
    path = os.path.join(S3_DIR, 'feature_classification_NSP12-NSP8.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['ca_x'] and row['ca_y'] and row['ca_z']:
                records.append({
                    'residue': row['residue'],
                    'chain': row['chain'],
                    'position': int(row['position']),
                    'aa': row['aa'],
                    'primary_feature': row['primary_feature'],
                    'contact_score': int(row['contact_score']),
                    'conservation': float(row['conservation']),
                    'composite': float(row['composite']),
                    'is_anchor': int(row['is_anchor']),
                    'x': float(row['ca_x']),
                    'y': float(row['ca_y']),
                    'z': float(row['ca_z']),
                })
    return records


# ── Run DBSCAN ───────────────────────────────────────────────
def run_dbscan(records):
    coords = np.array([[r['x'], r['y'], r['z']] for r in records])
    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES).fit(coords)
    return [int(l) for l in db.labels_]


# ── Analyze clusters ─────────────────────────────────────────
def analyze_clusters(records, labels):
    cluster_info = {}

    for label in sorted(set(labels)):
        members = [records[i] for i, l in enumerate(labels) if l == label]
        if not members:
            continue

        coords = np.array([[r['x'], r['y'], r['z']] for r in members])
        centroid = coords.mean(axis=0)

        radius = float(max(np.linalg.norm(c - centroid) for c in coords))
        mean_cons = float(np.mean([r['conservation'] for r in members]))
        mean_comp = float(np.mean([r['composite'] for r in members]))

        has_p = any(r['primary_feature'] == 'anchor_primary' for r in members)
        has_s = any(r['primary_feature'] == 'anchor_secondary' for r in members)

        if has_p and mean_cons >= 0.700:
            decision = 'PRIMARY PHARMACOPHORE'
        elif has_s:
            decision = 'ELECTROSTATIC ANCHOR'
        elif mean_cons >= 0.800 and mean_comp >= 0.500:
            decision = 'SECONDARY PHARMACOPHORE'
        elif mean_cons >= 0.600:
            decision = 'SUPPORTING'
        else:
            decision = 'DEPRIORITIZE'

        cluster_info[label] = {
            'label': int(label),
            'n_residues': len(members),
            'centroid_x': float(centroid[0]),
            'centroid_y': float(centroid[1]),
            'centroid_z': float(centroid[2]),
            'radius': round(radius, 2),
            'mean_cons': round(mean_cons, 3),
            'mean_comp': round(mean_comp, 3),
            'has_primary_anchor': has_p,
            'has_secondary_anchor': has_s,
            'decision': decision,
            'members': [r['residue'] for r in members],
        }

    return cluster_info


# ── Make figure ──────────────────────────────────────────────
def make_figure(records, labels, cluster_info):
    fig = plt.figure(figsize=(16.5, 6.5))
    gs  = fig.add_gridspec(1, 2, wspace=0.42,
                           left=0.07, right=0.97,
                           top=0.88, bottom=0.28)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # ── Panel A ──────────────────────────────────────────────
    plotted = set()

    for i, r in enumerate(records):
        lbl = labels[i]

        if lbl == -1:
            col, fill = NOISE_COL, NOISE_FILL
            label = 'Noise' if 'noise' not in plotted else None
            plotted.add('noise')
        else:
            idx = lbl % len(CLUSTER_COLORS)
            col, fill = CLUSTER_COLORS[idx], CLUSTER_FILLS[idx]
            label = f"C{lbl}" if lbl not in plotted else None
            plotted.add(lbl)

        is_p = r['primary_feature'] == 'anchor_primary'
        is_s = r['primary_feature'] == 'anchor_secondary'

        size = 180 if is_p else (110 if is_s else 40)
        marker = '*' if is_p else ('D' if is_s else 'o')

        ax1.scatter(r['x'], r['z'],
                    color=fill, edgecolors=col,
                    linewidths=1.0, s=size,
                    marker=marker, alpha=0.85,
                    label=label)

    # centroids
    for lbl, info in cluster_info.items():
        if lbl < 0:
            continue
        idx = lbl % len(CLUSTER_COLORS)
        ax1.scatter(info['centroid_x'], info['centroid_z'],
                    marker='+', color=CLUSTER_COLORS[idx],
                    s=180, linewidths=2.5)

    ax1.set_xlabel('X coordinate (Å)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Z coordinate (Å)', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=8, frameon=True)

    panel_title(ax1, 'A',
                'DBSCAN Spatial Clusters',
                f'eps={EPS}Å  |  ★ Primary  ◆ Secondary')

    # ── Panel B ──────────────────────────────────────────────
    cl_list = [l for l in sorted(cluster_info.keys()) if l >= 0]
    if -1 in cluster_info:
        cl_list.append(-1)

    x  = np.arange(len(cl_list))
    w  = 0.25
    os_ = 0.28

    for xi, lbl in enumerate(cl_list):
        info = cluster_info[lbl]
        col  = NOISE_COL if lbl == -1 else CLUSTER_COLORS[lbl % len(CLUSTER_COLORS)]
        fill = NOISE_FILL if lbl == -1 else CLUSTER_FILLS[lbl % len(CLUSTER_FILLS)]

        n = info['n_residues']
        ax2.bar(xi - os_, n, width=w,
                facecolor=fill, edgecolor=col)

        mc = info['mean_cons'] * 10
        ax2.bar(xi, mc, width=w,
                facecolor=COLORS['green_fill'],
                edgecolor=COLORS['green'])

        ax2.text(xi, -1.8, short_label(info['decision']),
                 ha='center', va='top', fontsize=7)

    ax2.set_xticks(x)
    ax2.set_xticklabels([f"C{l}" if l >= 0 else "Noise" for l in cl_list])

    for label in ax2.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(8)

    ax2.set_ylabel('Residue count / Mean conservation ×10',
                   fontsize=11, fontweight='bold')

    ax2.yaxis.grid(True, linestyle='--', linewidth=0.6, alpha=0.4)
    ax2.set_axisbelow(True)

    prism_axes(ax2)

    panel_title(ax2, 'B',
                'Cluster Summary',
                'Residue count + conservation')

    # ── Axis styling ─────────────────────────────────────────
    for ax in [ax1, ax2]:
        ax.spines['left'].set_linewidth(1.8)
        ax.spines['bottom'].set_linewidth(1.8)
        ax.tick_params(width=1.6, length=7)

    path = os.path.join(OUT_DIR,
        'Fig_Step04_DBSCAN_NSP12-NSP8.png')
    fig.savefig(path)
    plt.close()

    print(f"Figure saved to: {path}")


# ── Main ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Running Step 4 DBSCAN...')

    records = load_features()
    labels = run_dbscan(records)
    cluster_info = analyze_clusters(records, labels)

    make_figure(records, labels, cluster_info)

    print("Done.")

