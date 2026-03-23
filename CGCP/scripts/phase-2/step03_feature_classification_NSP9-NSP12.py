#!/usr/bin/env python3
"""
CGCP Phase 2 Step 3 - Feature Classification: NSP9-NSP12
Anchor: ARG733 (NSP12, cons=1.000)
NSP9: Chain G

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step03_feature_classification_NSP9-NSP12.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
from prism_style import (apply_prism, prism_axes,
                          set_xticklabels_vertical,
                          make_legend, panel_title,
                          COLORS)

apply_prism()

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
S2_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-02-contacts')
OUT_DIR = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-03-features')
os.makedirs(OUT_DIR, exist_ok=True)

FEATURE_MAP = {
    'PHE': 'aromatic',   'TYR': 'aromatic',
    'TRP': 'aromatic',   'HIS': 'aromatic',
    'LEU': 'hydrophobic','VAL': 'hydrophobic',
    'ILE': 'hydrophobic','ALA': 'hydrophobic',
    'MET': 'hydrophobic','PRO': 'hydrophobic',
    'GLY': 'hydrophobic',
    'ARG': 'charged_pos','LYS': 'charged_pos',
    'ASP': 'charged_neg','GLU': 'charged_neg',
    'SER': 'hbond_donor','THR': 'hbond_donor',
    'ASN': 'hbond_donor','GLN': 'hbond_donor',
    'CYS': 'hbond_donor',
}

ANCHOR_RES   = 733
ANCHOR_CHAIN = 'NSP12'

FEAT_COLORS = {
    'anchor':      COLORS['black'],
    'aromatic':    COLORS['red'],
    'hydrophobic': COLORS['gray'],
    'charged_pos': COLORS['blue'],
    'charged_neg': COLORS['green'],
    'hbond_donor': COLORS['hbond'],
    'unknown':     '#AAAAAA',
}
FEAT_FILLS = {
    'anchor':      '#555555',
    'aromatic':    COLORS['red_fill'],
    'hydrophobic': COLORS['gray_fill'],
    'charged_pos': COLORS['blue_fill'],
    'charged_neg': COLORS['green_fill'],
    'hbond_donor': COLORS['hbond_f'],
    'unknown':     '#DDDDDD',
}


def classify(aa, position, chain):
    if position == ANCHOR_RES and chain == ANCHOR_CHAIN:
        return 'anchor'
    return FEATURE_MAP.get(aa, 'unknown')


def load_contacts():
    path = os.path.join(S2_DIR,
        'contact_map_NSP9-NSP12.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            records.append(row)
    return records


def classify_all(contacts):
    classified = []
    for row in contacts:
        pos   = int(row['position'])
        chain = row['chain']
        aa    = row['aa']
        feat  = classify(aa, pos, chain)
        classified.append({
            'residue':         row['residue'],
            'chain':           chain,
            'position':        pos,
            'aa':              aa,
            'primary_feature': feat,
            'contact_score':   int(row['contact_score']),
            'n_partners':      int(row['n_partners']),
            'conservation':    float(row['conservation']),
            'composite':       float(row['composite']),
            'is_anchor':       int(row['is_anchor']),
            'ca_x': float(row['ca_x'])
                    if row['ca_x'] else None,
            'ca_y': float(row['ca_y'])
                    if row['ca_y'] else None,
            'ca_z': float(row['ca_z'])
                    if row['ca_z'] else None,
        })

    classified.sort(key=lambda x: (
        0 if x['primary_feature'] == 'anchor' else 1,
        -x['composite']))
    return classified


def save_outputs(classified):
    tsv_path = os.path.join(OUT_DIR,
        'feature_classification_NSP9-NSP12.tsv')
    fields = ['residue','chain','position','aa',
              'primary_feature','contact_score',
              'n_partners','conservation','composite',
              'is_anchor','ca_x','ca_y','ca_z']
    with open(tsv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(classified)
    print(f"  TSV: {tsv_path}")

    feat_counts = {}
    for r in classified:
        ft = r['primary_feature']
        feat_counts[ft] = feat_counts.get(ft, 0) + 1

    json_path = os.path.join(OUT_DIR,
        'feature_classification_NSP9-NSP12.json')
    with open(json_path, 'w') as f:
        json.dump({
            'interface':    'NSP9-NSP12',
            'anchor':       'ARG733',
            'n_total':      len(classified),
            'feature_counts': feat_counts,
            'records':      classified,
        }, f, indent=2)
    print(f"  JSON: {json_path}")
    return feat_counts


def make_figure(classified, feat_counts):
    fig = plt.figure(figsize=(13.0, 5.5))
    gs  = fig.add_gridspec(1, 2, wspace=0.45,
                           left=0.07, right=0.97,
                           top=0.88, bottom=0.26)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Panel A: feature distribution
    feat_order  = ['anchor','aromatic','hydrophobic',
                   'charged_pos','charged_neg',
                   'hbond_donor','unknown']
    feat_labels = ['Anchor','Aromatic','Hydrophobic',
                   'Charged+','Charged-',
                   'H-bond\ndonor','Unknown']
    counts = [feat_counts.get(f, 0) for f in feat_order]
    cols   = [FEAT_COLORS[f] for f in feat_order]
    fills  = [FEAT_FILLS[f]  for f in feat_order]

    x = np.arange(len(feat_order))
    for xi, val, fc, ec in zip(x, counts, fills, cols):
        if val == 0:
            continue
        ax1.bar(xi, val, width=0.50,
                facecolor=fc, edgecolor=ec,
                linewidth=1.2, zorder=3, clip_on=False)
        ax1.text(xi, val+0.15, str(val),
                 ha='center', va='bottom',
                 fontsize=9, fontweight='bold',
                 color=ec, clip_on=False)

    ax1.set_xticks(x)
    set_xticklabels_vertical(ax1, feat_labels, fontsize=9)
    ax1.set_ylabel('Number of residues',
                   fontsize=11, fontweight='bold')
    ymax = max(counts) + 5
    prism_axes(ax1, ymax=ymax,
               yticks=range(0, ymax+1, 5))
    ax1.set_xlim(-0.6, len(feat_order)-0.4)
    panel_title(ax1, 'A',
                'Feature distribution — NSP9-NSP12',
                f'total={sum(counts)} residues')

    # Panel B: composite by feature (scatter + mean)
    for feat in feat_order:
        members = [r for r in classified
                   if r['primary_feature'] == feat]
        if not members:
            continue
        col  = FEAT_COLORS[feat]
        fill = FEAT_FILLS[feat]
        xi   = feat_order.index(feat)
        comp_vals = [r['composite'] for r in members]
        jitter = np.random.uniform(-0.15, 0.15,
                                   len(members))
        ax2.scatter([xi+j for j in jitter],
                    comp_vals,
                    color=fill, edgecolors=col,
                    linewidths=0.8, s=35,
                    alpha=0.85, zorder=4)
        ax2.plot([xi-0.25, xi+0.25],
                 [np.mean(comp_vals)]*2,
                 color=col, linewidth=2.2, zorder=5)

    ax2.set_xticks(range(len(feat_order)))
    set_xticklabels_vertical(ax2, feat_labels, fontsize=9)
    ax2.set_ylabel('CGCP composite score',
                   fontsize=11, fontweight='bold')
    ax2.spines['left'].set_position(('outward', 6))
    ax2.spines['bottom'].set_position(('outward', 6))
    ax2.spines['left'].set_linewidth(1.2)
    ax2.spines['bottom'].set_linewidth(1.2)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylim(-0.05, 1.15)
    ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.tick_params(labelsize=9, width=1.2,
                    length=6, direction='out')
    panel_title(ax2, 'B',
                'Composite score by feature type',
                'dots = residues; line = mean')

    path = os.path.join(OUT_DIR,
        'Fig_Step03_FeatureClassification_NSP9-NSP12.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(classified, feat_counts):
    print('\n' + '='*60)
    print('STEP 3 FEATURE CLASSIFICATION — NSP9-NSP12')
    print('='*60)
    print(f"\n  Feature distribution:")
    for feat, count in sorted(feat_counts.items(),
                               key=lambda x: -x[1]):
        print(f"    {feat:<22} {count:>4} residues")

    print(f"\n  Top residues per feature:")
    shown = {}
    for r in classified:
        ft = r['primary_feature']
        if ft not in shown:
            shown[ft] = []
        if len(shown[ft]) < 3:
            shown[ft].append(r)

    for feat in ['anchor','aromatic','hydrophobic',
                 'charged_pos','charged_neg','hbond_donor']:
        members = shown.get(feat, [])
        if not members:
            continue
        print(f"\n  {feat}:")
        for r in members:
            print(f"    {r['aa']}{r['position']:<8} "
                  f"{r['chain']:<7} "
                  f"cons={r['conservation']:.3f} "
                  f"comp={r['composite']:.3f}")
    print('='*60)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 3 — NSP9-NSP12 Feature Classification')
    contacts    = load_contacts()
    classified  = classify_all(contacts)
    feat_counts = save_outputs(classified)
    make_figure(classified, feat_counts)
    print_summary(classified, feat_counts)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Step 4 — DBSCAN clustering")
