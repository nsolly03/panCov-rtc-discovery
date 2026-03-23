#!/usr/bin/env python3
"""
CGCP Phase 2 Step 3 - Feature Classification: NSP12-NSP8
Classifies interface residues by pharmacophore feature type.
Dual anchor: LEU387 (primary, hydrophobic) + LYS332 (secondary, electrostatic)

Feature types:
  anchor_primary      — LEU387 (hydrophobic core)
  anchor_secondary    — LYS332 (electrostatic salt bridge)
  aromatic            — PHE, TYR, TRP, HIS
  hydrophobic         — LEU, VAL, ILE, ALA, MET, PRO, GLY
  charged_pos         — ARG, LYS
  charged_neg         — ASP, GLU
  hbond_donor         — SER, THR, ASN, GLN, CYS

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step03_feature_classification_NSP12-NSP8.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
    'axes.grid':          False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
})

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
S2_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-02-contacts')
OUT_DIR = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-03-features')
os.makedirs(OUT_DIR, exist_ok=True)

# Feature classification rules
FEATURE_MAP = {
    # Aromatic
    'PHE': 'aromatic', 'TYR': 'aromatic',
    'TRP': 'aromatic', 'HIS': 'aromatic',
    # Hydrophobic
    'LEU': 'hydrophobic', 'VAL': 'hydrophobic',
    'ILE': 'hydrophobic', 'ALA': 'hydrophobic',
    'MET': 'hydrophobic', 'PRO': 'hydrophobic',
    'GLY': 'hydrophobic',
    # Charged positive
    'ARG': 'charged_pos', 'LYS': 'charged_pos',
    # Charged negative
    'ASP': 'charged_neg', 'GLU': 'charged_neg',
    # H-bond donor/acceptor
    'SER': 'hbond_donor', 'THR': 'hbond_donor',
    'ASN': 'hbond_donor', 'GLN': 'hbond_donor',
    'CYS': 'hbond_donor',
}

# Dual anchors
PRIMARY_ANCHOR   = (387, 'NSP12')   # LEU387
SECONDARY_ANCHOR = (332, 'NSP12')   # LYS332

# Feature colors (Prism aesthetic)
FEAT_COLORS = {
    'anchor_primary':     '#1A1A1A',
    'anchor_secondary':   '#D7191C',
    'aromatic':           '#E6550D',
    'hydrophobic':        '#636363',
    'charged_pos':        '#2166AC',
    'charged_neg':        '#4DAC26',
    'hbond_donor':        '#8B4513',
    'unknown':            '#AAAAAA',
}
FEAT_FILLS = {
    'anchor_primary':     '#555555',
    'anchor_secondary':   '#FDAE61',
    'aromatic':           '#FDD0A2',
    'hydrophobic':        '#CCCCCC',
    'charged_pos':        '#92C5DE',
    'charged_neg':        '#A6D96A',
    'hbond_donor':        '#DEB887',
    'unknown':            '#DDDDDD',
}


def classify_feature(aa, position, chain):
    if (position, chain) == PRIMARY_ANCHOR:
        return 'anchor_primary'
    if (position, chain) == SECONDARY_ANCHOR:
        return 'anchor_secondary'
    return FEATURE_MAP.get(aa, 'unknown')


def load_contacts():
    path = os.path.join(S2_DIR,
        'contact_map_NSP12-NSP8.tsv')
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
        feat  = classify_feature(aa, pos, chain)

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

    # Sort: primary anchor, secondary anchor, then by composite
    classified.sort(key=lambda x: (
        0 if x['primary_feature'] == 'anchor_primary' else
        1 if x['primary_feature'] == 'anchor_secondary' else
        2, -x['composite']
    ))
    return classified


def save_outputs(classified):
    tsv_path = os.path.join(OUT_DIR,
        'feature_classification_NSP12-NSP8.tsv')
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

    # Count by feature
    feat_counts = {}
    for r in classified:
        ft = r['primary_feature']
        feat_counts[ft] = feat_counts.get(ft, 0) + 1

    json_path = os.path.join(OUT_DIR,
        'feature_classification_NSP12-NSP8.json')
    with open(json_path, 'w') as f:
        json.dump({
            'interface':        'NSP12-NSP8',
            'anchor_primary':   'LEU387',
            'anchor_secondary': 'LYS332',
            'n_total':          len(classified),
            'feature_counts':   feat_counts,
            'records':          classified,
        }, f, indent=2)
    print(f"  JSON: {json_path}")
    return feat_counts


def make_figure(classified, feat_counts):
    fig, axes = plt.subplots(1, 2,
                             figsize=(13, 5.5))
    ax1, ax2 = axes

    # ── Panel A: Feature distribution bar chart ───────────────
    feat_order = ['anchor_primary', 'anchor_secondary',
                  'aromatic', 'hydrophobic',
                  'charged_pos', 'charged_neg',
                  'hbond_donor', 'unknown']
    feat_labels = ['Anchor\n(primary)', 'Anchor\n(secondary)',
                   'Aromatic', 'Hydrophobic',
                   'Charged+', 'Charged-',
                   'H-bond\ndonor', 'Unknown']
    counts = [feat_counts.get(f, 0) for f in feat_order]
    cols   = [FEAT_COLORS[f] for f in feat_order]
    fills  = [FEAT_FILLS[f] for f in feat_order]

    x = np.arange(len(feat_order))
    for xi, val, fc, ec in zip(x, counts, fills, cols):
        if val == 0:
            continue
        ax1.bar(xi, val, width=0.5,
                facecolor=fc, edgecolor=ec,
                linewidth=0.9, zorder=3,
                clip_on=False)
        ax1.text(xi, val+0.1, str(val),
                 ha='center', va='bottom',
                 fontsize=8, color=ec,
                 fontweight='bold', clip_on=False)

    ax1.set_xticks(x)
    ax1.set_xticklabels(feat_labels, fontsize=7.5)
    ax1.set_ylabel('Number of residues', fontsize=9)
    ymax = max(counts) + 5
    ax1.spines['left'].set_position(('outward', 5))
    ax1.spines['bottom'].set_position(('outward', 5))
    ax1.spines['left'].set_bounds(0, max(counts))
    ax1.set_ylim(-1, ymax)
    ax1.set_yticks(range(0, max(counts)+5, 5))
    ax1.set_xlim(-0.6, len(feat_order)-0.4)
    ax1.set_title(
        'A   Feature distribution — NSP12-NSP8\n'
        '(dual anchor: LEU387 primary + LYS332 secondary)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel B: Composite score by feature ───────────────────
    for feat in feat_order:
        members = [r for r in classified
                   if r['primary_feature'] == feat]
        if not members:
            continue
        col  = FEAT_COLORS[feat]
        fill = FEAT_FILLS[feat]
        xi_vals = [feat_order.index(feat)] * len(members)
        comp_vals = [r['composite'] for r in members]

        # Jitter
        jitter = np.random.uniform(-0.15, 0.15,
                                   len(members))
        ax2.scatter(
            [feat_order.index(feat) + j
             for j in jitter],
            comp_vals,
            color=fill, edgecolors=col,
            linewidths=0.7, s=30,
            alpha=0.8, zorder=4)

        # Mean bar
        mean_comp = np.mean(comp_vals)
        ax2.plot(
            [feat_order.index(feat)-0.25,
             feat_order.index(feat)+0.25],
            [mean_comp, mean_comp],
            color=col, linewidth=2.0, zorder=5)

    ax2.set_xticks(range(len(feat_order)))
    ax2.set_xticklabels(feat_labels, fontsize=7.5)
    ax2.set_ylabel('CGCP composite score', fontsize=9)
    ax2.spines['left'].set_position(('outward', 5))
    ax2.spines['bottom'].set_position(('outward', 5))
    ax2.set_ylim(-0.05, 1.15)
    ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_title(
        'B   Composite score by feature type\n'
        '(dots = residues; line = mean)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    fig.tight_layout(pad=2.0)
    path = os.path.join(OUT_DIR,
        'Fig_Step03_FeatureClassification_NSP12-NSP8.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(classified, feat_counts):
    print('\n' + '='*60)
    print('STEP 3 FEATURE CLASSIFICATION — NSP12-NSP8')
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

    for feat in ['anchor_primary', 'anchor_secondary',
                 'aromatic', 'hydrophobic',
                 'charged_pos', 'charged_neg',
                 'hbond_donor']:
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
    print('CGCP Phase 2 Step 3 — NSP12-NSP8 Feature Classification')
    contacts   = load_contacts()
    classified = classify_all(contacts)
    feat_counts = save_outputs(classified)
    make_figure(classified, feat_counts)
    print_summary(classified, feat_counts)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Step 4 — DBSCAN clustering")
