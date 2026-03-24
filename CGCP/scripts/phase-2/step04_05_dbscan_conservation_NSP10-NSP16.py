#!/usr/bin/env python3
"""
CGCP Phase 2 Steps 4+5 - DBSCAN + Conservation Overlay: NSP10-NSP16
Same method as NSP9-NSP12 steps 4+5.
Anchor: LYS93 (NSP10, B4346, polyprotein numbering)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step04_05_dbscan_conservation_NSP10-NSP16.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import DBSCAN
from prism_style import (apply_prism, prism_axes,
                          set_xticklabels_vertical,
                          make_legend, panel_title,
                          COLORS, CLUSTER_COLORS, CLUSTER_FILLS)

apply_prism()

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):  return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray):  return obj.tolist()
        if isinstance(obj, np.bool_):    return bool(obj)
        return super().default(obj)

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
S3_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP10-NSP16/step-03-features')
OUT4    = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP10-NSP16/step-04-clusters')
OUT5    = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP10-NSP16/step-05-conservation')
os.makedirs(OUT4, exist_ok=True)
os.makedirs(OUT5, exist_ok=True)

ANCHOR_RESNUM = 4346   # B4346 = LYS93 in NSP10
EPS           = 8.0
MIN_SAMPLES   = 2
NOISE_COL     = COLORS['gray']
NOISE_FILL    = COLORS['gray_fill']

TIER_COLORS = {
    'identical': '#D7191C',
    'high':      '#F46D43',
    'moderate':  '#FEE090',
    'variable':  '#AAAAAA',
}
TIER_FILLS = {
    'identical': '#FDAE61',
    'high':      '#FDD0A2',
    'moderate':  '#FFFFBF',
    'variable':  '#DDDDDD',
}


def get_tier(cons):
    if cons == 1.000:  return 'identical'
    if cons >= 0.800:  return 'high'
    if cons >= 0.600:  return 'moderate'
    return 'variable'


def load_features():
    path = os.path.join(S3_DIR,
        'feature_classification_NSP10-NSP16.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['ca_x'] and row['ca_y'] and row['ca_z']:
                records.append({
                    'residue':         row['residue'],
                    'chain':           row['chain'],
                    'position':        int(row['position']),
                    'prot_pos':        int(row['prot_pos']),
                    'aa':              row['aa'],
                    'primary_feature': row['primary_feature'],
                    'contact_score':   int(row['contact_score']),
                    'conservation':    float(row['conservation']),
                    'composite':       float(row['composite']),
                    'is_anchor':       int(row['is_anchor']),
                    'x': float(row['ca_x']),
                    'y': float(row['ca_y']),
                    'z': float(row['ca_z']),
                })
    return records


def run_dbscan(records):
    coords = np.array([[r['x'], r['y'], r['z']]
                       for r in records])
    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES,
                metric='euclidean').fit(coords)
    return [int(l) for l in db.labels_]


def analyze_clusters(records, labels):
    cluster_info = {}
    for label in sorted(set(labels)):
        members = [records[i] for i, l
                   in enumerate(labels) if l == label]
        if not members:
            continue
        coords    = np.array([[r['x'],r['y'],r['z']]
                               for r in members])
        centroid  = coords.mean(axis=0)
        radius    = float(max(
            np.linalg.norm(c - centroid) for c in coords))
        mean_cons = float(np.mean(
            [r['conservation'] for r in members]))
        mean_comp = float(np.mean(
            [r['composite'] for r in members]))
        has_anc   = bool(any(
            r['primary_feature'] == 'anchor'
            for r in members))

        if   has_anc and mean_cons >= 0.700:
            decision = 'PRIMARY PHARMACOPHORE'
        elif mean_cons >= 0.800 and mean_comp >= 0.500:
            decision = 'SECONDARY PHARMACOPHORE'
        elif mean_cons >= 0.600:
            decision = 'SUPPORTING'
        else:
            decision = 'DEPRIORITIZE'

        cluster_info[label] = {
            'label':      int(label),
            'n_residues': int(len(members)),
            'centroid_x': round(float(centroid[0]), 3),
            'centroid_y': round(float(centroid[1]), 3),
            'centroid_z': round(float(centroid[2]), 3),
            'radius':     round(radius, 2),
            'mean_cons':  round(mean_cons, 3),
            'mean_comp':  round(mean_comp, 3),
            'has_anchor': has_anc,
            'decision':   decision,
            'members':    [r['residue'] for r in members],
        }
    return cluster_info


def save_step4(records, labels, cluster_info):
    rows = [{
        'residue':         r['residue'],
        'chain':           r['chain'],
        'position':        int(r['position']),
        'prot_pos':        int(r['prot_pos']),
        'aa':              r['aa'],
        'primary_feature': r['primary_feature'],
        'conservation':    float(r['conservation']),
        'composite':       float(r['composite']),
        'cluster':         int(labels[i]),
        'x':               float(r['x']),
        'y':               float(r['y']),
        'z':               float(r['z']),
    } for i, r in enumerate(records)]

    tsv = os.path.join(OUT4, 'clusters_NSP10-NSP16.tsv')
    with open(tsv, 'w', newline='') as f:
        w = csv.DictWriter(f,
            fieldnames=list(rows[0].keys()),
            delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f"  Step4 TSV: {tsv}")

    json_out = {
        'interface':   'NSP10-NSP16',
        'anchor':      'LYS93(B4346)',
        'dbscan_eps':  float(EPS),
        'min_samples': int(MIN_SAMPLES),
        'n_clusters':  int(len([l for l in cluster_info if l >= 0])),
        'n_noise':     int(sum(1 for l in labels if l == -1)),
        'clusters':    {str(k): v for k,v in cluster_info.items()},
        'residues':    rows,
    }
    jpath = os.path.join(OUT4, 'clusters_NSP10-NSP16.json')
    with open(jpath, 'w') as f:
        json.dump(json_out, f, indent=2, cls=NpEncoder)
    print(f"  Step4 JSON: {jpath}")
    return rows


def save_step5(records, labels):
    result = [{**r, 'cons_tier': get_tier(r['conservation'])}
              for r in records]
    result.sort(key=lambda x: (
        0 if x['primary_feature']=='anchor' else 1,
        -x['conservation'], -x['composite']))

    rows_out = [{**r, 'cluster': int(labels[i])}
                for i, r in enumerate(records)]
    rows_out.sort(key=lambda x: (
        0 if x['primary_feature']=='anchor' else 1,
        -x['conservation'], -x['composite']))

    tsv = os.path.join(OUT5,
        'conservation_overlay_NSP10-NSP16.tsv')
    fields = ['residue','chain','position','prot_pos','aa',
              'primary_feature','cluster',
              'conservation','cons_tier','composite',
              'x','y','z']
    with open(tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(rows_out)
    print(f"  Step5 TSV: {tsv}")

    tier_stats = {}
    for tier in ['identical','high','moderate','variable']:
        members = [r for r in result if r['cons_tier']==tier]
        if members:
            tier_stats[tier] = {
                'n':         len(members),
                'mean_comp': round(float(np.mean(
                    [r['composite'] for r in members])), 3),
                'members':   [r['residue'] for r in members],
                'nsp16':     [r['residue'] for r in members
                              if r['chain']=='NSP16'],
                'nsp10':     [r['residue'] for r in members
                              if r['chain']=='NSP10'],
            }

    jpath = os.path.join(OUT5,
        'conservation_overlay_NSP10-NSP16.json')
    with open(jpath, 'w') as f:
        json.dump({
            'interface':  'NSP10-NSP16',
            'anchor':     'LYS93(B4346)',
            'tier_stats': tier_stats,
            'records':    rows_out,
        }, f, indent=2, cls=NpEncoder)
    print(f"  Step5 JSON: {jpath}")
    return result, tier_stats


def make_figures(records, labels, cluster_info,
                 result, tier_stats):
    # ── FIGURE 1: DBSCAN ─────────────────────────────────────
    fig1, axes1 = plt.subplots(1, 2, figsize=(13, 5.5))
    ax1, ax2 = axes1
    fig1.subplots_adjust(wspace=0.45, left=0.07,
                         right=0.97, top=0.88, bottom=0.22)

    plotted = set()
    for i, r in enumerate(records):
        lbl = labels[i]
        if lbl == -1:
            col, fill = NOISE_COL, NOISE_FILL
            lbl_label = ('Noise' if 'noise' not in plotted
                         else None)
            plotted.add('noise')
        else:
            idx  = lbl % len(CLUSTER_COLORS)
            col  = CLUSTER_COLORS[idx]
            fill = CLUSTER_FILLS[idx]
            ck   = f'C{lbl}'
            dec  = cluster_info[lbl]['decision']
            lbl_label = (f"C{lbl}: {dec[:18]}"
                         if ck not in plotted else None)
            plotted.add(ck)

        is_anc = r['primary_feature'] == 'anchor'
        size   = 220 if is_anc else 45
        marker = '*'  if is_anc else 'o'
        ax1.scatter(r['x'], r['z'], color=fill,
                    edgecolors=col, linewidths=1.0,
                    s=size, marker=marker,
                    alpha=0.85, zorder=4, label=lbl_label)
        if is_anc or r['composite'] >= 0.750:
            ax1.annotate(
                r['residue'].split('-')[1],
                (r['x'], r['z']),
                fontsize=7.5, color=col,
                fontweight='bold' if is_anc else 'normal',
                xytext=(4, 3),
                textcoords='offset points')

    for lbl, info in cluster_info.items():
        if lbl < 0: continue
        idx = lbl % len(CLUSTER_COLORS)
        ax1.scatter(info['centroid_x'], info['centroid_z'],
                    marker='+', color=CLUSTER_COLORS[idx],
                    s=180, linewidths=2.5, zorder=6)

    ax1.set_xlabel('X coordinate (Å)',
                   fontsize=11, fontweight='bold')
    ax1.set_ylabel('Z coordinate (Å)',
                   fontsize=11, fontweight='bold')
    ax1.spines['left'].set_position(('outward', 6))
    ax1.spines['bottom'].set_position(('outward', 6))
    ax1.spines['left'].set_linewidth(1.2)
    ax1.spines['bottom'].set_linewidth(1.2)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(labelsize=9, width=1.2,
                    length=6, direction='out')
    ax1.legend(fontsize=7, frameon=True,
               facecolor='white', edgecolor='#CCCCCC',
               loc='best', framealpha=0.9)
    panel_title(ax1, 'A', 'DBSCAN clusters (X-Z)',
                f'eps={EPS}Å | ★=LYS93 anchor')

    cl_list = [l for l in sorted(cluster_info.keys())
               if l >= 0]
    if -1 in cluster_info:
        cl_list.append(-1)
    x = np.arange(len(cl_list))

    for xi, lbl in enumerate(cl_list):
        info = cluster_info[lbl]
        col  = (NOISE_COL if lbl == -1
                else CLUSTER_COLORS[lbl % len(CLUSTER_COLORS)])
        fill = (NOISE_FILL if lbl == -1
                else CLUSTER_FILLS[lbl % len(CLUSTER_FILLS)])

        ax2.bar(xi-0.30, info['n_residues'], width=0.28,
                facecolor=fill, edgecolor=col,
                linewidth=1.2, zorder=3, clip_on=False)
        ax2.text(xi-0.30, info['n_residues']+0.3,
                 str(info['n_residues']),
                 ha='center', va='bottom',
                 fontsize=8, fontweight='bold',
                 color=col, clip_on=False)

        ax2.bar(xi, info['mean_cons']*10, width=0.28,
                facecolor=COLORS['green_fill'],
                edgecolor=COLORS['green'],
                linewidth=1.2, zorder=3, clip_on=False)
        ax2.text(xi, info['mean_cons']*10+0.3,
                 f"{info['mean_cons']:.2f}",
                 ha='center', va='bottom',
                 fontsize=8, fontweight='bold',
                 color=COLORS['green'], clip_on=False)

        ax2.text(xi, -1.8, info['decision'][:12],
                 ha='center', va='top',
                 fontsize=6, color=col, clip_on=False)

    ax2.set_xticks(x)
    set_xticklabels_vertical(
        ax2, [f"C{l}" if l >= 0 else "Noise"
              for l in cl_list], fontsize=9)
    ax2.set_ylabel('Count / Conservation ×10',
                   fontsize=11, fontweight='bold')
    ymax = max(cluster_info[l]['n_residues']
               for l in cl_list) + 5
    prism_axes(ax2, ymax=ymax,
               yticks=list(range(0, ymax, 5)))
    ax2.set_xlim(-0.6, len(cl_list)-0.4)
    ax2.tick_params(labelsize=9, width=1.2,
                    length=6, direction='out')
    make_legend(ax2, [
        (COLORS['gray_fill'],  COLORS['gray'],
         'Residue count'),
        (COLORS['green_fill'], COLORS['green'],
         'Mean conservation ×10'),
    ], loc='upper right', fontsize=8)
    panel_title(ax2, 'B', 'Cluster summary')

    p1 = os.path.join(OUT4,
        'Fig_Step04_DBSCAN_NSP10-NSP16.png')
    fig1.savefig(p1)
    plt.close()
    print(f"  Fig DBSCAN: {p1}")

    # ── FIGURE 2: Conservation overlay ───────────────────────
    fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5.5))
    ax_a, ax_b, ax_c = axes2
    fig2.subplots_adjust(wspace=0.48, left=0.06,
                         right=0.97, top=0.88, bottom=0.26)

    top = sorted(result, key=lambda x: -x['composite'])[:30]
    labels_a = [r['residue'].split('-')[1] for r in top]
    x = np.arange(len(top))

    for xi, r in enumerate(top):
        ec = TIER_COLORS[r['cons_tier']]
        fc = TIER_FILLS[r['cons_tier']]
        is_anc = r['primary_feature'] == 'anchor'
        lw = 2.0 if is_anc else 1.2
        ax_a.bar(xi, r['conservation'], width=0.40,
                 facecolor=fc, edgecolor=ec,
                 linewidth=lw, zorder=3, clip_on=False)
        if is_anc:
            ax_a.text(xi, r['conservation']+0.02, '★',
                      ha='center', va='bottom',
                      fontsize=10, color=ec, clip_on=False)

    for thresh, lbl, col in [
        (1.000, 'Identical', TIER_COLORS['identical']),
        (0.800, 'High',      TIER_COLORS['high']),
        (0.600, 'Moderate',  TIER_COLORS['moderate']),
    ]:
        ax_a.axhline(thresh, color=col, linewidth=0.8,
                     linestyle='--', alpha=0.6)
        ax_a.text(len(top)-0.5, thresh+0.01, lbl,
                  fontsize=7, color=col,
                  ha='right', va='bottom', clip_on=False)

    ax_a.set_xticks(x)
    set_xticklabels_vertical(ax_a, labels_a, fontsize=7)
    ax_a.set_ylabel('Conservation score',
                    fontsize=11, fontweight='bold')
    ax_a.set_xlim(-0.6, len(top)-0.4)
    prism_axes(ax_a, ymax=1.1,
               yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    make_legend(ax_a, [
        (TIER_FILLS['identical'], TIER_COLORS['identical'],
         'Identical (1.000)'),
        (TIER_FILLS['high'],      TIER_COLORS['high'],
         'High (≥0.800)'),
        (TIER_FILLS['moderate'],  TIER_COLORS['moderate'],
         'Moderate (≥0.600)'),
        (TIER_FILLS['variable'],  TIER_COLORS['variable'],
         'Variable (<0.600)'),
    ], loc='lower left', fontsize=7.5)
    panel_title(ax_a, 'A',
                'Conservation scores (top 30)',
                '★=LYS93 anchor')

    tier_order  = ['identical','high','moderate','variable']
    tier_labels = ['Identical','High','Moderate','Variable']
    tier_counts = [tier_stats.get(t, {}).get('n', 0)
                   for t in tier_order]
    x2 = np.arange(len(tier_order))
    for xi, val, t in zip(x2, tier_counts, tier_order):
        ec = TIER_COLORS[t]
        fc = TIER_FILLS[t]
        ax_b.bar(xi, val, width=0.50,
                 facecolor=fc, edgecolor=ec,
                 linewidth=1.2, zorder=3, clip_on=False)
        ax_b.text(xi, val+0.3, str(val),
                  ha='center', va='bottom',
                  fontsize=10, fontweight='bold',
                  color=ec, clip_on=False)
    ax_b.set_xticks(x2)
    set_xticklabels_vertical(ax_b, tier_labels, fontsize=9)
    ax_b.set_ylabel('Number of residues',
                    fontsize=11, fontweight='bold')
    ymax_b = max(tier_counts)+6
    prism_axes(ax_b, ymax=ymax_b,
               yticks=range(0, ymax_b, 5))
    ax_b.set_xlim(-0.6, len(tier_order)-0.4)
    panel_title(ax_b, 'B', 'Conservation tier distribution',
                f'total={sum(tier_counts)} residues')

    for r in result:
        ec = TIER_COLORS[r['cons_tier']]
        fc = TIER_FILLS[r['cons_tier']]
        is_anc = r['primary_feature'] == 'anchor'
        size   = 220 if is_anc else 35
        marker = '*'  if is_anc else 'o'
        ax_c.scatter(r['composite'], r['conservation'],
                     color=fc, edgecolors=ec,
                     linewidths=0.9, s=size,
                     marker=marker, alpha=0.85, zorder=4)
        if is_anc or (r['conservation'] == 1.000
                      and r['composite'] >= 0.700):
            ax_c.annotate(
                r['residue'].split('-')[1],
                (r['composite'], r['conservation']),
                fontsize=7, color=ec,
                xytext=(4, 2),
                textcoords='offset points')

    ax_c.axhline(0.800, color=TIER_COLORS['high'],
                 linewidth=0.8, linestyle='--', alpha=0.6)
    ax_c.axvline(0.600, color=COLORS['gray'],
                 linewidth=0.8, linestyle='--', alpha=0.6)
    ax_c.set_xlabel('CGCP composite score',
                    fontsize=11, fontweight='bold')
    ax_c.set_ylabel('Conservation score',
                    fontsize=11, fontweight='bold')
    ax_c.spines['left'].set_position(('outward', 6))
    ax_c.spines['bottom'].set_position(('outward', 6))
    ax_c.spines['left'].set_linewidth(1.2)
    ax_c.spines['bottom'].set_linewidth(1.2)
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)
    ax_c.set_xlim(-0.05, 1.10)
    ax_c.set_ylim(-0.05, 1.10)
    ax_c.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_c.tick_params(labelsize=9, width=1.2,
                     length=6, direction='out')
    panel_title(ax_c, 'C',
                'Composite vs conservation',
                'dashed = thresholds')

    p2 = os.path.join(OUT5,
        'Fig_Step05_ConservationOverlay_NSP10-NSP16.png')
    fig2.savefig(p2)
    plt.close()
    print(f"  Fig Conservation: {p2}")


def print_summary(cluster_info, labels, result, tier_stats):
    print('\n' + '='*65)
    print('STEPS 4+5: DBSCAN + CONSERVATION — NSP10-NSP16')
    print('='*65)

    n_clusters = len([l for l in cluster_info if l >= 0])
    n_noise    = sum(1 for l in labels if l == -1)
    print(f"\n  DBSCAN: {n_clusters} clusters | "
          f"{n_noise} noise")
    for lbl in sorted(cluster_info.keys()):
        info = cluster_info[lbl]
        name = f"Cluster {lbl}" if lbl >= 0 else "Noise"
        print(f"  {name}: {info['n_residues']} res | "
              f"cons={info['mean_cons']:.3f} | "
              f"comp={info['mean_comp']:.3f} | "
              f"{info['decision']}")
        mem = ', '.join(info['members'][:5])
        if len(info['members']) > 5:
            mem += f"... (+{len(info['members'])-5})"
        print(f"    → {mem}")

    print(f"\n  Conservation tiers:")
    for tier in ['identical','high','moderate','variable']:
        if tier in tier_stats:
            info = tier_stats[tier]
            print(f"    {tier:<12} n={info['n']:>3} "
                  f"mean_comp={info['mean_comp']:.3f}")

    candidates = [r for r in result
                  if r['conservation'] >= 0.800
                  and r['composite'] >= 0.500]
    candidates.sort(key=lambda x: -x['composite'])
    print(f"\n  Pharmacophore preview "
          f"(cons>=0.800, comp>=0.500): "
          f"{len(candidates)}")
    for r in candidates[:12]:
        sym = '★' if r['primary_feature'] == 'anchor' else ' '
        print(f"  {sym} {r['residue']:<24} "
              f"{r['chain']:<7} "
              f"cons={r['conservation']:.3f} "
              f"comp={r['composite']:.3f} "
              f"{r['cons_tier']}")
    print('='*65)


if __name__ == '__main__':
    print('CGCP Phase 2 Steps 4+5 — NSP10-NSP16')
    records      = load_features()
    labels       = run_dbscan(records)
    cluster_info = analyze_clusters(records, labels)
    rows         = save_step4(records, labels, cluster_info)
    result, tier_stats = save_step5(records, labels)
    make_figures(records, labels, cluster_info,
                 result, tier_stats)
    print_summary(cluster_info, labels, result, tier_stats)
    print(f"\nStep4: {OUT4}")
    print(f"Step5: {OUT5}")
    print("Next: Steps 6+7 — integrated assessment + pharmacophore")
