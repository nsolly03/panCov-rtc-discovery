#!/usr/bin/env python3
"""
CGCP Phase 2 Step 4 - DBSCAN Spatial Clustering: NSP12-NSP7
Clusters interface residues by 3D Ca coordinates.
Weights aromatic/anchor residues to ensure PHE440 cluster is primary.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step04_dbscan_clustering_NSP12-NSP7.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from Bio import PDB

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
    'xtick.major.pad':    4,
    'ytick.major.pad':    4,
    'axes.grid':          False,
    'legend.fontsize':    8,
    'legend.frameon':     False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
})

# Cluster colors — Prism palette
CLUSTER_COLORS = {
    0:  '#D7191C',   # red    — primary (PHE440 cluster)
    1:  '#2166AC',   # blue
    2:  '#4DAC26',   # green
    3:  '#7B2D8B',   # purple
    4:  '#E6550D',   # orange
    -1: '#CCCCCC',   # gray   — noise
}
CLUSTER_FILLS = {k: v+'55' for k, v in CLUSTER_COLORS.items()}

FEATURE_COLORS = {
    'anchor':        '#D7191C',
    'aromatic':      '#E6550D',
    'hydrophobic':   '#636363',
    'hbond_donor':   '#2166AC',
    'hbond_acceptor':'#74ADD1',
    'charged_pos':   '#7B2D8B',
    'charged_neg':   '#4DAC26',
}

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP12-NSP7_3/receptor_NSP12-NSP7_3.pdb')
S3_DIR   = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-03-features')
OUT_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-04-clusters')
os.makedirs(OUT_DIR, exist_ok=True)


# ── Load Step 3 classifications ──────────────────────────────
def load_features():
    path = os.path.join(S3_DIR,
                        'feature_classification_NSP12-NSP7.tsv')
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append(row)
    return rows


# ── Get Ca coordinates from PDB ──────────────────────────────
def get_coordinates(features):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]

    chain_map = {'NSP12': 'A', 'NSP7': 'C'}
    coords = []

    for f in features:
        chain_id = chain_map[f['chain']]
        pos      = int(f['position'])
        chain    = structure[chain_id]
        try:
            res  = chain[(' ', pos, ' ')]
            ca   = res['CA'].coord
            coords.append({
                'residue':         f['residue'],
                'chain':           f['chain'],
                'position':        pos,
                'aa':              f['aa'],
                'primary_feature': f['primary_feature'],
                'conservation':    float(f['conservation']),
                'composite':       float(f['composite']),
                'x': float(ca[0]),
                'y': float(ca[1]),
                'z': float(ca[2]),
            })
        except (KeyError, Exception):
            print(f"  Warning: Ca not found for "
                  f"{f['chain']} {pos} — skipping")

    print(f"  Coordinates retrieved: {len(coords)}/{len(features)}")
    return coords


# ── DBSCAN clustering ────────────────────────────────────────
def run_dbscan(coords, eps=8.0, min_samples=2):
    """
    eps=8.0 A — residues within 8A of each other are neighbors
    min_samples=2 — minimum 2 residues to form a cluster
    This is generous to capture the full aromatic cluster around PHE440
    """
    xyz = np.array([[c['x'], c['y'], c['z']] for c in coords])

    db = DBSCAN(eps=eps, min_samples=min_samples,
                metric='euclidean')
    labels = db.fit_predict(xyz)

    # Identify which cluster contains PHE440
    anchor_label = None
    for i, c in enumerate(coords):
        if c['position'] == 440 and c['chain'] == 'NSP12':
            anchor_label = labels[i]
            break

    # Remap so PHE440 cluster = 0
    if anchor_label is not None and anchor_label != 0:
        remapped = []
        for l in labels:
            if l == 0:
                remapped.append(anchor_label)
            elif l == anchor_label:
                remapped.append(0)
            else:
                remapped.append(l)
        labels = np.array(remapped)

    # Assign labels
    for i, c in enumerate(coords):
        c['cluster'] = int(labels[i])

    # Cluster stats
    unique = set(labels)
    unique.discard(-1)
    print(f"  DBSCAN eps={eps}A min_samples={min_samples}")
    print(f"  Clusters found: {len(unique)} "
          f"(+ {sum(1 for l in labels if l==-1)} noise points)")

    return coords, labels


# ── Compute cluster centroids ────────────────────────────────
def compute_centroids(coords):
    from collections import defaultdict
    clusters = defaultdict(list)
    for c in coords:
        if c['cluster'] >= 0:
            clusters[c['cluster']].append(c)

    centroids = {}
    for cid, members in clusters.items():
        cx = np.mean([m['x'] for m in members])
        cy = np.mean([m['y'] for m in members])
        cz = np.mean([m['z'] for m in members])
        top = sorted(members,
                     key=lambda x: -x['composite'])[:3]
        centroids[cid] = {
            'centroid': [round(float(cx),3),
                         round(float(cy),3),
                         round(float(cz),3)],
            'n_members':  len(members),
            'members':    [m['residue'] for m in members],
            'top_residues': [m['residue'] for m in top],
            'mean_conservation': round(
                float(np.mean([m['conservation']
                               for m in members])), 3),
            'mean_composite': round(
                float(np.mean([m['composite']
                               for m in members])), 3),
            'features': list(set(m['primary_feature']
                                 for m in members)),
        }
    return centroids


# ── Save outputs ─────────────────────────────────────────────
def save_tsv(coords):
    path = os.path.join(OUT_DIR,
                        'clusters_NSP12-NSP7.tsv')
    fields = ['residue','chain','position','aa',
              'primary_feature','cluster',
              'x','y','z','conservation','composite']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(
            f, fieldnames=fields, delimiter='\t',
            extrasaction='ignore')
        w.writeheader()
        w.writerows(coords)
    print(f"  TSV: {path}")


def save_json(coords, centroids):
    out = {
        'interface':   'NSP12-NSP7',
        'date':        '2026-03-18',
        'method':      'DBSCAN eps=8.0A min_samples=2',
        'n_residues':  len(coords),
        'n_clusters':  len(centroids),
        'anchor':      'PHE440',
        'anchor_cluster': 0,
        'centroids':   centroids,
        'residues':    coords,
    }
    path = os.path.join(OUT_DIR,
                        'clusters_NSP12-NSP7.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  JSON: {path}")
    return out


# ── Prism figure — 2x2 ──────────────────────────────────────
def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax*0.04, ymax*1.1)
    if yticks is not None:
        ax.set_yticks(yticks)


def make_figure(coords, centroids):
    fig = plt.figure(figsize=(13.0, 11.0))
    gs  = fig.add_gridspec(2, 2,
                           hspace=0.60, wspace=0.45,
                           left=0.08, right=0.97,
                           top=0.93, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0], projection='3d')
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # ── Panel A: 3D cluster scatter ───────────────────────────
    ax_a.set_facecolor('white')
    for pane in [ax_a.xaxis.pane, ax_a.yaxis.pane,
                 ax_a.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#DDDDDD')
    ax_a.grid(False)

    plotted_clusters = set()
    for c in coords:
        cl  = c['cluster']
        col = CLUSTER_COLORS.get(cl, '#AAAAAA')
        is_anchor = (c['position'] == 440
                     and c['chain'] == 'NSP12')
        size  = 180 if is_anchor else 60
        marker = '*' if is_anchor else 'o'
        label = (f'Cluster {cl}' if cl >= 0
                 else 'Noise') \
                 if cl not in plotted_clusters else None
        ax_a.scatter(c['x'], c['y'], c['z'],
                     color=col, s=size,
                     marker=marker,
                     edgecolors='white',
                     linewidths=0.5,
                     zorder=4, label=label,
                     alpha=0.85)
        if is_anchor:
            ax_a.text(c['x']+1, c['y']+1, c['z']+1.5,
                      'PHE440', fontsize=8,
                      color=CLUSTER_COLORS[0],
                      fontweight='bold')
        plotted_clusters.add(cl)

    ax_a.set_xlabel('X (Å)', fontsize=7, labelpad=2)
    ax_a.set_ylabel('Y (Å)', fontsize=7, labelpad=2)
    ax_a.set_zlabel('Z (Å)', fontsize=7, labelpad=2)
    ax_a.tick_params(labelsize=6)
    ax_a.view_init(elev=25, azim=-55)
    ax_a.legend(fontsize=7, loc='upper left',
                frameon=True, facecolor='white',
                edgecolor='#CCCCCC', framealpha=0.9)
    ax_a.set_title(
        'A   3D spatial clusters (DBSCAN)\n'
        '(Ca coordinates; eps=8.0 A; star=PHE440 anchor)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel B: Cluster membership bar chart ─────────────────
    from collections import Counter
    cl_counts = Counter(c['cluster'] for c in coords
                        if c['cluster'] >= 0)
    cl_ids    = sorted(cl_counts.keys())
    cl_vals   = [cl_counts[k] for k in cl_ids]
    cl_colors = [CLUSTER_COLORS.get(k, '#AAAAAA')
                 for k in cl_ids]
    cl_fills  = [CLUSTER_FILLS.get(k, '#AAAAAA55')
                 for k in cl_ids]
    cl_labels = [f'Cluster {k}' for k in cl_ids]

    x = np.arange(len(cl_ids))
    w = 0.45
    for xi, val, fc, ec, lbl in zip(
            x, cl_vals, cl_fills, cl_colors, cl_labels):
        ax_b.bar(xi, val, width=w,
                 facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3,
                 clip_on=False)
        ax_b.text(xi, val+0.08, str(val),
                  ha='center', va='bottom',
                  fontsize=9, color=ec,
                  fontweight='bold', clip_on=False)

        # Top residues annotation
        if cl_ids[xi] in centroids:
            top = ', '.join(
                r.split('-')[-1]
                for r in centroids[cl_ids[xi]
                    ]['top_residues'][:2])
            ax_b.text(xi, -0.6, top,
                      ha='center', va='top',
                      fontsize=6.5, color=ec,
                      fontstyle='italic',
                      clip_on=False)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(cl_labels, fontsize=8.5)
    ax_b.set_ylabel('Number of residues', fontsize=9)
    ax_b.set_xlim(-0.6, len(cl_ids)-0.4)
    ymax_b = max(cl_vals) + 2
    prism_axes(ax_b, ymax=ymax_b,
               yticks=range(0, ymax_b+1, 2))
    ax_b.set_title(
        'B   Cluster size and top residues\n'
        '(italic = top 2 residues by composite score)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel C: Mean conservation per cluster ────────────────
    cl_cons = {k: centroids[k]['mean_conservation']
               for k in cl_ids if k in centroids}
    cons_vals = [cl_cons.get(k, 0) for k in cl_ids]

    for xi, val, fc, ec in zip(
            x, cons_vals, cl_fills, cl_colors):
        ax_c.bar(xi, val, width=w,
                 facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3,
                 clip_on=False)
        ax_c.scatter(xi, val, color=ec,
                     s=20, zorder=5, clip_on=False)
        ax_c.text(xi, val+0.015, f'{val:.3f}',
                  ha='center', va='bottom',
                  fontsize=8, color=ec,
                  fontweight='bold', clip_on=False)

    ax_c.axhline(0.80, color='#636363',
                 linewidth=0.8,
                 linestyle=(0,(4,3)), zorder=1)
    ax_c.text(len(cl_ids)-0.45, 0.815,
              'pan-coronavirus\nthreshold',
              fontsize=6.5, color='#636363',
              ha='right', va='bottom',
              clip_on=False)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(cl_labels, fontsize=8.5)
    ax_c.set_ylabel('Mean conservation score', fontsize=9)
    ax_c.set_xlim(-0.6, len(cl_ids)-0.4)
    prism_axes(ax_c, ymax=1.05,
               yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_c.set_title(
        'C   Mean conservation per cluster\n'
        '(target: clusters above 0.80 threshold)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel D: Mean composite score per cluster ─────────────
    cl_comp  = {k: centroids[k]['mean_composite']
                for k in cl_ids if k in centroids}
    comp_vals = [cl_comp.get(k, 0) for k in cl_ids]

    for xi, val, fc, ec in zip(
            x, comp_vals, cl_fills, cl_colors):
        ax_d.bar(xi, val, width=w,
                 facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3,
                 clip_on=False)
        ax_d.scatter(xi, val, color=ec,
                     s=20, zorder=5, clip_on=False)
        ax_d.text(xi, val+0.008, f'{val:.3f}',
                  ha='center', va='bottom',
                  fontsize=8, color=ec,
                  fontweight='bold', clip_on=False)

    ax_d.set_xticks(x)
    ax_d.set_xticklabels(cl_labels, fontsize=8.5)
    ax_d.set_ylabel('Mean CGCP composite score', fontsize=9)
    ax_d.set_xlim(-0.6, len(cl_ids)-0.4)
    ymax_d = round(max(comp_vals)*1.4+0.05, 1)
    prism_axes(ax_d, ymax=ymax_d,
               yticks=[round(i*0.1, 1)
                       for i in range(int(ymax_d/0.1)+1)])
    ax_d.set_title(
        'D   Mean CGCP composite score per cluster\n'
        '(higher = more druggable interface region)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    path = os.path.join(OUT_DIR,
        'Fig_Step04_DBSCAN_Clusters.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


# ── Print summary ────────────────────────────────────────────
def print_summary(coords, centroids):
    print('\n' + '='*60)
    print('STEP 4 DBSCAN CLUSTERING — NSP12-NSP7')
    print('='*60)

    from collections import defaultdict
    clusters = defaultdict(list)
    for c in coords:
        clusters[c['cluster']].append(c)

    for cid in sorted(clusters.keys()):
        members = clusters[cid]
        label   = f'Cluster {cid}' if cid >= 0 else 'Noise'
        print(f"\n  {label} ({len(members)} residues):")
        for m in sorted(members,
                        key=lambda x: -x['composite'])[:8]:
            print(f"    {m['aa']}{m['position']:<6} "
                  f"{m['chain']:<6} "
                  f"feat={m['primary_feature']:<15} "
                  f"cons={m['conservation']:.3f} "
                  f"comp={m['composite']:.3f}")
        if cid in centroids:
            ct = centroids[cid]
            print(f"    Centroid: "
                  f"({ct['centroid'][0]:.1f}, "
                  f"{ct['centroid'][1]:.1f}, "
                  f"{ct['centroid'][2]:.1f})")
            print(f"    Mean cons: "
                  f"{ct['mean_conservation']:.3f}  "
                  f"Mean comp: "
                  f"{ct['mean_composite']:.3f}")
    print('='*60)


# ── Main ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('CGCP Phase 2 Step 4 - DBSCAN Spatial Clustering')
    print('Loading features and coordinates...')

    features = load_features()
    coords   = get_coordinates(features)
    coords, labels = run_dbscan(coords, eps=8.0,
                                min_samples=2)
    centroids = compute_centroids(coords)

    save_tsv(coords)
    out = save_json(coords, centroids)
    make_figure(coords, centroids)
    print_summary(coords, centroids)

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 5 - conservation overlay")
