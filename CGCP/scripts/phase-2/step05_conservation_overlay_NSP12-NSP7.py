#!/usr/bin/env python3
"""
CGCP Phase 2 Step 5 - Conservation Overlay: NSP12-NSP7
Overlays per-residue conservation onto DBSCAN clusters.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step05_conservation_overlay_NSP12-NSP7.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
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
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
})

CLUSTER_COLORS = {0:'#D7191C', 1:'#2166AC', 2:'#4DAC26', -1:'#CCCCCC'}

def cons_color(score):
    if score >= 1.000: return '#1A7D2E', '#A8D5B0'
    elif score >= 0.800: return '#4DAC26', '#C5E8B0'
    elif score >= 0.500: return '#E6A817', '#FAE0A0'
    else: return '#D7191C', '#FDAE61'

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VAL_DIR = os.path.join(BASE, '02-validation/NSP12-NSP7')
S4_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-04-clusters')
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-05-conservation')
os.makedirs(OUT_DIR, exist_ok=True)


def load_clusters():
    rows = []
    with open(os.path.join(S4_DIR, 'clusters_NSP12-NSP7.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append(row)
    return rows


def load_conservation():
    cons_nsp12, cons_nsp7 = {}, {}
    with open(os.path.join(VAL_DIR, 'conservation_NSP12.csv')) as f:
        for row in csv.DictReader(f):
            pos = int(row['position'])
            cons_nsp12[pos] = {
                'score': float(row['conservation']),
                'SARS2': row.get('SARS-CoV-2', '-'),
                'SARS1': row.get('SARS-CoV-1', '-'),
                'MERS':  row.get('MERS-CoV', '-'),
                'H229E': row.get('HCoV-229E', '-'),
                'HNL63': row.get('HCoV-NL63', '-'),
            }
    with open(os.path.join(VAL_DIR, 'conservation_NSP7.csv')) as f:
        for row in csv.DictReader(f):
            pos = int(row['position'])
            cons_nsp7[pos] = {
                'score': float(row['conservation']),
                'SARS2': row.get('SARS-CoV-2', '-'),
                'SARS1': row.get('SARS-CoV-1', '-'),
                'MERS':  row.get('MERS-CoV', '-'),
                'H229E': row.get('HCoV-229E', '-'),
                'HNL63': row.get('HCoV-NL63', '-'),
            }
    return cons_nsp12, cons_nsp7


def merge(clusters, cons_nsp12, cons_nsp7):
    merged = []
    for c in clusters:
        pos   = int(c['position'])
        chain = c['chain']
        cons  = cons_nsp12.get(pos, {}) if chain == 'NSP12' else cons_nsp7.get(pos, {})
        score = cons.get('score', float(c['conservation']))

        if score >= 1.000:   tier = 'identical'
        elif score >= 0.800: tier = 'pan-coronavirus'
        elif score >= 0.500: tier = 'moderate'
        else:                tier = 'variable'

        merged.append({
            **c,
            'conservation': round(score, 3),
            'cons_tier':    tier,
            'SARS2':  cons.get('SARS2', '-'),
            'SARS1':  cons.get('SARS1', '-'),
            'MERS':   cons.get('MERS',  '-'),
            'H229E':  cons.get('H229E', '-'),
            'HNL63':  cons.get('HNL63', '-'),
        })
    return merged


def save_tsv(merged):
    path = os.path.join(OUT_DIR, 'conservation_overlay_NSP12-NSP7.tsv')
    fields = ['residue','chain','position','aa','primary_feature',
              'cluster','conservation','cons_tier',
              'SARS2','SARS1','MERS','H229E','HNL63','composite']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(merged)
    print(f"  TSV: {path}")


def save_json(merged):
    from collections import defaultdict
    clusters = defaultdict(list)
    for m in merged:
        clusters[int(m['cluster'])].append(m)

    cluster_summary = {}
    for cid, members in clusters.items():
        scores = [float(m['conservation']) for m in members]
        tiers  = [m['cons_tier'] for m in members]
        n_ident   = tiers.count('identical')
        n_pancov  = n_ident + tiers.count('pan-coronavirus')
        cluster_summary[cid] = {
            'n_total':      len(members),
            'n_identical':  n_ident,
            'n_pan_cov':    n_pancov,
            'n_moderate':   tiers.count('moderate'),
            'n_variable':   tiers.count('variable'),
            'mean_cons':    round(float(np.mean(scores)), 3),
            'min_cons':     round(float(np.min(scores)),  3),
            'pan_cov_frac': round(n_pancov / len(members), 3),
        }

    out = {
        'interface':             'NSP12-NSP7',
        'date':                  '2026-03-18',
        'pharmacophore_cluster': 0,
        'cluster_summary':       cluster_summary,
        'residues':              merged,
    }
    path = os.path.join(OUT_DIR, 'conservation_overlay_NSP12-NSP7.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  JSON: {path}")
    return out


def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax*0.04, ymax*1.1)
    if yticks is not None:
        ax.set_yticks(yticks)


def make_figure(merged, out):
    from collections import defaultdict

    fig = plt.figure(figsize=(13.0, 11.0))
    gs  = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.45,
                           left=0.09, right=0.97,
                           top=0.93, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    cl0 = sorted([m for m in merged if int(m['cluster']) == 0],
                 key=lambda x: -float(x['composite']))
    labels_a = [m['aa'][:3]+str(m['position']) for m in cl0]
    x = np.arange(len(cl0))
    w = 0.42

    # ── Panel A: Conservation bars Cluster 0 ─────────────────
    for xi, m in enumerate(cl0):
        val = float(m['conservation'])
        ec, fc = cons_color(val)
        ax_a.bar(xi, val, width=w, facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3, clip_on=False)
        ax_a.scatter(xi, val, color=ec, s=18, zorder=5, clip_on=False)
        ax_a.text(xi, val+0.015, f'{val:.3f}', ha='center', va='bottom',
                  fontsize=6.5, color=ec, fontweight='bold', clip_on=False)

    for thresh, lbl, col in [
        (1.000, 'identical (5/5)', '#1A7D2E'),
        (0.800, 'pan-coronavirus', '#4DAC26'),
        (0.500, 'moderate',        '#E6A817'),
    ]:
        ax_a.axhline(thresh, color=col, linewidth=0.8,
                     linestyle=(0,(4,3)), zorder=1)
        ax_a.text(len(cl0)-0.45, thresh+0.015, lbl,
                  fontsize=6.5, color=col, ha='right',
                  va='bottom', clip_on=False)

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(labels_a, fontsize=7.5,
                         rotation=90, ha='center', va='top')
    ax_a.tick_params(axis='x', pad=2, length=3, width=0.75)
    ax_a.set_xlim(-0.6, len(cl0)-0.4)
    ax_a.set_ylabel('Conservation score (0-1)', fontsize=9)
    prism_axes(ax_a, ymax=1.05, yticks=[0,0.2,0.4,0.6,0.8,1.0])
    ax_a.set_title('A   Per-residue conservation — Cluster 0\n'
                   '(pharmacophore cluster; color = conservation tier)',
                   loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    # ── Panel B: Stacked tier counts per cluster ──────────────
    cl_ids = sorted([k for k in out['cluster_summary'].keys() if k >= 0])
    x2 = np.arange(len(cl_ids))
    tier_order  = ['identical','pan-coronavirus','moderate','variable']
    tier_colors = ['#1A7D2E','#4DAC26','#E6A817','#D7191C']
    tier_fills  = ['#A8D5B0','#C5E8B0','#FAE0A0','#FDAE61']
    bottoms = np.zeros(len(cl_ids))

    for tier, tc, tf in zip(tier_order, tier_colors, tier_fills):
        vals = []
        for cid in cl_ids:
            cs = out['cluster_summary'][cid]
            if tier == 'identical':
                vals.append(cs['n_identical'])
            elif tier == 'pan-coronavirus':
                vals.append(cs['n_pan_cov'] - cs['n_identical'])
            elif tier == 'moderate':
                vals.append(cs['n_moderate'])
            else:
                vals.append(cs['n_variable'])
        ax_b.bar(x2, vals, width=0.55, bottom=bottoms,
                 facecolor=tf, edgecolor=tc, linewidth=0.9,
                 zorder=3, label=tier)
        bottoms += np.array(vals, dtype=float)

    ax_b.set_xticks(x2)
    ax_b.set_xticklabels([f'Cluster {k}' for k in cl_ids], fontsize=8.5)
    ax_b.set_ylabel('Number of residues', fontsize=9)
    ax_b.set_xlim(-0.6, len(cl_ids)-0.4)
    ymax_b = int(max(bottoms)) + 2
    prism_axes(ax_b, ymax=ymax_b, yticks=range(0, ymax_b+1, 2))
    ax_b.legend(loc='upper right', fontsize=7.5, frameon=True,
                facecolor='white', edgecolor='#CCCCCC')
    ax_b.set_title('B   Conservation tier distribution by cluster\n'
                   '(stacked: identical / pan-cov / moderate / variable)',
                   loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    # ── Panel C: Species heatmap Cluster 0 ───────────────────
    species   = ['SARS2','SARS1','MERS','H229E','HNL63']
    sp_labels = ['SARS-CoV-2','SARS-CoV-1','MERS-CoV',
                 'HCoV-229E','HCoV-NL63']

    matrix = np.zeros((len(cl0), len(species)))
    for i, m in enumerate(cl0):
        ref = m.get('SARS2', '-')
        for j, sp in enumerate(species):
            aa = m.get(sp, '-')
            if not aa or aa == '-':  matrix[i,j] = 0.5
            elif aa == ref:          matrix[i,j] = 1.0
            else:                    matrix[i,j] = 0.0

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'cons_map',
        [(0.0,'#D7191C'),(0.5,'#FFFFFF'),(1.0,'#1A7D2E')])

    im = ax_c.imshow(matrix.T, aspect='auto', cmap=cmap,
                     vmin=0, vmax=1, interpolation='nearest')

    for i, m in enumerate(cl0):
        for j, sp in enumerate(species):
            aa = m.get(sp, '-')
            if aa and aa != '-':
                col = 'white' if matrix[i,j] in (0,1) else '#636363'
                ax_c.text(i, j, aa, ha='center', va='center',
                          fontsize=7, color=col, fontweight='bold')

    ax_c.set_xticks(range(len(cl0)))
    ax_c.set_xticklabels(labels_a, fontsize=7,
                         rotation=90, ha='center', va='top')
    ax_c.set_yticks(range(len(species)))
    ax_c.set_yticklabels(sp_labels, fontsize=8)
    ax_c.tick_params(axis='x', pad=2, length=0, width=0)
    ax_c.set_title('C   Species conservation heatmap — Cluster 0\n'
                   '(green=identical to SARS-CoV-2; red=different; white=gap)',
                   loc='left', fontsize=8.5, pad=6, linespacing=1.5)
    plt.colorbar(im, ax=ax_c, fraction=0.03, pad=0.02,
                 label='Conservation')

    # ── Panel D: Pan-cov fraction per cluster ────────────────
    pc_fracs  = [out['cluster_summary'][k]['pan_cov_frac'] for k in cl_ids]
    pc_colors = [CLUSTER_COLORS.get(k,'#AAAAAA') for k in cl_ids]
    pc_fills  = [c+'55' for c in pc_colors]

    for xi, val, fc, ec in zip(x2, pc_fracs, pc_fills, pc_colors):
        ax_d.bar(xi, val, width=w, facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3, clip_on=False)
        ax_d.scatter(xi, val, color=ec, s=20, zorder=5, clip_on=False)
        ax_d.text(xi, val+0.015, f'{val:.0%}',
                  ha='center', va='bottom', fontsize=9,
                  color=ec, fontweight='bold', clip_on=False)

    ax_d.axhline(0.80, color='#636363', linewidth=0.8,
                 linestyle=(0,(4,3)), zorder=1)
    ax_d.text(len(cl_ids)-0.45, 0.815, 'target >= 80%',
              fontsize=7, color='#636363', ha='right',
              va='bottom', clip_on=False)
    ax_d.set_xticks(x2)
    ax_d.set_xticklabels([f'Cluster {k}' for k in cl_ids], fontsize=8.5)
    ax_d.set_ylabel('Fraction pan-coronavirus conserved', fontsize=9)
    ax_d.set_xlim(-0.6, len(cl_ids)-0.4)
    prism_axes(ax_d, ymax=1.05, yticks=[0,0.2,0.4,0.6,0.8,1.0])
    ax_d.set_title('D   Pan-coronavirus fraction per cluster\n'
                   '(residues with conservation >= 0.80)',
                   loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    path = os.path.join(OUT_DIR, 'Fig_Step05_ConservationOverlay.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(merged, out):
    print('\n' + '='*60)
    print('STEP 5 CONSERVATION OVERLAY - NSP12-NSP7')
    print('='*60)
    for cid, cs in sorted(out['cluster_summary'].items()):
        print(f"\n  Cluster {cid}:")
        print(f"    Total: {cs['n_total']}  "
              f"Identical: {cs['n_identical']}  "
              f"Pan-cov: {cs['n_pan_cov']}  "
              f"Moderate: {cs['n_moderate']}  "
              f"Variable: {cs['n_variable']}")
        print(f"    Mean cons: {cs['mean_cons']:.3f}  "
              f"Pan-cov frac: {cs['pan_cov_frac']:.1%}")

    print(f"\n  Cluster 0 residue detail:")
    cl0 = sorted([m for m in merged if int(m['cluster'])==0],
                 key=lambda x: -float(x['conservation']))
    for m in cl0:
        sym = ('★' if m['cons_tier']=='identical' else
               '✓' if m['cons_tier']=='pan-coronavirus' else
               '~' if m['cons_tier']=='moderate' else '✗')
        print(f"    {sym} {m['aa']}{m['position']:<6} "
              f"{m['chain']:<6} "
              f"cons={float(m['conservation']):.3f} "
              f"({m['cons_tier']})")
    print('='*60)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 5 - Conservation Overlay')
    clusters = load_clusters()
    cons_nsp12, cons_nsp7 = load_conservation()
    merged = merge(clusters, cons_nsp12, cons_nsp7)
    save_tsv(merged)
    out = save_json(merged)
    make_figure(merged, out)
    print_summary(merged, out)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 6 - integrated assessment")
