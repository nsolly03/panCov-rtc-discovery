#!/usr/bin/env python3
"""
CGCP Phase 2 Step 3 - Feature Classification: NSP12-NSP7
Assigns each interface residue a biophysical feature type
based on chemistry, interaction data, and 3D context.

Feature types:
  - aromatic        : PHE, TYR, TRP (pi-stacking capable)
  - hydrophobic     : ALA, VAL, ILE, LEU, MET, PRO, GLY
  - hbond_donor     : SER, THR, TYR, ASN, GLN, LYS, ARG, HIS
  - hbond_acceptor  : ASP, GLU, SER, THR, ASN, GLN
  - charged_pos     : LYS, ARG, HIS
  - charged_neg     : ASP, GLU
  - anchor          : top-ranked residue with cons=1.0 (PHE440)

Note: residues can have multiple features (e.g. TYR = aromatic + hbond_donor)
Primary feature assigned by interaction type from Step 2b data.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step03_feature_classification_NSP12-NSP7.py
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
    'xtick.labelsize':    8,
    'ytick.labelsize':    8,
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

# Feature colors — Prism palette
FC = {
    'anchor':        '#D7191C',   # red
    'aromatic':      '#E6550D',   # orange
    'hydrophobic':   '#636363',   # gray
    'hbond_donor':   '#2166AC',   # blue
    'hbond_acceptor':'#74ADD1',   # light blue
    'charged_pos':   '#7B2D8B',   # purple
    'charged_neg':   '#4DAC26',   # green
}
FC_FILL = {k: v + '55' for k, v in FC.items()}

FEATURE_LABEL = {
    'anchor':        'Anchor (aromatic hydrophobic core)',
    'aromatic':      'Aromatic (pi-stacking)',
    'hydrophobic':   'Hydrophobic (van der Waals)',
    'hbond_donor':   'H-bond donor',
    'hbond_acceptor':'H-bond acceptor',
    'charged_pos':   'Charged+ (SB donor)',
    'charged_neg':   'Charged- (SB acceptor)',
}

# ── Chemistry rules ──────────────────────────────────────────
def get_primary_feature(aa, sb_loss, hb_loss, hy_loss, position, chain):
    """
    Assign primary feature based on:
    1. If anchor position (PHE440 chain A) -> anchor
    2. If salt bridge loss > 0 -> charged
    3. If H-bond loss > 0 -> hbond
    4. If aromatic AA -> aromatic
    5. Otherwise -> hydrophobic
    """
    aa = aa.upper()

    # Rule 1: anchor
    if position == 440 and chain == 'NSP12':
        return 'anchor'

    # Rule 2: salt bridge dominant
    if sb_loss > 0:
        if aa in ('LYS', 'ARG', 'HIS'):
            return 'charged_pos'
        if aa in ('ASP', 'GLU'):
            return 'charged_neg'

    # Rule 3: H-bond dominant
    if hb_loss > 0:
        if aa in ('PHE', 'TYR', 'TRP'):
            return 'aromatic'
        if aa in ('SER', 'THR', 'CYS', 'ASN', 'GLN'):
            return 'hbond_donor'
        if aa in ('ASP', 'GLU'):
            return 'hbond_acceptor'
        if aa in ('LYS', 'ARG', 'HIS'):
            return 'hbond_donor'

    # Rule 4: aromatic
    if aa in ('PHE', 'TYR', 'TRP'):
        return 'aromatic'

    # Rule 5: charged without SB
    if aa in ('LYS', 'ARG', 'HIS'):
        return 'charged_pos'
    if aa in ('ASP', 'GLU'):
        return 'charged_neg'

    # Rule 6: polar
    if aa in ('SER', 'THR', 'CYS', 'ASN', 'GLN'):
        return 'hbond_donor'

    # Default: hydrophobic
    return 'hydrophobic'


def get_secondary_features(aa):
    """All applicable features for a residue (for multi-feature annotation)."""
    aa = aa.upper()
    feats = []
    if aa in ('PHE', 'TYR', 'TRP'):
        feats.append('aromatic')
    if aa in ('ALA','VAL','ILE','LEU','MET','PRO','GLY','PHE','TRP'):
        feats.append('hydrophobic')
    if aa in ('SER','THR','TYR','ASN','GLN','LYS','ARG','HIS'):
        feats.append('hbond_donor')
    if aa in ('ASP','GLU','SER','THR','ASN','GLN','TYR'):
        feats.append('hbond_acceptor')
    if aa in ('LYS','ARG','HIS'):
        feats.append('charged_pos')
    if aa in ('ASP','GLU'):
        feats.append('charged_neg')
    return feats if feats else ['hydrophobic']


# ── Load Step 2 data ─────────────────────────────────────────
BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VAL_DIR = os.path.join(BASE, '02-validation/NSP12-NSP7')
S2_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-02-contacts')
OUT_DIR = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-03-features')
os.makedirs(OUT_DIR, exist_ok=True)


def load_contacts():
    path = os.path.join(S2_DIR, 'contact_map_NSP12-NSP7.tsv')
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append(row)
    return rows


def load_csv():
    rows = []
    with open(os.path.join(VAL_DIR,
              'composite_ranking_NSP12-NSP7_3.csv')) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return {r['position']: r for r in rows}


# ── Classify all residues ────────────────────────────────────
def classify(contacts, csv_data):
    classified = []
    for c in contacts:
        pos   = c['position']
        chain = c['chain']
        aa    = c['aa'].upper() if c['aa'] else 'UNK'
        # Skip residues not present in structure
        if not aa or aa == '???':
            continue

        # Get interaction losses from CSV
        csv_row = csv_data.get(pos, {})
        sb_loss = float(csv_row.get('sb_loss', 0))
        hb_loss = float(csv_row.get('hb_loss', 0))
        hy_loss = float(csv_row.get('hy_loss', 0))

        primary = get_primary_feature(
            aa, sb_loss, hb_loss, hy_loss,
            int(pos), chain
        )
        secondary = get_secondary_features(aa)

        classified.append({
            'residue':          c['residue'],
            'chain':            chain,
            'position':         int(pos),
            'aa':               aa,
            'primary_feature':  primary,
            'secondary_features': ','.join(secondary),
            'sb_loss':          int(sb_loss),
            'hb_loss':          int(hb_loss),
            'hy_loss':          int(hy_loss),
            'bsa':              float(c['bsa']),
            'contacts':         int(c['contacts']),
            'conservation':     float(c['conservation']),
            'composite':        float(c['composite']),
            'is_hotspot':       c['is_hotspot'],
        })

    # Filter out residues with unknown AA (not in this PDB)
    classified = [c for c in classified if c['aa'] != '???']
    classified.sort(key=lambda x: -x['composite'])
    return classified


# ── Save outputs ─────────────────────────────────────────────
def save_tsv(classified):
    path = os.path.join(OUT_DIR,
                        'feature_classification_NSP12-NSP7.tsv')
    fields = ['residue','chain','position','aa',
              'primary_feature','secondary_features',
              'sb_loss','hb_loss','hy_loss',
              'bsa','contacts','conservation','composite',
              'is_hotspot']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(classified)
    print(f"  TSV: {path}")


def save_json(classified):
    # Summary by feature type
    from collections import Counter
    feat_counts = Counter(c['primary_feature']
                          for c in classified)
    nsp12 = [c for c in classified if c['chain'] == 'NSP12']
    nsp7  = [c for c in classified if c['chain'] == 'NSP7']

    out = {
        'interface':    'NSP12-NSP7',
        'date':         '2026-03-18',
        'n_total':      len(classified),
        'feature_counts': dict(feat_counts),
        'nsp12_features': {
            c['aa']+str(c['position']): c['primary_feature']
            for c in nsp12
        },
        'nsp7_features': {
            c['aa']+str(c['position']): c['primary_feature']
            for c in nsp7
        },
        'anchor':       'PHE440',
        'classified':   classified,
    }
    path = os.path.join(OUT_DIR,
                        'feature_classification_NSP12-NSP7.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  JSON: {path}")
    return out


# ── Prism figure — 2 panels ──────────────────────────────────
def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax * 0.04, ymax * 1.08)
    if yticks is not None:
        ax.set_yticks(yticks)


def make_figure(classified, summary):
    fig = plt.figure(figsize=(13.0, 10.0))
    gs  = fig.add_gridspec(2, 2,
                           hspace=0.65, wspace=0.45,
                           left=0.09, right=0.97,
                           top=0.93, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])   # NSP12 feature bars
    ax_b = fig.add_subplot(gs[0, 1])   # NSP7 feature bars
    ax_c = fig.add_subplot(gs[1, 0])   # Feature count donut
    ax_d = fig.add_subplot(gs[1, 1])   # Conservation by feature

    nsp12 = [c for c in classified
             if c['chain'] == 'NSP12'][:10]
    nsp7  = [c for c in classified
             if c['chain'] == 'NSP7'][:8]

    # ── Panels A & B: composite bars colored by feature ──────
    for ax, group, panel, title in [
        (ax_a, nsp12, 'A', 'NSP12 (chain A)'),
        (ax_b, nsp7,  'B', 'NSP7 (chain C)'),
    ]:
        labels = [c['aa'][:3]+str(c['position'])
                  for c in group]
        vals   = [c['composite'] for c in group]
        feats  = [c['primary_feature'] for c in group]
        x      = np.arange(len(group))
        w      = 0.42

        for xi, val, feat in zip(x, vals, feats):
            ax.bar(xi, val, width=w,
                   facecolor=FC_FILL[feat],
                   edgecolor=FC[feat],
                   linewidth=0.9, zorder=3,
                   clip_on=False)
            ax.scatter(xi, val, color=FC[feat],
                       s=16, zorder=5, clip_on=False)
            ax.text(xi, val + 0.012, f'{val:.3f}',
                    ha='center', va='bottom',
                    fontsize=6.5, color=FC[feat],
                    clip_on=False)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.5,
                           rotation=90, ha='center',
                           va='top')
        ax.tick_params(axis='x', pad=2,
                       length=3, width=0.75)
        ax.set_xlim(-0.6, len(group) - 0.4)
        ax.set_ylabel('CGCP composite score', fontsize=9)

        ymax   = round(max(vals) * 1.35 + 0.1, 1)
        step   = 0.2
        yticks = [round(i*step, 1)
                  for i in range(int(ymax/step)+1)]
        prism_axes(ax, ymax=max(yticks), yticks=yticks)
        ax.set_title(
            f'{panel}   Feature classification — {title}\n'
            f'(color = primary biophysical feature)',
            loc='left', fontsize=8.5, pad=6,
            linespacing=1.5)

    # ── Panel C: Feature count bar chart ─────────────────────
    feat_counts = summary['feature_counts']
    fc_labels   = list(feat_counts.keys())
    fc_vals     = [feat_counts[k] for k in fc_labels]
    fc_colors   = [FC[k] for k in fc_labels]
    fc_fills    = [FC_FILL[k] for k in fc_labels]

    y = np.arange(len(fc_labels))
    ax_c.barh(y, fc_vals, height=0.5,
              facecolor=fc_fills,
              edgecolor=fc_colors,
              linewidth=0.9, zorder=3)
    for yi, val, col in zip(y, fc_vals, fc_colors):
        ax_c.text(val + 0.1, yi, str(val),
                  va='center', ha='left',
                  fontsize=8.5, color=col,
                  fontweight='bold', clip_on=False)

    short_labels = [FEATURE_LABEL[k] for k in fc_labels]
    ax_c.set_yticks(y)
    ax_c.set_yticklabels(short_labels, fontsize=8)
    ax_c.set_xlabel('Number of residues', fontsize=9)
    ax_c.set_xlim(0, max(fc_vals) + 2)
    ax_c.spines['left'].set_position(('outward', 5))
    ax_c.spines['bottom'].set_position(('outward', 5))
    ax_c.spines['left'].set_bounds(0, len(fc_labels)-1)
    ax_c.set_title(
        'C   Feature type distribution\n'
        '(all 29 interface residues, NSP12 + NSP7)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel D: Conservation by feature type ────────────────
    feat_types = list(FC.keys())
    feat_cons  = {f: [] for f in feat_types}
    for c in classified:
        feat_cons[c['primary_feature']].append(
            c['conservation'])

    positions, medians, q1s, q3s = [], [], [], []
    labels_d, colors_d, fills_d  = [], [], []
    xi = 0
    for f in feat_types:
        vals = feat_cons[f]
        if not vals:
            continue
        positions.append(xi)
        medians.append(np.median(vals))
        q1s.append(np.percentile(vals, 25))
        q3s.append(np.percentile(vals, 75))
        labels_d.append(f.replace('_', '\n'))
        colors_d.append(FC[f])
        fills_d.append(FC_FILL[f])

        # Individual points
        for v in vals:
            ax_d.scatter(xi, v,
                         color=FC[f], s=22,
                         alpha=0.7, zorder=4,
                         edgecolors='white',
                         linewidths=0.4)
        xi += 1

    # IQR bars
    for xi2, med, q1, q3, fc, ff in zip(
            positions, medians, q1s, q3s,
            colors_d, fills_d):
        ax_d.bar(xi2, q3-q1, width=0.4,
                 bottom=q1,
                 facecolor=ff, edgecolor=fc,
                 linewidth=0.9, zorder=2,
                 clip_on=False)
        ax_d.plot([xi2-0.2, xi2+0.2],
                  [med, med],
                  color=fc, linewidth=1.5,
                  zorder=5)

    ax_d.set_xticks(positions)
    ax_d.set_xticklabels(labels_d, fontsize=7.5)
    ax_d.set_ylabel('Conservation score (0-1)', fontsize=9)
    ax_d.set_xlim(-0.6, len(positions) - 0.4)
    ax_d.axhline(0.8, color='#636363', linewidth=0.8,
                 linestyle=(0,(4,3)), zorder=1)
    ax_d.text(len(positions)-0.5, 0.82,
              'pan-coronavirus\nthreshold',
              fontsize=6.5, color='#636363',
              ha='right', va='bottom')
    prism_axes(ax_d, ymax=1.05,
               yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_d.set_title(
        'D   Conservation score by feature type\n'
        '(bar = IQR; line = median; dots = individual residues)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    path = os.path.join(OUT_DIR,
        'Fig_Step03_FeatureClassification.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


# ── Print summary ────────────────────────────────────────────
def print_summary(classified, summary):
    print('\n' + '='*60)
    print('STEP 3 FEATURE CLASSIFICATION — NSP12-NSP7')
    print('='*60)
    print(f"\n  Feature distribution:")
    for feat, n in summary['feature_counts'].items():
        print(f"    {feat:<20} {n} residues")
    print(f"\n  NSP12 residues:")
    for res, feat in summary['nsp12_features'].items():
        print(f"    {res:<12} {feat}")
    print(f"\n  NSP7 residues:")
    for res, feat in summary['nsp7_features'].items():
        print(f"    {res:<12} {feat}")
    print('='*60)


# ── Main ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('CGCP Phase 2 Step 3 - Feature Classification')

    contacts = load_contacts()
    csv_data = load_csv()
    classified = classify(contacts, csv_data)

    save_tsv(classified)
    summary = save_json(classified)
    make_figure(classified, summary)
    print_summary(classified, summary)

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 4 - DBSCAN spatial clustering")
