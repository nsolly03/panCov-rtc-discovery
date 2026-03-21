#!/usr/bin/env python3
"""
CGCP Phase 2 Step 6 - Integrated Assessment: NSP12-NSP7
Combines all evidence from Steps 1-5 into one decision table.
Produces the CGCP pharmacophore candidate list.

Evidence integrated:
  Step 1: Structure verified (7/7 PASS)
  Step 2: Contact map (29 residues, PHE440 anchor)
  Step 3: Feature classification
  Step 4: DBSCAN clusters (3 clusters)
  Step 5: Conservation overlay (6 NSP12 residues cons=1.000)

Decision criteria for pharmacophore inclusion:
  MUST: cluster == 0 (primary cluster)
  MUST: chain == NSP12 (NSP12 side of interface only)
  MUST: conservation >= 0.689 (moderate or above)
  SHOULD: composite >= 0.40
  BONUS: conservation == 1.000

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step06_integrated_assessment_NSP12-NSP7.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
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

P = {
    'include':   '#1A7D2E', 'include_f':  '#A8D5B0',
    'secondary': '#2166AC', 'secondary_f':'#92C5DE',
    'exclude':   '#D7191C', 'exclude_f':  '#FDAE61',
    'anchor':    '#D7191C', 'anchor_f':   '#FDAE61',
    'gray':      '#636363', 'text':       '#1A1A1A',
}

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
S2_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-02-contacts')
S3_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-03-features')
S4_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-04-clusters')
S5_DIR  = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-05-conservation')
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-06-assessment')
os.makedirs(OUT_DIR, exist_ok=True)


# ── Load all step outputs ────────────────────────────────────
def load_all():
    # Step 2: contact map
    s2 = {}
    with open(os.path.join(S2_DIR,
              'contact_map_NSP12-NSP7.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            s2[row['residue']] = row

    # Step 3: features
    s3 = {}
    with open(os.path.join(S3_DIR,
              'feature_classification_NSP12-NSP7.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            s3[row['residue']] = row

    # Step 5: conservation overlay (has cluster info too)
    s5 = []
    with open(os.path.join(S5_DIR,
              'conservation_overlay_NSP12-NSP7.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            s5.append(row)

    return s2, s3, s5


# ── Score each residue against inclusion criteria ────────────
def assess(s2, s3, s5):
    records = []
    for row in s5:
        res    = row['residue']
        chain  = row['chain']
        pos    = int(row['position'])
        aa     = row['aa']
        cl     = int(row['cluster'])
        cons   = float(row['conservation'])
        tier   = row['cons_tier']
        feat   = row['primary_feature']
        comp   = float(row['composite'])

        # Evidence scores (0 or 1 each)
        e1 = 1 if cl == 0 else 0               # primary cluster
        e2 = 1 if chain == 'NSP12' else 0      # NSP12 side
        e3 = 1 if cons >= 0.689 else 0         # moderate+ conservation
        e4 = 1 if comp >= 0.40 else 0          # composite threshold
        e5 = 1 if cons == 1.000 else 0         # identical (bonus)
        e6 = 1 if feat in ('anchor','aromatic',
                            'hydrophobic') else 0  # druggable feature

        total = e1 + e2 + e3 + e4 + e5 + e6

        # Special rule: high-scoring noise residues
        # TYR420 and GLU431 in noise cluster (-1) due to
        # spatial isolation but biologically critical
        is_high_noise = (cl == -1 and chain == 'NSP12'
                         and cons == 1.000 and comp >= 0.60)

        # Decision
        if pos == 440 and chain == 'NSP12':
            decision = 'ANCHOR'
        elif e1 and e2 and e3 and e4 and e5:
            decision = 'INCLUDE'
        elif is_high_noise:
            decision = 'INCLUDE'
        elif e1 and e2 and e3:
            decision = 'SECONDARY'
        elif cl == -1 and chain == 'NSP12' and cons >= 0.689:
            decision = 'SECONDARY'
        else:
            decision = 'EXCLUDE'

        records.append({
            'residue':         res,
            'chain':           chain,
            'position':        pos,
            'aa':              aa,
            'primary_feature': feat,
            'cluster':         cl,
            'conservation':    cons,
            'cons_tier':       tier,
            'composite':       comp,
            'e1_cluster0':     e1,
            'e2_nsp12':        e2,
            'e3_cons_mod':     e3,
            'e4_composite':    e4,
            'e5_identical':    e5,
            'e6_druggable':    e6,
            'evidence_score':  total,
            'decision':        decision,
        })

    records.sort(key=lambda x: (
        0 if x['decision'] == 'ANCHOR' else
        1 if x['decision'] == 'INCLUDE' else
        2 if x['decision'] == 'SECONDARY' else 3,
        -x['evidence_score']
    ))
    return records


# ── Save outputs ─────────────────────────────────────────────
def save_tsv(records):
    path = os.path.join(OUT_DIR,
        'integrated_assessment_NSP12-NSP7.tsv')
    fields = ['residue','chain','position','aa',
              'primary_feature','cluster',
              'conservation','cons_tier','composite',
              'e1_cluster0','e2_nsp12','e3_cons_mod',
              'e4_composite','e5_identical','e6_druggable',
              'evidence_score','decision']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(records)
    print(f"  TSV: {path}")


def save_json(records):
    anchor    = [r for r in records
                 if r['decision'] == 'ANCHOR']
    included  = [r for r in records
                 if r['decision'] == 'INCLUDE']
    secondary = [r for r in records
                 if r['decision'] == 'SECONDARY']
    excluded  = [r for r in records
                 if r['decision'] == 'EXCLUDE']

    out = {
        'interface':      'NSP12-NSP7',
        'date':           '2026-03-18',
        'step':           'Phase 2 Step 6 - Integrated Assessment',
        'criteria': {
            'E1': 'Primary cluster (Cluster 0)',
            'E2': 'NSP12 chain (interface anchor side)',
            'E3': 'Conservation >= 0.689 (moderate+)',
            'E4': 'Composite score >= 0.40',
            'E5': 'Conservation = 1.000 (identical)',
            'E6': 'Druggable feature (anchor/aromatic/hydrophobic)',
        },
        'decisions': {
            'ANCHOR':    len(anchor),
            'INCLUDE':   len(included),
            'SECONDARY': len(secondary),
            'EXCLUDE':   len(excluded),
        },
        'pharmacophore_candidates': [
            r['residue'] for r in anchor + included
        ],
        'secondary_candidates': [
            r['residue'] for r in secondary
        ],
        'anchor':    anchor,
        'included':  included,
        'secondary': secondary,
        'excluded':  excluded,
        'all':       records,
    }

    path = os.path.join(OUT_DIR,
        'integrated_assessment_NSP12-NSP7.json')
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


def decision_color(d):
    if d == 'ANCHOR':    return P['anchor'],    P['anchor_f']
    if d == 'INCLUDE':   return P['include'],   P['include_f']
    if d == 'SECONDARY': return P['secondary'], P['secondary_f']
    return P['exclude'], P['exclude_f']


def make_figure(records, out):
    fig = plt.figure(figsize=(14.0, 11.0))
    gs  = fig.add_gridspec(2, 2,
                           hspace=0.62, wspace=0.45,
                           left=0.09, right=0.97,
                           top=0.93, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, :])   # full-width top
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # ── Panel A: Evidence score per residue ───────────────────
    # Show all residues sorted by decision then evidence score
    top = [r for r in records
           if r['decision'] in ('ANCHOR','INCLUDE',
                                'SECONDARY')]
    top = sorted(top, key=lambda x: -x['evidence_score'])

    labels = [r['aa'][:3]+str(r['position'])
              for r in top]
    evals  = [r['evidence_score'] for r in top]
    x      = np.arange(len(top))
    w      = 0.42

    for xi, r in enumerate(top):
        ec, fc = decision_color(r['decision'])
        ax_a.bar(xi, r['evidence_score'],
                 width=w, facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3, clip_on=False)
        ax_a.text(xi, r['evidence_score']+0.06,
                  str(r['evidence_score']),
                  ha='center', va='bottom',
                  fontsize=7.5, color=ec,
                  fontweight='bold', clip_on=False)

        # Decision badge below bar
        badge = ('A' if r['decision']=='ANCHOR' else
                 'I' if r['decision']=='INCLUDE' else 'S')
        ax_a.text(xi, -0.35, badge,
                  ha='center', va='top',
                  fontsize=7, color=ec,
                  fontweight='bold', clip_on=False)

    # Threshold line — minimum INCLUDE
    ax_a.axhline(5, color=P['include'],
                 linewidth=0.8, linestyle=(0,(4,3)),
                 zorder=1)
    ax_a.text(len(top)-0.45, 5.08,
              'INCLUDE threshold (5/6)',
              fontsize=7, color=P['include'],
              ha='right', va='bottom', clip_on=False)

    ax_a.axhline(3, color=P['secondary'],
                 linewidth=0.8, linestyle=(0,(4,3)),
                 zorder=1)
    ax_a.text(len(top)-0.45, 3.08,
              'SECONDARY threshold (3/6)',
              fontsize=7, color=P['secondary'],
              ha='right', va='bottom', clip_on=False)

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(labels, fontsize=8,
                         rotation=90, ha='center',
                         va='top')
    ax_a.tick_params(axis='x', pad=2,
                     length=3, width=0.75)
    ax_a.set_xlim(-0.6, len(top)-0.4)
    ax_a.set_ylabel('Evidence score (max 6)', fontsize=9)
    prism_axes(ax_a, ymax=7,
               yticks=[0,1,2,3,4,5,6])
    ax_a.set_title(
        'A   Integrated evidence score per residue\n'
        '(A=anchor; I=include; S=secondary; '
        'criteria: cluster0 + NSP12 + cons + composite + identical + druggable)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    # Legend
    handles = [
        mpatches.Patch(facecolor=P['anchor_f'],
                       edgecolor=P['anchor'],
                       linewidth=0.9,
                       label='ANCHOR — PHE440'),
        mpatches.Patch(facecolor=P['include_f'],
                       edgecolor=P['include'],
                       linewidth=0.9,
                       label='INCLUDE — pharmacophore core'),
        mpatches.Patch(facecolor=P['secondary_f'],
                       edgecolor=P['secondary'],
                       linewidth=0.9,
                       label='SECONDARY — supporting residues'),
    ]
    ax_a.legend(handles=handles, loc='upper right',
                fontsize=7.5, frameon=True,
                facecolor='white', edgecolor='#CCCCCC')

    # ── Panel B: Decision count bar chart ─────────────────────
    dec_order  = ['ANCHOR','INCLUDE','SECONDARY','EXCLUDE']
    dec_counts = [out['decisions'].get(d, 0)
                  for d in dec_order]
    dec_colors = [P['anchor'], P['include'],
                  P['secondary'], P['exclude']]
    dec_fills  = [P['anchor_f'], P['include_f'],
                  P['secondary_f'], P['exclude_f']]

    x2 = np.arange(len(dec_order))
    for xi, val, fc, ec in zip(x2, dec_counts,
                                dec_fills, dec_colors):
        ax_b.bar(xi, val, width=0.5,
                 facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3,
                 clip_on=False)
        ax_b.text(xi, val+0.08, str(val),
                  ha='center', va='bottom',
                  fontsize=10, color=ec,
                  fontweight='bold', clip_on=False)

    ax_b.set_xticks(x2)
    ax_b.set_xticklabels(dec_order, fontsize=8.5)
    ax_b.set_ylabel('Number of residues', fontsize=9)
    ax_b.set_xlim(-0.6, len(dec_order)-0.4)
    ymax_b = max(dec_counts) + 3
    prism_axes(ax_b, ymax=ymax_b,
               yticks=range(0, ymax_b+1, 2))
    ax_b.set_title(
        'B   Decision summary\n'
        '(pharmacophore candidates = ANCHOR + INCLUDE)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    # ── Panel C: Composite vs conservation scatter ────────────
    for r in records:
        ec, fc = decision_color(r['decision'])
        size = 120 if r['decision'] == 'ANCHOR' else 50
        marker = '*' if r['decision'] == 'ANCHOR' else 'o'
        ax_c.scatter(r['composite'],
                     r['conservation'],
                     color=fc, edgecolors=ec,
                     linewidths=0.8, s=size,
                     marker=marker, zorder=4,
                     alpha=0.85)
        if r['decision'] in ('ANCHOR','INCLUDE'):
            ax_c.annotate(
                r['aa'][:3]+str(r['position']),
                (r['composite'], r['conservation']),
                fontsize=6.5, color=ec,
                xytext=(4, 2),
                textcoords='offset points',
            )

    # Decision boundaries
    ax_c.axhline(0.689, color=P['gray'],
                 linewidth=0.8,
                 linestyle=(0,(4,3)), zorder=1)
    ax_c.axvline(0.40, color=P['gray'],
                 linewidth=0.8,
                 linestyle=(0,(4,3)), zorder=1)
    ax_c.text(0.41, 0.02, 'composite\nthreshold',
              fontsize=6.5, color=P['gray'],
              va='bottom')
    ax_c.text(0.01, 0.70, 'conservation\nthreshold',
              fontsize=6.5, color=P['gray'],
              va='bottom')

    ax_c.set_xlabel('CGCP composite score', fontsize=9)
    ax_c.set_ylabel('Conservation score', fontsize=9)
    ax_c.set_xlim(-0.05, 1.1)
    ax_c.spines['left'].set_position(('outward', 5))
    ax_c.spines['bottom'].set_position(('outward', 5))
    ax_c.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_c.set_ylim(-0.05, 1.1)
    ax_c.set_title(
        'C   Composite score vs conservation\n'
        '(dashed = decision thresholds; star = PHE440)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    path = os.path.join(OUT_DIR,
        'Fig_Step06_IntegratedAssessment.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


# ── Print decision table ─────────────────────────────────────
def print_summary(records, out):
    print('\n' + '='*65)
    print('STEP 6 INTEGRATED ASSESSMENT — NSP12-NSP7')
    print('='*65)

    print(f"\n  {'Residue':<18} {'Chain':<7} {'Feat':<16} "
          f"{'Cons':>6} {'Comp':>6} {'Score':>6} {'Decision'}")
    print(f"  {'-'*65}")

    for r in records:
        sym = ('★' if r['decision']=='ANCHOR' else
               '✓' if r['decision']=='INCLUDE' else
               '·' if r['decision']=='SECONDARY' else ' ')
        print(f"  {sym} {r['aa']}{r['position']:<14} "
              f"{r['chain']:<7} "
              f"{r['primary_feature']:<16} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f} "
              f"{r['evidence_score']:>4}/6  "
              f"{r['decision']}")

    print(f"\n  PHARMACOPHORE CANDIDATES "
          f"(ANCHOR + INCLUDE):")
    for res in out['pharmacophore_candidates']:
        print(f"    → {res}")

    print(f"\n  SECONDARY CANDIDATES:")
    for res in out['secondary_candidates']:
        print(f"    · {res}")

    print(f"\n  Summary: "
          f"ANCHOR={out['decisions']['ANCHOR']}  "
          f"INCLUDE={out['decisions']['INCLUDE']}  "
          f"SECONDARY={out['decisions']['SECONDARY']}  "
          f"EXCLUDE={out['decisions']['EXCLUDE']}")
    print('='*65)


# ── Main ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('CGCP Phase 2 Step 6 - Integrated Assessment')
    print('Loading all step outputs...')

    s2, s3, s5 = load_all()
    records = assess(s2, s3, s5)

    save_tsv(records)
    out = save_json(records)
    make_figure(records, out)
    print_summary(records, out)

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 7 - pharmacophore hypothesis")
