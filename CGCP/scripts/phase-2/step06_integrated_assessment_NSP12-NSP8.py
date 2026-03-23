#!/usr/bin/env python3
"""
CGCP Phase 2 Step 6 - Integrated Assessment: NSP12-NSP8
Combines all evidence from Steps 1-5 into one decision table.
Dual anchor: LEU387 (primary) + LYS332 (secondary electrostatic)

Decision criteria:
  ANCHOR_PRIMARY:   LEU387 — highest composite, cons=1.000
  ANCHOR_SECONDARY: LYS332 — salt bridge, cons=0.600
  INCLUDE:  cons >= 0.800 AND composite >= 0.600
  SECONDARY: cons >= 0.600 AND composite >= 0.400
  EXCLUDE:  everything else

Evidence scores (6 criteria):
  E1: In cluster 0 (all residues pass — single cluster)
  E2: NSP12 chain (interface anchor side)
  E3: conservation >= 0.800
  E4: composite >= 0.600
  E5: conservation = 1.000 (identical)
  E6: druggable feature (anchor/aromatic/hydrophobic/charged)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step06_integrated_assessment_NSP12-NSP8.py
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
S5_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-05-conservation')
OUT_DIR = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-06-assessment')
os.makedirs(OUT_DIR, exist_ok=True)

DRUGGABLE_FEATURES = {
    'anchor_primary', 'anchor_secondary',
    'aromatic', 'hydrophobic',
    'charged_pos', 'charged_neg',
}

P = {
    'anchor_p':  '#1A1A1A', 'anchor_p_f':  '#555555',
    'anchor_s':  '#D7191C', 'anchor_s_f':  '#FDAE61',
    'include':   '#1A7D2E', 'include_f':   '#A8D5B0',
    'secondary': '#2166AC', 'secondary_f': '#92C5DE',
    'exclude':   '#636363', 'exclude_f':   '#CCCCCC',
}


def load_conservation():
    path = os.path.join(S5_DIR,
        'conservation_overlay_NSP12-NSP8.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            records.append(row)
    return records


def assess(records):
    assessed = []
    for row in records:
        res   = row['residue']
        chain = row['chain']
        pos   = int(row['position'])
        aa    = row['aa']
        feat  = row['primary_feature']
        cl    = int(row['cluster'])
        cons  = float(row['conservation'])
        tier  = row['cons_tier']
        comp  = float(row['composite'])

        # Evidence scores
        e1 = 1  # all residues in cluster 0
        e2 = 1 if chain == 'NSP12' else 0
        e3 = 1 if cons >= 0.800 else 0
        e4 = 1 if comp >= 0.600 else 0
        e5 = 1 if cons == 1.000 else 0
        e6 = 1 if feat in DRUGGABLE_FEATURES else 0

        total = e1 + e2 + e3 + e4 + e5 + e6

        # Decision
        if pos == 387 and chain == 'NSP12':
            decision = 'ANCHOR_PRIMARY'
        elif pos == 332 and chain == 'NSP12':
            decision = 'ANCHOR_SECONDARY'
        elif cons >= 0.800 and comp >= 0.600:
            decision = 'INCLUDE'
        elif cons >= 0.600 and comp >= 0.400:
            decision = 'SECONDARY'
        else:
            decision = 'EXCLUDE'

        assessed.append({
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
            'e3_cons_high':    e3,
            'e4_composite':    e4,
            'e5_identical':    e5,
            'e6_druggable':    e6,
            'evidence_score':  total,
            'decision':        decision,
        })

    assessed.sort(key=lambda x: (
        0 if x['decision'] == 'ANCHOR_PRIMARY'   else
        1 if x['decision'] == 'ANCHOR_SECONDARY' else
        2 if x['decision'] == 'INCLUDE'          else
        3 if x['decision'] == 'SECONDARY'        else 4,
        -x['evidence_score'],
        -x['composite'],
    ))
    return assessed


def save_outputs(assessed):
    tsv_path = os.path.join(OUT_DIR,
        'integrated_assessment_NSP12-NSP8.tsv')
    fields = ['residue','chain','position','aa',
              'primary_feature','cluster',
              'conservation','cons_tier','composite',
              'e1_cluster0','e2_nsp12','e3_cons_high',
              'e4_composite','e5_identical','e6_druggable',
              'evidence_score','decision']
    with open(tsv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(assessed)
    print(f"  TSV: {tsv_path}")

    dec_counts = {}
    for r in assessed:
        d = r['decision']
        dec_counts[d] = dec_counts.get(d, 0) + 1

    pharmacophore = [r['residue'] for r in assessed
                     if r['decision'] in
                     ('ANCHOR_PRIMARY','ANCHOR_SECONDARY',
                      'INCLUDE')]
    secondary = [r['residue'] for r in assessed
                 if r['decision'] == 'SECONDARY']

    out = {
        'interface':      'NSP12-NSP8',
        'date':           '2026-03-23',
        'step':           'Phase 2 Step 6',
        'anchor_primary':   'LEU387',
        'anchor_secondary': 'LYS332',
        'criteria': {
            'E1': 'Cluster 0 (all pass — single cluster)',
            'E2': 'NSP12 chain',
            'E3': 'Conservation >= 0.800',
            'E4': 'Composite >= 0.600',
            'E5': 'Conservation = 1.000 (identical)',
            'E6': 'Druggable feature',
        },
        'decisions':              dec_counts,
        'pharmacophore_candidates': pharmacophore,
        'secondary_candidates':   secondary,
        'all':                    assessed,
    }
    json_path = os.path.join(OUT_DIR,
        'integrated_assessment_NSP12-NSP8.json')
    with open(json_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  JSON: {json_path}")
    return out


def dec_style(d):
    if d == 'ANCHOR_PRIMARY':   return P['anchor_p'],   P['anchor_p_f']
    if d == 'ANCHOR_SECONDARY': return P['anchor_s'],   P['anchor_s_f']
    if d == 'INCLUDE':          return P['include'],     P['include_f']
    if d == 'SECONDARY':        return P['secondary'],   P['secondary_f']
    return P['exclude'], P['exclude_f']


def make_figure(assessed, out):
    fig = plt.figure(figsize=(16.0, 6.5))
    gs  = fig.add_gridspec(1, 3, wspace=0.48,
                           left=0.06, right=0.97,
                           top=0.88, bottom=0.26)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    # ── Panel A: Evidence score — top 35 residues ─────────────
    top = [r for r in assessed
           if r['decision'] != 'EXCLUDE'][:35]
    labels_a = [f"{r['aa'][:3]}{r['position']}"
                for r in top]
    x = np.arange(len(top))

    for xi, r in enumerate(top):
        ec, fc = dec_style(r['decision'])
        ax_a.bar(xi, r['evidence_score'],
                 width=0.40, facecolor=fc,
                 edgecolor=ec, linewidth=1.2,
                 zorder=3, clip_on=False)
        ax_a.text(xi, r['evidence_score'] + 0.06,
                  str(r['evidence_score']),
                  ha='center', va='bottom',
                  fontsize=7, fontweight='bold',
                  color=ec, clip_on=False)
        badge = ('A1' if r['decision']=='ANCHOR_PRIMARY'
                 else 'A2' if r['decision']=='ANCHOR_SECONDARY'
                 else 'I'  if r['decision']=='INCLUDE'
                 else 'S')
        ax_a.text(xi, -0.4, badge,
                  ha='center', va='top',
                  fontsize=6.5, color=ec,
                  clip_on=False)

    # Threshold lines
    ax_a.axhline(5, color=P['include'], linewidth=0.8,
                 linestyle='--', alpha=0.7)
    ax_a.axhline(3, color=P['secondary'], linewidth=0.8,
                 linestyle='--', alpha=0.7)
    ax_a.text(len(top)-0.5, 5.08,
              'INCLUDE (5/6)', fontsize=7,
              color=P['include'], ha='right', clip_on=False)
    ax_a.text(len(top)-0.5, 3.08,
              'SECONDARY (3/6)', fontsize=7,
              color=P['secondary'], ha='right', clip_on=False)

    ax_a.set_xticks(x)
    set_xticklabels_vertical(ax_a, labels_a, fontsize=7)
    ax_a.set_ylabel('Evidence score (max 6)',
                    fontsize=11, fontweight='bold')
    ax_a.set_xlim(-0.6, len(top) - 0.4)
    prism_axes(ax_a, ymax=7, yticks=[0,1,2,3,4,5,6])

    make_legend(ax_a, [
        (P['anchor_p_f'],  P['anchor_p'],  'ANCHOR primary'),
        (P['anchor_s_f'],  P['anchor_s'],  'ANCHOR secondary'),
        (P['include_f'],   P['include'],   'INCLUDE'),
        (P['secondary_f'], P['secondary'], 'SECONDARY'),
    ], loc='upper right', fontsize=7.5)
    panel_title(ax_a, 'A',
                'Integrated evidence score (non-excluded)',
                'A1=primary anchor  A2=secondary  I=include  S=secondary')

    # ── Panel B: Decision counts ──────────────────────────────
    dec_order  = ['ANCHOR_PRIMARY','ANCHOR_SECONDARY',
                  'INCLUDE','SECONDARY','EXCLUDE']
    dec_labels = ['Anchor\n(primary)','Anchor\n(secondary)',
                  'Include','Secondary','Exclude']
    dec_counts = [out['decisions'].get(d, 0)
                  for d in dec_order]

    x2 = np.arange(len(dec_order))
    for xi, val, d in zip(x2, dec_counts, dec_order):
        ec, fc = dec_style(d)
        ax_b.bar(xi, val, width=0.50,
                 facecolor=fc, edgecolor=ec,
                 linewidth=1.2, zorder=3, clip_on=False)
        ax_b.text(xi, val + 0.3, str(val),
                  ha='center', va='bottom',
                  fontsize=10, fontweight='bold',
                  color=ec, clip_on=False)

    ax_b.set_xticks(x2)
    set_xticklabels_vertical(ax_b, dec_labels, fontsize=9)
    ax_b.set_ylabel('Number of residues',
                    fontsize=11, fontweight='bold')
    ymax_b = max(dec_counts) + 8
    prism_axes(ax_b, ymax=ymax_b,
               yticks=range(0, ymax_b, 10))
    ax_b.set_xlim(-0.6, len(dec_order) - 0.4)
    panel_title(ax_b, 'B',
                'Decision summary',
                'pharmacophore = ANCHOR + INCLUDE')

    # ── Panel C: Composite vs conservation scatter ────────────
    for r in assessed:
        ec, fc = dec_style(r['decision'])
        is_p = r['decision'] == 'ANCHOR_PRIMARY'
        is_s = r['decision'] == 'ANCHOR_SECONDARY'
        size   = 220 if is_p else (130 if is_s else 35)
        marker = '*'  if is_p else ('D' if is_s else 'o')
        ax_c.scatter(r['composite'], r['conservation'],
                     color=fc, edgecolors=ec,
                     linewidths=0.9, s=size,
                     marker=marker, alpha=0.85, zorder=4)
        if r['decision'] in ('ANCHOR_PRIMARY',
                              'ANCHOR_SECONDARY',
                              'INCLUDE'):
            ax_c.annotate(
                f"{r['aa'][:3]}{r['position']}",
                (r['composite'], r['conservation']),
                fontsize=6.5, color=ec,
                xytext=(4, 2),
                textcoords='offset points')

    ax_c.axhline(0.800, color=P['include'],
                 linewidth=0.8, linestyle='--', alpha=0.7)
    ax_c.axvline(0.600, color=P['include'],
                 linewidth=0.8, linestyle='--', alpha=0.7)
    ax_c.text(0.01, 0.81, 'cons threshold',
              fontsize=7, color=P['include'])
    ax_c.text(0.61, 0.02, 'comp\nthreshold',
              fontsize=7, color=P['include'], va='bottom')

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
    ax_c.tick_params(axis='both', labelsize=9,
                     width=1.2, length=6, direction='out')
    panel_title(ax_c, 'C',
                'Composite vs conservation',
                '★=primary  ◆=secondary  dashed=thresholds')

    path = os.path.join(OUT_DIR,
        'Fig_Step06_IntegratedAssessment_NSP12-NSP8.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(assessed, out):
    print('\n' + '='*65)
    print('STEP 6 INTEGRATED ASSESSMENT — NSP12-NSP8')
    print('='*65)

    print(f"\n  {'Residue':<22} {'Chain':<7} {'Feat':<18} "
          f"{'Cons':>6} {'Comp':>6} {'Score':>6} {'Decision'}")
    print(f"  {'-'*75}")

    for r in assessed:
        if r['decision'] == 'EXCLUDE':
            continue
        sym = ('★' if r['decision']=='ANCHOR_PRIMARY'
               else '◆' if r['decision']=='ANCHOR_SECONDARY'
               else '✓' if r['decision']=='INCLUDE'
               else '·')
        print(f"  {sym} {r['aa']}{r['position']:<17} "
              f"{r['chain']:<7} "
              f"{r['primary_feature']:<18} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f} "
              f"{r['evidence_score']:>4}/6  "
              f"{r['decision']}")

    print(f"\n  PHARMACOPHORE CANDIDATES:")
    for res in out['pharmacophore_candidates']:
        print(f"    → {res}")

    print(f"\n  SECONDARY CANDIDATES ({len(out['secondary_candidates'])}):")
    for res in out['secondary_candidates'][:10]:
        print(f"    · {res}")
    if len(out['secondary_candidates']) > 10:
        print(f"    ... (+{len(out['secondary_candidates'])-10} more)")

    print(f"\n  Summary: "
          + " | ".join(f"{k}={v}"
                       for k, v in out['decisions'].items()))
    print('='*65)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 6 — NSP12-NSP8 Integrated Assessment')
    records  = load_conservation()
    assessed = assess(records)
    out      = save_outputs(assessed)
    make_figure(assessed, out)
    print_summary(assessed, out)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 7 — pharmacophore hypothesis")
