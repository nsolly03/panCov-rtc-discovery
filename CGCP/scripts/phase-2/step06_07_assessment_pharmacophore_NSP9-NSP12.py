#!/usr/bin/env python3
"""
CGCP Phase 2 Steps 6+7 - Assessment + Pharmacophore: NSP9-NSP12
Anchor: ARG733 (NSP12, cons=1.000)

Pharmacophoric elements:
  E1 — H-bond/hydrophobic core (NSP9 N-terminus)
       ASN2, LEU4, ASN96, LEU97, ASN1, ASN95, LEU103, LEU112
  E2 — Electrostatic anchor (NSP12)
       ARG733(anchor), ARG735, VAL233
  E3 — Supporting contacts
       GLU3(NSP9), ARG99(NSP9), ASP221(NSP12), ASP291(NSP12)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step06_07_assessment_pharmacophore_NSP9-NSP12.py
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
from Bio import PDB
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

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP9-NSP12_5/receptor_NSP9-NSP12_5.pdb')
S5_DIR   = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-05-conservation')
OUT6     = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-06-assessment')
OUT7     = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-07-pharmacophore')
os.makedirs(OUT6, exist_ok=True)
os.makedirs(OUT7, exist_ok=True)

DRUGGABLE = {'anchor','aromatic','hydrophobic',
             'charged_pos','charged_neg','hbond_donor'}

P = {
    'anchor':    '#1A1A1A', 'anchor_f':    '#555555',
    'include':   '#1A7D2E', 'include_f':   '#A8D5B0',
    'secondary': '#2166AC', 'secondary_f': '#92C5DE',
    'exclude':   '#636363', 'exclude_f':   '#CCCCCC',
}

EL_COLORS = {'E1':'#D7191C','E2':'#2166AC','E3':'#4DAC26'}
EL_FILLS  = {'E1':'#FDAE61','E2':'#92C5DE','E3':'#A6D96A'}
EL_NAMES  = {
    'E1': 'E1 — H-bond/hydrophobic core (NSP9)',
    'E2': 'E2 — Electrostatic anchor (NSP12)',
    'E3': 'E3 — Supporting contacts',
}

# Pharmacophore element assignments
# (chain, position): element
ELEMENT_MAP = {
    ('NSP9',  2): 'E1', ('NSP9',  4): 'E1',
    ('NSP9', 96): 'E1', ('NSP9', 97): 'E1',
    ('NSP9',  1): 'E1', ('NSP9', 95): 'E1',
    ('NSP9',103): 'E1', ('NSP9',112): 'E1',
    ('NSP12',733): 'E2',('NSP12',735): 'E2',
    ('NSP12',233): 'E2',
    ('NSP9',  3): 'E3', ('NSP9', 99): 'E3',
    ('NSP12',221): 'E3',('NSP12',291): 'E3',
}


def load_conservation():
    path = os.path.join(S5_DIR,
        'conservation_overlay_NSP9-NSP12.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            records.append(row)
    return records


# ── STEP 6: Integrated Assessment ───────────────────────────
def assess(records):
    assessed = []
    for row in records:
        res   = row['residue']
        chain = row['chain']
        pos   = int(row['position'])
        aa    = row['aa']
        feat  = row['primary_feature']
        cons  = float(row['conservation'])
        tier  = row['cons_tier']
        comp  = float(row['composite'])

        e1 = 1  # all in cluster 0
        e2 = 1 if chain == 'NSP12' else 0
        e3 = 1 if cons >= 0.800 else 0
        e4 = 1 if comp >= 0.600 else 0
        e5 = 1 if cons == 1.000 else 0
        e6 = 1 if feat in DRUGGABLE else 0
        total = e1+e2+e3+e4+e5+e6

        if pos == 733 and chain == 'NSP12':
            decision = 'ANCHOR'
        elif cons >= 0.800 and comp >= 0.600:
            decision = 'INCLUDE'
        elif cons >= 0.600 and comp >= 0.400:
            decision = 'SECONDARY'
        else:
            decision = 'EXCLUDE'

        assessed.append({
            'residue':        res,
            'chain':          chain,
            'position':       pos,
            'aa':             aa,
            'primary_feature':feat,
            'conservation':   cons,
            'cons_tier':      tier,
            'composite':      comp,
            'e1_cluster0':    e1,
            'e2_nsp12':       e2,
            'e3_cons_high':   e3,
            'e4_composite':   e4,
            'e5_identical':   e5,
            'e6_druggable':   e6,
            'evidence_score': total,
            'decision':       decision,
        })

    assessed.sort(key=lambda x: (
        0 if x['decision']=='ANCHOR'    else
        1 if x['decision']=='INCLUDE'   else
        2 if x['decision']=='SECONDARY' else 3,
        -x['evidence_score'], -x['composite']))
    return assessed


def save_step6(assessed):
    tsv = os.path.join(OUT6,
        'integrated_assessment_NSP9-NSP12.tsv')
    fields = ['residue','chain','position','aa',
              'primary_feature','conservation',
              'cons_tier','composite',
              'e1_cluster0','e2_nsp12','e3_cons_high',
              'e4_composite','e5_identical','e6_druggable',
              'evidence_score','decision']
    with open(tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(assessed)
    print(f"  Step6 TSV: {tsv}")

    dec_counts = {}
    for r in assessed:
        d = r['decision']
        dec_counts[d] = dec_counts.get(d,0)+1

    pharmacophore = [r['residue'] for r in assessed
                     if r['decision'] in ('ANCHOR','INCLUDE')]
    secondary     = [r['residue'] for r in assessed
                     if r['decision']=='SECONDARY']

    out = {
        'interface':              'NSP9-NSP12',
        'anchor':                 'ARG733',
        'decisions':              dec_counts,
        'pharmacophore_candidates': pharmacophore,
        'secondary_candidates':   secondary,
        'all':                    assessed,
    }
    jpath = os.path.join(OUT6,
        'integrated_assessment_NSP9-NSP12.json')
    with open(jpath, 'w') as f:
        json.dump(out, f, indent=2, cls=NpEncoder)
    print(f"  Step6 JSON: {jpath}")
    return out


# ── STEP 7: Pharmacophore Hypothesis ────────────────────────
def get_coords():
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]
    chain_a   = structure['A']
    chain_g   = structure['G']
    coords    = {}

    for (chain_name, pos), el in ELEMENT_MAP.items():
        chain = chain_a if chain_name=='NSP12' else chain_g
        try:
            res = chain[(' ', pos, ' ')]
            ca  = res['CA'].coord
            coords[(chain_name, pos)] = {
                'x': round(float(ca[0]),3),
                'y': round(float(ca[1]),3),
                'z': round(float(ca[2]),3),
                'aa': res.resname,
            }
        except KeyError:
            print(f"  Warning: {chain_name}-{pos} not found")
    return coords


def compute_centroids(coords):
    elements = {}
    for (chain_name, pos), el in ELEMENT_MAP.items():
        if el not in elements:
            elements[el] = []
        if (chain_name, pos) in coords:
            elements[el].append(coords[(chain_name, pos)])

    centroids = {}
    for el, pts in elements.items():
        cx = float(np.mean([p['x'] for p in pts]))
        cy = float(np.mean([p['y'] for p in pts]))
        cz = float(np.mean([p['z'] for p in pts]))
        radius = float(max(
            np.sqrt((p['x']-cx)**2 +
                    (p['y']-cy)**2 +
                    (p['z']-cz)**2)
            for p in pts))
        centroids[el] = {
            'x': round(cx,3), 'y': round(cy,3),
            'z': round(cz,3), 'radius': round(radius,2),
        }
    return centroids


def build_pharmacophore(coords, centroids, assessed):
    # Get composite + conservation from assessment
    assess_map = {
        (r['chain'], r['position']): r
        for r in assessed}

    features = []
    for (chain_name, pos), el in ELEMENT_MAP.items():
        if (chain_name, pos) not in coords:
            continue
        c    = coords[(chain_name, pos)]
        amap = assess_map.get((chain_name, pos), {})
        cons = amap.get('conservation', 0.750)
        comp = amap.get('composite', 0.500)
        feat = amap.get('primary_feature', 'unknown')

        pharm_type = ('HYDROPHOBIC' if feat=='hydrophobic'
                      else 'AROMATIC'     if feat=='aromatic'
                      else 'CHARGED_POS'  if feat=='charged_pos'
                      else 'CHARGED_NEG'  if feat=='charged_neg'
                      else 'HBOND'        if feat=='hbond_donor'
                      else 'ANCHOR'       if feat=='anchor'
                      else 'UNKNOWN')
        tol = (1.0 if pharm_type in ('CHARGED_POS',
                                      'CHARGED_NEG','ANCHOR')
               else 1.5 if pharm_type=='AROMATIC'
               else 2.0)

        role = ('ANCHOR' if pos==733 and
                chain_name=='NSP12' else 'INCLUDE')

        features.append({
            'id':           f"{chain_name}-{c['aa']}{pos}",
            'chain':        chain_name,
            'position':     pos,
            'aa':           c['aa'],
            'element':      el,
            'role':         role,
            'pharm_type':   pharm_type,
            'tolerance_A':  tol,
            'conservation': float(cons),
            'composite':    float(comp),
            'x':            c['x'],
            'y':            c['y'],
            'z':            c['z'],
        })

    features.sort(key=lambda x: (
        x['element'],
        0 if x['role']=='ANCHOR' else 1,
        -x['composite']))

    el_summary = {}
    for el in ['E1','E2','E3']:
        members = [f for f in features if f['element']==el]
        el_summary[el] = {
            'name':       EL_NAMES[el],
            'n_residues': len(members),
            'centroid':   centroids.get(el, {}),
            'members':    [f['id'] for f in members],
        }

    return {
        'name':       'CGCP-NSP9-NSP12-v1',
        'interface':  'NSP9-NSP12',
        'date':       '2026-03-23',
        'method':     'CGCP pipeline v1.0',
        'n_features': len(features),
        'elements':   el_summary,
        'features':   features,
        'screening_notes': (
            'Compounds targeting CGCP-NSP9-NSP12-v1 should: '
            '(1) present H-bond donor/acceptor groups to engage '
            'NSP9 N-terminus (ASN1/ASN2/ASN96 — E1); '
            '(2) carry positively charged group to engage '
            'ARG733/ARG735 electrostatic anchor (E2); '
            '(3) fit within 15A radius of E1 centroid. '
            'NSP9 N-terminal asparagine cluster is unique '
            'druggability feature of this interface.'
        ),
    }


def save_step7(pharmacophore):
    jpath = os.path.join(OUT7,
        'pharmacophore_NSP9-NSP12.json')
    with open(jpath, 'w') as f:
        json.dump(pharmacophore, f, indent=2, cls=NpEncoder)
    print(f"  Step7 JSON: {jpath}")

    tsv = os.path.join(OUT7,
        'pharmacophore_NSP9-NSP12.tsv')
    fields = ['id','chain','position','aa','element',
              'role','pharm_type','tolerance_A',
              'conservation','composite','x','y','z']
    with open(tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(pharmacophore['features'])
    print(f"  Step7 TSV: {tsv}")


# ── Figures ──────────────────────────────────────────────────
def dec_style(d):
    if d == 'ANCHOR':    return P['anchor'],    P['anchor_f']
    if d == 'INCLUDE':   return P['include'],   P['include_f']
    if d == 'SECONDARY': return P['secondary'], P['secondary_f']
    return P['exclude'], P['exclude_f']


def make_figures(assessed, out6, pharmacophore, centroids):
    features = pharmacophore['features']

    # ── Figure 1: Assessment ─────────────────────────────────
    fig1, axes1 = plt.subplots(1, 2, figsize=(14, 5.5))
    ax1, ax2 = axes1
    fig1.subplots_adjust(wspace=0.45, left=0.07,
                         right=0.97, top=0.88,
                         bottom=0.26)

    top = [r for r in assessed
           if r['decision'] != 'EXCLUDE'][:30]
    labels_a = [f"{r['aa'][:3]}{r['position']}"
                for r in top]
    x = np.arange(len(top))

    for xi, r in enumerate(top):
        ec, fc = dec_style(r['decision'])
        ax1.bar(xi, r['evidence_score'], width=0.40,
                facecolor=fc, edgecolor=ec,
                linewidth=1.2, zorder=3, clip_on=False)
        ax1.text(xi, r['evidence_score']+0.07,
                 str(r['evidence_score']),
                 ha='center', va='bottom',
                 fontsize=7, fontweight='bold',
                 color=ec, clip_on=False)
        badge = ('A' if r['decision']=='ANCHOR'
                 else 'I' if r['decision']=='INCLUDE'
                 else 'S')
        ax1.text(xi, -0.4, badge, ha='center', va='top',
                 fontsize=6.5, color=ec, clip_on=False)

    ax1.axhline(4, color=P['include'], linewidth=0.8,
                linestyle='--', alpha=0.7)
    ax1.axhline(3, color=P['secondary'], linewidth=0.8,
                linestyle='--', alpha=0.7)
    ax1.text(len(top)-0.5, 4.08, 'INCLUDE (4/6)',
             fontsize=7, color=P['include'],
             ha='right', clip_on=False)
    ax1.text(len(top)-0.5, 3.08, 'SECONDARY (3/6)',
             fontsize=7, color=P['secondary'],
             ha='right', clip_on=False)

    ax1.set_xticks(x)
    set_xticklabels_vertical(ax1, labels_a, fontsize=7)
    ax1.set_ylabel('Evidence score (max 6)',
                   fontsize=11, fontweight='bold')
    ax1.set_xlim(-0.6, len(top)-0.4)
    prism_axes(ax1, ymax=7, yticks=[0,1,2,3,4,5,6])
    make_legend(ax1, [
        (P['anchor_f'],    P['anchor'],    'ANCHOR'),
        (P['include_f'],   P['include'],   'INCLUDE'),
        (P['secondary_f'], P['secondary'], 'SECONDARY'),
    ], loc='upper right', fontsize=7.5)
    panel_title(ax1, 'A',
                'Integrated evidence scores (non-excluded)',
                'A=anchor  I=include  S=secondary')

    dec_order  = ['ANCHOR','INCLUDE','SECONDARY','EXCLUDE']
    dec_labels = ['Anchor','Include','Secondary','Exclude']
    dec_counts = [out6['decisions'].get(d,0)
                  for d in dec_order]
    x2 = np.arange(len(dec_order))
    for xi, val, d in zip(x2, dec_counts, dec_order):
        ec, fc = dec_style(d)
        ax2.bar(xi, val, width=0.50,
                facecolor=fc, edgecolor=ec,
                linewidth=1.2, zorder=3, clip_on=False)
        ax2.text(xi, val+0.3, str(val),
                 ha='center', va='bottom',
                 fontsize=10, fontweight='bold',
                 color=ec, clip_on=False)
    ax2.set_xticks(x2)
    set_xticklabels_vertical(ax2, dec_labels, fontsize=9)
    ax2.set_ylabel('Number of residues',
                   fontsize=11, fontweight='bold')
    ymax_b = max(dec_counts)+6
    prism_axes(ax2, ymax=ymax_b,
               yticks=range(0, ymax_b, 5))
    ax2.set_xlim(-0.6, len(dec_order)-0.4)
    panel_title(ax2, 'B', 'Decision summary',
                'pharmacophore = ANCHOR + INCLUDE')

    p1 = os.path.join(OUT6,
        'Fig_Step06_Assessment_NSP9-NSP12.png')
    fig1.savefig(p1)
    plt.close()
    print(f"  Fig Assessment: {p1}")

    # ── Figure 2: Pharmacophore ───────────────────────────────
    fig2 = plt.figure(figsize=(16.0, 6.0))
    gs   = fig2.add_gridspec(1, 3, wspace=0.48,
                              left=0.06, right=0.97,
                              top=0.88, bottom=0.26)
    ax_a = fig2.add_subplot(gs[0,0])
    ax_b = fig2.add_subplot(gs[0,1])
    ax_c = fig2.add_subplot(gs[0,2], polar=True)

    # Panel A: 2D X-Z spatial map
    for el_id, ct in centroids.items():
        col  = EL_COLORS[el_id]
        fill = EL_FILLS[el_id]
        circ = Circle((ct['x'], ct['z']),
                      radius=ct['radius'],
                      facecolor=fill, edgecolor=col,
                      linewidth=1.2, alpha=0.20, zorder=1)
        ax_a.add_patch(circ)
        ax_a.scatter(ct['x'], ct['z'], marker='+',
                     color=col, s=180, linewidths=2.5,
                     zorder=5)
        ax_a.text(ct['x'],
                  ct['z']+ct['radius']+1.0,
                  el_id, ha='center', va='bottom',
                  fontsize=10, color=col,
                  fontweight='bold')

    for f in features:
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role'] == 'ANCHOR'
        size   = 220 if is_anc else 50
        marker = '*'  if is_anc else 'o'
        ax_a.scatter(f['x'], f['z'], color=fill,
                     edgecolors=col, linewidths=1.0,
                     s=size, marker=marker,
                     alpha=0.88, zorder=4)
        if is_anc or f['composite'] >= 0.750:
            ax_a.text(f['x']+1.5, f['z']+0.8,
                      f['id'].split('-')[-1],
                      fontsize=7.5, color=col,
                      fontweight='bold' if is_anc else 'normal',
                      va='center')

    ax_a.set_xlabel('X coordinate (Å)',
                    fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Z coordinate (Å)',
                    fontsize=11, fontweight='bold')
    ax_a.spines['left'].set_position(('outward', 6))
    ax_a.spines['bottom'].set_position(('outward', 6))
    ax_a.spines['left'].set_linewidth(1.2)
    ax_a.spines['bottom'].set_linewidth(1.2)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.tick_params(labelsize=9, width=1.2,
                     length=6, direction='out')
    make_legend(ax_a, [
        (EL_FILLS[el], EL_COLORS[el], EL_NAMES[el])
        for el in ['E1','E2','E3']
    ], loc='best', fontsize=7.5)
    panel_title(ax_a, 'A',
                '2D pharmacophore map (X–Z)',
                '★=ARG733 anchor')

    # Panel B: composite scores
    feats_s = sorted(features, key=lambda x: -x['composite'])
    labels_b = [f['id'].split('-')[-1] for f in feats_s]
    x3 = np.arange(len(feats_s))
    for xi, f in enumerate(feats_s):
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role']=='ANCHOR'
        ax_b.bar(xi, f['composite'], width=0.38,
                 facecolor=fill, edgecolor=col,
                 linewidth=2.0 if is_anc else 1.2,
                 zorder=3, clip_on=False)
        ax_b.errorbar(xi, f['composite'],
                      yerr=f['tolerance_A']/10.0,
                      color=col, linewidth=0.9,
                      capsize=3, capthick=0.9, zorder=6)

    ax_b.set_xticks(x3)
    set_xticklabels_vertical(ax_b, labels_b, fontsize=7)
    ax_b.set_ylabel('CGCP composite score',
                    fontsize=11, fontweight='bold')
    ax_b.set_xlim(-0.6, len(feats_s)-0.4)
    prism_axes(ax_b, ymax=1.25,
               yticks=[0,0.2,0.4,0.6,0.8,1.0])
    ax_b.tick_params(labelsize=9, width=1.2,
                     length=6, direction='out')
    make_legend(ax_b, [
        (EL_FILLS[el], EL_COLORS[el], EL_NAMES[el])
        for el in ['E1','E2','E3']
    ], loc='upper right', fontsize=7.5)
    panel_title(ax_b, 'B',
                'Composite score per feature',
                'error bars = tolerance ÷ 10')

    # Panel C: polar wheel
    ax_c.set_facecolor('white')
    n      = len(features)
    angles = np.linspace(0,2*np.pi,n,endpoint=False)
    theta  = np.linspace(0,2*np.pi,100)
    for r_val, col in [(1.00,'#D7191C'),
                       (0.75,'#2166AC'),
                       (0.50,'#4DAC26')]:
        ax_c.plot(theta,[r_val]*100,color=col,
                  linewidth=0.7,linestyle='--',alpha=0.35)

    for i, f in enumerate(features):
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role']=='ANCHOR'
        ax_c.bar(angles[i], f['composite'],
                 width=2*np.pi/n*0.70,
                 facecolor=fill, edgecolor=col,
                 linewidth=1.0, alpha=0.88, zorder=3)
        ax_c.scatter(angles[i], f['composite'],
                     color=col,
                     s=80 if is_anc else 25,
                     marker='*' if is_anc else 'o',
                     zorder=5, edgecolors='white',
                     linewidths=0.5)
        rot = np.degrees(angles[i])
        if angles[i]>np.pi: rot-=180
        ax_c.text(angles[i], f['composite']+0.14,
                  f['id'].split('-')[-1],
                  ha='center', va='center',
                  fontsize=6.5, color=col,
                  fontweight='bold' if is_anc else 'normal',
                  rotation=rot)

    for r_val, lbl in [(0.25,'0.25'),(0.50,'0.50'),
                       (0.75,'0.75'),(1.00,'1.00')]:
        ax_c.text(0,r_val,lbl,ha='center',va='bottom',
                  fontsize=6.5,color='#AAAAAA')
    ax_c.set_ylim(0,1.45)
    ax_c.set_yticks([]); ax_c.set_xticks([])
    ax_c.spines['polar'].set_visible(False)
    make_legend(ax_c, [
        (EL_FILLS[el], EL_COLORS[el], EL_NAMES[el])
        for el in ['E1','E2','E3']
    ], loc='lower center', fontsize=7.5)
    panel_title(ax_c, 'C',
                'Feature wheel\n(radius = composite; ★=anchor)',
                fontsize=9)

    p2 = os.path.join(OUT7,
        'Fig_Step07_Pharmacophore_NSP9-NSP12.png')
    fig2.savefig(p2)
    plt.close()
    print(f"  Fig Pharmacophore: {p2}")


def print_summary(assessed, out6, pharmacophore, centroids):
    print('\n' + '='*65)
    print('STEPS 6+7 — NSP9-NSP12')
    print('='*65)

    print(f"\n  STEP 6 — Decision table (non-excluded):")
    print(f"  {'Residue':<20} {'Chain':<7} "
          f"{'Cons':>6} {'Comp':>6} {'Score':>6} Decision")
    print(f"  {'-'*60}")
    for r in assessed:
        if r['decision'] == 'EXCLUDE': continue
        sym = ('★' if r['decision']=='ANCHOR'
               else '✓' if r['decision']=='INCLUDE'
               else '·')
        print(f"  {sym} {r['aa']}{r['position']:<14} "
              f"{r['chain']:<7} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f} "
              f"{r['evidence_score']:>4}/6  "
              f"{r['decision']}")

    print(f"\n  Pharmacophore: "
          f"{len(out6['pharmacophore_candidates'])} candidates")
    for res in out6['pharmacophore_candidates']:
        print(f"    → {res}")

    print(f"\n  STEP 7 — CGCP-NSP9-NSP12-v1:")
    for el_id in ['E1','E2','E3']:
        el  = pharmacophore['elements'][el_id]
        ct  = el['centroid']
        print(f"\n  {el_id}: {el['name']}")
        print(f"    Residues: {el['n_residues']} | "
              f"Centroid: ({ct.get('x',0):.2f}, "
              f"{ct.get('y',0):.2f}, "
              f"{ct.get('z',0):.2f}) | "
              f"Radius: {ct.get('radius',0):.2f}Å")
        print(f"    Members: {', '.join(el['members'])}")
    print('='*65)


if __name__ == '__main__':
    print('CGCP Phase 2 Steps 6+7 — NSP9-NSP12')

    # Step 6
    records  = load_conservation()
    assessed = assess(records)
    out6     = save_step6(assessed)

    # Step 7
    coords      = get_coords()
    centroids   = compute_centroids(coords)
    pharmacophore = build_pharmacophore(
        coords, centroids, assessed)
    save_step7(pharmacophore)

    # Figures
    make_figures(assessed, out6, pharmacophore, centroids)
    print_summary(assessed, out6, pharmacophore, centroids)

    print(f"\nStep6 outputs: {OUT6}")
    print(f"Step7 outputs: {OUT7}")
    print("Next: Step 8 — ChimeraX visualization")
