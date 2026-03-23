#!/usr/bin/env python3
"""
CGCP Phase 2 Step 7 - Pharmacophore Hypothesis: NSP12-NSP8
Formalizes the 21 pharmacophore candidates into named elements
with 3D coordinates, tolerances, and interaction types.

Three pharmacophoric elements:
  E1 — Hydrophobic/aromatic core
       LEU387(anchor), PHE368, TYR273, VAL330, LEU389, LEU329
       Feature: HYDROPHOBIC + AROMATIC
       Interaction: van der Waals, pi-stacking

  E2 — Electrostatic cluster
       ARG80(NSP8), LYS127(NSP8), ARG392, LYS332(anchor)
       Feature: POSITIVE CHARGE + NEG CHARGE
       Interaction: salt bridge network

  E3 — H-bond zone
       ASN386, SER518, THR324, ASN118(NSP8)
       Feature: H-BOND DONOR/ACCEPTOR
       Interaction: hydrogen bonding

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step07_pharmacophore_hypothesis_NSP12-NSP8.py
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

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP12-NSP8_4/receptor_NSP12-NSP8_4.pdb')
S6_DIR   = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-06-assessment')
OUT_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-07-pharmacophore')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Pharmacophore definition ─────────────────────────────────
# Based on Step 6 decisions
# (chain, resnum): {element, role, pharm_type, tolerance, interaction, cons, composite}
PHARMACOPHORE_DEF = {
    # E1 — Hydrophobic/aromatic core (NSP12)
    ('A', 387): {'element':'E1','role':'ANCHOR_PRIMARY',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':1.5,
                 'interaction':'vdW core — primary anchor',
                 'cons':1.000,'composite':1.000},
    ('A', 368): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'AROMATIC','tolerance_A':1.5,
                 'interaction':'pi-stacking',
                 'cons':1.000,'composite':0.881},
    ('A', 273): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'AROMATIC','tolerance_A':1.5,
                 'interaction':'pi-stacking + H-bond',
                 'cons':1.000,'composite':0.770},
    ('A', 330): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic',
                 'cons':1.000,'composite':0.867},
    ('A', 389): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic',
                 'cons':0.800,'composite':0.881},
    ('A', 329): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic',
                 'cons':1.000,'composite':0.606},
    # E2 — Electrostatic cluster
    ('B', 80):  {'element':'E2','role':'INCLUDE',
                 'pharm_type':'CHARGED_POS','tolerance_A':1.0,
                 'interaction':'salt bridge network (NSP8)',
                 'cons':1.000,'composite':0.931},
    ('B', 127): {'element':'E2','role':'INCLUDE',
                 'pharm_type':'CHARGED_POS','tolerance_A':1.0,
                 'interaction':'salt bridge (NSP8)',
                 'cons':1.000,'composite':0.781},
    ('A', 392): {'element':'E2','role':'INCLUDE',
                 'pharm_type':'CHARGED_POS','tolerance_A':1.0,
                 'interaction':'charged contact',
                 'cons':1.000,'composite':0.639},
    ('A', 332): {'element':'E2','role':'ANCHOR_SECONDARY',
                 'pharm_type':'CHARGED_POS','tolerance_A':1.0,
                 'interaction':'salt bridge LYS332-ASP99',
                 'cons':0.600,'composite':0.508},
    # E3 — H-bond zone
    ('A', 386): {'element':'E3','role':'INCLUDE',
                 'pharm_type':'HBOND','tolerance_A':1.5,
                 'interaction':'H-bond donor/acceptor',
                 'cons':1.000,'composite':0.727},
    ('A', 518): {'element':'E3','role':'INCLUDE',
                 'pharm_type':'HBOND','tolerance_A':1.5,
                 'interaction':'H-bond donor',
                 'cons':1.000,'composite':0.628},
    ('A', 324): {'element':'E3','role':'INCLUDE',
                 'pharm_type':'HBOND','tolerance_A':1.5,
                 'interaction':'H-bond donor',
                 'cons':1.000,'composite':0.606},
    ('B', 118): {'element':'E3','role':'INCLUDE',
                 'pharm_type':'HBOND','tolerance_A':1.5,
                 'interaction':'H-bond acceptor (NSP8)',
                 'cons':0.800,'composite':0.743},
    # Additional INCLUDE residues
    ('B', 91):  {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic (NSP8)',
                 'cons':1.000,'composite':0.672},
    ('A', 371): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic',
                 'cons':0.800,'composite':0.645},
    ('B', 119): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic (NSP8)',
                 'cons':1.000,'composite':0.617},
    ('B', 149): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'AROMATIC','tolerance_A':1.5,
                 'interaction':'pi-stacking (NSP8)',
                 'cons':1.000,'composite':0.617},
    ('B', 116): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic (NSP8)',
                 'cons':0.800,'composite':0.792},
    ('B', 87):  {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic (NSP8)',
                 'cons':0.800,'composite':0.754},
    ('B', 115): {'element':'E1','role':'INCLUDE',
                 'pharm_type':'HYDROPHOBIC','tolerance_A':2.0,
                 'interaction':'vdW hydrophobic (NSP8)',
                 'cons':0.800,'composite':0.732},
}

EL_COLORS = {
    'E1': '#D7191C', 'E2': '#2166AC', 'E3': '#4DAC26',
}
EL_FILLS  = {
    'E1': '#FDAE61', 'E2': '#92C5DE', 'E3': '#A6D96A',
}
EL_NAMES  = {
    'E1': 'E1 — Hydrophobic/aromatic core',
    'E2': 'E2 — Electrostatic cluster',
    'E3': 'E3 — H-bond zone',
}
FEAT_COLORS = {
    'HYDROPHOBIC':  '#636363',
    'AROMATIC':     '#D7191C',
    'CHARGED_POS':  '#2166AC',
    'CHARGED_NEG':  '#4DAC26',
    'HBOND':        '#8B4513',
}
FEAT_FILLS = {k: v + '55' for k, v in FEAT_COLORS.items()}


def get_coords():
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]
    coords    = {}

    chain_map = {'A': structure['A'], 'B': structure['B']}
    for (chain_id, pos), defn in PHARMACOPHORE_DEF.items():
        chain = chain_map[chain_id]
        try:
            res = chain[(' ', pos, ' ')]
            ca  = res['CA'].coord
            coords[(chain_id, pos)] = {
                'x': round(float(ca[0]), 3),
                'y': round(float(ca[1]), 3),
                'z': round(float(ca[2]), 3),
            }
        except KeyError:
            print(f"  Warning: chain {chain_id} res {pos} not found")

    return coords


def compute_centroids(coords):
    elements = {}
    for (chain_id, pos), defn in PHARMACOPHORE_DEF.items():
        el = defn['element']
        if el not in elements:
            elements[el] = []
        if (chain_id, pos) in coords:
            elements[el].append(coords[(chain_id, pos)])

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
            'x': round(cx, 3), 'y': round(cy, 3),
            'z': round(cz, 3), 'radius': round(radius, 2),
        }
    return centroids


def build_pharmacophore(coords, centroids):
    features = []
    for (chain_id, pos), defn in PHARMACOPHORE_DEF.items():
        if (chain_id, pos) not in coords:
            continue
        chain_name = 'NSP12' if chain_id == 'A' else 'NSP8'
        aa         = ''
        # Get AA from PDB
        try:
            parser = PDB.PDBParser(QUIET=True)
            s = parser.get_structure('cx', PDB_FILE)[0]
            c = s['A'] if chain_id == 'A' else s['B']
            aa = c[(' ', pos, ' ')].resname
        except Exception:
            pass

        features.append({
            'id':           f"{chain_name}-{aa}{pos}",
            'chain':        chain_name,
            'chain_id':     chain_id,
            'position':     pos,
            'aa':           aa,
            'element':      defn['element'],
            'role':         defn['role'],
            'pharm_type':   defn['pharm_type'],
            'tolerance_A':  defn['tolerance_A'],
            'interaction':  defn['interaction'],
            'conservation': defn['cons'],
            'composite':    defn['composite'],
            'x':            coords[(chain_id, pos)]['x'],
            'y':            coords[(chain_id, pos)]['y'],
            'z':            coords[(chain_id, pos)]['z'],
        })

    features.sort(key=lambda x: (
        x['element'],
        0 if 'ANCHOR' in x['role'] else 1,
        -x['composite'],
    ))

    el_summary = {}
    for el in ['E1','E2','E3']:
        members = [f for f in features if f['element']==el]
        el_summary[el] = {
            'name':        EL_NAMES[el],
            'n_residues':  len(members),
            'centroid':    centroids.get(el, {}),
            'members':     [f['id'] for f in members],
        }

    pharmacophore = {
        'name':        'CGCP-NSP12-NSP8-v1',
        'interface':   'NSP12-NSP8',
        'date':        '2026-03-23',
        'method':      'CGCP pipeline v1.0',
        'reference_pdb': 'NSP12-NSP8_4',
        'n_features':  len(features),
        'n_elements':  3,
        'elements':    el_summary,
        'features':    features,
        'screening_notes': (
            'Compounds targeting CGCP-NSP12-NSP8-v1 should: '
            '(1) engage E1 hydrophobic/aromatic core '
            '(LEU387/PHE368/TYR273) with aromatic/hydrophobic scaffold; '
            '(2) carry positively charged or H-bond donor group '
            'to engage E2 electrostatic cluster (ARG80/LYS127 NSP8); '
            '(3) optionally present H-bond acceptor for E3 zone '
            '(ASN386/SER518). '
            'Unlike NSP12-NSP7, this interface requires both '
            'hydrophobic AND electrostatic features.'
        ),
    }
    return pharmacophore


def save_outputs(pharmacophore):
    json_path = os.path.join(OUT_DIR,
        'pharmacophore_NSP12-NSP8.json')
    with open(json_path, 'w') as f:
        json.dump(pharmacophore, f, indent=2)
    print(f"  JSON: {json_path}")

    tsv_path = os.path.join(OUT_DIR,
        'pharmacophore_NSP12-NSP8.tsv')
    fields = ['id','chain','position','aa','element',
              'role','pharm_type','tolerance_A',
              'interaction','conservation','composite',
              'x','y','z']
    with open(tsv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(pharmacophore['features'])
    print(f"  TSV: {tsv_path}")


def make_figure(pharmacophore, centroids):
    features = pharmacophore['features']
    fig = plt.figure(figsize=(16.0, 6.5))
    gs  = fig.add_gridspec(1, 3, wspace=0.48,
                           left=0.06, right=0.97,
                           top=0.88, bottom=0.26)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2], polar=True)

    # ── Panel A: 2D spatial map (X-Z) ────────────────────────
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
        ax_a.text(ct['x'], ct['z'] + ct['radius'] + 0.8,
                  el_id, ha='center', va='bottom',
                  fontsize=10, color=col,
                  fontweight='bold')

    for f in features:
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_pa = f['role'] == 'ANCHOR_PRIMARY'
        is_sa = f['role'] == 'ANCHOR_SECONDARY'
        size   = 220 if is_pa else (130 if is_sa else 50)
        marker = '*'  if is_pa else ('D' if is_sa else 'o')
        ax_a.scatter(f['x'], f['z'], color=fill,
                     edgecolors=col, linewidths=1.0,
                     s=size, marker=marker,
                     alpha=0.88, zorder=4)
        if is_pa or is_sa or f['composite'] >= 0.800:
            ox = 1.5
            if f['x'] < centroids.get(
                    f['element'],{}).get('x', f['x']):
                ox = -3.5
            ax_a.text(f['x'] + ox, f['z'] + 0.8,
                      f['id'].split('-')[-1],
                      fontsize=7.5, color=col,
                      fontweight='bold' if is_pa else 'normal',
                      va='center')

    ax_a.set_xlabel('X coordinate (Å)', fontsize=11,
                    fontweight='bold')
    ax_a.set_ylabel('Z coordinate (Å)', fontsize=11,
                    fontweight='bold')
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
    ] + [('#FFFFFF','#1A1A1A','★ = primary anchor'),
         ('#FFFFFF','#636363','◆ = secondary anchor')],
    loc='lower right', fontsize=7.5)
    panel_title(ax_a, 'A',
                '2D spatial projection (X–Z plane)',
                'circle = element zone; cross = centroid')

    # ── Panel B: Composite score per feature ─────────────────
    feats_s = sorted(features, key=lambda x: -x['composite'])
    labels_b = [f['id'].split('-')[-1] for f in feats_s]
    vals_b   = [f['composite'] for f in feats_s]
    cols_b   = [EL_COLORS[f['element']] for f in feats_s]
    fills_b  = [EL_FILLS[f['element']] for f in feats_s]

    x2 = np.arange(len(feats_s))
    for xi, val, fc, ec, f in zip(x2, vals_b,
                                   fills_b, cols_b, feats_s):
        is_pa = f['role'] == 'ANCHOR_PRIMARY'
        is_sa = f['role'] == 'ANCHOR_SECONDARY'
        lw    = 2.0 if is_pa else (1.5 if is_sa else 1.2)
        ax_b.bar(xi, val, width=0.38,
                 facecolor=fc, edgecolor=ec,
                 linewidth=lw, zorder=3, clip_on=False)
        ax_b.scatter(xi, val, color=ec, s=20,
                     zorder=5, clip_on=False)
        # Tolerance error bar
        ax_b.errorbar(xi, val,
                      yerr=f['tolerance_A'] / 10.0,
                      color=ec, linewidth=0.9,
                      capsize=3, capthick=0.9,
                      zorder=6)

    ax_b.set_xticks(x2)
    set_xticklabels_vertical(ax_b, labels_b, fontsize=6.5)
    ax_b.set_ylabel('CGCP composite score',
                    fontsize=11, fontweight='bold')
    ax_b.set_xlim(-0.6, len(feats_s) - 0.4)
    prism_axes(ax_b, ymax=1.25,
               yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_b.tick_params(labelsize=9, width=1.2,
                     length=6, direction='out')
    make_legend(ax_b, [
        (EL_FILLS[el], EL_COLORS[el], EL_NAMES[el])
        for el in ['E1','E2','E3']
    ], loc='upper right', fontsize=7.5)
    panel_title(ax_b, 'B',
                'Composite score per feature',
                'error bars = spatial tolerance ÷ 10')

    # ── Panel C: Feature wheel (polar) ───────────────────────
    ax_c.set_facecolor('white')
    n      = len(features)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    theta  = np.linspace(0, 2*np.pi, 100)

    for r_val, col in [(1.00,'#D7191C'),
                       (0.75,'#2166AC'),
                       (0.50,'#4DAC26')]:
        ax_c.plot(theta, [r_val]*100,
                  color=col, linewidth=0.7,
                  linestyle='--', alpha=0.35)

    for i, f in enumerate(features):
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_pa = f['role'] == 'ANCHOR_PRIMARY'
        is_sa = f['role'] == 'ANCHOR_SECONDARY'
        ax_c.bar(angles[i], f['composite'],
                 width=2*np.pi/n*0.70,
                 facecolor=fill, edgecolor=col,
                 linewidth=1.0, alpha=0.88, zorder=3)
        ax_c.scatter(angles[i], f['composite'],
                     color=col,
                     s=80 if is_pa else (50 if is_sa else 25),
                     marker='*' if is_pa else ('D' if is_sa else 'o'),
                     zorder=5, edgecolors='white',
                     linewidths=0.5)
        label_r = f['composite'] + 0.14
        rot = np.degrees(angles[i])
        if angles[i] > np.pi:
            rot -= 180
        ax_c.text(angles[i], label_r,
                  f['id'].split('-')[-1],
                  ha='center', va='center',
                  fontsize=6.5, color=col,
                  fontweight='bold' if is_pa else 'normal',
                  rotation=rot)

    for r_val, lbl in [(0.25,'0.25'),(0.50,'0.50'),
                       (0.75,'0.75'),(1.00,'1.00')]:
        ax_c.text(0, r_val, lbl, ha='center',
                  va='bottom', fontsize=6.5,
                  color='#AAAAAA')

    ax_c.set_ylim(0, 1.45)
    ax_c.set_yticks([])
    ax_c.set_xticks([])
    ax_c.spines['polar'].set_visible(False)
    make_legend(ax_c, [
        (EL_FILLS[el], EL_COLORS[el], EL_NAMES[el])
        for el in ['E1','E2','E3']
    ], loc='lower center',
    fontsize=7.5)
    panel_title(ax_c, 'C',
                'Pharmacophore feature wheel',
                'radius = composite score; ★=primary anchor',
                fontsize=9)

    path = os.path.join(OUT_DIR,
        'Fig_Step07_PharmacophorHypothesis_NSP12-NSP8.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(pharmacophore, centroids):
    print('\n' + '='*65)
    print(f"STEP 7 PHARMACOPHORE HYPOTHESIS — NSP12-NSP8")
    print(f"Name: {pharmacophore['name']}")
    print('='*65)

    for el_id in ['E1','E2','E3']:
        el_def  = pharmacophore['elements'][el_id]
        ct      = el_def['centroid']
        members = [f for f in pharmacophore['features']
                   if f['element'] == el_id]
        print(f"\n  {el_id}: {el_def['name']}")
        print(f"    Residues: {el_def['n_residues']}")
        print(f"    Centroid: ({ct.get('x',0):.2f}, "
              f"{ct.get('y',0):.2f}, "
              f"{ct.get('z',0):.2f})")
        print(f"    Radius:   {ct.get('radius',0):.2f} Å")
        print(f"    Members:")
        for f in members:
            sym = ('★' if f['role']=='ANCHOR_PRIMARY'
                   else '◆' if f['role']=='ANCHOR_SECONDARY'
                   else ' ')
            print(f"      {sym} {f['id']:<22} "
                  f"{f['pharm_type']:<14} "
                  f"tol={f['tolerance_A']}Å  "
                  f"cons={f['conservation']:.3f}  "
                  f"comp={f['composite']:.3f}")

    print(f"\n  Screening notes:")
    print(f"  {pharmacophore['screening_notes']}")
    print('='*65)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 7 — NSP12-NSP8 Pharmacophore Hypothesis')
    coords      = get_coords()
    centroids   = compute_centroids(coords)
    pharmacophore = build_pharmacophore(coords, centroids)
    save_outputs(pharmacophore)
    make_figure(pharmacophore, centroids)
    print_summary(pharmacophore, centroids)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 8 — controls + ChimeraX figures")
