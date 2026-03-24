#!/usr/bin/env python3
"""
CGCP Phase 2 Step 7 fix — NSP13-Helicase and NSP15
Corrects element assignment for homodimer interfaces where
all residues are at moderate conservation (0.750).

For homodimers:
  E1 = anchor + top same-chain residues by composite
  E2 = partner chain residues
  E3 = remaining same-chain residues

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step07_fix_homodimer_elements.py
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
                          COLORS)

apply_prism()

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):  return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray):  return obj.tolist()
        if isinstance(obj, np.bool_):    return bool(obj)
        return super().default(obj)

BASE = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VS   = os.path.join(BASE, '03-virtual-screening')

EL_COLORS = {'E1':'#D7191C','E2':'#2166AC','E3':'#4DAC26'}
EL_FILLS  = {'E1':'#FDAE61','E2':'#92C5DE','E3':'#A6D96A'}

HOMODIMER_CONFIGS = {
    'NSP13-Helicase': {
        'pdb':          f'{VS}/NSP13-Helicase_7/receptor_NSP13-Helicase_7.pdb',
        'chain1_id':    'A', 'name1': 'NSP13a',
        'chain2_id':    'E', 'name2': 'NSP13b',
        'anchor_chain': 'NSP13a', 'anchor_res': 414,
        # E1 = anchor + top 4 same-chain by composite
        'e1_chain1_top': 4,
        'el_names': {
            'E1': 'E1 — LYS414 salt bridge + ILE480/HIS482 core',
            'E2': 'E2 — NSP13b symmetric contacts',
            'E3': 'E3 — NSP13a extended contacts',
        },
    },
    'NSP15': {
        'pdb':          f'{VS}/NSP15_9/receptor_NSP15_9.pdb',
        'chain1_id':    'A', 'name1': 'NSP15a',
        'chain2_id':    'B', 'name2': 'NSP15b',
        'anchor_chain': 'NSP15a', 'anchor_res': 40,
        # E1 = anchor + top 5 same-chain by composite
        'e1_chain1_top': 5,
        'el_names': {
            'E1': 'E1 — ASP40 + GLU171/ALA172 homodimer core',
            'E2': 'E2 — NSP15b symmetric contacts',
            'E3': 'E3 — NSP15a extended contacts',
        },
    },
}


def load_assessed(iface):
    path = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-06-assessment/'
        f'integrated_assessment_{iface}.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            records.append(row)
    return records


def get_coords(assessed, cfg):
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', cfg['pdb'])[0]
    chain1    = structure[cfg['chain1_id']]
    chain2    = structure[cfg['chain2_id']]

    coords = {}
    for r in assessed:
        chain_obj = (chain1 if r['chain']==cfg['name1']
                     else chain2)
        pos = int(r['position'])
        try:
            res = chain_obj[(' ', pos, ' ')]
            ca  = res['CA'].coord
            coords[(r['chain'], pos)] = {
                'x': round(float(ca[0]),3),
                'y': round(float(ca[1]),3),
                'z': round(float(ca[2]),3),
                'aa': res.resname,
            }
        except KeyError:
            pass
    return coords


def assign_homodimer_elements(assessed, coords, cfg):
    """
    For homodimer interfaces:
    E1 = anchor + top N residues from chain1 (same chain as anchor)
    E2 = all chain2 (partner chain) INCLUDE residues
    E3 = remaining chain1 INCLUDE residues
    """
    features = []
    include_res = [r for r in assessed
                   if r['decision'] in ('ANCHOR','INCLUDE')]

    # Separate by chain
    chain1_res = [r for r in include_res
                  if r['chain'] == cfg['name1']]
    chain2_res = [r for r in include_res
                  if r['chain'] == cfg['name2']]

    # Sort chain1 by composite desc (anchor first)
    chain1_res.sort(key=lambda x: (
        0 if x['primary_feature']=='anchor' else 1,
        -float(x['composite'])))

    # E1 = anchor + top N chain1
    n_e1   = cfg['e1_chain1_top'] + 1  # +1 for anchor
    e1_res = chain1_res[:n_e1]
    e3_res = chain1_res[n_e1:]

    # Assign elements
    def make_feature(r, el):
        key = (r['chain'], int(r['position']))
        if key not in coords:
            return None
        c    = coords[key]
        feat = r['primary_feature']
        pharm_type = (
            'ANCHOR'      if feat=='anchor'      else
            'AROMATIC'    if feat=='aromatic'     else
            'HYDROPHOBIC' if feat=='hydrophobic'  else
            'CHARGED_POS' if feat=='charged_pos'  else
            'CHARGED_NEG' if feat=='charged_neg'  else
            'HBOND'       if feat=='hbond_donor'  else
            'UNKNOWN')
        tol = (1.0 if pharm_type in ('ANCHOR','CHARGED_POS','CHARGED_NEG')
               else 1.5 if pharm_type=='AROMATIC'
               else 2.0)
        return {
            'id':           r['residue'],
            'chain':        r['chain'],
            'position':     int(r['position']),
            'aa':           c['aa'],
            'element':      el,
            'role':         r['decision'],
            'pharm_type':   pharm_type,
            'tolerance_A':  tol,
            'conservation': float(r['conservation']),
            'composite':    float(r['composite']),
            'x':            c['x'],
            'y':            c['y'],
            'z':            c['z'],
        }

    for r in e1_res:
        f = make_feature(r, 'E1')
        if f: features.append(f)

    for r in chain2_res:
        f = make_feature(r, 'E2')
        if f: features.append(f)

    for r in e3_res:
        f = make_feature(r, 'E3')
        if f: features.append(f)

    features.sort(key=lambda x:(
        x['element'],
        0 if x['role']=='ANCHOR' else 1,
        -x['composite']))
    return features


def compute_centroids(features):
    centroids = {}
    for el in ['E1','E2','E3']:
        pts = [f for f in features if f['element']==el]
        if not pts:
            continue
        cx = float(np.mean([p['x'] for p in pts]))
        cy = float(np.mean([p['y'] for p in pts]))
        cz = float(np.mean([p['z'] for p in pts]))
        radius = float(max(
            np.sqrt((p['x']-cx)**2+(p['y']-cy)**2+(p['z']-cz)**2)
            for p in pts))
        centroids[el] = {
            'x': round(cx,3), 'y': round(cy,3),
            'z': round(cz,3), 'radius': round(radius,2),
        }
    return centroids


def save_pharmacophore(iface, features, centroids, cfg):
    out7 = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-07-pharmacophore')
    os.makedirs(out7, exist_ok=True)

    el_summary = {}
    for el in ['E1','E2','E3']:
        members = [f for f in features if f['element']==el]
        el_summary[el] = {
            'name':       cfg['el_names'][el],
            'n_residues': len(members),
            'centroid':   centroids.get(el,{}),
            'members':    [f['id'] for f in members],
        }

    pharmacophore = {
        'name':       f'CGCP-{iface}-v1',
        'interface':  iface,
        'date':       '2026-03-24',
        'method':     'CGCP pipeline v1.0',
        'note':       'Homodimer — E1 top chain1 residues, E2 chain2, E3 remaining',
        'n_features': len(features),
        'elements':   el_summary,
        'features':   features,
    }

    jpath = os.path.join(out7, f'pharmacophore_{iface}.json')
    with open(jpath, 'w') as f:
        json.dump(pharmacophore, f, indent=2, cls=NpEncoder)
    print(f"  JSON: {jpath}")

    tsv = os.path.join(out7, f'pharmacophore_{iface}.tsv')
    fields = ['id','chain','position','aa','element',
              'role','pharm_type','tolerance_A',
              'conservation','composite','x','y','z']
    with open(tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(features)
    print(f"  TSV: {tsv}")

    return pharmacophore, out7


def make_figure(iface, features, centroids,
                pharmacophore, out7, cfg):
    fig = plt.figure(figsize=(14, 5.5))
    gs  = fig.add_gridspec(1, 2, wspace=0.45,
                            left=0.07, right=0.97,
                            top=0.88, bottom=0.26)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # Panel A: 2D X-Z
    for el_id, ct in centroids.items():
        col  = EL_COLORS[el_id]
        fill = EL_FILLS[el_id]
        circ = Circle((ct['x'], ct['z']), radius=ct['radius'],
                      facecolor=fill, edgecolor=col,
                      linewidth=1.2, alpha=0.20, zorder=1)
        ax_a.add_patch(circ)
        ax_a.scatter(ct['x'], ct['z'], marker='+',
                     color=col, s=180, linewidths=2.5, zorder=5)
        ax_a.text(ct['x'], ct['z']+ct['radius']+1.0,
                  el_id, ha='center', va='bottom',
                  fontsize=10, color=col, fontweight='bold')

    for f in features:
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role'] == 'ANCHOR'
        ax_a.scatter(f['x'], f['z'], color=fill,
                     edgecolors=col, linewidths=1.0,
                     s=220 if is_anc else 50,
                     marker='*' if is_anc else 'o',
                     alpha=0.88, zorder=4)
        if is_anc or f['composite'] >= 0.700:
            ax_a.text(f['x']+1.0, f['z']+0.8,
                      f['id'].split('-')[1],
                      fontsize=7.5, color=col,
                      fontweight='bold' if is_anc else 'normal',
                      va='center')

    ax_a.set_xlabel('X coordinate (Å)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Z coordinate (Å)', fontsize=11, fontweight='bold')
    ax_a.spines['left'].set_position(('outward', 6))
    ax_a.spines['bottom'].set_position(('outward', 6))
    ax_a.spines['left'].set_linewidth(1.2)
    ax_a.spines['bottom'].set_linewidth(1.2)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.tick_params(labelsize=9, width=1.2,
                     length=6, direction='out')
    active_els = [el for el in ['E1','E2','E3']
                  if el in centroids]
    make_legend(ax_a, [
        (EL_FILLS[el], EL_COLORS[el],
         pharmacophore['elements'][el]['name'])
        for el in active_els
    ], loc='best', fontsize=7.5)
    panel_title(ax_a, 'A',
                f'2D pharmacophore map — {iface} (fixed)',
                '★=anchor')

    # Panel B: composite scores per element
    feats_s = sorted(features, key=lambda x: -x['composite'])
    x3 = np.arange(len(feats_s))
    for xi, f in enumerate(feats_s):
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role'] == 'ANCHOR'
        ax_b.bar(xi, f['composite'], width=0.38,
                 facecolor=fill, edgecolor=col,
                 linewidth=2.0 if is_anc else 1.2,
                 zorder=3, clip_on=False)
        ax_b.errorbar(xi, f['composite'],
                      yerr=f['tolerance_A']/10.0,
                      color=col, linewidth=0.9,
                      capsize=3, capthick=0.9, zorder=6)

    ax_b.set_xticks(x3)
    set_xticklabels_vertical(
        ax_b, [f['id'].split('-')[1] for f in feats_s],
        fontsize=7)
    ax_b.set_ylabel('CGCP composite score',
                    fontsize=11, fontweight='bold')
    ax_b.set_xlim(-0.6, len(feats_s)-0.4)
    prism_axes(ax_b, ymax=1.25, yticks=[0,0.2,0.4,0.6,0.8,1.0])
    ax_b.tick_params(labelsize=9, width=1.2,
                     length=6, direction='out')
    make_legend(ax_b, [
        (EL_FILLS[el], EL_COLORS[el],
         pharmacophore['elements'][el]['name'])
        for el in active_els
    ], loc='upper right', fontsize=7.5)
    panel_title(ax_b, 'B', 'Composite score per feature')

    path = os.path.join(out7,
        f'Fig_Step07_Pharmacophore_{iface}_fixed.png')
    fig.savefig(path); plt.close()
    print(f"  Figure: {path}")


def print_summary(iface, features, centroids,
                  pharmacophore, cfg):
    print(f"\n{'='*60}")
    print(f"CGCP-{iface}-v1 (FIXED element assignment)")
    print(f"{'='*60}")
    for el_id in ['E1','E2','E3']:
        el = pharmacophore['elements'][el_id]
        ct = el['centroid']
        if el['n_residues'] == 0:
            continue
        print(f"\n  {el_id}: {el['name']}")
        print(f"    {el['n_residues']} residues | "
              f"centroid ({ct.get('x',0):.2f}, "
              f"{ct.get('y',0):.2f}, "
              f"{ct.get('z',0):.2f}) | "
              f"radius {ct.get('radius',0):.2f}Å")
        for m in el['members']:
            f = next(x for x in features if x['id']==m)
            sym = '★' if f['role']=='ANCHOR' else ' '
            print(f"    {sym} {m:<26} "
                  f"cons={f['conservation']:.3f} "
                  f"comp={f['composite']:.3f}")


if __name__ == '__main__':
    print('CGCP Step 7 Fix — Homodimer element assignment')
    print('Interfaces: NSP13-Helicase, NSP15')

    for iface, cfg in HOMODIMER_CONFIGS.items():
        print(f"\nProcessing {iface}...")
        assessed  = load_assessed(iface)
        coords    = get_coords(assessed, cfg)
        features  = assign_homodimer_elements(
            assessed, coords, cfg)
        centroids = compute_centroids(features)
        pharmacophore, out7 = save_pharmacophore(
            iface, features, centroids, cfg)
        make_figure(iface, features, centroids,
                    pharmacophore, out7, cfg)
        print_summary(iface, features, centroids,
                      pharmacophore, cfg)

    print('\nDone — regenerate ChimeraX Step 8 for both interfaces')
