#!/usr/bin/env python3
"""
CGCP Phase 2 Step 1 - Structural Verification: NSP9-NSP12
Anchor: ARG733 (NSP12, Chain A, cons=1.000)
NSP9: Chain G

7 checks:
  C1: Chains A (NSP12) + G (NSP9) present
  C2: Chain sizes reasonable
  C3: ARG733 adjacent to NSP9 (< 15A)
  C4: Backbone completeness
  C5: ARG733 present and correct
  C6: Positive control — ARG733 to NSP9 interface
  C7: Negative control — ARG733 far from catalytic triad

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step01_structural_verification_NSP9-NSP12.py
"""

import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from Bio import PDB
from prism_style import apply_prism, prism_axes, COLORS

apply_prism()

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP9-NSP12_5/receptor_NSP9-NSP12_5.pdb')
OUT_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-01-structure')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'visuals'), exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'pymol-sessions'), exist_ok=True)

ANCHOR_RES   = 733
ANCHOR_AA    = 'ARG'
ANCHOR_CHAIN = 'A'
NSP9_CHAIN   = 'G'
NEG_CTRL     = [618, 759, 760]


def get_ca(chain, resnum):
    try:
        return chain[(' ', resnum, ' ')]['CA'].coord
    except KeyError:
        return None


def dist(a, b):
    return float(np.linalg.norm(
        np.array(a) - np.array(b)))


def run_checks(pdb_file):
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', pdb_file)[0]
    results   = {}

    # C1: chains A + G present
    chains = [c.id for c in structure.get_chains()]
    c1 = 'A' in chains and NSP9_CHAIN in chains
    results['C1'] = {
        'name':   'Chains A+G present',
        'pass':   c1,
        'detail': f"Chains found: {chains}",
    }

    chain_a = structure['A']
    chain_g = structure[NSP9_CHAIN]

    # C2: chain sizes
    na = sum(1 for r in chain_a.get_residues()
             if PDB.is_aa(r))
    ng = sum(1 for r in chain_g.get_residues()
             if PDB.is_aa(r))
    c2 = na > 500 and ng > 50
    results['C2'] = {
        'name':   'Chain sizes reasonable',
        'pass':   c2,
        'detail': f"NSP12={na} res, NSP9={ng} res",
    }

    # C3: ARG733 adjacent to NSP9
    ca_anchor = get_ca(chain_a, ANCHOR_RES)
    if ca_anchor is not None:
        min_dist = min(
            dist(ca_anchor, r['CA'].coord)
            for r in chain_g.get_residues()
            if PDB.is_aa(r) and 'CA' in r
        )
        c3 = min_dist < 15.0
        results['C3'] = {
            'name':   'Anchor adjacent to NSP9',
            'pass':   c3,
            'detail': f"ARG733 min dist to NSP9: {min_dist:.2f}Å",
        }
    else:
        results['C3'] = {
            'name':   'Anchor adjacent to NSP9',
            'pass':   False,
            'detail': 'ARG733 CA not found',
        }

    # C4: backbone completeness
    missing = []
    for r in chain_a.get_residues():
        if not PDB.is_aa(r):
            continue
        if 'CA' not in r or 'N' not in r or 'C' not in r:
            missing.append(r.id[1])
    c4 = len(missing) == 0
    results['C4'] = {
        'name':   'Backbone complete (Chain A)',
        'pass':   c4,
        'detail': f"{len(missing)} missing backbone atoms"
                  + (f": {missing[:5]}" if missing else ""),
    }

    # C5: ARG733 present and correct
    try:
        res733 = chain_a[(' ', ANCHOR_RES, ' ')]
        aa733  = res733.resname
        c5     = aa733 == ANCHOR_AA
        ca733  = res733['CA'].coord
        results['C5'] = {
            'name':   'ARG733 present and correct',
            'pass':   c5,
            'detail': f"Res 733 = {aa733} at "
                      f"({ca733[0]:.3f}, {ca733[1]:.3f}, "
                      f"{ca733[2]:.3f})",
        }
    except KeyError:
        results['C5'] = {
            'name':   'ARG733 present and correct',
            'pass':   False,
            'detail': 'Residue 733 not found',
        }

    # C6: positive control
    if ca_anchor is not None:
        nsp9_res = [r for r in chain_g.get_residues()
                    if PDB.is_aa(r) and 'CA' in r]
        close = [r for r in nsp9_res
                 if dist(ca_anchor, r['CA'].coord) < 15.0]
        c6 = len(close) >= 3
        results['C6'] = {
            'name':   'Positive ctrl: ARG733 near NSP9',
            'pass':   c6,
            'detail': f"{len(close)} NSP9 residues within 15Å",
        }
    else:
        results['C6'] = {
            'name':   'Positive ctrl: ARG733 near NSP9',
            'pass':   False,
            'detail': 'ARG733 not found',
        }

    # C7: negative control
    if ca_anchor is not None:
        triad_dists = []
        for pos in NEG_CTRL:
            ca_neg = get_ca(chain_a, pos)
            if ca_neg is not None:
                triad_dists.append(
                    (pos, dist(ca_anchor, ca_neg)))
        c7 = all(d > 20.0 for _, d in triad_dists)
        detail = " | ".join(
            [f"res{p}={d:.1f}Å" for p, d in triad_dists])
        results['C7'] = {
            'name':   'Neg ctrl: ARG733 far from catalytic triad',
            'pass':   c7,
            'detail': detail,
        }
    else:
        results['C7'] = {
            'name':   'Neg ctrl: ARG733 far from catalytic triad',
            'pass':   False,
            'detail': 'ARG733 not found',
        }

    return results


def save_report(results):
    n_pass  = sum(1 for r in results.values() if r['pass'])
    n_total = len(results)

    path = os.path.join(OUT_DIR, 'audit_report.md')
    with open(path, 'w') as f:
        f.write('# CGCP Step 1 Structural Audit — NSP9-NSP12\n\n')
        f.write(f'**Result: {n_pass}/{n_total} PASS**\n\n')
        f.write('| Check | Name | Result | Detail |\n')
        f.write('|-------|------|--------|--------|\n')
        for ck, r in results.items():
            sym = '✅' if r['pass'] else '❌'
            f.write(f"| {ck} | {r['name']} | "
                    f"{sym} {'PASS' if r['pass'] else 'FAIL'} | "
                    f"{r['detail']} |\n")
    print(f"  Report: {path}")

    path_json = os.path.join(OUT_DIR,
        'structural_verification.json')
    with open(path_json, 'w') as f:
        json.dump({
            'interface': 'NSP9-NSP12',
            'pdb':       PDB_FILE,
            'anchor':    f'{ANCHOR_AA}{ANCHOR_RES}',
            'n_pass':    n_pass,
            'n_total':   n_total,
            'checks':    results,
        }, f, indent=2)
    print(f"  JSON: {path_json}")
    return n_pass, n_total


def make_figure(results):
    checks = list(results.keys())
    colors = [COLORS['green'] if results[c]['pass']
              else COLORS['red'] for c in checks]
    fills  = [COLORS['green_fill'] if results[c]['pass']
              else COLORS['anchor_f'] for c in checks]

    fig, ax = plt.subplots(figsize=(10, 5.0))
    x = np.arange(len(checks))

    for xi, ck in enumerate(checks):
        r  = results[ck]
        ec = COLORS['green'] if r['pass'] else COLORS['red']
        fc = COLORS['green_fill'] if r['pass'] \
             else COLORS['anchor_f']
        val = 1 if r['pass'] else 0
        ax.bar(xi, val, width=0.50,
               facecolor=fc, edgecolor=ec,
               linewidth=1.2, zorder=3, clip_on=False)
        sym = 'PASS' if r['pass'] else 'FAIL'
        ax.text(xi, val + 0.03, sym,
                ha='center', va='bottom',
                fontsize=9, fontweight='bold',
                color=ec, clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{c}\n{results[c]['name'][:22]}"
         for c in checks],
        fontsize=8, rotation=0, ha='center')
    ax.tick_params(axis='x', pad=4,
                   length=6, width=1.2)
    ax.set_ylabel('Check result (1=PASS)',
                  fontsize=11, fontweight='bold')
    ax.set_ylim(-0.1, 1.6)
    ax.set_xlim(-0.6, len(checks) - 0.4)
    prism_axes(ax, ymax=1.2, yticks=[0, 1])

    n_pass = sum(1 for r in results.values() if r['pass'])
    ax.set_title(
        f'NSP9-NSP12 Structural Audit — {n_pass}/{len(checks)} PASS\n'
        f'Anchor: ARG733 (cons=1.000) | PDB: NSP9-NSP12_5',
        loc='left', fontsize=9.5, pad=8,
        linespacing=1.5)

    path = os.path.join(OUT_DIR,
        'Fig_Step01_StructuralAudit_NSP9-NSP12.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(results, n_pass, n_total):
    print('\n' + '='*60)
    print('STEP 1 STRUCTURAL VERIFICATION — NSP9-NSP12')
    print('='*60)
    for ck, r in results.items():
        sym = '✅' if r['pass'] else '❌'
        print(f"  {sym} {ck}: {r['name']}")
        print(f"       {r['detail']}")
    print(f"\n  RESULT: {n_pass}/{n_total} PASS")
    if n_pass == n_total:
        print("  → Structure APPROVED for CGCP pipeline")
    else:
        print("  → Structure needs review")
    print('='*60)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 1 — NSP9-NSP12 Structural Verification')
    results         = run_checks(PDB_FILE)
    n_pass, n_total = save_report(results)
    make_figure(results)
    print_summary(results, n_pass, n_total)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Step 2 — contact mapping")
