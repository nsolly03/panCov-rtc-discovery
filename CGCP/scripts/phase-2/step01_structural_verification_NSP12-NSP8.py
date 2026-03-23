#!/usr/bin/env python3
"""
CGCP Phase 2 Step 1 - Structural Verification: NSP12-NSP8
Verifies structure integrity for CGCP pipeline.
Anchor: LYS332 (NSP12, cons=1.000)

7 checks:
  C1: Chains A (NSP12) + B (NSP8) present
  C2: Chain sizes reasonable
  C3: LYS332 adjacent to NSP8 (< 15A)
  C4: Backbone completeness (missing residues)
  C5: LYS332 present and correct
  C6: Positive control — LYS332 to NSP8 interface
  C7: Negative control — LYS332 to RdRp catalytic triad (> 20A)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step01_structural_verification_NSP12-NSP8.py
"""

import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
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
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'axes.grid':          False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
})

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP12-NSP8_4/receptor_NSP12-NSP8_4.pdb')
OUT_DIR = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-01-structure')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'visuals'), exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'pymol-sessions'), exist_ok=True)

ANCHOR_RES  = 332
ANCHOR_AA   = 'LYS'
ANCHOR_CHAIN = 'A'

# RdRp catalytic triad for negative control
NEG_CTRL = [618, 759, 760]  # ASP618, SER759, ASP760


def get_ca(chain, resnum):
    try:
        return chain[(' ', resnum, ' ')]['CA'].coord
    except KeyError:
        return None


def dist(a, b):
    return float(np.linalg.norm(np.array(a) - np.array(b)))


def run_checks(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', pdb_file)[0]

    results = {}

    # C1: chains present
    chains = [c.id for c in structure.get_chains()]
    c1 = 'A' in chains and 'B' in chains
    results['C1'] = {
        'name':    'Chains A+B present',
        'pass':    c1,
        'detail':  f"Chains found: {chains}",
    }

    chain_a = structure['A']
    chain_b = structure['B']

    # C2: chain sizes
    na = sum(1 for _ in chain_a.get_residues()
             if PDB.is_aa(_))
    nb = sum(1 for _ in chain_b.get_residues()
             if PDB.is_aa(_))
    c2 = na > 500 and nb > 50
    results['C2'] = {
        'name':   'Chain sizes reasonable',
        'pass':   c2,
        'detail': f"NSP12={na} res, NSP8={nb} res",
    }

    # C3: LYS332 adjacent to NSP8
    ca_anchor = get_ca(chain_a, ANCHOR_RES)
    if ca_anchor is not None:
        min_dist = min(
            dist(ca_anchor, r['CA'].coord)
            for r in chain_b.get_residues()
            if PDB.is_aa(r) and 'CA' in r
        )
        c3 = min_dist < 15.0
        results['C3'] = {
            'name':   'Anchor adjacent to NSP8',
            'pass':   c3,
            'detail': f"LYS332 min dist to NSP8: {min_dist:.2f}A",
        }
    else:
        results['C3'] = {
            'name':   'Anchor adjacent to NSP8',
            'pass':   False,
            'detail': 'LYS332 CA not found',
        }

    # C4: backbone completeness
    missing = []
    res_list = [r for r in chain_a.get_residues()
                if PDB.is_aa(r)]
    for r in res_list:
        if 'CA' not in r or 'N' not in r or 'C' not in r:
            missing.append(r.id[1])
    c4 = len(missing) == 0
    results['C4'] = {
        'name':   'Backbone complete (Chain A)',
        'pass':   c4,
        'detail': f"{len(missing)} missing backbone atoms"
                  + (f": {missing[:5]}" if missing else ""),
    }

    # C5: LYS332 present and correct AA
    try:
        res332 = chain_a[(' ', ANCHOR_RES, ' ')]
        aa332  = res332.resname
        c5 = aa332 == ANCHOR_AA
        ca332 = res332['CA'].coord
        results['C5'] = {
            'name':   'LYS332 present and correct',
            'pass':   c5,
            'detail': f"Res 332 = {aa332} at "
                      f"({ca332[0]:.3f}, {ca332[1]:.3f}, {ca332[2]:.3f})",
        }
    except KeyError:
        results['C5'] = {
            'name':   'LYS332 present and correct',
            'pass':   False,
            'detail': 'Residue 332 not found',
        }

    # C6: positive control — LYS332 to NSP8 interface residues
    if ca_anchor is not None:
        nsp8_res = [r for r in chain_b.get_residues()
                    if PDB.is_aa(r) and 'CA' in r]
        close = [r for r in nsp8_res
                 if dist(ca_anchor, r['CA'].coord) < 15.0]
        c6 = len(close) >= 3
        results['C6'] = {
            'name':   'Positive ctrl: LYS332 near NSP8',
            'pass':   c6,
            'detail': f"{len(close)} NSP8 residues within 15A of LYS332",
        }
    else:
        results['C6'] = {
            'name':   'Positive ctrl: LYS332 near NSP8',
            'pass':   False,
            'detail': 'LYS332 not found',
        }

    # C7: negative control — LYS332 far from catalytic triad
    if ca_anchor is not None:
        triad_dists = []
        for pos in NEG_CTRL:
            ca_neg = get_ca(chain_a, pos)
            if ca_neg is not None:
                d = dist(ca_anchor, ca_neg)
                triad_dists.append((pos, d))
        c7 = all(d > 20.0 for _, d in triad_dists)
        detail = " | ".join(
            [f"res{p}={d:.1f}A" for p, d in triad_dists])
        results['C7'] = {
            'name':   'Neg ctrl: LYS332 far from catalytic triad',
            'pass':   c7,
            'detail': detail,
        }
    else:
        results['C7'] = {
            'name':   'Neg ctrl: LYS332 far from catalytic triad',
            'pass':   False,
            'detail': 'LYS332 not found',
        }

    return results, structure


def save_report(results):
    n_pass = sum(1 for r in results.values() if r['pass'])
    n_total = len(results)

    path = os.path.join(OUT_DIR, 'audit_report.md')
    with open(path, 'w') as f:
        f.write('# CGCP Step 1 Structural Audit — NSP12-NSP8\n\n')
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
            'interface':   'NSP12-NSP8',
            'pdb':         PDB_FILE,
            'anchor':      f'{ANCHOR_AA}{ANCHOR_RES}',
            'n_pass':      n_pass,
            'n_total':     n_total,
            'checks':      results,
        }, f, indent=2)
    print(f"  JSON:   {path_json}")
    return n_pass, n_total


def make_figure(results):
    checks = list(results.keys())
    colors = ['#1A7D2E' if results[c]['pass']
              else '#D7191C' for c in checks]
    fills  = ['#A8D5B0' if results[c]['pass']
              else '#FDAE61' for c in checks]

    fig, ax = plt.subplots(figsize=(9, 4.5))

    x = np.arange(len(checks))
    vals = [1 if results[c]['pass'] else 0 for c in checks]

    for xi, val, fc, ec in zip(x, vals, fills, colors):
        ax.bar(xi, val, width=0.5,
               facecolor=fc, edgecolor=ec,
               linewidth=0.9, zorder=3,
               clip_on=False)
        sym = 'PASS' if val else 'FAIL'
        ax.text(xi, val + 0.03, sym,
                ha='center', va='bottom',
                fontsize=8, color=ec,
                fontweight='bold', clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{c}\n{results[c]['name'][:20]}"
         for c in checks],
        fontsize=7.5, ha='center')
    ax.set_ylabel('Check result (1=PASS)', fontsize=9)
    ax.set_ylim(-0.1, 1.5)
    ax.set_xlim(-0.6, len(checks) - 0.4)
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.set_yticks([0, 1])

    n_pass = sum(vals)
    ax.set_title(
        f'NSP12-NSP8 Structural Audit — {n_pass}/{len(checks)} PASS\n'
        f'Anchor: LYS332 (cons=1.000) | PDB: NSP12-NSP8_4',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    path = os.path.join(OUT_DIR,
        'Fig_Step01_StructuralAudit_NSP12-NSP8.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(results, n_pass, n_total):
    print('\n' + '='*60)
    print('STEP 1 STRUCTURAL VERIFICATION — NSP12-NSP8')
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
    print('CGCP Phase 2 Step 1 — NSP12-NSP8 Structural Verification')
    results, structure = run_checks(PDB_FILE)
    n_pass, n_total    = save_report(results)
    make_figure(results)
    print_summary(results, n_pass, n_total)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Step 2 — contact mapping")
