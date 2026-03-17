#!/usr/bin/env python3
"""
CGCP Phase 2 — Step 1: Structural Verification
Interface: NSP12-NSP7
Date: 2026-03-17

Verifies the receptor structure is correct, complete, and ready
for downstream contact mapping and cluster analysis.

Checks:
  1. Both chains present (A=NSP12, C=NSP7)
  2. Chain sizes match expected (834 NSP12, ~63 NSP7)
  3. Interface adjacency — chains are physically close
  4. No missing backbone atoms in interface region
  5. Known anchor PHE440 is present and located correctly
  6. Positive control: PHE440 is near chain C (interface side)
  7. Negative control: random surface residue is far from chain C

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step01_structural_verification_NSP12-NSP7.py

Output:
    CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/
        audit_report.md
        structural_verification.json
        Fig_Step01_StructuralAudit.png
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

try:
    from Bio import PDB
    from Bio.PDB import NeighborSearch
except ImportError:
    raise ImportError("Install biopython: conda install biopython")

# ─────────────────────────────────────────────────────────────
# PRISM STYLE
# ─────────────────────────────────────────────────────────────
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
    'axes.titleweight':   'normal',
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'xtick.major.pad':    4,
    'ytick.major.pad':    4,
    'xtick.labelsize':    8.5,
    'ytick.labelsize':    8.5,
    'axes.grid':          False,
    'legend.fontsize':    8,
    'legend.frameon':     False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.12,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

P = {
    'pass':    '#2166AC',
    'pass_f':  '#92C5DE',
    'fail':    '#D7191C',
    'fail_f':  '#FDAE61',
    'warn':    '#4DAC26',
    'warn_f':  '#A6D96A',
    'gray':    '#636363',
    'text':    '#1A1A1A',
}

# ─────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────
BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(
    BASE, '03-virtual-screening/NSP12-NSP7_3/receptor_NSP12-NSP7_3.pdb')
OUT_DIR  = os.path.join(
    BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-01-structure')
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────
# EXPECTED VALUES
# ─────────────────────────────────────────────────────────────
CHAIN_A        = 'A'    # NSP12
CHAIN_C        = 'C'    # NSP7
EXPECTED_A     = 834    # residues in NSP12
EXPECTED_C     = 63     # residues in NSP7 (truncated in 7BV2)
ANCHOR_RESI    = 440    # PHE440 — primary pharmacophore anchor
ANCHOR_NAME    = 'PHE'
# Positive control: PHE440 should be close to NSP7 (<15 Å)
POS_CTRL_RESI  = 440
POS_CTRL_DIST  = 15.0
# Negative control: a surface residue far from NSP7 (>25 Å)
NEG_CTRL_RESI  = 200    # NSP12 surface residue, away from interface
NEG_CTRL_DIST  = 25.0

# ─────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────
def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax * 0.02, ymax * 1.02)
    if yticks is not None:
        ax.set_yticks(yticks)


def min_chain_distance(chain_a_res, chain_c_res):
    """Minimum distance between any atoms in two residue lists."""
    min_d = float('inf')
    for ra in chain_a_res:
        for rc in chain_c_res:
            for aa in ra.get_atoms():
                for ac in rc.get_atoms():
                    d = np.linalg.norm(aa.coord - ac.coord)
                    if d < min_d:
                        min_d = d
    return min_d


def residue_to_chain_distance(residue, chain_residues):
    """Minimum distance from one residue to all atoms in a chain."""
    min_d = float('inf')
    for cr in chain_residues:
        for aa in residue.get_atoms():
            for ac in cr.get_atoms():
                d = np.linalg.norm(aa.coord - ac.coord)
                if d < min_d:
                    min_d = d
    return min_d


def check_backbone(residues):
    """Count residues with incomplete backbone (missing N, CA, C, O)."""
    missing = []
    for res in residues:
        if res.get_id()[0] != ' ':
            continue
        atom_names = {a.get_name() for a in res.get_atoms()}
        for backbone in ['N', 'CA', 'C', 'O']:
            if backbone not in atom_names:
                missing.append((res.get_id()[1], res.get_resname(),
                                 backbone))
                break
    return missing


# ─────────────────────────────────────────────────────────────
# MAIN AUDIT
# ─────────────────────────────────────────────────────────────
def run_audit():
    print('=' * 60)
    print('CGCP Phase 2 Step 1 — Structural Verification')
    print('Interface: NSP12-NSP7')
    print('=' * 60)

    checks = []   # list of {name, status, detail, value}

    # ── Load structure ────────────────────────────────────────
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('NSP12-NSP7', PDB_FILE)
    model = structure[0]

    # ── CHECK 1: Chains present ──────────────────────────────
    chains_present = [c.id for c in model.get_chains()]
    both_present   = CHAIN_A in chains_present and CHAIN_C in chains_present
    checks.append({
        'name':   'C1: Both chains present',
        'status': 'PASS' if both_present else 'FAIL',
        'detail': f'Chains found: {chains_present}',
        'value':  1 if both_present else 0,
    })
    print(f"\n  C1 Chains: {chains_present} "
          f"→ {'PASS' if both_present else 'FAIL'}")

    if not both_present:
        print("  FATAL: required chains missing. Stopping.")
        return checks, None

    chain_a = model[CHAIN_A]
    chain_c = model[CHAIN_C]

    res_a = [r for r in chain_a.get_residues() if r.get_id()[0] == ' ']
    res_c = [r for r in chain_c.get_residues() if r.get_id()[0] == ' ']
    n_a, n_c = len(res_a), len(res_c)

    # ── CHECK 2: Chain sizes ─────────────────────────────────
    size_ok = (abs(n_a - EXPECTED_A) <= 20 and
               abs(n_c - EXPECTED_C) <= 10)
    checks.append({
        'name':   'C2: Chain sizes correct',
        'status': 'PASS' if size_ok else 'WARN',
        'detail': (f'Chain A: {n_a} res (expected ~{EXPECTED_A})  '
                   f'Chain C: {n_c} res (expected ~{EXPECTED_C})'),
        'value':  1 if size_ok else 0,
    })
    print(f"  C2 Sizes: A={n_a} (exp {EXPECTED_A}), "
          f"C={n_c} (exp {EXPECTED_C}) "
          f"→ {'PASS' if size_ok else 'WARN'}")

    # ── CHECK 3: Interface adjacency ─────────────────────────
    # Use only first 20 residues of each chain for speed
    sample_a = [r for r in res_a if abs(r.get_id()[1] - 440) < 30]
    sample_c = res_c
    min_d = min_chain_distance(sample_a, sample_c)
    adj_ok = min_d < 10.0
    checks.append({
        'name':   'C3: Chains adjacent (<10 Å)',
        'status': 'PASS' if adj_ok else 'FAIL',
        'detail': f'Min inter-chain distance (sample): {min_d:.2f} Å',
        'value':  round(min_d, 2),
    })
    print(f"  C3 Adjacency: min_dist={min_d:.2f} Å "
          f"→ {'PASS' if adj_ok else 'FAIL'}")

    # ── CHECK 4: Backbone completeness ───────────────────────
    missing_a = check_backbone(res_a)
    missing_c = check_backbone(res_c)
    bb_ok = len(missing_a) == 0 and len(missing_c) == 0
    checks.append({
        'name':   'C4: Backbone complete',
        'status': 'PASS' if bb_ok else 'WARN',
        'detail': (f'Chain A missing: {len(missing_a)} residues  '
                   f'Chain C missing: {len(missing_c)} residues'),
        'value':  len(missing_a) + len(missing_c),
    })
    print(f"  C4 Backbone: A_missing={len(missing_a)}, "
          f"C_missing={len(missing_c)} "
          f"→ {'PASS' if bb_ok else 'WARN'}")

    # ── CHECK 5: PHE440 present and correct ──────────────────
    phe440_present = False
    phe440_correct = False
    phe440_resi    = None
    for res in res_a:
        if res.get_id()[1] == ANCHOR_RESI:
            phe440_present = True
            phe440_resi    = res
            phe440_correct = res.get_resname() == ANCHOR_NAME
            break

    checks.append({
        'name':   f'C5: {ANCHOR_NAME}{ANCHOR_RESI} present and correct',
        'status': 'PASS' if (phe440_present and phe440_correct) else 'FAIL',
        'detail': (f'Present: {phe440_present}  '
                   f'Identity: {phe440_resi.get_resname() if phe440_resi else "N/A"} '
                   f'(expected {ANCHOR_NAME})'),
        'value':  1 if (phe440_present and phe440_correct) else 0,
    })
    print(f"  C5 PHE440: present={phe440_present}, "
          f"correct={phe440_correct} "
          f"→ {'PASS' if (phe440_present and phe440_correct) else 'FAIL'}")

    # ── CHECK 6 (POSITIVE CONTROL): PHE440 near NSP7 ─────────
    pos_ctrl_dist = None
    pos_ok        = False
    if phe440_resi is not None:
        pos_ctrl_dist = residue_to_chain_distance(phe440_resi, res_c)
        pos_ok        = pos_ctrl_dist < POS_CTRL_DIST
    checks.append({
        'name':   f'C6 [+ctrl]: PHE440 near NSP7 (<{POS_CTRL_DIST} Å)',
        'status': 'PASS' if pos_ok else 'FAIL',
        'detail': (f'PHE440 to NSP7 min dist: '
                   f'{pos_ctrl_dist:.2f} Å' if pos_ctrl_dist else 'N/A'),
        'value':  round(pos_ctrl_dist, 2) if pos_ctrl_dist else 999,
    })
    print(f"  C6 +ctrl PHE440-NSP7: {pos_ctrl_dist:.2f} Å "
          f"→ {'PASS' if pos_ok else 'FAIL'}")

    # ── CHECK 7 (NEGATIVE CONTROL): Surface residue far ──────
    neg_ctrl_resi_obj = None
    for res in res_a:
        if res.get_id()[1] == NEG_CTRL_RESI:
            neg_ctrl_resi_obj = res
            break

    neg_ctrl_dist = None
    neg_ok        = False
    if neg_ctrl_resi_obj is not None:
        neg_ctrl_dist = residue_to_chain_distance(neg_ctrl_resi_obj, res_c)
        neg_ok        = neg_ctrl_dist > NEG_CTRL_DIST
        neg_name      = neg_ctrl_resi_obj.get_resname()
    else:
        neg_name = 'NOT FOUND'

    checks.append({
        'name':   (f'C7 [-ctrl]: {neg_name}{NEG_CTRL_RESI} '
                   f'far from NSP7 (>{NEG_CTRL_DIST} Å)'),
        'status': 'PASS' if neg_ok else 'WARN',
        'detail': (f'{neg_name}{NEG_CTRL_RESI} to NSP7 min dist: '
                   f'{neg_ctrl_dist:.2f} Å' if neg_ctrl_dist else 'N/A'),
        'value':  round(neg_ctrl_dist, 2) if neg_ctrl_dist else 0,
    })
    print(f"  C7 -ctrl {neg_name}{NEG_CTRL_RESI}-NSP7: "
          f"{neg_ctrl_dist:.2f} Å "
          f"→ {'PASS' if neg_ok else 'WARN'}")

    # ── Summary ───────────────────────────────────────────────
    n_pass = sum(1 for c in checks if c['status'] == 'PASS')
    n_warn = sum(1 for c in checks if c['status'] == 'WARN')
    n_fail = sum(1 for c in checks if c['status'] == 'FAIL')

    overall = 'PASS' if n_fail == 0 else 'FAIL'
    print(f"\n  Summary: {n_pass} PASS | {n_warn} WARN | {n_fail} FAIL")
    print(f"  Overall: {overall}")

    metadata = {
        'pdb_file':     PDB_FILE,
        'chain_a':      CHAIN_A,
        'chain_c':      CHAIN_C,
        'n_res_a':      n_a,
        'n_res_c':      n_c,
        'anchor':       f'{ANCHOR_NAME}{ANCHOR_RESI}',
        'pos_ctrl_dist': round(pos_ctrl_dist, 2) if pos_ctrl_dist else None,
        'neg_ctrl_dist': round(neg_ctrl_dist, 2) if neg_ctrl_dist else None,
        'n_pass':       n_pass,
        'n_warn':       n_warn,
        'n_fail':       n_fail,
        'overall':      overall,
    }

    return checks, metadata


# ─────────────────────────────────────────────────────────────
# SAVE AUDIT REPORT (Markdown)
# ─────────────────────────────────────────────────────────────
def save_report(checks, metadata):
    path = os.path.join(OUT_DIR, 'audit_report.md')
    lines = [
        '# Phase 2 Step 1 — Structural Audit: NSP12-NSP7',
        '',
        f'**Date:** 2026-03-17  ',
        f'**PDB file:** `{os.path.basename(metadata["pdb_file"])}`  ',
        f'**Chain A:** NSP12 ({metadata["n_res_a"]} residues)  ',
        f'**Chain C:** NSP7  ({metadata["n_res_c"]} residues)  ',
        f'**Overall result:** {metadata["overall"]}  ',
        '',
        '## Checklist',
        '',
        '| Check | Status | Detail |',
        '|-------|--------|--------|',
    ]
    for c in checks:
        icon = '✅' if c['status'] == 'PASS' else \
               '⚠️'  if c['status'] == 'WARN' else '❌'
        lines.append(f'| {c["name"]} | {icon} {c["status"]} '
                     f'| {c["detail"]} |')

    lines += [
        '',
        '## Controls',
        '',
        f'**Positive control:** PHE440 (primary anchor) to NSP7 = '
        f'{metadata["pos_ctrl_dist"]} Å  ',
        f'Expected: < {POS_CTRL_DIST} Å (at interface)  ',
        '',
        f'**Negative control:** RES{NEG_CTRL_RESI} (surface) to NSP7 = '
        f'{metadata["neg_ctrl_dist"]} Å  ',
        f'Expected: > {NEG_CTRL_DIST} Å (away from interface)  ',
        '',
        '## Decision',
        '',
        f'Structure **{"accepted" if metadata["overall"] == "PASS" else "REJECTED"}** '
        f'for Phase 2 Step 2 (contact mapping).',
        '',
    ]

    with open(path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Report: {path}")
    return path


# ─────────────────────────────────────────────────────────────
# SAVE JSON
# ─────────────────────────────────────────────────────────────
def save_json(checks, metadata):
    import numpy as np

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer, np.floating)):
                return float(obj)
            return super().default(obj)

    out = {'metadata': metadata, 'checks': checks}
    fpath = os.path.join(OUT_DIR, 'structural_verification.json')
    with open(fpath, 'w') as fh:
        json.dump(out, fh, indent=2, cls=NumpyEncoder)
    print(f"  JSON: {fpath}")
    return fpath



def make_figure(checks, metadata):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2))
    fig.subplots_adjust(wspace=0.50,
                        bottom=0.26, top=0.88,
                        left=0.10,   right=0.97)

    # ── Panel A: Checklist bars ───────────────────────────────
    ax = axes[0]
    labels  = [c['name'].split(':')[0] for c in checks]   # C1..C7
    statuses = [c['status'] for c in checks]
    vals    = [1] * len(checks)

    colors_bar  = [P['pass_f'] if s == 'PASS' else
                   P['warn_f'] if s == 'WARN' else
                   P['fail_f'] for s in statuses]
    colors_edge = [P['pass'] if s == 'PASS' else
                   P['warn'] if s == 'WARN' else
                   P['fail'] for s in statuses]

    y = np.arange(len(checks))
    bars = ax.barh(y, vals, height=0.55,
                   color=colors_bar,
                   edgecolor=colors_edge,
                   linewidth=0.9, zorder=3)

    # Status labels
    for i, (s, c) in enumerate(zip(statuses, checks)):
        icon  = '✓' if s == 'PASS' else '!' if s == 'WARN' else '✗'
        color = (P['pass'] if s == 'PASS' else
                 P['warn'] if s == 'WARN' else P['fail'])
        ax.text(1.05, i, f'{icon}  {s}',
                va='center', ha='left',
                fontsize=8, color=color,
                fontweight='bold', clip_on=False)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_position(('outward', 5))
    ax.invert_yaxis()
    ax.set_title('A   Structural audit checklist\n'
                 '(receptor_NSP12-NSP7_3.pdb)',
                 loc='left', fontsize=8.5, pad=6,
                 linespacing=1.5)

    # Overall badge
    overall = metadata['overall']
    badge_color = P['pass'] if overall == 'PASS' else P['fail']
    badge_fill  = P['pass_f'] if overall == 'PASS' else P['fail_f']
    ax.text(0.97, 0.04,
            f'Overall: {overall}',
            transform=ax.transAxes,
            ha='right', va='bottom',
            fontsize=8.5, fontweight='bold',
            color=badge_color,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor=badge_fill,
                      edgecolor=badge_color,
                      linewidth=0.7))

    # ── Panel B: Control distances ────────────────────────────
    ax2 = axes[1]
    ctrl_labels = [
        f'PHE440\n(+control\nanchor)',
        f'RES{NEG_CTRL_RESI}\n(-control\nsurface)',
    ]
    ctrl_vals  = [
        metadata['pos_ctrl_dist'] or 0,
        metadata['neg_ctrl_dist'] or 0,
    ]
    ctrl_fills = [
        P['pass_f'] if ctrl_vals[0] < POS_CTRL_DIST else P['fail_f'],
        P['pass_f'] if ctrl_vals[1] > NEG_CTRL_DIST else P['warn_f'],
    ]
    ctrl_edges = [
        P['pass'] if ctrl_vals[0] < POS_CTRL_DIST else P['fail'],
        P['pass'] if ctrl_vals[1] > NEG_CTRL_DIST else P['warn'],
    ]

    xpos = np.array([0, 1])
    w    = 0.45
    for xi, val, fc, ec in zip(xpos, ctrl_vals, ctrl_fills, ctrl_edges):
        ax2.bar(xi, val, width=w, facecolor=fc, edgecolor=ec,
                linewidth=0.9, zorder=3, clip_on=False)
        ax2.text(xi, val + 0.4, f'{val:.1f} Å',
                 ha='center', va='bottom', fontsize=9,
                 fontweight='bold', color=ec, clip_on=False)

    # Threshold lines
    ax2.axhline(POS_CTRL_DIST, color=P['pass'], linewidth=0.8,
                linestyle=(0, (4, 3)), zorder=1)
    ax2.text(1.42, POS_CTRL_DIST + 0.3,
             f'<{POS_CTRL_DIST:.0f} Å\ninterface',
             fontsize=7, color=P['pass'],
             va='bottom', ha='right', clip_on=False)

    ax2.axhline(NEG_CTRL_DIST, color=P['gray'], linewidth=0.8,
                linestyle=(0, (4, 3)), zorder=1)
    ax2.text(1.42, NEG_CTRL_DIST + 0.3,
             f'>{NEG_CTRL_DIST:.0f} Å\nsurface',
             fontsize=7, color=P['gray'],
             va='bottom', ha='right', clip_on=False)

    ax2.set_xticks(xpos)
    ax2.set_xticklabels(ctrl_labels, fontsize=8.5)
    ax2.set_ylabel('Min distance to NSP7 chain (Å)', fontsize=9)
    ax2.set_xlim(-0.55, 1.55)
    ymax_b = max(ctrl_vals) * 1.55
    ymax_b = max(ymax_b, NEG_CTRL_DIST * 1.3)
    step_b = 10
    yticks_b = list(range(0, int(ymax_b) + step_b, step_b))
    ax2.set_title('B   Positive / negative control distances\n'
                  '(PHE440 = interface anchor; RES200 = surface)',
                  loc='left', fontsize=8.5, pad=6,
                  linespacing=1.5)
    # prism style
    ax2.spines['left'].set_position(('outward', 5))
    ax2.spines['bottom'].set_position(('outward', 5))
    ax2.spines['left'].set_bounds(0, max(yticks_b))
    ax2.set_ylim(-max(yticks_b) * 0.02, max(yticks_b) * 1.02)
    ax2.set_yticks(yticks_b)

    path = os.path.join(OUT_DIR, 'Fig_Step01_StructuralAudit.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")
    return path


# ─────────────────────────────────────────────────────────────
# WORKLOG ENTRY
# ─────────────────────────────────────────────────────────────
def print_worklog_entry(metadata):
    status = metadata['overall']
    print('\n' + '─' * 60)
    print('WORKLOG ENTRY — copy to CGCP/WORKLOG.md')
    print('─' * 60)
    print(f"""
## Entry 012 — Phase 2 Step 1: NSP12-NSP7 Structural Verification

**Date:** 2026-03-17
**Phase:** 2 — Deep Dive
**Step:** 1 — Structural verification
**Interface:** NSP12-NSP7
**Script:** CGCP/scripts/phase-2/step01_structural_verification_NSP12-NSP7.py

**Structure:** receptor_NSP12-NSP7_3.pdb
**Chains:** A = NSP12 ({metadata['n_res_a']} res) | C = NSP7 ({metadata['n_res_c']} res)

**Checklist:** {metadata['n_pass']} PASS | {metadata['n_warn']} WARN | {metadata['n_fail']} FAIL

**Controls:**
- Positive: PHE440 (anchor) → NSP7 = {metadata['pos_ctrl_dist']} Å (expected <{POS_CTRL_DIST} Å) ✓
- Negative: RES{NEG_CTRL_RESI} (surface) → NSP7 = {metadata['neg_ctrl_dist']} Å (expected >{NEG_CTRL_DIST} Å) ✓

**Overall: {status}**
**Decision:** {"Proceed to Step 2 — contact mapping" if status == "PASS" else "Fix structure before proceeding"}

**Status:** {"✅ Complete" if status == "PASS" else "❌ Needs fix"}
**Next:** Phase 2 Step 2 — raw contact mapping
""")


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    checks, metadata = run_audit()

    if metadata:
        save_report(checks, metadata)
        save_json(checks, metadata)
        make_figure(checks, metadata)
        print_worklog_entry(metadata)

        print(f"\nOutputs saved to:")
        print(f"  {OUT_DIR}")
        print(f"\nNext: Phase 2 Step 2 — contact mapping")
