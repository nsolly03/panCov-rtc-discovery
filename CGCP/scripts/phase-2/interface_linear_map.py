#!/usr/bin/env python3
"""
Interface Linear Contact Maps — All 9 RTC interfaces
Two protein bars per interface, coloured ticks at contact residues,
arcs connecting contact pairs.  Style: panel A from Gao et al. 2020.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/interface_linear_map.py
"""

import os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from Bio import PDB

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

BASE   = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDBREF = f'{BASE}/00-reference/pdb_structures'
OUT    = os.path.join(BASE, 'CGCP/02-deep-dive/interface_contacts')
os.makedirs(OUT, exist_ok=True)

CUTOFF = 8.0   # Å Ca-Ca

# ── Protein lengths (literature / UniProt P0DTD1) ─────────────
LENGTHS = {
    'NSP7':   83,  'NSP8':  198, 'NSP9':  113,
    'NSP10': 139,  'NSP12': 932, 'NSP13': 601,
    'NSP14': 527,  'NSP15': 346, 'NSP16': 298,
}

# ── Domain colours per protein (for the bar fill) ─────────────
PROT_COLOR = {
    'NSP7':  '#B8A0D8',  'NSP8':  '#C8945A',  'NSP9':  '#98C8A0',
    'NSP10': '#E8C870',  'NSP12': '#E8A0A0',  'NSP13': '#7DA8C8',
    'NSP14': '#C0A0E0',  'NSP15': '#F0A080',  'NSP16': '#A0C8A0',
    'NSP13b\n(homodimer)': '#5D88A8',
    'NSP15b\n(homodimer)': '#D08060',
}

TICK_COLOR_A = '#1A5276'   # contact ticks chain A
TICK_COLOR_B = '#922B21'   # contact ticks chain B
ARC_COLOR    = '#AAAAAA'

# ── Interface registry ────────────────────────────────────────
INTERFACES = [
    dict(key='NSP12–NSP7',
         pdb=f'{PDBREF}/7C2K.pdb',
         chainA='A', nameA='NSP12',
         chainB='C', nameB='NSP7'),
    dict(key='NSP12–NSP8',
         pdb=f'{PDBREF}/7C2K.pdb',
         chainA='A', nameA='NSP12',
         chainB='B', nameB='NSP8'),
    dict(key='NSP9–NSP12',
         pdb=f'{PDBREF}/8SQK.pdb',
         chainA='G', nameA='NSP9',
         chainB='A', nameB='NSP12'),
    dict(key='NSP10–NSP16',
         pdb=f'{PDBREF}/6W4H.pdb',
         chainA='B', nameA='NSP10',
         chainB='A', nameB='NSP16'),
    dict(key='NSP7–NSP8 (AF3)',
         pdb=f'{BASE}/03-virtual-screening/NSP7-NSP8_6/receptor_NSP7-NSP8_ModeB_AF3_6.pdb',
         chainA='B', nameA='NSP7',
         chainB='A', nameB='NSP8'),
    dict(key='NSP10–NSP14',
         pdb=f'{PDBREF}/7DIY.pdb',
         chainA='A', nameA='NSP10',
         chainB='B', nameB='NSP14'),
    dict(key='NSP13–Helicase',
         pdb=f'{PDBREF}/7NIO.pdb',
         chainA='A', nameA='NSP13a',
         chainB='E', nameB='NSP13b'),
    dict(key='NSP12–NSP13',
         pdb=f'{PDBREF}/7RDY.pdb',
         chainA='A', nameA='NSP12',
         chainB='E', nameB='NSP13'),
    dict(key='NSP15 homodimer',
         pdb=f'{PDBREF}/NSP15/9HH5.pdb',
         chainA='A', nameA='NSP15a',
         chainB='B', nameB='NSP15b'),
]

parser = PDB.PDBParser(QUIET=True)

def get_contacts(pdb, cA, cB, cutoff=8.0):
    if not os.path.exists(pdb):
        print(f"  MISSING {pdb}"); return {}, {}, []
    model = parser.get_structure('s', pdb)[0]
    def ca(cid):
        try:    ch = model[cid]
        except: return {}
        return {r.get_id()[1]: (r.resname, r['CA'].coord)
                for r in ch if PDB.is_aa(r) and 'CA' in r}
    mA, mB = ca(cA), ca(cB)
    pairs = [(ra, rb) for ra,(_,va) in mA.items()
                       for rb,(_,vb) in mB.items()
                       if np.linalg.norm(va-vb) <= cutoff]
    resA = {p[0] for p in pairs}
    resB = {p[1] for p in pairs}
    # actual resnum range for scaling
    all_a = sorted(mA.keys()); all_b = sorted(mB.keys())
    return (mA, mB, pairs,
            min(all_a), max(all_a),
            min(all_b), max(all_b))

# ── Draw one interface panel ──────────────────────────────────
def draw_panel(ifc, ax):
    key  = ifc['key']
    nA, nB = ifc['nameA'], ifc['nameB']
    # strip homodimer tag for lookup
    nA_key = nA.replace('a','').replace('b','').replace('\n(homodimer)','')
    nB_key = nB.replace('a','').replace('b','').replace('\n(homodimer)','')

    result = get_contacts(ifc['pdb'], ifc['chainA'], ifc['chainB'], CUTOFF)
    if len(result) == 3:
        ax.axis('off'); return

    mA, mB, pairs, minA, maxA, minB, maxB = result

    resA_contact = sorted({p[0] for p in pairs})
    resB_contact = sorted({p[1] for p in pairs})

    # Full protein lengths from literature
    LA = LENGTHS.get(nA_key, maxA - minA + 1)
    LB = LENGTHS.get(nB_key, maxB - minB + 1)

    # For proteins where PDB starts at residue >1 (e.g. NSP12 from res 31)
    # we keep offset for correct positioning
    offA = minA - 1   # first residue in PDB structure
    offB = minB - 1

    # Scale everything to [0, 1] using full protein length
    def norm_a(r): return (r - 1) / LA
    def norm_b(r): return (r - 1) / LB

    cA = PROT_COLOR.get(nA_key, '#AAAAAA')
    cB = PROT_COLOR.get(nB_key, '#AAAAAA')

    BAR_H  = 0.22
    Y_A    = 0.72
    Y_B    = 0.28
    Y_MID  = (Y_A + Y_B) / 2

    ax.set_xlim(-0.05, 1.08)
    ax.set_ylim(0.0, 1.05)
    ax.set_facecolor('white')
    ax.axis('off')

    # ── Protein bars ─────────────────────────────────────────
    for (y, col, lname, L, off) in [
            (Y_A, cA, nA, LA, offA),
            (Y_B, cB, nB, LB, offB)]:

        # Main bar — from first residue in structure to end
        x_start = (off) / L
        x_end   = 1.0
        w       = x_end - x_start

        # Grey placeholder for unresolved region if any
        if off > 0:
            ax.add_patch(plt.Rectangle(
                (0, y - BAR_H/2), x_start, BAR_H,
                color='#E8E8E8', zorder=2, linewidth=0))

        ax.add_patch(plt.Rectangle(
            (x_start, y - BAR_H/2), w, BAR_H,
            color=col, zorder=2, linewidth=0, alpha=0.85))

        # border
        ax.add_patch(plt.Rectangle(
            (0, y - BAR_H/2), 1.0, BAR_H,
            fill=False, edgecolor='#CCCCCC', linewidth=0.8, zorder=3))

        # protein label left
        ax.text(-0.02, y, f"{lname}  (1–{L})",
                ha='right', va='center',
                fontsize=9, fontweight='bold', color='#222222')

        # position ticks
        step = 50 if L <= 250 else 100 if L <= 600 else 200
        for t in range(0, L+1, step):
            xp = t / L
            ax.plot([xp, xp], [y - BAR_H/2 - 0.025,
                                y - BAR_H/2],
                    color='#AAAAAA', lw=0.7, zorder=4)
            ax.text(xp, y - BAR_H/2 - 0.04, str(t),
                    ha='center', va='top', fontsize=5.5,
                    color='#AAAAAA')

    # ── Contact ticks ─────────────────────────────────────────
    TICK_H = 0.16

    for rn in resA_contact:
        xp = norm_a(rn)
        ax.add_patch(plt.Rectangle(
            (xp - 0.003, Y_A - BAR_H/2 - TICK_H),
            0.006, TICK_H,
            color=TICK_COLOR_A, zorder=5, linewidth=0))

    for rn in resB_contact:
        xp = norm_b(rn)
        ax.add_patch(plt.Rectangle(
            (xp - 0.003, Y_B + BAR_H/2),
            0.006, TICK_H,
            color=TICK_COLOR_B, zorder=5, linewidth=0))

    # ── Contact arcs ──────────────────────────────────────────
    drawn = set()
    for rn_a, rn_b in pairs:
        xa = norm_a(rn_a)
        xb = norm_b(rn_b)
        k  = (round(xa, 3), round(xb, 3))
        if k in drawn: continue
        drawn.add(k)
        ya_tip = Y_A - BAR_H/2 - TICK_H
        yb_tip = Y_B + BAR_H/2 + TICK_H
        dx = xb - xa
        rad = np.clip(dx * 0.5, -0.45, 0.45)
        ax.annotate('',
            xy=(xb, yb_tip),
            xytext=(xa, ya_tip),
            arrowprops=dict(
                arrowstyle='-',
                color=ARC_COLOR, lw=0.7, alpha=0.7,
                connectionstyle=f'arc3,rad={rad:.2f}'),
            zorder=1)

    # ── Title + stats ─────────────────────────────────────────
    ax.set_title(
        f"{key}  interface  ({CUTOFF:.0f} Å Cα–Cα)\n"
        f"{len(resA_contact)} {nA} residues  +  "
        f"{len(resB_contact)} {nB} residues  |  "
        f"{len(pairs)} contact pairs",
        fontsize=9.5, fontweight='bold',
        color='#111111', pad=6)

    # ── Legend patches ────────────────────────────────────────
    leg = [
        mpatches.Patch(color=TICK_COLOR_A,
                       label=f'{nA} contact residues'),
        mpatches.Patch(color=TICK_COLOR_B,
                       label=f'{nB} contact residues'),
    ]
    ax.legend(handles=leg, loc='upper right',
              fontsize=7.5, frameon=False)


# ── Individual figures ────────────────────────────────────────
print(f"Drawing interface linear maps  (cutoff = {CUTOFF} Å)\n")
for ifc in INTERFACES:
    k = ifc['key']
    print(f"  {k} ...", end=' ', flush=True)
    fig, ax = plt.subplots(figsize=(14, 3.2))
    fig.patch.set_facecolor('white')
    draw_panel(ifc, ax)
    plt.tight_layout(pad=1.5)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'linear_{k.replace(" ","_")}.{ext}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("saved")

# ── Combined 3x3 ──────────────────────────────────────────────
print("\nGenerating combined figure ...", end=' ', flush=True)
fig, axes = plt.subplots(3, 3, figsize=(24, 13))
fig.patch.set_facecolor('white')
fig.suptitle(
    'Interface Contact Maps — All 9 RTC Binary Interfaces\n'
    f'Blue ticks = top chain contacts  |  Red ticks = bottom chain contacts  |  '
    f'Cα–Cα ≤ {CUTOFF:.0f} Å',
    fontsize=13, fontweight='bold', y=1.01)
fig.subplots_adjust(hspace=0.75, wspace=0.35,
                    left=0.09, right=0.99,
                    top=0.95, bottom=0.04)

for ax, ifc in zip(axes.flat, INTERFACES):
    draw_panel(ifc, ax)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT, f'Fig_LinearContactMaps.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("saved")
print(f"\nAll outputs → {OUT}")
