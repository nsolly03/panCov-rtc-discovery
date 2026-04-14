#!/usr/bin/env python3
"""
Interface Contact Diagrams — All 9 RTC interfaces
Uses ORIGINAL PDB structures, 8 A Ca-Ca cutoff (= ~5 A heavy atom).
Residues coloured by chemical feature type.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/interface_contact_diagrams.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from Bio import PDB

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'figure.facecolor': 'white',
    'pdf.fonttype': 42,
    'ps.fonttype':  42,
})

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDBREF  = f'{BASE}/00-reference/pdb_structures'
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/interface_contacts')
os.makedirs(OUT_DIR, exist_ok=True)

CUTOFF = 8.0   # Å Ca-Ca  (= ~5 A heavy atom distance)

# ── Chemical colours ──────────────────────────────────────────
AROMATIC    = {'PHE','TYR','TRP','HIS'}
CHARGED_POS = {'ARG','LYS'}
CHARGED_NEG = {'ASP','GLU'}
HBOND       = {'SER','THR','ASN','GLN','CYS'}
HYDROPHOBIC = {'ALA','VAL','LEU','ILE','MET','PRO','GLY'}

COL = {
    'aromatic':    '#7B5EA7',  # purple
    'charged':     '#C08B2A',  # gold
    'hbond':       '#2E8B8B',  # teal
    'hydro_A':     '#3D6A9E',  # blue  — top chain
    'hydro_B':     '#8B2525',  # red   — bottom chain
    'other':       '#888888',
}

def get_color(resname, role):
    if resname in AROMATIC:                        return COL['aromatic']
    if resname in CHARGED_POS|CHARGED_NEG:        return COL['charged']
    if resname in HBOND:                           return COL['hbond']
    if resname in HYDROPHOBIC:
        return COL['hydro_A'] if role == 'A' else COL['hydro_B']
    return COL['other']

AA3 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}
def fmt(resname, resnum):
    return f"{AA3.get(resname, resname[:1])}{resnum}"

# ── Interface definitions ─────────────────────────────────────
INTERFACES = [
    dict(key='NSP12-NSP7',
         pdb=f'{PDBREF}/7C2K.pdb',   chainA='A', chainB='C',
         labelA='NSP12',             labelB='NSP7'),
    dict(key='NSP12-NSP8',
         pdb=f'{PDBREF}/7C2K.pdb',   chainA='A', chainB='B',
         labelA='NSP12',             labelB='NSP8'),
    dict(key='NSP9-NSP12',
         pdb=f'{PDBREF}/8SQK.pdb',   chainA='A', chainB='G',
         labelA='NSP12',             labelB='NSP9'),
    dict(key='NSP10-NSP16',
         pdb=f'{PDBREF}/6W4H.pdb',   chainA='A', chainB='B',
         labelA='NSP16',             labelB='NSP10'),
    dict(key='NSP7-NSP8',
         pdb=f'{BASE}/03-virtual-screening/NSP7-NSP8_6/receptor_NSP7-NSP8_ModeB_AF3_6.pdb',
         chainA='A', chainB='B',
         labelA='NSP8',              labelB='NSP7'),
    dict(key='NSP10-NSP14',
         pdb=f'{PDBREF}/7DIY.pdb',   chainA='A', chainB='B',
         labelA='NSP10',             labelB='NSP14'),
    dict(key='NSP13-Helicase',
         pdb=f'{PDBREF}/7NIO.pdb',   chainA='A', chainB='E',
         labelA='NSP13a',            labelB='NSP13b'),
    dict(key='NSP12-NSP13',
         pdb=f'{PDBREF}/7RDY.pdb',   chainA='A', chainB='E',
         labelA='NSP12',             labelB='NSP13'),
    dict(key='NSP15',
         pdb=f'{PDBREF}/NSP15/9HH5.pdb', chainA='A', chainB='B',
         labelA='NSP15a',            labelB='NSP15b'),
]

# ── Load contacts ─────────────────────────────────────────────
def load_interface(ifc):
    parser = PDB.PDBParser(QUIET=True)
    path, ca_id, cb_id = ifc['pdb'], ifc['chainA'], ifc['chainB']

    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        return [], [], []

    model = parser.get_structure('s', path)[0]

    def ca_map(cid):
        try:    ch = model[cid]
        except: return {}
        return {r.get_id()[1]: (r.resname, r['CA'].coord)
                for r in ch if PDB.is_aa(r) and 'CA' in r}

    cA = ca_map(ca_id)
    cB = ca_map(cb_id)

    pairs = [(ra, rb)
             for ra, (_, va) in cA.items()
             for rb, (_, vb) in cB.items()
             if np.linalg.norm(va - vb) <= CUTOFF]

    res_a = [(rn, cA[rn][0]) for rn in sorted({p[0] for p in pairs})]
    res_b = [(rn, cB[rn][0]) for rn in sorted({p[1] for p in pairs})]
    return res_a, res_b, pairs

# ── Draw one panel ────────────────────────────────────────────
def draw_panel(ifc, ax, show_legend=True):
    res_a, res_b, pairs = load_interface(ifc)
    lA, lB = ifc['labelA'], ifc['labelB']

    if not res_a:
        ax.set_facecolor('#F5F0EC')
        ax.text(0.5, 0.5, f"No data — {ifc['key']}",
                ha='center', va='center',
                transform=ax.transAxes, color='red', fontsize=9)
        ax.axis('off')
        return

    n_a, n_b = len(res_a), len(res_b)
    W = float(max(n_a, n_b) + 2)

    ax.set_xlim(-0.6, W - 0.4)
    ax.set_ylim(-0.55, 3.8)
    ax.set_facecolor('#F5F0EC')
    ax.axis('off')

    R = 0.30
    LABEL_PAD = -0.8   # x offset of chain label from first circle

    def draw_row(residues, y, role):
        n   = len(residues)
        x0  = (W - n) / 2.0
        pos = {}
        for i, (rn, rname) in enumerate(residues):
            x   = x0 + i
            col = get_color(rname, role)
            # white halo
            ax.add_patch(plt.Circle((x, y), R + 0.055,
                                    color='white', zorder=3, linewidth=0))
            # coloured circle
            ax.add_patch(plt.Circle((x, y), R,
                                    color=col, zorder=4, linewidth=0))
            lbl = fmt(rname, rn)
            fs  = 5.2 if len(lbl) >= 5 else 6.0 if len(lbl) == 4 else 7.0
            ax.text(x, y, lbl, ha='center', va='center',
                    fontsize=fs, fontweight='bold',
                    color='white', zorder=5)
            pos[rn] = (x, y)
        # chain label
        x_first = (W - n) / 2.0
        ax.text(x_first + LABEL_PAD, y, role == 'A' and lA or lB,
                va='center', ha='right',
                fontsize=9, fontweight='bold', color='#222222')
        return pos

    pos_a = draw_row(res_a, y=2.7, role='A')
    pos_b = draw_row(res_b, y=0.7, role='B')

    # Contact arcs — deduplicated by x positions
    drawn = set()
    for rn_a, rn_b in pairs:
        if rn_a not in pos_a or rn_b not in pos_b:
            continue
        x1, y1 = pos_a[rn_a]
        x2, y2 = pos_b[rn_b]
        key_arc = (round(x1, 2), round(x2, 2))
        if key_arc in drawn:
            continue
        drawn.add(key_arc)
        dx  = x2 - x1
        rad = np.clip(dx / W * 0.4 + 0.05, -0.35, 0.35)
        ax.annotate('',
            xy=(x2, y2 + R + 0.04),
            xytext=(x1, y1 - R - 0.04),
            arrowprops=dict(
                arrowstyle='-',
                color='#BBBBBB', lw=0.75, alpha=0.85,
                connectionstyle=f'arc3,rad={rad:.2f}'),
            zorder=2)

    # Title + stats
    ax.set_title(
        f"{ifc['key'].replace('-', '–')}  interface  ({CUTOFF:.0f} Å Ca–Ca)",
        fontsize=9.5, fontweight='bold', pad=5, color='#111111')
    ax.text((W - 1) / 2, 0.08,
            f"{n_a} {lA}  +  {n_b} {lB}  residues   |   {len(pairs)} contacts",
            ha='center', fontsize=7, color='#888888', style='italic')

    if show_legend:
        leg = [
            mpatches.Patch(color=COL['aromatic'],
                           label='Aromatic  (PHE · TYR · TRP · HIS)'),
            mpatches.Patch(color=COL['hydro_A'],
                           label=f'Hydrophobic  {lA}'),
            mpatches.Patch(color=COL['hydro_B'],
                           label=f'Hydrophobic  {lB}'),
            mpatches.Patch(color=COL['charged'],
                           label='Charged  (ARG · LYS · ASP · GLU)'),
            mpatches.Patch(color=COL['hbond'],
                           label='Polar / H-bond  (SER · THR · ASN · GLN · CYS)'),
        ]
        ax.legend(handles=leg, loc='lower center', ncol=3,
                  fontsize=7.5, frameon=False,
                  bbox_to_anchor=(0.5, -0.32))

# ── Individual figures ────────────────────────────────────────
print(f"Generating contact diagrams  (cutoff = {CUTOFF} A Ca-Ca)\n")
for ifc in INTERFACES:
    k = ifc['key']
    print(f"  {k} ...", end=' ', flush=True)
    fig, ax = plt.subplots(figsize=(18, 5.8))
    fig.patch.set_facecolor('white')
    draw_panel(ifc, ax, show_legend=True)
    fig.tight_layout(rect=[0, 0.12, 1, 1])
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT_DIR, f'contacts_{k}.{ext}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("saved")

# ── Combined 3x3 ──────────────────────────────────────────────
print("\nGenerating combined 3x3 figure ...", end=' ', flush=True)
fig = plt.figure(figsize=(30, 16))
fig.patch.set_facecolor('white')
gs  = GridSpec(3, 3, figure=fig,
               hspace=0.60, wspace=0.12,
               left=0.02, right=0.99,
               top=0.93,  bottom=0.08)
for i, ifc in enumerate(INTERFACES):
    ax = fig.add_subplot(gs[i // 3, i % 3])
    draw_panel(ifc, ax, show_legend=False)

fig.suptitle(
    'RTC Interface Contact Maps — All 9 Binary Interfaces\n'
    f'Ca–Ca distance ≤ {CUTOFF:.0f} Å  |  Coloured by chemical feature type',
    fontsize=14, fontweight='bold', y=0.99)

fig.legend(handles=[
    mpatches.Patch(color=COL['aromatic'], label='Aromatic (PHE/TYR/TRP/HIS)'),
    mpatches.Patch(color=COL['hydro_A'],  label='Hydrophobic — top chain'),
    mpatches.Patch(color=COL['hydro_B'],  label='Hydrophobic — bottom chain'),
    mpatches.Patch(color=COL['charged'],  label='Charged (ARG/LYS/ASP/GLU)'),
    mpatches.Patch(color=COL['hbond'],    label='Polar/H-bond (SER/THR/ASN/GLN/CYS)'),
], loc='lower center', ncol=5, fontsize=10.5,
   frameon=False, bbox_to_anchor=(0.5, 0.005))

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR, f'Fig_AllInterfaceContacts.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("saved")
print(f"\nAll figures saved to:\n  {OUT_DIR}")
