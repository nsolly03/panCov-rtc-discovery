"""
Script 11_3: 3D Visualization Notebook — NSP12-NSP7 (nglview)
=============================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

8 interactive views using nglview:
  1. Full complex — all conserved hotspots
  2. Salt bridge zoomed — GLU431–LYS1
  3. Hotspots colored by composite score
  4. Hotspots colored by BSA burial depth
  5. Docking box visualization
  6. Structural overlay — 7BV2 + 6NUR + 7C2K + AF3
  7. Full hotspot ranking table
  8. Hydrophobic core — PHE440/PHE442/PHE415/TYR420 ★

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/11_3D_visualization_NSP12-NSP7_3.py
  jupyter notebook notebooks/NSP12-NSP7_3D_3.ipynb
"""

import json
from pathlib import Path

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP7"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP7"
NB_DIR   = PROJECT / "notebooks"
NB_DIR.mkdir(exist_ok=True)

HOTSPOTS_NSP12  = [440,412,442,443,420,843,409,
                    40,33,41,37,413,415,14,23]
HOTSPOTS_NSP7   = [40,14,33,41,37,11,23,5,15,29,12,4,1]
PRIMARY_NSP12   = {431}
PRIMARY_NSP7    = {1}
HYDROPHOBIC_CORE= {412,413,415,420,440,442,843}
NSP12_OFFSET    = 0
NSP7_OFFSET     = 1

CELLS = []

def code(src):
    CELLS.append({"cell_type":"code",
                  "execution_count":None,
                  "metadata":{},
                  "outputs":[],
                  "source": src.strip()})

def md(src):
    CELLS.append({"cell_type":"markdown",
                  "metadata":{},
                  "source": src.strip()})


def build():
    # ── Header ─────────────────────────────────────────────
    md("""# NSP12–NSP7 Interface: 3D Structural Visualization
**Project:** panCov-rtc-discovery
**Complex:** NSP12–NSP7 (RdRp core — thumb subdomain)
**Primary PDB:** 7BV2_NSP12-NSP7 (SARS-CoV-2, 2.90 Å)
**AF3 model:** iptm=0.81, ranking=0.84

## Key findings:
- **PHE440(NSP12) primary pharmacophore** — hydrophobic anchor, cons=1.000
- **Hydrophobic core:** PHE412/413/415/420/440/442/PHE843 — all conserved
- **Salt bridge:** GLU431(NSP12)–LYS1(NSP7) — SARS-CoV-1/2 only
- **fpocket druggability = 0.961** (6NUR) — highest in project
- **AF3 validation F1=0.951** — best in project

## Views:
1. Full complex  2. Salt bridge  3. Composite score
4. BSA burial  5. Docking box  6. Structural overlay
7. Ranking table  8. Hydrophobic core ★
""")

    # ── Setup ──────────────────────────────────────────────
    code(f"""import nglview as nv
import pandas as pd
import json
import numpy as np
from pathlib import Path
from IPython.display import display

PROJECT  = Path.home() / 'projects' / 'rtc-pan-coronavirus'
PDB_DIR  = PROJECT / '00-reference' / 'pdb_structures'
AF3_DIR  = PROJECT / '01-alphafold3' / 'NSP12-NSP7'
RES_DIR  = PROJECT / '02-validation' / 'NSP12-NSP7'

p_7bv2 = str(PDB_DIR / '7BV2_NSP12-NSP7.pdb')
p_6nur = str(PDB_DIR / '6NUR_NSP12-NSP7.pdb')
p_7c2k = str(PDB_DIR / '7C2K_NSP12-NSP7.pdb')
p_af3  = str(AF3_DIR / 'NSP12_NSP7_best_model.pdb')

df_rank = pd.read_csv(RES_DIR / 'composite_ranking_NSP12-NSP7_3.csv')
pocket  = json.load(open(RES_DIR / 'pocket_analysis_3.json'))
bsa_raw = json.load(open(RES_DIR / 'bsa_alascan_NSP12-NSP7_3.json'))['bsa']
box     = pocket['docking_box']

NSP12_OFFSET    = {NSP12_OFFSET}
NSP7_OFFSET     = {NSP7_OFFSET}
CHAIN_NSP12     = 'A'
CHAIN_NSP7      = 'C'
HOTSPOTS_NSP12  = {HOTSPOTS_NSP12}
HOTSPOTS_NSP7   = {HOTSPOTS_NSP7}
PRIMARY_NSP12   = {set(PRIMARY_NSP12)}
PRIMARY_NSP7    = {set(PRIMARY_NSP7)}
HYDROPHOBIC_CORE= {set(HYDROPHOBIC_CORE)}

def g12(local_list):
    return [p + NSP12_OFFSET for p in local_list]
def g7(local_list):
    return [p + NSP7_OFFSET for p in local_list]
def sel_resi(chain, resi_list):
    return ' or '.join(':' + str(chain) + ' and ' + str(r)
                       for r in resi_list)

print('Setup complete.')
print(f'Ranking: {{len(df_rank)}} residues')
print(f'Top pharmacophore: {{df_rank.iloc[0].residue}} '
      f'score={{df_rank.iloc[0].composite:.4f}}')
""")

    # ── View 1: Full complex ───────────────────────────────
    md("## View 1 — Full Complex: All Conserved Hotspots (7BV2)")
    code(f"""v1 = nv.show_file(p_7bv2, ext='pdb')
v1.clear_representations()

# Cartoon backbone
v1.add_representation('cartoon', selection=':A',
                      color='#AED6F1', opacity=0.85)
v1.add_representation('cartoon', selection=':C',
                      color='#A9DFBF', opacity=0.85)

# NSP12 hotspots — blue licorice
nsp12_sel = sel_resi('A', g12(HOTSPOTS_NSP12))
v1.add_representation('licorice', selection=nsp12_sel,
                      color='#2C5F8A', radius=0.3)

# NSP7 hotspots — green licorice
nsp7_sel = sel_resi('C', g7(HOTSPOTS_NSP7))
v1.add_representation('licorice', selection=nsp7_sel,
                      color='#27AE60', radius=0.3)

# Hydrophobic core — orange ball+stick
hcore_sel = sel_resi('A', g12(list(HYDROPHOBIC_CORE)))
v1.add_representation('ball+stick', selection=hcore_sel,
                      color='#E67E22', radius=0.4)

# Primary salt bridge — red ball+stick
prim12_sel = sel_resi('A', g12(list(PRIMARY_NSP12)))
prim7_sel  = sel_resi('C', g7(list(PRIMARY_NSP7)))
v1.add_representation('ball+stick', selection=prim12_sel,
                      color='#C0392B', radius=0.45)
v1.add_representation('ball+stick', selection=prim7_sel,
                      color='#C0392B', radius=0.45)

v1.center()
v1
""")

    # ── View 2: Salt bridge ────────────────────────────────
    md("""## View 2 — Salt Bridge Zoomed
**GLU431(NSP12) — LYS1(NSP7): 4.24 Å**
Crystallographically confirmed in 7BV2. AF3 predicted at 4.26 Å.
Note: SARS-CoV-1/2 specific — absent in MERS/HCoV (N-terminal truncation).
""")
    code(f"""v2 = nv.show_file(p_7bv2, ext='pdb')
v2.clear_representations()
v2.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v2.add_representation('cartoon', selection=':C',
                      color='#D5D8DC', opacity=0.3)

# GLU431 NSP12 — red
glu431 = sel_resi('A', g12([431]))
v2.add_representation('ball+stick', selection=glu431,
                      color='#C0392B', radius=0.55)
v2.add_representation('spacefill',  selection=glu431,
                      color='#C0392B', radius=0.35)

# LYS1 NSP7 — blue
lys1 = sel_resi('C', g7([1]))
v2.add_representation('ball+stick', selection=lys1,
                      color='#2980B9', radius=0.55)
v2.add_representation('spacefill',  selection=lys1,
                      color='#2980B9', radius=0.35)

# Context residues
ctx12 = sel_resi('A', g12([420,429,430,431,432,440]))
ctx7  = sel_resi('C', g7([1,2,3,4,5]))
v2.add_representation('licorice', selection=ctx12,
                      color='#E67E22', radius=0.25)
v2.add_representation('licorice', selection=ctx7,
                      color='#27AE60', radius=0.25)

v2.center(glu431)
v2
""")

    # ── View 3: Composite score ────────────────────────────
    md("## View 3 — Hotspots Colored by Composite Drug Target Score")
    code(f"""v3 = nv.show_file(p_7bv2, ext='pdb')
v3.clear_representations()
v3.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v3.add_representation('cartoon', selection=':C',
                      color='#D5D8DC', opacity=0.3)

for _, row in df_rank.iterrows():
    s = row['composite']
    if s >= 0.7:
        color = '#C0392B'
    elif s >= 0.4:
        color = '#E67E22'
    elif s >= 0.15:
        color = '#2C5F8A'
    else:
        color = '#BDC3C7'

    if row['chain'] == 'NSP12':
        gpos = int(row['position']) + NSP12_OFFSET
        ch   = 'A'
    else:
        gpos = int(row['position']) + NSP7_OFFSET
        ch   = 'C'
    v3.add_representation('ball+stick',
        selection=':' + ch + ' and ' + str(gpos),
        color=color, radius=0.35)

v3.center()
v3
""")

    # ── View 4: BSA ────────────────────────────────────────
    md("## View 4 — Hotspots Colored by BSA Burial Depth")
    code(f"""v4 = nv.show_file(p_7bv2, ext='pdb')
v4.clear_representations()
v4.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v4.add_representation('cartoon', selection=':C',
                      color='#D5D8DC', opacity=0.3)

bsa_map = {{}}
for key, val in bsa_raw.items():
    ch, pos = key.split(',')
    bsa_map[(ch.strip(), int(pos))] = float(val)

max_bsa = max(bsa_map.values()) if bsa_map else 1

for (chain, pos), bsa_val in bsa_map.items():
    norm = bsa_val / max_bsa
    if norm >= 0.6:
        color = '#148F77'
    elif norm >= 0.3:
        color = '#2980B9'
    elif norm >= 0.1:
        color = '#F39C12'
    else:
        continue
    if chain == 'A':
        gpos = pos + NSP12_OFFSET
    else:
        gpos = pos + NSP7_OFFSET
    v4.add_representation('ball+stick',
        selection=':' + chain + ' and ' + str(gpos),
        color=color, radius=0.32)

v4.center()
v4
""")

    # ── View 5: Docking box ────────────────────────────────
    md("## View 5 — Docking Box (centered on hydrophobic core + salt bridge)")
    code(f"""v5 = nv.show_file(p_7bv2, ext='pdb')
v5.clear_representations()
v5.add_representation('cartoon', selection=':A',
                      color='#AED6F1', opacity=0.7)
v5.add_representation('cartoon', selection=':C',
                      color='#A9DFBF', opacity=0.7)

# Hydrophobic core + primary SB
hcore_sel = sel_resi('A', g12(list(HYDROPHOBIC_CORE)))
prim12    = sel_resi('A', g12(list(PRIMARY_NSP12)))
prim7     = sel_resi('C', g7(list(PRIMARY_NSP7)))
v5.add_representation('ball+stick', selection=hcore_sel,
                      color='#E67E22', radius=0.4)
v5.add_representation('ball+stick', selection=prim12,
                      color='#C0392B', radius=0.45)
v5.add_representation('ball+stick', selection=prim7,
                      color='#C0392B', radius=0.45)

# Docking box
cx,cy,cz = box['center_x'],box['center_y'],box['center_z']
sx,sy,sz = box['size_x']/2,box['size_y']/2,box['size_z']/2
corners = [
    [cx-sx,cy-sy,cz-sz],[cx+sx,cy-sy,cz-sz],
    [cx+sx,cy+sy,cz-sz],[cx-sx,cy+sy,cz-sz],
    [cx-sx,cy-sy,cz+sz],[cx+sx,cy-sy,cz+sz],
    [cx+sx,cy+sy,cz+sz],[cx-sx,cy+sy,cz+sz],
]
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),
         (6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
for i,j in edges:
    v5.shape.add_cylinder(
        corners[i], corners[j], [0.9,0.5,0.1], 0.3)

v5.center()
v5
""")

    # ── View 6: Structural overlay ─────────────────────────
    md("""## View 6 — Structural Overlay: 7BV2 + 6NUR + 7C2K + AF3
Confirms interface conservation across all structures.""")
    code(f"""v6 = nv.show_file(p_7bv2, ext='pdb')
v6.clear_representations()
v6.add_component(p_6nur, ext='pdb')
v6.add_component(p_7c2k, ext='pdb')
v6.add_component(p_af3,  ext='pdb')

# Style each component
v6.add_representation('cartoon', color='#2C5F8A',
                      opacity=0.85, component=0)
for i, color in enumerate(['#C0392B','#27AE60','#F39C12'],1):
    v6.add_representation('cartoon', color=color,
                          opacity=0.45, component=i)

# Salt bridge + hydrophobic core on 7BV2
glu431 = sel_resi('A', g12([431]))
lys1   = sel_resi('C', g7([1]))
hcore  = sel_resi('A', g12(list(HYDROPHOBIC_CORE)))
v6.add_representation('ball+stick', selection=glu431,
                      color='#C0392B', radius=0.5, component=0)
v6.add_representation('ball+stick', selection=lys1,
                      color='#C0392B', radius=0.5, component=0)
v6.add_representation('ball+stick', selection=hcore,
                      color='#E67E22', radius=0.4, component=0)

v6.center()
print('Overlay: BLUE=7BV2 | RED=6NUR | GREEN=7C2K | GOLD=AF3')
v6
""")

    # ── View 7: Ranking table ──────────────────────────────
    md("## View 7 — Complete Drug Target Ranking Table")
    code("""display_cols = ['residue','chain','bsa','contact_score',
                'conservation','total_loss','composite',
                'primary_sb','hcore']
df_show = df_rank[display_cols].copy()
df_show.columns = ['Residue','Chain','BSA(A2)','Contact',
                   'Conservation','AlaLoss','Score',
                   'PrimarySB','HydroCore']

def highlight(row):
    if row['PrimarySB']:
        return ['background-color:#FADBD8']*len(row)
    elif row['HydroCore']:
        return ['background-color:#FDEBD0']*len(row)
    elif row['Chain'] == 'NSP7':
        return ['background-color:#D6EAF8']*len(row)
    return ['']*len(row)

styled = (df_show.style
          .apply(highlight, axis=1)
          .format({'Score':'{:.4f}','BSA(A2)':'{:.1f}',
                   'Conservation':'{:.3f}'})
          .set_caption(
              'NSP12-NSP7 Drug Target Ranking | '
              'Red=primary SB | Orange=hydrophobic core | '
              'Blue=NSP7'))
display(styled)
print(f'Top pharmacophore: {df_rank.iloc[0].residue} '
      f'(PHE440 — hydrophobic anchor)')
print(f'Primary SB: {list(df_rank[df_rank.primary_sb==True].residue)}')
""")

    # ── View 8: Hydrophobic core ★ ─────────────────────────
    md("""## View 8 ★ — Hydrophobic Core: PHE440/PHE442/PHE415/TYR420

**Key finding:** NSP12-NSP7 interface is hydrophobic-dominated.
The aromatic cluster (PHE440, PHE442, PHE415, TYR420, PHE843) forms
the primary binding pocket — all conserved 1.000 across 5 coronaviruses.

| Color | Meaning |
|-------|---------|
| 🟠 Orange | Hydrophobic core residues |
| 🔴 Red | PHE440 — primary pharmacophore |
| 🟣 Purple | GLU431 — salt bridge anchor |
| 🟢 Green | LYS1(NSP7) — salt bridge partner |

**Drug design:** Target the aromatic cluster for pan-coronavirus coverage.
GLU431–LYS1 salt bridge adds selectivity for SARS-CoV-1/2.
""")
    code(f"""v8 = nv.show_file(p_7bv2, ext='pdb')
v8.clear_representations()

# Full complex — very muted
v8.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.2)
v8.add_representation('cartoon', selection=':C',
                      color='#D5D8DC', opacity=0.2)

# Hydrophobic core neighbourhood cartoon — orange
hcore_nb = g12([409,410,411,412,413,414,415,
                419,420,421,438,439,440,441,
                442,443,841,842,843,844])
hcore_nb_sel = sel_resi('A', hcore_nb)
v8.add_representation('cartoon', selection=hcore_nb_sel,
                      color='#E67E22', opacity=0.9)

# All hydrophobic core — orange ball+stick
hcore_sel = sel_resi('A', g12(list(HYDROPHOBIC_CORE)))
v8.add_representation('ball+stick', selection=hcore_sel,
                      color='#E67E22', radius=0.4)

# PHE440 — red (primary pharmacophore)
phe440 = sel_resi('A', g12([440]))
v8.add_representation('ball+stick', selection=phe440,
                      color='#C0392B', radius=0.55)
v8.add_representation('spacefill',  selection=phe440,
                      color='#C0392B', radius=0.35)

# GLU431 — purple (salt bridge NSP12)
glu431 = sel_resi('A', g12([431]))
v8.add_representation('ball+stick', selection=glu431,
                      color='#7D3C98', radius=0.5)

# LYS1 NSP7 — green (salt bridge partner)
lys1 = sel_resi('C', g7([1]))
v8.add_representation('ball+stick', selection=lys1,
                      color='#27AE60', radius=0.5)
v8.add_representation('spacefill',  selection=lys1,
                      color='#27AE60', radius=0.3)

v8.center(hcore_nb_sel)
print('Hydrophobic core visualization:')
print('  RED    = PHE440 — primary pharmacophore (score=1.000)')
print('  ORANGE = Hydrophobic core (PHE412/413/415/420/442/843)')
print('  PURPLE = GLU431 (NSP12) — salt bridge anchor')
print('  GREEN  = LYS1 (NSP7)  — salt bridge partner')
print()
print('All hydrophobic core residues: cons=1.000 (pan-coronavirus)')
print('GLU431-LYS1 salt bridge: SARS-CoV-1/2 specific only')
v8
""")

    # ── Save notebook ──────────────────────────────────────
    nb = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language":     "python",
                "name":         "python3"
            },
            "language_info": {
                "name":    "python",
                "version": "3.10.0"
            }
        },
        "cells": CELLS
    }
    out = NB_DIR / "NSP12-NSP7_3D_3.ipynb"
    import json as _j
    _j.dump(nb, open(out,"w"), indent=2)
    print(f"  Saved: {out}")
    return out


def main():
    print("\n" + "="*57)
    print("  Script 11_3: 3D Notebook (nglview) — NSP12-NSP7")
    print("  8 views including hydrophobic core ★")
    print("="*57)
    out = build()
    print(f"\n  Launch: jupyter notebook {out}")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
