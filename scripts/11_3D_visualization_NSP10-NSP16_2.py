"""
Script 11_2: 3D Visualization Notebook — NSP10-NSP16 (nglview)
==============================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

8 interactive views using nglview (JupyterLab native):
  1. Full complex — all conserved hotspots
  2. Salt bridges zoomed — 3 confirmed pairs
  3. Hotspots colored by composite score
  4. Hotspots colored by BSA burial depth
  5. Docking box visualization
  6. Structural overlay — 6W4H + 6WVN + 6WKQ + AF3
  7. Full hotspot ranking table
  8. Zn1 zinc finger context ★
"""

import json
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP10-NSP16"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
NB_DIR  = PROJECT / "notebooks"
NB_DIR.mkdir(exist_ok=True)

NSP10_OFFSET   = 4270
NSP10_SHIFT    = 17
NSP16_OFFSET   = 6797
HOTSPOTS_NSP10 = [40,42,43,44,45,71,76,78,80,93,94,95,96]
HOTSPOTS_NSP16 = [40,41,44,76,83,87,102,104,106,244,247]
PRIMARY_NSP10  = {80, 93, 95}
PRIMARY_NSP16  = {102, 106}
ZN1_COORD      = {74, 77, 83, 90}
ZN1_NEIGHBOURS = {71,74,76,77,78,80,83,90,93}

def g10(local_list):
    return [p + NSP10_OFFSET - NSP10_SHIFT
            for p in local_list]

def g16(local_list):
    return [p + NSP16_OFFSET for p in local_list]

def sel10(local_list):
    return " or ".join(f":{chr(66)} and {r}"
                       for r in g10(local_list))

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
    # ── Header ────────────────────────────────────────────
    md("""# NSP10–NSP16 Interface: 3D Structural Visualization
**Project:** panCov-rtc-discovery  
**Complex:** NSP10–NSP16 (2'-O-methyltransferase activation)  
**Primary PDB:** 6W4H (SARS-CoV-2, 1.80 Å)  

## Key findings:
- **3 confirmed salt bridges:** ASP106–LYS93, ASP106–LYS95, ASP102–HIS80
- **NSP16 interface: 11/11 hotspots conserved = 1.000**
- **HIS80 inside Zn1 zinc finger loop** — dual mechanism potential
- **fpocket druggability = 0.546** — highest in project
- **Water bridges at HIS80 and LYS93** — displacement pharmacophores

## Views:
1. Full complex  2. Salt bridges  3. Composite score
4. BSA burial  5. Docking box  6. Structural overlay
7. Ranking table  8. Zn1 finger ★
""")

    # ── Setup cell ────────────────────────────────────────
    code(f"""import nglview as nv
import pandas as pd
import json
import numpy as np
from pathlib import Path

PROJECT = Path.home() / 'projects' / 'rtc-pan-coronavirus'
PDB_DIR = PROJECT / '00-reference' / 'pdb_structures'
AF3_DIR = PROJECT / '01-alphafold3' / 'NSP10-NSP16'
RES_DIR = PROJECT / '02-validation' / 'NSP10-NSP16'

# Paths
p_6w4h = str(PDB_DIR / '6W4H.pdb')
p_6wvn = str(PDB_DIR / '6WVN.pdb')
p_6wkq = str(PDB_DIR / '6WKQ.pdb')
p_af3  = str(AF3_DIR / 'NSP10_NSP16_best_model.pdb')

df_rank = pd.read_csv(RES_DIR / 'composite_ranking_NSP10-NSP16_2.csv')
pocket  = json.load(open(RES_DIR / 'pocket_analysis_2.json'))
bsa_raw = json.load(open(RES_DIR / 'bsa_alascan_NSP10-NSP16_2.json'))['bsa']
box     = pocket['docking_box']

# Residue position constants
NSP10_OFFSET = {NSP10_OFFSET}
NSP10_SHIFT  = {NSP10_SHIFT}
NSP16_OFFSET = {NSP16_OFFSET}

def g10(local_list):
    return [p + NSP10_OFFSET - NSP10_SHIFT for p in local_list]
def g16(local_list):
    return [p + NSP16_OFFSET for p in local_list]
def sel_resi(chain, resi_list):
    return ' or '.join(':' + str(chain) + ' and ' + str(r) for r in resi_list)

HOTSPOTS_NSP10 = [40,42,43,44,45,71,76,78,80,93,94,95,96]
HOTSPOTS_NSP16 = [40,41,44,76,83,87,102,104,106,244,247]
PRIMARY_NSP10  = {{80,93,95}}
PRIMARY_NSP16  = {{102,106}}
ZN1_COORD      = {{74,77,83,90}}
ZN1_NEIGHBOURS = {{71,74,76,77,78,80,83,90,93}}

print('Setup complete.')
print(f'Ranking: {{len(df_rank)}} residues')
print(f'Box center: ({{box[\"center_x\"]}}, {{box[\"center_y\"]}}, {{box[\"center_z\"]}})')
""")

    # ── View 1: Full complex ───────────────────────────────
    md("## View 1 — Full Complex: All Conserved Hotspots (6W4H)")
    code(f"""v1 = nv.show_file(p_6w4h, ext='pdb')
v1.clear_representations()

# Cartoon backbone
v1.add_representation('cartoon', selection=':A', color='#AED6F1', opacity=0.85)
v1.add_representation('cartoon', selection=':B', color='#A9DFBF', opacity=0.85)

# NSP16 hotspots — blue licorice
nsp16_sel = sel_resi('A', g16(HOTSPOTS_NSP16))
v1.add_representation('licorice', selection=nsp16_sel, color='#2C5F8A', radius=0.3)

# NSP10 hotspots — green licorice
nsp10_sel = sel_resi('B', g10(HOTSPOTS_NSP10))
v1.add_representation('licorice', selection=nsp10_sel, color='#27AE60', radius=0.3)

# Primary salt bridge residues — red ball+stick
prim10_sel = sel_resi('B', g10(list(PRIMARY_NSP10)))
prim16_sel = sel_resi('A', g16(list(PRIMARY_NSP16)))
v1.add_representation('ball+stick', selection=prim10_sel, color='#C0392B', radius=0.4)
v1.add_representation('ball+stick', selection=prim16_sel, color='#C0392B', radius=0.4)

v1.center()
v1
""")

    # ── View 2: Salt bridges ───────────────────────────────
    md("""## View 2 — Salt Bridges Zoomed
**3 confirmed pairs:**
- ASP106(NSP16) — LYS93(NSP10): 3.90 Å [all 4 structures]
- ASP106(NSP16) — LYS95(NSP10): 3.90 Å [all 4 structures]
- ASP102(NSP16) — HIS80(NSP10): 3.33 Å [6WVN + 6WKQ + AF3]
""")
    sb16 = g16([102, 106])
    sb10 = g10([80, 93, 95])
    code(f"""v2 = nv.show_file(p_6w4h, ext='pdb')
v2.clear_representations()

# Muted cartoon
v2.add_representation('cartoon', selection=':A', color='#D5D8DC', opacity=0.3)
v2.add_representation('cartoon', selection=':B', color='#D5D8DC', opacity=0.3)

# NSP16 salt bridge residues — red
sb16_sel = sel_resi('A', {sb16})
v2.add_representation('ball+stick', selection=sb16_sel, color='#C0392B', radius=0.5)
v2.add_representation('spacefill',  selection=sb16_sel, color='#C0392B', radius=0.3)

# NSP10 salt bridge residues — blue
sb10_sel = sel_resi('B', {sb10})
v2.add_representation('ball+stick', selection=sb10_sel, color='#2980B9', radius=0.5)
v2.add_representation('spacefill',  selection=sb10_sel, color='#2980B9', radius=0.3)

# Zoom to interface
v2.center(sb16_sel)
v2
""")

    # ── View 3: Composite score ────────────────────────────
    md("## View 3 — Hotspots Colored by Composite Drug Target Score")
    code(f"""v3 = nv.show_file(p_6w4h, ext='pdb')
v3.clear_representations()
v3.add_representation('cartoon', selection=':A', color='#D5D8DC', opacity=0.3)
v3.add_representation('cartoon', selection=':B', color='#D5D8DC', opacity=0.3)

# Color tiers by composite score
high   = df_rank[df_rank['composite'] >= 0.5]
mid    = df_rank[(df_rank['composite'] >= 0.15) & (df_rank['composite'] < 0.5)]
low    = df_rank[(df_rank['composite'] >= 0.05) & (df_rank['composite'] < 0.15)]

for tier, color in [(high,'#C0392B'),(mid,'#E67E22'),(low,'#2C5F8A')]:
    for _, row in tier.iterrows():
        if row['chain'] == 'NSP10':
            gpos = int(row['position']) + NSP10_OFFSET - NSP10_SHIFT
            ch = 'B'
        else:
            gpos = int(row['position']) + NSP16_OFFSET
            ch = 'A'
        sel = f':{{ch}} and {{gpos}}'
        v3.add_representation('ball+stick', selection=sel,
                              color=color, radius=0.4)

v3.center()
v3
""")

    # ── View 4: BSA burial ─────────────────────────────────
    md("## View 4 — Hotspots Colored by BSA Burial Depth")
    code(f"""v4 = nv.show_file(p_6w4h, ext='pdb')
v4.clear_representations()
v4.add_representation('cartoon', selection=':A', color='#D5D8DC', opacity=0.3)
v4.add_representation('cartoon', selection=':B', color='#D5D8DC', opacity=0.3)

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
    if chain == 'B':
        gpos = pos + NSP10_OFFSET - NSP10_SHIFT
    else:
        gpos = pos + NSP16_OFFSET
    v4.add_representation('ball+stick',
                          selection=f':{{chain}} and {{gpos}}',
                          color=color, radius=0.35)

v4.center()
v4
""")

    # ── View 5: Docking box ────────────────────────────────
    md("## View 5 — Docking Box (primary hotspot-centered)")
    prim10_g = g10(list(PRIMARY_NSP10))
    prim16_g = g16(list(PRIMARY_NSP16))
    code(f"""v5 = nv.show_file(p_6w4h, ext='pdb')
v5.clear_representations()
v5.add_representation('cartoon', selection=':A', color='#AED6F1', opacity=0.7)
v5.add_representation('cartoon', selection=':B', color='#A9DFBF', opacity=0.7)

# Primary residues
v5.add_representation('ball+stick',
    selection=sel_resi('B', {prim10_g}), color='#C0392B', radius=0.4)
v5.add_representation('ball+stick',
    selection=sel_resi('A', {prim16_g}), color='#C0392B', radius=0.4)

# Docking box as shape
cx,cy,cz = box['center_x'],box['center_y'],box['center_z']
sx,sy,sz = box['size_x']/2,box['size_y']/2,box['size_z']/2

shape = v5.shape
corners = [
    [cx-sx,cy-sy,cz-sz],[cx+sx,cy-sy,cz-sz],
    [cx+sx,cy+sy,cz-sz],[cx-sx,cy+sy,cz-sz],
    [cx-sx,cy-sy,cz+sz],[cx+sx,cy-sy,cz+sz],
    [cx+sx,cy+sy,cz+sz],[cx-sx,cy+sy,cz+sz],
]
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),
         (6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
for i,j in edges:
    shape.add_cylinder(corners[i], corners[j],
                       [0.9,0.5,0.1], 0.3)

v5.center()
v5
""")

    # ── View 6: Structural overlay ─────────────────────────
    md("""## View 6 — Structural Overlay: 6W4H + 6WVN + 6WKQ + AF3
Confirms salt bridge conservation across species and prediction methods.""")
    code("""v6 = nv.NGLWidget()

# Load all 4 structures
c1 = v6.add_component(nv.FileStructure(p_6w4h, ext='pdb'))
c2 = v6.add_component(nv.FileStructure(p_6wvn, ext='pdb'))
c3 = v6.add_component(nv.FileStructure(p_6wkq, ext='pdb'))
c4 = v6.add_component(nv.FileStructure(p_af3,  ext='pdb'))

c1.clear_representations()
c2.clear_representations()
c3.clear_representations()
c4.clear_representations()

# 6W4H — blue (primary)
c1.add_representation('cartoon', color='#2C5F8A', opacity=0.8)
# 6WVN — red (SARS-CoV-1)
c2.add_representation('cartoon', color='#C0392B', opacity=0.5)
# 6WKQ — green
c3.add_representation('cartoon', color='#27AE60', opacity=0.5)
# AF3  — gold
c4.add_representation('cartoon', color='#F39C12', opacity=0.5)

# Salt bridge residues on 6W4H
sb16 = sel_resi('A', g16([102,106]))
sb10 = sel_resi('B', g10([80,93,95]))
c1.add_representation('ball+stick', selection=sb16, color='#C0392B', radius=0.5)
c1.add_representation('ball+stick', selection=sb10, color='#2C5F8A', radius=0.5)

v6.center()
v6
""")

    # ── View 7: Ranking table ──────────────────────────────
    md("## View 7 — Complete Drug Target Ranking Table")
    code("""from IPython.display import display

display_cols = ['residue','chain','bsa','contact_score',
                'conservation','total_loss','composite',
                'primary_sb','zn1_region']
df_show = df_rank[display_cols].copy()
df_show.columns = ['Residue','Chain','BSA(Å²)','Contact',
                   'Conservation','AlaLoss','Score',
                   'PrimarySB','Zn1Region']

def highlight(row):
    if row['PrimarySB']:
        return ['background-color:#FADBD8']*len(row)
    elif row['Zn1Region']:
        return ['background-color:#FDEBD0']*len(row)
    elif row['Chain'] == 'NSP16':
        return ['background-color:#D6EAF8']*len(row)
    return ['']*len(row)

styled = (df_show.style
          .apply(highlight, axis=1)
          .format({'Score':'{:.4f}','BSA(Å²)':'{:.1f}',
                   'Conservation':'{:.3f}'})
          .set_caption(
              'NSP10-NSP16 Drug Target Ranking | '
              'Red=primary SB | Orange=Zn1 | Blue=NSP16'))
display(styled)
print(f'Primary SB residues: '
      f'{list(df_rank[df_rank.primary_sb==True].residue)}')
print(f'Zn1 region: '
      f'{list(df_rank[df_rank.zn1_region==True].residue)}')
""")

    # ── View 8: Zn1 finger ────────────────────────────────
    md("""## View 8 ★ — Zn1 Zinc Finger: Coordination + Interface Context

**HIS80** (primary salt bridge partner of ASP102) sits **inside the Zn1 finger loop**.

| Color | Meaning |
|-------|---------|
| 🟡 Gold | Zn1 coordinators: CYS74, CYS77, HIS83, CYS90 |
| 🔴 Red | HIS80 — salt bridge anchor inside Zn1 loop |
| 🟣 Purple | ASP102 (NSP16) — salt bridge partner |
| 🔵 Cyan sphere | Zn atom |

**Drug implication:** Targeting HIS80–ASP102 = PPI disruption + potential Zn perturbation
""")
    zn1_g   = g10(list(ZN1_COORD))
    zn1nb_g = g10(list(ZN1_NEIGHBOURS))
    his80_g = g10([80])[0]
    asp102_g= g16([102])[0]
    code(f"""v8 = nv.show_file(p_6w4h, ext='pdb')
v8.clear_representations()

# Full complex — very muted
v8.add_representation('cartoon', selection=':A', color='#D5D8DC', opacity=0.2)
v8.add_representation('cartoon', selection=':B', color='#D5D8DC', opacity=0.2)

# Zn1 neighbourhood cartoon — gold
zn1nb_sel = sel_resi('B', {zn1nb_g})
v8.add_representation('cartoon', selection=zn1nb_sel,
                      color='#F39C12', opacity=0.9)

# Zn1 coordinators — gold ball+stick
zn1_sel = sel_resi('B', {zn1_g})
v8.add_representation('ball+stick', selection=zn1_sel,
                      color='#F39C12', radius=0.4)

# HIS80 — red (salt bridge + Zn1 context)
v8.add_representation('ball+stick', selection=':B and {his80_g}',
                      color='#C0392B', radius=0.55)
v8.add_representation('spacefill',  selection=':B and {his80_g}',
                      color='#C0392B', radius=0.35)

# ASP102 NSP16 — purple
v8.add_representation('ball+stick', selection=':A and {asp102_g}',
                      color='#7D3C98', radius=0.55)
v8.add_representation('spacefill',  selection=':A and {asp102_g}',
                      color='#7D3C98', radius=0.35)

# Zn atom — cyan spacefill
v8.add_representation('spacefill', selection=':B and ZN',
                      color='#00BCD4', radius=1.2)

# Zoom to Zn1 region
v8.center(zn1nb_sel)
v8
""")

    # ── Save notebook ─────────────────────────────────────
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
    out = NB_DIR / "NSP10-NSP16_3D_2.ipynb"
    import json as _j
    _j.dump(nb, open(out,"w"), indent=2)
    print(f"Saved: {out}")
    return out


def main():
    print("\n" + "="*55)
    print("  Script 11_2: 3D Notebook (nglview) — NSP10-NSP16")
    print("="*55)
    out = build()
    print(f"\n  Launch: jupyter lab {out}")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
