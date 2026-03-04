"""
Script 11_4: 3D Visualization Notebook — NSP12-NSP8 (nglview)
=============================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

8 interactive views:
  1. Full complex — all conserved hotspots
  2. Salt bridges zoomed — ASP523-ARG80 + LYS332-ASP99
  3. Hotspots colored by composite score
  4. Hotspots colored by BSA burial depth
  5. Docking box visualization
  6. Structural overlay — 7BV2 + 6NUR + 7C2K + AF3
  7. Full ranking table
  8. LYS332 primary pharmacophore ★

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/11_3D_visualization_NSP12-NSP8_4.py
  jupyter notebook notebooks/NSP12-NSP8_3D_4.ipynb
"""

import json
from pathlib import Path

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP8"
NB_DIR   = PROJECT / "notebooks"
NB_DIR.mkdir(exist_ok=True)

HOTSPOTS_NSP12 = [387,129,389,271,330,131,380,523,
                   91,87,332,95,117,517,99]
HOTSPOTS_NSP8  = [117,129,80,115,131,112,91,87,
                   116,95,98,83,128,90,121]
PRIMARY_NSP12  = {523,332,517}
PRIMARY_NSP8   = {80,99,79}
NSP12_OFFSET   = 0
NSP8_OFFSET    = 77

CELLS = []

def code(src):
    CELLS.append({"cell_type":"code",
                  "execution_count":None,
                  "metadata":{},"outputs":[],
                  "source":src.strip()})

def md(src):
    CELLS.append({"cell_type":"markdown",
                  "metadata":{},
                  "source":src.strip()})


def build():
    md("""# NSP12–NSP8 Interface: 3D Structural Visualization
**Project:** panCov-rtc-discovery
**Complex:** NSP12–NSP8 (RdRp core — fingers/palm subdomain)
**Primary PDB:** 7BV2_NSP12-NSP8 (SARS-CoV-2, 2.90 Å)
**AF3 model:** iptm=0.85, ranking=0.86

## Key findings:
- **LYS332(NSP12) primary pharmacophore** — SB loss=30, cons=1.000 pan-coronavirus
- **3 salt bridges confirmed:** ASP523-ARG80 (all 4), LYS332-ASP99 (3 structs), ASP517-LYS79 (3 structs)
- **Large hydrophobic anchors:** LEU387/LEU389 (BSA 155/118 Å²)
- **fpocket druggability = 0.874**
- **AF3 validation F1=0.934**

## Views:
1. Full complex  2. Salt bridges  3. Composite score
4. BSA burial    5. Docking box   6. Structural overlay
7. Ranking table 8. LYS332 pharmacophore ★
""")

    code(f"""import nglview as nv
import pandas as pd
import json
import numpy as np
from pathlib import Path
from IPython.display import display

PROJECT  = Path.home() / 'projects' / 'rtc-pan-coronavirus'
PDB_DIR  = PROJECT / '00-reference' / 'pdb_structures'
AF3_DIR  = PROJECT / '01-alphafold3' / 'NSP12-NSP8'
RES_DIR  = PROJECT / '02-validation' / 'NSP12-NSP8'

p_7bv2 = str(PDB_DIR / '7BV2_NSP12-NSP8.pdb')
p_6nur = str(PDB_DIR / '6NUR_NSP12-NSP8.pdb')
p_7c2k = str(PDB_DIR / '7C2K_NSP12-NSP8.pdb')
p_af3  = str(AF3_DIR / 'NSP12_NSP8_best_model.pdb')

df_rank = pd.read_csv(RES_DIR / 'composite_ranking_NSP12-NSP8_4.csv')
pocket  = json.load(open(RES_DIR / 'pocket_analysis_4.json'))
bsa_raw = json.load(open(RES_DIR / 'bsa_alascan_NSP12-NSP8_4.json'))['bsa']
box     = pocket['docking_box']

NSP12_OFFSET   = {NSP12_OFFSET}
NSP8_OFFSET    = {NSP8_OFFSET}
CHAIN_NSP12    = 'A'
CHAIN_NSP8     = 'B'
HOTSPOTS_NSP12 = {HOTSPOTS_NSP12}
HOTSPOTS_NSP8  = {HOTSPOTS_NSP8}
PRIMARY_NSP12  = {set(PRIMARY_NSP12)}
PRIMARY_NSP8   = {set(PRIMARY_NSP8)}

def g12(local_list):
    return [p + NSP12_OFFSET for p in local_list]
def g8(local_list):
    return [p + NSP8_OFFSET for p in local_list]
def sel(chain, resi_list):
    return ' or '.join(':' + chain + ' and ' + str(r)
                       for r in resi_list)

print('Setup complete.')
print(f'Top pharmacophore: {{df_rank.iloc[0].residue}} '
      f'score={{df_rank.iloc[0].composite:.4f}}')
print(f'Primary SB residues NSP12: {set(PRIMARY_NSP12)}')
print(f'Primary SB residues NSP8:  {set(PRIMARY_NSP8)}')
""")

    # View 1: Full complex
    md("## View 1 — Full Complex: All Conserved Hotspots (7BV2)")
    code(f"""v1 = nv.show_file(p_7bv2, ext='pdb')
v1.clear_representations()
v1.add_representation('cartoon', selection=':A',
                      color='#AED6F1', opacity=0.85)
v1.add_representation('cartoon', selection=':B',
                      color='#A9DFBF', opacity=0.85)

# NSP12 hotspots — blue
v1.add_representation('licorice',
    selection=sel('A', g12(HOTSPOTS_NSP12)),
    color='#2C5F8A', radius=0.3)
# NSP8 hotspots — green
v1.add_representation('licorice',
    selection=sel('B', g8(HOTSPOTS_NSP8)),
    color='#27AE60', radius=0.3)
# Primary SB NSP12 — red
v1.add_representation('ball+stick',
    selection=sel('A', g12(list(PRIMARY_NSP12))),
    color='#C0392B', radius=0.45)
# Primary SB NSP8 — red
v1.add_representation('ball+stick',
    selection=sel('B', g8(list(PRIMARY_NSP8))),
    color='#C0392B', radius=0.45)

v1.center()
v1
""")

    # View 2: Salt bridges
    md("""## View 2 — Salt Bridges Zoomed
**ASP523(NSP12)–ARG80(NSP8): 3.58 Å** — present in ALL 4 structures
**LYS332(NSP12)–ASP99(NSP8): 3.72 Å** — 7BV2+7C2K+AF3, cons=1.000 ★
**ASP517(NSP12)–LYS79(NSP8): 4.66 Å** — 6NUR+7C2K+AF3
""")
    code(f"""v2 = nv.show_file(p_7bv2, ext='pdb')
v2.clear_representations()
v2.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.25)
v2.add_representation('cartoon', selection=':B',
                      color='#D5D8DC', opacity=0.25)

# ASP523 — red
v2.add_representation('ball+stick',
    selection=sel('A', g12([523])),
    color='#C0392B', radius=0.55)
# ARG80 — orange
v2.add_representation('ball+stick',
    selection=sel('B', g8([80])),
    color='#E67E22', radius=0.55)

# LYS332 — dark red (primary pharmacophore)
v2.add_representation('ball+stick',
    selection=sel('A', g12([332])),
    color='#922B21', radius=0.6)
v2.add_representation('spacefill',
    selection=sel('A', g12([332])),
    color='#922B21', radius=0.4)
# ASP99 — blue
v2.add_representation('ball+stick',
    selection=sel('B', g8([99])),
    color='#2980B9', radius=0.55)

# ASP517 — purple
v2.add_representation('ball+stick',
    selection=sel('A', g12([517])),
    color='#7D3C98', radius=0.5)
# LYS79 — teal
v2.add_representation('ball+stick',
    selection=sel('B', g8([79])),
    color='#148F77', radius=0.5)

# Context
v2.add_representation('licorice',
    selection=sel('A', g12([330,331,332,333,
                             521,522,523,524,
                             515,516,517,518])),
    color='#AED6F1', radius=0.2)
v2.add_representation('licorice',
    selection=sel('B', g8([78,79,80,81,98,99,100])),
    color='#A9DFBF', radius=0.2)

v2.center(sel('A', g12([332])))
print('RED    = ASP523 (NSP12) | ORANGE = ARG80 (NSP8)')
print('DKRED  = LYS332 (NSP12) | BLUE   = ASP99 (NSP8)')
print('PURPLE = ASP517 (NSP12) | TEAL   = LYS79 (NSP8)')
v2
""")

    # View 3: Composite score
    md("## View 3 — Hotspots Colored by Composite Drug Target Score")
    code(f"""v3 = nv.show_file(p_7bv2, ext='pdb')
v3.clear_representations()
v3.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v3.add_representation('cartoon', selection=':B',
                      color='#D5D8DC', opacity=0.3)

for _, row in df_rank.iterrows():
    s = row['composite']
    color = ('#C0392B' if s >= 0.7 else
             '#E67E22' if s >= 0.4 else
             '#2C5F8A' if s >= 0.15 else '#BDC3C7')
    ch   = 'A' if row['chain']=='NSP12' else 'B'
    gpos = (int(row['position']) + NSP12_OFFSET
            if ch=='A'
            else int(row['position']) + NSP8_OFFSET)
    v3.add_representation('ball+stick',
        selection=':' + ch + ' and ' + str(gpos),
        color=color, radius=0.35)

v3.center()
print('RED=score≥0.7 | ORANGE=0.4-0.7 | '
      'BLUE=0.15-0.4 | GREY<0.15')
v3
""")

    # View 4: BSA
    md("## View 4 — Hotspots Colored by BSA Burial Depth")
    code(f"""v4 = nv.show_file(p_7bv2, ext='pdb')
v4.clear_representations()
v4.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v4.add_representation('cartoon', selection=':B',
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
    gpos = (pos + NSP12_OFFSET if chain=='A'
            else pos + NSP8_OFFSET)
    v4.add_representation('ball+stick',
        selection=':' + chain + ' and ' + str(gpos),
        color=color, radius=0.32)

v4.center()
print('TEAL=BSA≥60% | BLUE=30-60% | GOLD=10-30%')
print(f'Max BSA: {{max_bsa:.1f}} Å² (NSP12-387/LEU387)')
v4
""")

    # View 5: Docking box
    md("## View 5 — Docking Box (centered on primary SB cluster)")
    code(f"""v5 = nv.show_file(p_7bv2, ext='pdb')
v5.clear_representations()
v5.add_representation('cartoon', selection=':A',
                      color='#AED6F1', opacity=0.7)
v5.add_representation('cartoon', selection=':B',
                      color='#A9DFBF', opacity=0.7)

prim12 = sel('A', g12(list(PRIMARY_NSP12)))
prim8  = sel('B', g8(list(PRIMARY_NSP8)))
v5.add_representation('ball+stick', selection=prim12,
                      color='#C0392B', radius=0.45)
v5.add_representation('ball+stick', selection=prim8,
                      color='#C0392B', radius=0.45)

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
print(f'Box: {{box["size_x"]:.1f}} x {{box["size_y"]:.1f}} '
      f'x {{box["size_z"]:.1f}} A')
v5
""")

    # View 6: Overlay
    md("## View 6 — Structural Overlay: 7BV2 + 6NUR + 7C2K + AF3")
    code(f"""v6 = nv.show_file(p_7bv2, ext='pdb')
v6.clear_representations()
v6.add_component(p_6nur, ext='pdb')
v6.add_component(p_7c2k, ext='pdb')
v6.add_component(p_af3,  ext='pdb')

for i, color in enumerate(
        ['#2C5F8A','#C0392B','#27AE60','#F39C12']):
    v6.add_representation('cartoon', color=color,
                          opacity=0.5 if i>0 else 0.9,
                          component=i)

# Primary SB on 7BV2
v6.add_representation('ball+stick',
    selection=sel('A', g12(list(PRIMARY_NSP12))),
    color='#C0392B', radius=0.5, component=0)
v6.add_representation('ball+stick',
    selection=sel('B', g8(list(PRIMARY_NSP8))),
    color='#C0392B', radius=0.5, component=0)

v6.center()
print('BLUE=7BV2 | RED=6NUR | GREEN=7C2K | GOLD=AF3')
v6
""")

    # View 7: Ranking table
    md("## View 7 — Complete Drug Target Ranking Table")
    code("""display_cols = ['residue','chain','bsa',
                'contact_score','conservation',
                'total_loss','composite','primary_sb']
df_show = df_rank[display_cols].copy()
df_show.columns = ['Residue','Chain','BSA(A2)',
                   'Contact','Conservation',
                   'AlaLoss','Score','PrimarySB']

def highlight(row):
    if row['PrimarySB']:
        return ['background-color:#FADBD8']*len(row)
    elif row['Chain']=='NSP8':
        return ['background-color:#D6EAF8']*len(row)
    return ['']*len(row)

styled = (df_show.style
    .apply(highlight, axis=1)
    .format({'Score':'{:.4f}','BSA(A2)':'{:.1f}',
             'Conservation':'{:.3f}'})
    .set_caption(
        'NSP12-NSP8 Drug Target Ranking | '
        'Red=primary SB | Blue=NSP8'))
display(styled)
print(f'Top: {df_rank.iloc[0].residue} '
      f'(LYS332 — pan-coronavirus SB, score=1.000)')
""")

    # View 8: LYS332 pharmacophore ★
    md("""## View 8 ★ — LYS332 Primary Pharmacophore

**LYS332(NSP12)** is the dominant drug target:
- Alanine scanning loss = 30 (highest in complex)
- Conservation = 1.000 across all 5 coronaviruses
- Confirmed in 7BV2 (3.89 Å), 7C2K (4.94 Å), AF3 (3.72 Å)
- BSA = 95.4 Å²

| Color | Residue | Role |
|-------|---------|------|
| Dark Red | LYS332 (NSP12) | Primary pharmacophore ★ |
| Blue | ASP99 (NSP8) | Salt bridge partner (cons=1.000) |
| Red | ASP523 (NSP12) | Secondary SB (all 4 structs) |
| Orange | ARG80 (NSP8) | Secondary SB partner |
| Purple | ASP517 (NSP12) | Tertiary SB |
| Teal | LYS79 (NSP8) | Tertiary SB partner |
""")
    code(f"""v8 = nv.show_file(p_7bv2, ext='pdb')
v8.clear_representations()
v8.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.2)
v8.add_representation('cartoon', selection=':B',
                      color='#D5D8DC', opacity=0.2)

# SB neighbourhood cartoon
nb12 = g12([328,329,330,331,332,333,334,
            387,388,389,390,
            519,520,521,522,523,524,
            515,516,517,518])
nb8  = g8([77,78,79,80,81,82,
           97,98,99,100,101])
v8.add_representation('cartoon',
    selection=sel('A', nb12),
    color='#AED6F1', opacity=0.8)
v8.add_representation('cartoon',
    selection=sel('B', nb8),
    color='#A9DFBF', opacity=0.8)

# LYS332 — dark red primary
lys332 = sel('A', g12([332]))
v8.add_representation('ball+stick', selection=lys332,
                      color='#922B21', radius=0.6)
v8.add_representation('spacefill',  selection=lys332,
                      color='#922B21', radius=0.45)

# ASP99 — blue partner
v8.add_representation('ball+stick',
    selection=sel('B', g8([99])),
    color='#2980B9', radius=0.55)

# ASP523 — red secondary
v8.add_representation('ball+stick',
    selection=sel('A', g12([523])),
    color='#C0392B', radius=0.5)

# ARG80 — orange
v8.add_representation('ball+stick',
    selection=sel('B', g8([80])),
    color='#E67E22', radius=0.5)

# ASP517 + LYS79 — tertiary
v8.add_representation('ball+stick',
    selection=sel('A', g12([517])),
    color='#7D3C98', radius=0.45)
v8.add_representation('ball+stick',
    selection=sel('B', g8([79])),
    color='#148F77', radius=0.45)

# LEU387/389 hydrophobic anchors — gold
v8.add_representation('ball+stick',
    selection=sel('A', g12([387,389])),
    color='#F39C12', radius=0.4)

v8.center(lys332)
print('DKRED  = LYS332 (primary pharmacophore, score=1.000)')
print('BLUE   = ASP99  (pan-coronavirus SB partner)')
print('RED    = ASP523 (all 4 structures)')
print('ORANGE = ARG80  (all 4 structures)')
print('PURPLE = ASP517 | TEAL = LYS79')
print('GOLD   = LEU387/389 (hydrophobic anchors, BSA>115A2)')
v8
""")

    nb = {
        "nbformat":4, "nbformat_minor":5,
        "metadata":{
            "kernelspec":{
                "display_name":"Python 3",
                "language":"python",
                "name":"python3"},
            "language_info":{
                "name":"python","version":"3.10.0"}},
        "cells":CELLS}

    out = NB_DIR / "NSP12-NSP8_3D_4.ipynb"
    import json as _j
    _j.dump(nb, open(out,"w"), indent=2)
    print(f"  Saved: {out}")
    return out


def main():
    print("\n" + "="*57)
    print("  Script 11_4: 3D Notebook — NSP12-NSP8")
    print("  8 views including LYS332 pharmacophore ★")
    print("="*57)
    out = build()
    print(f"\n  Launch: jupyter notebook {out}")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
