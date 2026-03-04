"""
Script 11_5: 3D Visualization Notebook — NSP9-NSP12
====================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

8 interactive views:
  1. Full complex — all conserved hotspots
  2. ARG733 primary pharmacophore ★ zoomed
  3. AF3 salt bridges — ASP740/GLU744–LYS36
  4. Hotspots colored by composite score
  5. Hotspots colored by BSA burial
  6. Docking box visualization
  7. Structural overlay — 8SQK + AF3
  8. Full ranking table

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/11_3D_visualization_NSP9-NSP12_5.py
  jupyter notebook notebooks/NSP9-NSP12_3D_5.ipynb
"""

import json
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP9-NSP12"
RES_DIR = PROJECT / "02-validation" / "NSP9-NSP12"
NB_DIR  = PROJECT / "notebooks"
NB_DIR.mkdir(exist_ok=True)

HOTSPOTS_NSP12 = [38,1,3,4,96,733,202,103,
                   221,233,291,2,223]
HOTSPOTS_NSP9  = [38,1,3,4,96,103,2,97]
PRIMARY_NSP12  = {740,744}
PRIMARY_NSP9   = {36}
PAN_COV_NSP12  = {733,202,221}
PAN_COV_NSP9   = {97,103}

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
    md("""# NSP9–NSP12 Interface: 3D Structural Visualization
**Project:** panCov-rtc-discovery
**Complex:** NSP9–NSP12 (NiRAN domain interface — novel target)
**Primary PDB:** 8SQK (SARS-CoV-2, Chain A=NSP12, Chain G=NSP9)
**AF3 model:** iptm=0.77, ranking=0.79, PAE=3.52 Å

## Key findings:
- **ARG733(NSP12) primary pharmacophore** — NiRAN domain, cons=1.000 pan-coronavirus
- **H-bond + hydrophobic interface** — no crystallographic salt bridges
- **AF3 predicts ASP740/GLU744–LYS36 SBs** — SARS-CoV-1/2 specific only
- **Pan-coronavirus anchors:** ARG733, ASP221, VAL202 (NSP12) + LEU97, LEU103 (NSP9)
- **fpocket druggability = 0.895** — highest single-structure score in project
- **AF3 validation F1=0.837**

## Views:
1. Full complex  2. ARG733 pharmacophore ★  3. AF3 salt bridges
4. Composite score  5. BSA burial  6. Docking box
7. Structural overlay  8. Ranking table
""")

    code(f"""import nglview as nv
import pandas as pd
import json
import numpy as np
from pathlib import Path
from IPython.display import display

PROJECT = Path.home() / 'projects' / 'rtc-pan-coronavirus'
PDB_DIR = PROJECT / '00-reference' / 'pdb_structures'
AF3_DIR = PROJECT / '01-alphafold3' / 'NSP9-NSP12'
RES_DIR = PROJECT / '02-validation' / 'NSP9-NSP12'

p_8sqk = str(PDB_DIR / '8SQK_NSP9-NSP12.pdb')
p_af3  = str(AF3_DIR / 'NSP9_NSP12_best_model.pdb')

df_rank = pd.read_csv(RES_DIR / 'composite_ranking_NSP9-NSP12_5.csv')
pocket  = json.load(open(RES_DIR / 'pocket_analysis_5.json'))
bsa_raw = json.load(open(RES_DIR / 'bsa_alascan_NSP9-NSP12_5.json'))['bsa']
box     = pocket['docking_box']

CHAIN_NSP12    = 'A'
CHAIN_NSP9     = 'G'
HOTSPOTS_NSP12 = {HOTSPOTS_NSP12}
HOTSPOTS_NSP9  = {HOTSPOTS_NSP9}
PRIMARY_NSP12  = {set(PRIMARY_NSP12)}
PRIMARY_NSP9   = {set(PRIMARY_NSP9)}
PAN_COV_NSP12  = {set(PAN_COV_NSP12)}
PAN_COV_NSP9   = {set(PAN_COV_NSP9)}

def sel(chain, resi_list):
    return ' or '.join(f':{{chain}} and {{r}}'
                       for r in resi_list)

print('Setup complete.')
print(f'Top pharmacophore: {{df_rank.iloc[0].residue_aa}} '
      f'score={{df_rank.iloc[0].composite:.4f}}')
print(f'Primary pan-cov NSP12: ARG733, ASP221, VAL202')
print(f'Primary pan-cov NSP9:  LEU97, LEU103')
print(f'AF3-only SB (SARS-specific): ASP740/GLU744–LYS36')
""")

    # View 1: Full complex
    md("## View 1 — Full Complex: All Hotspots (8SQK)")
    code("""v1 = nv.show_file(p_8sqk, ext='pdb')
v1.clear_representations()
v1.add_representation('cartoon', selection=':A',
                      color='#AED6F1', opacity=0.85)
v1.add_representation('cartoon', selection=':G',
                      color='#A9DFBF', opacity=0.85)
# NSP12 hotspots — blue
v1.add_representation('licorice',
    selection=sel('A', HOTSPOTS_NSP12),
    color='#2C5F8A', radius=0.3)
# NSP9 hotspots — green
v1.add_representation('licorice',
    selection=sel('G', HOTSPOTS_NSP9),
    color='#27AE60', radius=0.3)
# Pan-cov NSP12 — gold
v1.add_representation('ball+stick',
    selection=sel('A', PAN_COV_NSP12),
    color='#F39C12', radius=0.45)
# Pan-cov NSP9 — gold
v1.add_representation('ball+stick',
    selection=sel('G', PAN_COV_NSP9),
    color='#F39C12', radius=0.45)
v1.center()
v1
""")

    # View 2: ARG733 primary pharmacophore
    md("""## View 2 ★ — ARG733 Primary Pharmacophore
**ARG733(NSP12)** — NiRAN domain anchor:
- Conservation = 1.000 (all 5 coronaviruses)
- Composite score = 1.000 (top ranked)
- H-bond dominant contact — druggability 0.895

| Color | Residue | Role |
|-------|---------|------|
| Dark Red | ARG733 (NSP12) | Primary pharmacophore ★ |
| Gold | ASP221 (NSP12) | Pan-cov anchor ◆ |
| Gold | VAL202 (NSP12) | Pan-cov anchor ◆ |
| Teal | LEU97 (NSP9) | Pan-cov hydrophobic ◆ |
| Teal | LEU103 (NSP9) | Pan-cov hydrophobic ◆ |
| Green | GLU3 (NSP9) | H-bond anchor |
""")
    code("""v2 = nv.show_file(p_8sqk, ext='pdb')
v2.clear_representations()
v2.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.2)
v2.add_representation('cartoon', selection=':G',
                      color='#D5D8DC', opacity=0.2)
# Context neighbourhood
nb12 = [730,731,732,733,734,735,
        200,201,202,203,219,220,221,222]
nb9  = [1,2,3,4,95,96,97,98,102,103,104]
v2.add_representation('cartoon',
    selection=sel('A', nb12), color='#AED6F1', opacity=0.8)
v2.add_representation('cartoon',
    selection=sel('G', nb9),  color='#A9DFBF', opacity=0.8)

# ARG733 — dark red primary
arg733 = sel('A', [733])
v2.add_representation('ball+stick', selection=arg733,
                      color='#922B21', radius=0.6)
v2.add_representation('spacefill',  selection=arg733,
                      color='#922B21', radius=0.45)
# Pan-cov NSP12 — gold
for pos in [221, 202]:
    v2.add_representation('ball+stick',
        selection=sel('A',[pos]),
        color='#F39C12', radius=0.45)
# Pan-cov NSP9 — teal
for pos in [97, 103]:
    v2.add_representation('ball+stick',
        selection=sel('G',[pos]),
        color='#148F77', radius=0.45)
# GLU3 — green H-bond anchor
v2.add_representation('ball+stick',
    selection=sel('G',[3]),
    color='#27AE60', radius=0.45)
# ASN1, ASN2 — grey context
v2.add_representation('licorice',
    selection=sel('G',[1,2]),
    color='#AED6F1', radius=0.25)

v2.center(arg733)
print('DKRED = ARG733 (primary pharmacophore, score=1.000)')
print('GOLD  = ASP221/VAL202 (pan-cov NSP12)')
print('TEAL  = LEU97/LEU103 (pan-cov NSP9)')
print('GREEN = GLU3 (H-bond anchor)')
v2
""")

    # View 3: AF3 salt bridges
    md("""## View 3 — AF3-Predicted Salt Bridges (SARS-selective)
**ASP740(NSP12)–LYS36(NSP9): 4.92 Å** — AF3 only, cons=0.172
**GLU744(NSP12)–LYS36(NSP9): 3.73 Å** — AF3 only, cons=0.345
⚠️ Charge reversal in MERS (K→Q/G) — NOT pan-coronavirus
""")
    code("""v3 = nv.show_file(p_af3, ext='pdb')
v3.clear_representations()
v3.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.2)
v3.add_representation('cartoon', selection=':B',
                      color='#D5D8DC', opacity=0.2)
# Context
v3.add_representation('cartoon',
    selection=sel('A',[737,738,739,740,741,742,743,744,745]),
    color='#AED6F1', opacity=0.8)
v3.add_representation('cartoon',
    selection=sel('B',[33,34,35,36,37,38]),
    color='#A9DFBF', opacity=0.8)
# ASP740 — red
v3.add_representation('ball+stick',
    selection=sel('A',[740]),
    color='#C0392B', radius=0.55)
# GLU744 — purple
v3.add_representation('ball+stick',
    selection=sel('A',[744]),
    color='#7D3C98', radius=0.55)
# LYS36 — orange (AF3 chain B)
v3.add_representation('ball+stick',
    selection=sel('B',[36]),
    color='#E67E22', radius=0.6)
v3.add_representation('spacefill',
    selection=sel('B',[36]),
    color='#E67E22', radius=0.4)
v3.center(sel('B',[36]))
print('RED    = ASP740 (NSP12) cons=0.172 — SARS only')
print('PURPLE = GLU744 (NSP12) cons=0.345 — SARS only')
print('ORANGE = LYS36  (NSP9)  cons=0.250 — SARS only')
print('WARNING: charge reversal in MERS/HCoV')
v3
""")

    # View 4: Composite score
    md("## View 4 — Hotspots Colored by Composite Score")
    code("""v4 = nv.show_file(p_8sqk, ext='pdb')
v4.clear_representations()
v4.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v4.add_representation('cartoon', selection=':G',
                      color='#D5D8DC', opacity=0.3)
for _, row in df_rank.iterrows():
    s = row['composite']
    color = ('#922B21' if s >= 0.8 else
             '#E67E22' if s >= 0.5 else
             '#2C5F8A' if s >= 0.2 else '#BDC3C7')
    ch   = 'A' if row['chain']=='NSP12' else 'G'
    gpos = int(row['position'])
    v4.add_representation('ball+stick',
        selection=f':{ch} and {gpos}',
        color=color, radius=0.35)
v4.center()
print('DKRED=score>=0.8 | ORANGE=0.5-0.8 | '
      'BLUE=0.2-0.5 | GREY<0.2')
v4
""")

    # View 5: BSA
    md("## View 5 — Hotspots Colored by BSA Burial Depth")
    code("""v5 = nv.show_file(p_8sqk, ext='pdb')
v5.clear_representations()
v5.add_representation('cartoon', selection=':A',
                      color='#D5D8DC', opacity=0.3)
v5.add_representation('cartoon', selection=':G',
                      color='#D5D8DC', opacity=0.3)
bsa_map = {}
for key, val in bsa_raw.items():
    ch, pos = key.split(',')
    bsa_map[(ch.strip(), int(pos))] = float(val)
max_bsa = max(bsa_map.values()) if bsa_map else 1
for (chain, pos), bsa_val in bsa_map.items():
    norm = bsa_val / max_bsa
    if norm >= 0.6:   color = '#148F77'
    elif norm >= 0.3: color = '#2980B9'
    elif norm >= 0.1: color = '#F39C12'
    else: continue
    pdb_ch = 'A' if chain=='A' else 'G'
    v5.add_representation('ball+stick',
        selection=f':{pdb_ch} and {pos}',
        color=color, radius=0.32)
v5.center()
print(f'TEAL=BSA>=60% | BLUE=30-60% | GOLD=10-30%')
print(f'Max BSA: {max_bsa:.1f} A2 (NSP9-ASN2)')
v5
""")

    # View 6: Docking box
    md("## View 6 — Docking Box (druggability=0.895)")
    code("""v6 = nv.show_file(p_8sqk, ext='pdb')
v6.clear_representations()
v6.add_representation('cartoon', selection=':A',
                      color='#AED6F1', opacity=0.7)
v6.add_representation('cartoon', selection=':G',
                      color='#A9DFBF', opacity=0.7)
v6.add_representation('ball+stick',
    selection=sel('A', list(PAN_COV_NSP12)),
    color='#F39C12', radius=0.45)
v6.add_representation('ball+stick',
    selection=sel('G', list(PAN_COV_NSP9)),
    color='#F39C12', radius=0.45)
v6.add_representation('ball+stick',
    selection=sel('A', [733]),
    color='#922B21', radius=0.5)
cx,cy,cz = box['center_x'],box['center_y'],box['center_z']
sx,sy,sz = box['size_x']/2,box['size_y']/2,box['size_z']/2
corners = [
    [cx-sx,cy-sy,cz-sz],[cx+sx,cy-sy,cz-sz],
    [cx+sx,cy+sy,cz-sz],[cx-sx,cy+sy,cz-sz],
    [cx-sx,cy-sy,cz+sz],[cx+sx,cy-sy,cz+sz],
    [cx+sx,cy+sy,cz+sz],[cx-sx,cy+sy,cz+sz],
]
for i,j in [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),
             (6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]:
    v6.shape.add_cylinder(
        corners[i],corners[j],[0.9,0.5,0.1],0.3)
v6.center()
print(f'Box: {box["size_x"]:.1f} x {box["size_y"]:.1f} '
      f'x {box["size_z"]:.1f} A | druggability=0.895')
v6
""")

    # View 7: Overlay
    md("## View 7 — Structural Overlay: 8SQK + AF3")
    code("""v7 = nv.show_file(p_8sqk, ext='pdb')
v7.clear_representations()
v7.add_component(p_af3, ext='pdb')
for i,color in enumerate(['#2C5F8A','#C0392B']):
    v7.add_representation('cartoon', color=color,
                          opacity=0.6 if i>0 else 0.9,
                          component=i)
# Pan-cov on 8SQK
v7.add_representation('ball+stick',
    selection=sel('A', list(PAN_COV_NSP12)),
    color='#F39C12', radius=0.45, component=0)
v7.add_representation('ball+stick',
    selection=sel('G', list(PAN_COV_NSP9)),
    color='#F39C12', radius=0.45, component=0)
v7.center()
print('BLUE=8SQK crystal | RED=AF3 model')
print('GOLD=pan-coronavirus anchors')
v7
""")

    # View 8: Ranking table
    md("## View 8 — Complete Drug Target Ranking Table")
    code("""display_cols = ['residue_aa','chain','bsa',
                'total_loss','conservation','composite',
                'primary_sb','pan_cov']
df_show = df_rank[display_cols].copy()
df_show.columns = ['Residue','Chain','BSA(A2)',
                   'AlaLoss','Conservation',
                   'Score','PrimarySB','PanCov']
def highlight(row):
    if row['PrimarySB']:
        return ['background-color:#FADBD8']*len(row)
    elif row['PanCov']:
        return ['background-color:#FEF9E7']*len(row)
    elif row['Chain']=='NSP9':
        return ['background-color:#D6EAF8']*len(row)
    return ['']*len(row)
styled = (df_show.style
    .apply(highlight, axis=1)
    .format({'Score':'{:.4f}','BSA(A2)':'{:.1f}',
             'Conservation':'{:.3f}'})
    .set_caption(
        'NSP9-NSP12 Drug Target Ranking | '
        'Red=AF3 SB | Gold=pan-cov | Blue=NSP9'))
display(styled)
top = df_rank.iloc[0]
print(f'Primary pharmacophore: {top.residue_aa} '
      f'(NiRAN domain, score=1.000)')
print(f'Druggability: 0.895 — highest in project')
""")

    nb = {
        "nbformat":4,"nbformat_minor":5,
        "metadata":{
            "kernelspec":{
                "display_name":"Python 3",
                "language":"python",
                "name":"python3"},
            "language_info":{
                "name":"python","version":"3.10.0"}},
        "cells":CELLS}

    out = NB_DIR / "NSP9-NSP12_3D_5.ipynb"
    json.dump(nb, open(out,"w"), indent=2)
    print(f"  Saved: {out}")
    return out


def main():
    print("\n" + "="*57)
    print("  Script 11_5: 3D Notebook — NSP9-NSP12")
    print("  8 views incl. ARG733 NiRAN pharmacophore ★")
    print("="*57)
    out = build()
    print(f"\n  Launch: jupyter notebook {out}")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
