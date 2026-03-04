"""
Script 11_6: 3D Visualization Notebook — NSP7-NSP8
===================================================
8 views:
  1. Full complex Mode B (AF3)
  2. PHE92 primary pharmacophore zoomed
  3. Mode A vs Mode B interface comparison
  4. Hotspots by composite score
  5. Hotspots by BSA burial
  6. Mode B docking box
  7. Structural overlay 7BV2 + 6NUR + AF3
  8. Full ranking table
"""
import json
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
NB_DIR  = PROJECT / "notebooks"
NB_DIR.mkdir(exist_ok=True)

HOTSPOTS_NSP8_B = [84,87,89,90,91,92,94,95,96,98,102,
                    103,106,107,110,111,116,119,120,150,190]
HOTSPOTS_NSP7_B = [2,6,12,13,16,19,28,35,49,50,52,53,
                    54,56,57,58,59,60,66,68,69,71,74,75,76]
MODE_A_NSP8     = [163,178,179,180]
MODE_A_NSP7     = [24,26,27]
PAN_COV_NSP8    = [87,91,92,94,98,110,111,116,120,150,190]

CELLS = []
def code(src): CELLS.append({"cell_type":"code",
    "execution_count":None,"metadata":{},"outputs":[],"source":src.strip()})
def md(src):   CELLS.append({"cell_type":"markdown",
    "metadata":{},"source":src.strip()})

md("""# NSP7–NSP8 Interface: 3D Structural Visualization
**Complex:** NSP7–NSP8 (RdRp subunit assembly)
**Primary structure:** AF3 (iptm=0.87, PAE=0.90 Å) — Mode B
**Secondary:** 7BV2/6NUR crystal — Mode A (undruggable)

## Dual binding mode discovery:
| Mode | Source | Interface | Contacts | Druggability |
|------|--------|-----------|----------|-------------|
| A | Crystal (7BV2/6NUR) | NSP8 C-term (163-180) / NSP7 loop | 2-4 | 0.276 ❌ |
| B | AF3 (iptm=0.87) | NSP8 N-helix (80-120) / NSP7 body | 45 | 0.531 ✅ |

## Primary pharmacophore: NSP8-PHE92 (hydrophobic core, cons=1.000)
## Views: 1.Full  2.PHE92★  3.Dual mode  4.Score  5.BSA  6.Box  7.Overlay  8.Table
""")

code(f"""import nglview as nv
import pandas as pd
import json
import numpy as np
from pathlib import Path
from IPython.display import display

PROJECT = Path.home() / 'projects' / 'rtc-pan-coronavirus'
PDB_DIR = PROJECT / '00-reference' / 'pdb_structures'
AF3_DIR = PROJECT / '01-alphafold3' / 'NSP7-NSP8'
RES_DIR = PROJECT / '02-validation' / 'NSP7-NSP8'

p_7bv2 = str(PDB_DIR / '7BV2_NSP7-NSP8.pdb')
p_6nur = str(PDB_DIR / '6NUR_NSP7-NSP8.pdb')
p_af3  = str(AF3_DIR / 'NSP7_NSP8_best_model.pdb')

df_rank = pd.read_csv(RES_DIR / 'composite_ranking_NSP7-NSP8_6.csv')
pocket  = json.load(open(RES_DIR / 'pocket_analysis_6.json'))
box_B   = pocket['docking_box_B']

# Chain IDs
CH8_PDB='B'; CH7_PDB='C'   # crystal
CH8_AF3='A'; CH7_AF3='B'   # AF3

HOT8_B = {HOTSPOTS_NSP8_B}
HOT7_B = {HOTSPOTS_NSP7_B}
HOT8_A = {MODE_A_NSP8}
HOT7_A = {MODE_A_NSP7}
PANCOV = {PAN_COV_NSP8}

def sel(ch, rlist):
    return ' or '.join(f':{{ch}} and {{r}}' for r in rlist)

print('NSP7-NSP8 dual binding mode visualization')
print(f'Primary pharmacophore: {{df_rank.iloc[0].residue}} '
      f'score={{df_rank.iloc[0].composite:.4f}}')
print(f'Mode B fpocket: 0.531 | Mode A: 0.276 (undruggable)')
""")

# View 1: Full complex Mode B
md("## View 1 — Full Complex Mode B (AF3 physiological interface)")
code("""v1 = nv.show_file(p_af3, ext='pdb')
v1.clear_representations()
v1.add_representation('cartoon', selection=f':{CH8_AF3}',
                      color='#AED6F1', opacity=0.85)
v1.add_representation('cartoon', selection=f':{CH7_AF3}',
                      color='#A9DFBF', opacity=0.85)
v1.add_representation('licorice',
    selection=sel(CH8_AF3, HOT8_B),
    color='#2C5F8A', radius=0.3)
v1.add_representation('licorice',
    selection=sel(CH7_AF3, HOT7_B),
    color='#27AE60', radius=0.3)
# Pan-cov NSP8 — gold
v1.add_representation('ball+stick',
    selection=sel(CH8_AF3, PANCOV),
    color='#F39C12', radius=0.4)
v1.center()
v1
""")

# View 2: PHE92 pharmacophore
md("""## View 2 ★ — PHE92 Primary Pharmacophore
**NSP8-PHE92** — hydrophobic aromatic core anchor:
- Conservation = 1.000 (all 5 coronaviruses)
- BSA = 102.2 Å², composite score = 1.000
- Supported by LEU91, LEU98, MET87, MET94 (all cons=1.000)
- ARG190–GLU50 SB in vicinity (Mode B only)

| Residue | Role | Cons |
|---------|------|------|
| PHE92 (NSP8) | Aromatic anchor ★ | 1.000 |
| LEU98 (NSP8) | Hydrophobic core | 1.000 |
| LEU91 (NSP8) | Hydrophobic core | 1.000 |
| ARG111 (NSP8) | Charged anchor | 1.000 |
| ARG190 (NSP8) | SB anchor | 1.000 |
| GLU50  (NSP7) | SB partner | 0.410 |
""")
code("""v2 = nv.show_file(p_af3, ext='pdb')
v2.clear_representations()
v2.add_representation('cartoon', selection=f':{CH8_AF3}',
                      color='#D5D8DC', opacity=0.2)
v2.add_representation('cartoon', selection=f':{CH7_AF3}',
                      color='#D5D8DC', opacity=0.2)
# Context
nb8 = list(range(88,100)) + [87,110,111,150,190]
nb7 = [49,50,51,52,53,54,55,56]
v2.add_representation('cartoon',
    selection=sel(CH8_AF3, nb8), color='#AED6F1', opacity=0.8)
v2.add_representation('cartoon',
    selection=sel(CH7_AF3, nb7), color='#A9DFBF', opacity=0.8)
# PHE92 — dark red primary
v2.add_representation('ball+stick',
    selection=f':{CH8_AF3} and 92',
    color='#922B21', radius=0.6)
v2.add_representation('spacefill',
    selection=f':{CH8_AF3} and 92',
    color='#922B21', radius=0.4)
# Hydrophobic core — blue
for pos in [91,98,87,94]:
    v2.add_representation('ball+stick',
        selection=f':{CH8_AF3} and {pos}',
        color='#2C5F8A', radius=0.45)
# ARG111, ARG190 — teal
for pos in [111,190]:
    v2.add_representation('ball+stick',
        selection=f':{CH8_AF3} and {pos}',
        color='#148F77', radius=0.45)
# GLU50 NSP7 — orange SB partner
v2.add_representation('ball+stick',
    selection=f':{CH7_AF3} and 50',
    color='#E67E22', radius=0.45)
v2.center(f':{CH8_AF3} and 92')
print('DKRED=PHE92★ | BLUE=LEU91/LEU98/MET87/MET94 | '
      'TEAL=ARG111/ARG190 | ORANGE=GLU50(NSP7 SB)')
v2
""")

# View 3: Dual mode comparison
md("""## View 3 — Dual Mode Interface Comparison
**BLUE (7BV2 crystal)** = Mode A interface (C-terminal, undruggable)
**RED (AF3)** = Mode B interface (N-terminal helix, physiological)
""")
code("""v3 = nv.show_file(p_7bv2, ext='pdb')
v3.clear_representations()
v3.add_component(p_af3, ext='pdb')
# Crystal — grey cartoon
v3.add_representation('cartoon', color='#D5D8DC',
                      opacity=0.5, component=0)
# AF3 — light cartoon
v3.add_representation('cartoon', color='#D5D8DC',
                      opacity=0.3, component=1)
# Mode A residues — purple (crystal)
v3.add_representation('ball+stick',
    selection=sel(CH8_PDB, HOT8_A),
    color='#7D3C98', radius=0.5, component=0)
v3.add_representation('ball+stick',
    selection=sel(CH7_PDB, HOT7_A),
    color='#7D3C98', radius=0.5, component=0)
# Mode B residues — red/green (AF3)
v3.add_representation('ball+stick',
    selection=sel(CH8_AF3, list(PANCOV)),
    color='#C0392B', radius=0.45, component=1)
v3.add_representation('ball+stick',
    selection=sel(CH7_AF3, [50,56,49]),
    color='#27AE60', radius=0.45, component=1)
v3.center()
print('PURPLE=Mode A crystal (undruggable, fpocket=0.276)')
print('RED=Mode B NSP8 pan-cov core (fpocket=0.531)')
print('GREEN=Mode B NSP7 contacts')
v3
""")

# View 4: Score
md("## View 4 — Hotspots Colored by Composite Score (Mode B)")
code("""v4 = nv.show_file(p_af3, ext='pdb')
v4.clear_representations()
v4.add_representation('cartoon', selection=f':{CH8_AF3}',
                      color='#D5D8DC', opacity=0.3)
v4.add_representation('cartoon', selection=f':{CH7_AF3}',
                      color='#D5D8DC', opacity=0.3)
mode_b = df_rank[df_rank['mode']=='B']
for _,row in mode_b.iterrows():
    s  = row['composite']
    ch = CH8_AF3 if row['chain']=='NSP8' else CH7_AF3
    color = ('#922B21' if s>=0.8 else
             '#E67E22' if s>=0.5 else
             '#2C5F8A' if s>=0.2 else '#BDC3C7')
    v4.add_representation('ball+stick',
        selection=f':{ch} and {int(row["position"])}',
        color=color, radius=0.35)
v4.center()
print('DKRED≥0.8 | ORANGE 0.5-0.8 | BLUE 0.2-0.5 | GREY<0.2')
v4
""")

# View 5: BSA
md("## View 5 — Hotspots Colored by BSA Burial")
code("""v5 = nv.show_file(p_af3, ext='pdb')
v5.clear_representations()
v5.add_representation('cartoon', selection=f':{CH8_AF3}',
                      color='#D5D8DC', opacity=0.3)
v5.add_representation('cartoon', selection=f':{CH7_AF3}',
                      color='#D5D8DC', opacity=0.3)
mode_b = df_rank[df_rank['mode']=='B']
max_bsa = mode_b['bsa'].max()
for _,row in mode_b.iterrows():
    norm = row['bsa']/max_bsa if max_bsa>0 else 0
    if norm>=0.5:   color='#148F77'
    elif norm>=0.25: color='#2980B9'
    elif norm>=0.1:  color='#F39C12'
    else: continue
    ch = CH8_AF3 if row['chain']=='NSP8' else CH7_AF3
    v5.add_representation('ball+stick',
        selection=f':{ch} and {int(row["position"])}',
        color=color, radius=0.32)
v5.center()
print(f'TEAL≥50% | BLUE 25-50% | GOLD 10-25% | max_bsa={max_bsa:.1f} A2')
v5
""")

# View 6: Docking box
md("## View 6 — Mode B Docking Box (druggability=0.531)")
code("""v6 = nv.show_file(p_af3, ext='pdb')
v6.clear_representations()
v6.add_representation('cartoon', selection=f':{CH8_AF3}',
                      color='#AED6F1', opacity=0.7)
v6.add_representation('cartoon', selection=f':{CH7_AF3}',
                      color='#A9DFBF', opacity=0.7)
v6.add_representation('ball+stick',
    selection=sel(CH8_AF3, list(PANCOV)),
    color='#F39C12', radius=0.45)
# PHE92 primary
v6.add_representation('spacefill',
    selection=f':{CH8_AF3} and 92',
    color='#922B21', radius=0.5)
cx,cy,cz = box_B['center_x'],box_B['center_y'],box_B['center_z']
sx,sy,sz = box_B['size_x']/2,box_B['size_y']/2,box_B['size_z']/2
corners = [[cx+dx,cy+dy,cz+dz]
           for dx in [-sx,sx] for dy in [-sy,sy]
           for dz in [-sz,sz]]
for i,j in [(0,1),(0,2),(0,4),(1,3),(1,5),(2,3),
             (2,6),(3,7),(4,5),(4,6),(5,7),(6,7)]:
    v6.shape.add_cylinder(corners[i],corners[j],
                          [0.9,0.5,0.1],0.3)
v6.center()
print(f'Box: {box_B["size_x"]:.1f} x {box_B["size_y"]:.1f} '
      f'x {box_B["size_z"]:.1f} A | druggability=0.531')
v6
""")

# View 7: Overlay
md("## View 7 — Structural Overlay: 7BV2 + 6NUR + AF3")
code("""v7 = nv.show_file(p_7bv2, ext='pdb')
v7.clear_representations()
v7.add_component(p_6nur, ext='pdb')
v7.add_component(p_af3,  ext='pdb')
for i,color,label in [(0,'#2C5F8A','7BV2'),
                       (1,'#27AE60','6NUR'),
                       (2,'#C0392B','AF3')]:
    v7.add_representation('cartoon', color=color,
                          opacity=0.7, component=i)
# PHE92 on AF3
v7.add_representation('ball+stick',
    selection=f':{CH8_AF3} and 92',
    color='#922B21', radius=0.6, component=2)
# Mode A on crystal
v7.add_representation('ball+stick',
    selection=sel(CH8_PDB, HOT8_A),
    color='#7D3C98', radius=0.45, component=0)
v7.center()
print('BLUE=7BV2 | GREEN=6NUR | RED=AF3')
print('DKRED=PHE92 | PURPLE=Mode A crystal residues')
v7
""")

# View 8: Ranking table
md("## View 8 — Complete Dual Mode Ranking Table")
code("""display_cols = ['residue','chain','mode','bsa',
                'total_loss','conservation','composite','primary_sb']
df_show = df_rank[display_cols].copy()
df_show.columns = ['Residue','Chain','Mode','BSA',
                   'AlaLoss','Cons','Score','PrimSB']
def hl(row):
    if row['PrimSB']:     return ['background:#FADBD8']*len(row)
    if row['Mode']=='B':
        if row['Chain']=='NSP8': return ['background:#D6EAF8']*len(row)
        return ['background:#D5F5E3']*len(row)
    return ['background:#F5EEF8']*len(row)
styled = (df_show.style
    .apply(hl, axis=1)
    .format({'Score':'{:.4f}','BSA':'{:.1f}','Cons':'{:.3f}'})
    .set_caption('NSP7-NSP8 Dual Mode Ranking | '
                 'Red=SB | Blue=NSP8 B | Green=NSP7 B | Purple=Mode A'))
display(styled)
top = df_rank.iloc[0]
print(f'Primary pharmacophore: {top.residue} '
      f'[Mode {top.mode}] score={top.composite:.4f}')
print(f'PHE92 cons=1.000 | hydrophobic aromatic core anchor')
print(f'Dual mode: Mode A undruggable, Mode B physiological')
""")

nb = {"nbformat":4,"nbformat_minor":5,
      "metadata":{"kernelspec":{"display_name":"Python 3",
                                 "language":"python",
                                 "name":"python3"},
                   "language_info":{"name":"python",
                                     "version":"3.10.0"}},
      "cells":CELLS}
out = NB_DIR / "NSP7-NSP8_3D_6.ipynb"
json.dump(nb, open(out,"w"), indent=2)
print(f"\n  Saved: {out}")
print(f"  Launch: jupyter notebook {out}\n")
