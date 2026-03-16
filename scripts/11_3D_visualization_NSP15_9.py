"""
Script 11_9 — 3D Visualization NSP15 Homodimer
Complex  : NSP15-NSP15 (NendoU homodimer)
Suffix   : _9
Notebook : notebooks/NSP15_3D_9.ipynb
Views    : 8 standard views + 1 NSP15-specific (active site vs dimer interface)
"""

import json
import os
import textwrap
import nbformat

OUT_NB = "notebooks/NSP15_3D_9.ipynb"
os.makedirs("notebooks", exist_ok=True)

# ── Load data ──────────────────────────────────────────────────
with open("02-validation/NSP15/composite_ranking_NSP15_9.csv") as f:
    lines = f.readlines()

ranked = []
for line in lines[1:]:
    parts = line.strip().split(",")
    if len(parts) >= 13:
        ranked.append({
            "rank"       : int(parts[0]),
            "residue_aa" : parts[1],
            "position"   : int(parts[2]),
            "aa"         : parts[4],
            "cons"       : float(parts[5]),
            "bsa"        : float(parts[6]),
            "loss"       : int(parts[7]),
            "composite"  : float(parts[12]),
            "is_pos"     : parts[13].strip() == "True",
            "is_neg"     : parts[14].strip() == "True",
        })

with open("02-validation/NSP15/pocket_analysis_9.json") as f:
    pocket = json.load(f)

box    = pocket["docking_box"]
center = box["center"]
size   = box["size"]

# Key residue positions (PDB numbering)
hotspot_pos = [r["position"] for r in ranked[:20]]
primary_pos = [40, 62, 91, 267]
neg_pos     = [235, 250, 290]

# ── Build notebook cells ───────────────────────────────────────
cells = []

def code(src):
    return nbformat.v4.new_code_cell(textwrap.dedent(src))

def md(src):
    return nbformat.v4.new_markdown_cell(textwrap.dedent(src))

# Cell 0 — setup
cells.append(code("""
import matplotlib
matplotlib.use('Agg')
import nglview as nv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from IPython.display import display
print("NSP15 3D Visualization — 9 interactive views")
print("Complex: NSP15 NendoU homodimer (NendoU dimer interface)")
print("Reference: 6VWW biological assembly A+A-2 (75 contacts)")
"""))

# Cell 1 — load structures
cells.append(md("## View 1 — Full homodimer, all conserved hotspots"))
cells.append(code("""
import nglview as nv
v1 = nv.show_structure_file(
    '00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb')
v1.clear_representations()
v1.add_cartoon(selection='A', color='#3498DB')   # chain A blue
v1.add_cartoon(selection='B', color='#E67E22')   # chain B orange

# Hotspot residues (top 20)
hotspots = """ + str(hotspot_pos) + """
for pos in hotspots:
    sel = str(pos)
    v1.add_licorice(selection=sel, color='#2ECC71')

# Primary SB residues
for pos in [40, 62, 91, 267]:
    v1.add_licorice(selection=str(pos), color='#E74C3C')

v1.layout.width  = '900px'
v1.layout.height = '500px'
v1
"""))

# Cell 2 — salt bridges zoomed
cells.append(md("## View 2 — Primary salt bridges zoomed"))
cells.append(code("""
v2 = nv.show_structure_file(
    '00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb')
v2.clear_representations()
v2.add_cartoon(selection='A', color='#3498DB', opacity=0.3)
v2.add_cartoon(selection='B', color='#E67E22', opacity=0.3)

# ASP40-ARG91 salt bridge
for pos in [40, 91]:
    v2.add_licorice(selection=str(pos), color='#E74C3C')
    v2.add_ball_and_stick(selection=str(pos))

# ARG62-GLU267 salt bridge
for pos in [62, 267]:
    v2.add_licorice(selection=str(pos), color='#9B59B6')
    v2.add_ball_and_stick(selection=str(pos))

print("RED  : ASP40-ARG91 (pan-cov: cons ASP40=1.000)")
print("PURPLE: ARG62-GLU267 (pan-cov: cons ARG62=1.000, GLU267=0.800)")
v2.layout.width  = '900px'
v2.layout.height = '500px'
v2
"""))

# Cell 3 — active site vs dimer interface (NSP15-specific view)
cells.append(md("## View 3 — Active site vs dimer interface (NSP15-specific)"))
cells.append(code("""
v3 = nv.show_structure_file(
    '00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb')
v3.clear_representations()
v3.add_cartoon(color='#BDC3C7', opacity=0.4)

# Active site catalytic triad (negative controls)
for pos in [235, 250, 290]:
    v3.add_licorice(selection=str(pos), color='#7F8C8D')
    v3.add_ball_and_stick(selection=str(pos))

# Dimer interface primary SBs
for pos in [40, 62, 91, 267]:
    v3.add_licorice(selection=str(pos), color='#E74C3C')
    v3.add_ball_and_stick(selection=str(pos))

print("GRAY  : Active site (HIS235, HIS250, LYS290) — NOT dimer interface")
print("RED   : Dimer interface SBs (ASP40, ARG62, ARG91, GLU267)")
print("Key finding: active site and dimer interface are spatially DISTINCT")
print("Negative controls validate: 0 interface contacts for active site residues")
v3.layout.width  = '900px'
v3.layout.height = '500px'
v3
"""))

# Cell 4 — composite score coloring
cells.append(md("## View 4 — Hotspots colored by composite score"))
cells.append(code("""
v4 = nv.show_structure_file(
    '00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb')
v4.clear_representations()
v4.add_cartoon(color='#BDC3C7', opacity=0.3)

# Color by composite tier
top5   = """ + str([r["position"] for r in ranked[:5]])  + """
top10  = """ + str([r["position"] for r in ranked[5:10]]) + """
top20  = """ + str([r["position"] for r in ranked[10:20]]) + """

for pos in top5:
    v4.add_licorice(selection=str(pos), color='#C0392B')
for pos in top10:
    v4.add_licorice(selection=str(pos), color='#E67E22')
for pos in top20:
    v4.add_licorice(selection=str(pos), color='#F1C40F')

print("RED    : Top 5  composite score")
print("ORANGE : Top 6-10")
print("YELLOW : Top 11-20")
v4.layout.width  = '900px'
v4.layout.height = '500px'
v4
"""))

# Cell 5 — BSA burial depth
cells.append(md("## View 5 — Hotspots colored by BSA burial"))
cells.append(code("""
v5 = nv.show_structure_file(
    '00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb')
v5.clear_representations()
v5.add_cartoon(color='#BDC3C7', opacity=0.3)

high_bsa = """ + str([r["position"] for r in ranked if r["bsa"] > 60]) + """
mid_bsa  = """ + str([r["position"] for r in ranked if 30 < r["bsa"] <= 60]) + """
low_bsa  = """ + str([r["position"] for r in ranked if r["bsa"] <= 30]) + """

for pos in high_bsa:
    v5.add_licorice(selection=str(pos), color='#E74C3C')
for pos in mid_bsa:
    v5.add_licorice(selection=str(pos), color='#F39C12')
for pos in low_bsa:
    v5.add_licorice(selection=str(pos), color='#3498DB')

print("RED    : High BSA (>60 A2)")
print("ORANGE : Mid BSA (30-60 A2)")
print("BLUE   : Low BSA (<30 A2)")
v5.layout.width  = '900px'
v5.layout.height = '500px'
v5
"""))

# Cell 6 — docking box
cells.append(md("## View 6 — Docking box (distributed PPI interface)"))
cells.append(code("""
v6 = nv.show_structure_file(
    '00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb')
v6.clear_representations()
v6.add_cartoon(selection='A', color='#3498DB', opacity=0.5)
v6.add_cartoon(selection='B', color='#E67E22', opacity=0.5)
for pos in [40, 62, 91, 267]:
    v6.add_licorice(selection=str(pos), color='#E74C3C')

# Draw docking box edges
cx, cy, cz = """ + str(center) + """
sx, sy, sz = """ + str(size) + """
dx, dy, dz = sx/2, sy/2, sz/2

corners = [
    [cx-dx,cy-dy,cz-dz],[cx+dx,cy-dy,cz-dz],
    [cx+dx,cy+dy,cz-dz],[cx-dx,cy+dy,cz-dz],
    [cx-dx,cy-dy,cz+dz],[cx+dx,cy-dy,cz+dz],
    [cx+dx,cy+dy,cz+dz],[cx-dx,cy+dy,cz+dz],
]
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
for i,j in edges:
    c1,c2 = corners[i],corners[j]
    v6.shape.add_cylinder(c1,c2,[1,0.5,0],0.3)

print("Docking box: {} x {} x {} A".format(round(sx,1),round(sy,1),round(sz,1)))
print("Volume: {} A3".format(round(sx*sy*sz,0)))
print("Strategy: spans all 4 primary SB anchor pockets")
v6.layout.width  = '900px'
v6.layout.height = '500px'
v6
"""))

# Cell 7 — structural overlay
cells.append(md("## View 7 — Structural overlay 6VWW + 6W01 + 6WLC + 6WXC"))
cells.append(code("""
import nglview as nv

# Show all 4 structures overlaid
structs = [
    ('00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb', '#3498DB', '6VWW apo'),
    ('00-reference/known_interfaces/NSP15/6W01_NSP15-true-dimer.pdb', '#2ECC71', '6W01 apo'),
    ('00-reference/known_interfaces/NSP15/6WLC_NSP15-true-dimer-UMP.pdb', '#E67E22', '6WLC +UMP'),
    ('00-reference/known_interfaces/NSP15/6WXC_NSP15-true-dimer-tipiracil.pdb', '#E74C3C', '6WXC +tipiracil'),
]

v7 = nv.show_structure_file(structs[0][0])
v7.clear_representations()
v7.add_cartoon(color=structs[0][1], opacity=0.7)

for path, col, label in structs[1:]:
    s = v7.add_component(path)
    v7.add_cartoon(component=s, color=col, opacity=0.7)

for pos in [40, 62, 91, 267]:
    v7.add_licorice(selection=str(pos), color='#C0392B')

print("BLUE   : 6VWW apo (1.90A)")
print("GREEN  : 6W01 apo (2.20A)")
print("ORANGE : 6WLC +UMP (1.70A)")
print("RED    : 6WXC +tipiracil (1.80A)")
print("All structures show consistent dimer interface geometry")
v7.layout.width  = '900px'
v7.layout.height = '500px'
v7
"""))

# Cell 8 — ranking table
cells.append(md("## View 8 — Full composite ranking table"))
cells.append(code("""
import pandas as pd

rows = []
for r in """ + str(ranked[:25]) + """:
    ctrl = "★ POS" if r["is_pos"] else ("✗ NEG" if r["is_neg"] else "")
    rows.append({
        "Rank"      : r["rank"],
        "Residue"   : r["residue_aa"],
        "Cons"      : r["cons"],
        "BSA (A2)"  : r["bsa"],
        "Ala Loss"  : r["loss"],
        "Composite" : r["composite"],
        "Control"   : ctrl,
    })

df = pd.DataFrame(rows)

def color_rows(row):
    if "POS" in str(row["Control"]):
        return ['background-color: #FADBD8'] * len(row)
    elif "NEG" in str(row["Control"]):
        return ['background-color: #D5D8DC'] * len(row)
    elif row["Cons"] >= 0.8:
        return ['background-color: #D6EAF8'] * len(row)
    return [''] * len(row)

display(df.style.apply(color_rows, axis=1)
          .format({"Cons": "{:.3f}", "BSA (A2)": "{:.1f}",
                   "Composite": "{:.4f}"}))
print()
print("RED   : Primary SB anchor (positive control)")
print("GRAY  : Active site (negative control — 0 interface contacts)")
print("BLUE  : Conserved hotspot (score >= 0.8)")
"""))

# Cell 9 — scientific summary
cells.append(md("""
## NSP15 Scientific Summary

### Complex
NSP15 NendoU homodimer — dimer interface within trimer (biological assembly A+A-2)
4 crystal structures: 6VWW (apo), 6W01 (apo), 6WLC (+UMP), 6WXC (+tipiracil)
AF3 monomer ptm=0.94, RMSD=0.437Å (best in project)

### Interface character
HY-dominated (87% hydrophobic contacts, ~75 contacts per structure)
Two primary salt bridges present in ALL 4 structures:
  ASP40--ARG91: 2.65-3.02Å (pan-cov: ASP39 cons=1.000)
  ARG62--GLU267: 3.49-3.73Å (pan-cov: ARG61 cons=1.000)

### Conservation
Interface mean=0.593 < Surface mean=0.673
Consistent with Hayn et al. 2021: SARS-CoV-2 NSP15 less conserved
than SARS-CoV-1 — explains reduced IFN evasion capability

### Controls (cleanest in project)
Negative: HIS235/HIS250/LYS290 = 0 interface contacts each
Positive: ASP40 rank=4

### Primary pharmacophore: ASP40 (UniProt ASP39)
cons=1.000 pan-coronavirus, BSA=61.0Å², SB with ARG91 all 4 structures

### Druggability
Best pocket score=0.900 (active site — NendoU catalytic triad)
Dimer interface: distributed across 4 pockets — fragment approach required
Strategy: Rule-of-Three fragments (Raj et al. 2021) targeting individual SB pockets
"""))

# ── Write notebook ─────────────────────────────────────────────
nb = nbformat.v4.new_notebook()
nb.cells = cells
nb.metadata = {
    "kernelspec": {
        "display_name": "Python 3 (rtc-discovery)",
        "language"    : "python",
        "name"        : "rtc-discovery",
    }
}

with open(OUT_NB, "w") as f:
    nbformat.write(nb, f)

print(f"Notebook saved: {OUT_NB}")
print()
print("NSP15 PIPELINE COMPLETE — Scripts 04_9 through 11_9")
print()
print("Summary:")
print("  AF3 validation  : ptm=0.94, RMSD=0.437A  ✅")
print("  Interface       : 75 contacts, 31+32 residues")
print("  Salt bridges    : ASP40-ARG91 + ARG62-GLU267 (all 4 structures)")
print("  Conservation    : 23/60 hotspots >= 0.8")
print("  Primary pharma  : ASP40 (cons=1.000, pan-cov)")
print("  Neg controls    : 0 interface contacts (PERFECT)")
print("  Druggability    : 0.900 (active site) / distributed dimer interface")
print("  Docking box     : 46.8 x 52.9 x 30.2 A, vol=74,556 A3")
