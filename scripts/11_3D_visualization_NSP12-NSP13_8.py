#!/usr/bin/env python3
"""
Script 11_3D_visualization_NSP12-NSP13_8.py
Generates NSP12-NSP13_3D_8.ipynb for Jupyter Notebook.
"""

import os, nbformat
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell
from textwrap import dedent

BASE   = os.path.expanduser("~/projects/rtc-pan-coronavirus")
NB_DIR = os.path.join(BASE, "notebooks")
os.makedirs(NB_DIR, exist_ok=True)

cells = []

# ── Cell 0: Header ─────────────────────────────────────────────────────────
cells.append(new_markdown_cell(dedent("""\
# NSP12–NSP13 Interface — 3D Visualization
## Pan-Coronavirus RTC Inhibitor Discovery Pipeline
### Final complex (_8) — Completes all 8 PPI targets

| Property | Value |
|---|---|
| Interface size | 5 NSP12 + 7 NSP13 residues (smallest in project) |
| AF3 status | iptm=0.20 FAIL — transient/context-dependent |
| Primary pharmacophore | TYR93(NSP13)–MET902(NSP12) aromatic-hydrophobic |
| Salt bridge | ASP901(NSP12)–LYS94(NSP13) 3.95 Å |
| LYS94 conservation | 1.000 — K in all 5 coronaviruses |
| MET902 selectivity | SARS-CoV-1/2 only (M→S in MERS/229E/NL63) |
| Drug strategy A | pan-cov: LYS94–ASP901/E charged SB |
| Drug strategy B | SARS-selective: MET902 hydrophobic groove |
""")))

# ── Cell 1: Setup ──────────────────────────────────────────────────────────
cells.append(new_code_cell(dedent("""\
import nglview as nv
import json, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from IPython.display import display

BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = os.path.join(BASE, "00-reference", "pdb_structures")
AF3_DIR = os.path.join(BASE, "01-alphafold3", "NSP12-NSP13")
VAL_DIR = os.path.join(BASE, "02-validation", "NSP12-NSP13")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(RES_DIR, exist_ok=True)

PDB_7RDY = os.path.join(PDB_DIR, "7RDY.pdb")
PDB_6XEZ = os.path.join(PDB_DIR, "6XEZ.pdb")
PDB_7CXM = os.path.join(PDB_DIR, "7CXM.pdb")

NSP12_HOTSPOTS = [900, 901, 902, 903, 904]
NSP13_HOTSPOTS = [90, 91, 92, 93, 94, 95, 96]
PRIMARY_NSP12  = [901, 902]
PRIMARY_NSP13  = [93, 94]
SB_NSP12       = [901]
SB_NSP13       = [94]

DOCKING_BOX = {
    "center": [148.232, 152.667, 158.950],
    "size"  : [33.185,  31.428,  37.767],
}

with open(os.path.join(VAL_DIR, "bsa_alascan_NSP12-NSP13_8.json")) as f:
    bsa_data = json.load(f)
df_rank = pd.DataFrame(bsa_data["composite_ranking"]) \\
            .sort_values("composite", ascending=False) \\
            .reset_index(drop=True)
df_rank["rank"] = df_rank.index + 1

print("Setup complete")
for label, path in [("7RDY", PDB_7RDY), ("6XEZ", PDB_6XEZ), ("7CXM", PDB_7CXM)]:
    print(f"  {label}: {os.path.exists(path)}")
print("Top 3 hotspots:")
print(df_rank[["rank","nsp","res_name","res_num","composite"]].head(3).to_string(index=False))
""")))

# ── Cell 2: View 1 — Full heterodimer ─────────────────────────────────────
cells.append(new_markdown_cell("## View 1 — Full heterodimer: all 12 hotspots"))
cells.append(new_code_cell(dedent("""\
view1 = nv.show_file(PDB_7RDY)
view1.clear_representations()
view1.add_representation("cartoon", selection=":A", color="steelblue",    opacity=0.7)
view1.add_representation("cartoon", selection=":E", color="mediumorchid", opacity=0.7)

for rn in NSP12_HOTSPOTS:
    color = "#d62728" if rn in SB_NSP12 else "#ff7f0e" if rn in PRIMARY_NSP12 else "#1f77b4"
    view1.add_representation("ball+stick", selection=":A and " + str(rn),
                             color=color, radius=0.35)

for rn in NSP13_HOTSPOTS:
    color = "#d62728" if rn in SB_NSP13 else "#9467bd" if rn in PRIMARY_NSP13 else "#e377c2"
    view1.add_representation("ball+stick", selection=":E and " + str(rn),
                             color=color, radius=0.35)

view1.layout.width  = "900px"
view1.layout.height = "500px"
display(view1)
print("Red=SB anchor | Orange=hydrophobic primary | Purple=TYR93 #1 | Pink=NSP13 other")
""")))

# ── Cell 3: View 2 — TYR93-MET902 primary pharmacophore ───────────────────
cells.append(new_markdown_cell(
    "## View 2 ★ TYR93(NSP13)–MET902(NSP12): Primary pharmacophore\n"
    "TYR93: #1 composite (0.9486), AlaLoss=29, pan-cov conserved (1.000)  \n"
    "MET902: #2 composite (0.7785), BSA=73.6 Å², SARS-CoV-1/2 selective"
))
cells.append(new_code_cell(dedent("""\
view2 = nv.show_file(PDB_7RDY)
view2.clear_representations()
view2.add_representation("cartoon", selection=":A", color="steelblue",    opacity=0.25)
view2.add_representation("cartoon", selection=":E", color="mediumorchid", opacity=0.25)

# Primary pair — large radius
view2.add_representation("ball+stick", selection=":A and 902", color="#ff7f0e", radius=0.6)
view2.add_representation("ball+stick", selection=":E and 93",  color="#9467bd", radius=0.6)

# Supporting context — thin licorice
for rn in [900, 901, 903, 904]:
    view2.add_representation("licorice", selection=":A and " + str(rn),
                             color="#1f77b4", radius=0.15)
for rn in [90, 91, 92, 94, 95, 96]:
    view2.add_representation("licorice", selection=":E and " + str(rn),
                             color="#e377c2", radius=0.15)

view2.layout.width  = "900px"
view2.layout.height = "500px"
display(view2)
print("MET902 (orange, NSP12) + TYR93 (purple, NSP13) = primary pharmacophore")
""")))

# ── Cell 4: View 3 — ASP901-LYS94 salt bridge ─────────────────────────────
cells.append(new_markdown_cell(
    "## View 3 ★ ASP901–LYS94 salt bridge: Pan-cov charged anchor\n"
    "Distance: 3.95 Å [7RDY] | 4.12 Å [7CXM] | absent in 6XEZ (resolution artifact)  \n"
    "LYS94 conservation = 1.000 (K in all 5 species)"
))
cells.append(new_code_cell(dedent("""\
view3 = nv.show_file(PDB_7RDY)
view3.clear_representations()
view3.add_representation("cartoon", selection=":A", color="steelblue",    opacity=0.25)
view3.add_representation("cartoon", selection=":E", color="mediumorchid", opacity=0.25)

# SB pair
view3.add_representation("ball+stick", selection=":A and 901", color="#d62728", radius=0.6)
view3.add_representation("ball+stick", selection=":E and 94",  color="#d62728", radius=0.6)

# Context
view3.add_representation("licorice", selection=":A and 902", color="#ff7f0e", radius=0.2)
view3.add_representation("licorice", selection=":E and 93",  color="#9467bd", radius=0.2)

view3.layout.width  = "900px"
view3.layout.height = "500px"
display(view3)
print("ASP901 (red, NSP12) -- LYS94 (red, NSP13): 3.95 A [7RDY]")
print("MET902 (orange) + TYR93 (purple) shown for spatial context")
""")))

# ── Cell 5: View 4 — Composite score coloring ─────────────────────────────
cells.append(new_markdown_cell("## View 4 — All hotspots colored by composite score"))
cells.append(new_code_cell(dedent("""\
view4 = nv.show_file(PDB_7RDY)
view4.clear_representations()
view4.add_representation("cartoon", selection=":A", color="lightgray", opacity=0.35)
view4.add_representation("cartoon", selection=":E", color="lightgray", opacity=0.35)

def score_to_color(score, max_score=0.9486):
    r = score / max_score
    if r >= 0.70: return "#d62728"
    if r >= 0.50: return "#ff7f0e"
    if r >= 0.30: return "#2ca02c"
    return "#1f77b4"

for _, row in df_rank.iterrows():
    chain = ":A" if row["nsp"] == "NSP12" else ":E"
    color = score_to_color(row["composite"])
    view4.add_representation("ball+stick",
                             selection=chain + " and " + str(row["res_num"]),
                             color=color, radius=0.4)

view4.layout.width  = "900px"
view4.layout.height = "500px"
display(view4)
print("Red>=0.70 | Orange 0.50-0.70 | Green 0.30-0.50 | Blue<0.30")
print("Top tier (red): TYR93, MET902")
print("Mid tier (orange): SER904, LEU92, ASN95, ASP901")
""")))

# ── Cell 6: View 5 — BSA coloring ─────────────────────────────────────────
cells.append(new_markdown_cell("## View 5 — Hotspots colored by BSA burial depth"))
cells.append(new_code_cell(dedent("""\
view5 = nv.show_file(PDB_7RDY)
view5.clear_representations()
view5.add_representation("cartoon", selection=":A", color="lightgray", opacity=0.35)
view5.add_representation("cartoon", selection=":E", color="lightgray", opacity=0.35)

bsa_all = {
    "A": {900: 5.4,  901: 60.2, 902: 73.6, 903: 26.8, 904: 51.4},
    "E": {90:  22.5, 91:  1.8,  92:  44.7, 93:  64.1, 94:  19.4,
          95:  34.7, 96:  19.0},
}
max_bsa = 73.6

def bsa_to_color(bsa):
    r = bsa / max_bsa
    if r >= 0.70: return "#d62728"
    if r >= 0.50: return "#ff7f0e"
    if r >= 0.25: return "#2ca02c"
    return "#1f77b4"

for chain_id, rn_dict in bsa_all.items():
    for rn, bsa_val in rn_dict.items():
        view5.add_representation("ball+stick",
                                 selection=":" + chain_id + " and " + str(rn),
                                 color=bsa_to_color(bsa_val), radius=0.4)

view5.layout.width  = "900px"
view5.layout.height = "500px"
display(view5)
print("Red>=52.5 A2 | Orange 36.8-52.5 | Green 18.4-36.8 | Blue<18.4")
print("MET902=73.6 > TYR93=64.1 > ASP901=60.2 > SER904=51.4 > LEU92=44.7")
""")))

# ── Cell 7: View 6 — Docking box ──────────────────────────────────────────
cells.append(new_markdown_cell(
    "## View 6 — Docking box (cylinder wireframe)\n"
    "Center: (148.232, 152.667, 158.950) | Size: 33.2 × 31.4 × 37.8 Å | Volume: 39,389 Å³"
))
cells.append(new_code_cell(dedent("""\
view6 = nv.show_file(PDB_7RDY)
view6.clear_representations()
view6.add_representation("cartoon", selection=":A", color="steelblue",    opacity=0.5)
view6.add_representation("cartoon", selection=":E", color="mediumorchid", opacity=0.5)

for rn in NSP12_HOTSPOTS:
    color = "#d62728" if rn in SB_NSP12 else "#ff7f0e"
    view6.add_representation("ball+stick", selection=":A and " + str(rn),
                             color=color, radius=0.3)
for rn in NSP13_HOTSPOTS:
    color = "#d62728" if rn in SB_NSP13 else "#9467bd" if rn in PRIMARY_NSP13 else "#e377c2"
    view6.add_representation("ball+stick", selection=":E and " + str(rn),
                             color=color, radius=0.3)

cx, cy, cz = DOCKING_BOX["center"]
sx, sy, sz = [s/2 for s in DOCKING_BOX["size"]]
corners = [
    [cx-sx, cy-sy, cz-sz], [cx+sx, cy-sy, cz-sz],
    [cx+sx, cy+sy, cz-sz], [cx-sx, cy+sy, cz-sz],
    [cx-sx, cy-sy, cz+sz], [cx+sx, cy-sy, cz+sz],
    [cx+sx, cy+sy, cz+sz], [cx-sx, cy+sy, cz+sz],
]
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
shape = view6.shape
for i, j in edges:
    shape.add_cylinder(corners[i], corners[j], [0.2, 0.8, 0.2], 0.3)

view6.camera = "orthographic"
view6.layout.width  = "900px"
view6.layout.height = "500px"
display(view6)
print("Docking box: center =", DOCKING_BOX["center"])
print("             size   =", DOCKING_BOX["size"])
print("             volume = 39,389 A3")
""")))

# ── Cell 8: View 7 — Structural overlay ───────────────────────────────────
cells.append(new_markdown_cell(
    "## View 7 — Structural overlay: 7RDY + 6XEZ + 7CXM\n"
    "Blue/orchid = 7RDY (3.10 Å, primary) | Green = 6XEZ (3.50 Å) | Orange = 7CXM (3.20 Å)"
))
cells.append(new_code_cell(dedent("""\
view7 = nv.show_file(PDB_7RDY)
view7.clear_representations()

view7.add_structure(nv.FileStructure(PDB_6XEZ))
view7.add_structure(nv.FileStructure(PDB_7CXM))

# 7RDY — primary
view7.add_representation("cartoon", selection=":A", color="steelblue",    opacity=0.8, component=0)
view7.add_representation("cartoon", selection=":E", color="mediumorchid", opacity=0.8, component=0)

# 6XEZ — secondary (green)
view7.add_representation("cartoon", selection=":A", color="#2ca02c", opacity=0.35, component=1)
view7.add_representation("cartoon", selection=":E", color="#2ca02c", opacity=0.35, component=1)

# 7CXM — tertiary (orange)
view7.add_representation("cartoon", selection=":A", color="#ff7f0e", opacity=0.35, component=2)
view7.add_representation("cartoon", selection=":E", color="#ff7f0e", opacity=0.35, component=2)

# Primary hotspots on 7RDY
for rn, color in [(901,"#d62728"),(902,"#ff7f0e")]:
    view7.add_representation("ball+stick", selection=":A and " + str(rn),
                             color=color, radius=0.4, component=0)
for rn, color in [(93,"#9467bd"),(94,"#d62728")]:
    view7.add_representation("ball+stick", selection=":E and " + str(rn),
                             color=color, radius=0.4, component=0)

view7.layout.width  = "900px"
view7.layout.height = "500px"
display(view7)
print("3 crystal structures overlaid — confirms consistent NSP12 C-tail / NSP13 N-term interface")
""")))

# ── Cell 9: View 8 — Ranking table ────────────────────────────────────────
cells.append(new_markdown_cell("## View 8 — Full composite ranking table"))
cells.append(new_code_cell(dedent("""\
fig, ax = plt.subplots(figsize=(14, 5))
ax.axis("off")
fig.patch.set_facecolor("white")

ROLE_MAP = {
    ("NSP13", 93): "Aromatic core #1 ★",
    ("NSP12", 902): "Hydrophobic ★ SARS-sel",
    ("NSP12", 901): "SB anchor",
    ("NSP13", 94):  "SB pan-cov ★",
}

table_data = []
for _, row in df_rank.iterrows():
    role = ROLE_MAP.get((row["nsp"], row["res_num"]), "")
    table_data.append([
        int(row["rank"]),
        row["nsp"],
        row["res_name"] + str(row["res_num"]),
        f"{row['bsa']:.1f}",
        int(row["ala_loss"]),
        f"{row['conservation']:.3f}",
        f"{row['composite']:.4f}",
        role,
    ])

cols = ["Rank", "NSP", "Residue", "BSA (A2)", "AlaLoss", "Cons", "Composite", "Role"]
tbl  = ax.table(cellText=table_data, colLabels=cols, loc="center", cellLoc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.7)

# Header style
for j in range(len(cols)):
    tbl[(0,j)].set_facecolor("#2c3e50")
    tbl[(0,j)].set_text_props(color="white", fontweight="bold")

# Highlight primary rows
primary_idx = [i+1 for i, row in enumerate(df_rank.itertuples())
               if (row.nsp == "NSP12" and row.res_num in {901, 902}) or
                  (row.nsp == "NSP13" and row.res_num in {93, 94})]
for ri in primary_idx:
    for j in range(len(cols)):
        tbl[(ri,j)].set_facecolor("#ffeeba")

ax.set_title(
    "NSP12-NSP13 Interface: Full Composite Hotspot Ranking\\n"
    "Primary: TYR93-MET902 aromatic-hydrophobic core + ASP901-LYS94 SB",
    fontsize=11, fontweight="bold", pad=12, y=0.95)

out_table = os.path.join(RES_DIR, "Fig8_NSP12-NSP13_ranking_table_8.png")
fig.savefig(out_table, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print("Ranking table saved:", out_table)
display(df_rank[["rank","nsp","res_name","res_num","bsa","ala_loss","conservation","composite"]])
""")))

# ── Cell 10: Project completion summary ───────────────────────────────────
cells.append(new_markdown_cell(dedent("""\
## 🎉 PROJECT COMPLETE — All 8 PPI Interfaces Analyzed

| # | Complex | Primary Pharmacophore | Druggability | Selectivity |
|---|---------|----------------------|--------------|-------------|
| 1 | NSP10-NSP14 | HIS80–ASP126 SB | 0.000 | pan-cov |
| 2 | NSP10-NSP16 | LYS93–ASP106 SB + Zn1 | 0.546 | pan-cov |
| 3 | NSP12-NSP7  | PHE440 aromatic core | 0.961 | pan-cov |
| 4 | NSP12-NSP8  | LYS332–ASP99 SB | 0.874 | pan-cov |
| 5 | NSP9-NSP12  | ARG733 NiRAN domain | 0.895 | pan-cov |
| 6 | NSP7-NSP8   | PHE92 hydrophobic core | 0.531 | pan-cov |
| 7 | NSP13-Helicase | ILE480 + LYS414 dual SB | 0.001 | SARS-selective |
| **8** | **NSP12-NSP13** | **TYR93–MET902 + ASP901–LYS94** | **0.000** | **dual** |

---
**Next phase:** Virtual screening with AutoDock Vina / VirtualFlow  
**Repository:** https://github.com/nsolly03/panCov-rtc-discovery
""")))

# ── Write notebook with proper kernel metadata ─────────────────────────────
nb = new_notebook(cells=cells)
nb.metadata.update({
    "kernelspec": {
        "display_name": "Python 3 (rtc-discovery)",
        "language": "python",
        "name": "python3"
    },
    "language_info": {
        "name": "python",
        "version": "3.10.0"
    }
})

nb_path = os.path.join(NB_DIR, "NSP12-NSP13_3D_8.ipynb")
with open(nb_path, "w") as f:
    nbformat.write(nb, f)

print(f"Notebook written : {nb_path}")
print(f"Cells            : {len(cells)}")
print("\nTo open:")
print("  jupyter notebook notebooks/NSP12-NSP13_3D_8.ipynb --ip=0.0.0.0 --port=8889 --no-browser")
print("  Then open in Windows browser: http://127.0.0.1:8889/notebooks/NSP12-NSP13_3D_8.ipynb?token=...")
