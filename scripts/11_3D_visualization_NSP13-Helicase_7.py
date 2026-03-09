#!/usr/bin/env python3
"""
Script 11_3D_visualization_NSP13-Helicase_7.py

Generates NSP13-Helicase_3D_7.ipynb — interactive 3D visualization
notebook using nglview.

8 views:
  1. Full homodimer — all 17 conserved hotspots
  2. Dual salt bridge zoomed — LYS414-ASP580 + LYS414-ASP583
  3. ILE480 primary hydrophobic pharmacophore ★
  4. Hotspots colored by composite score
  5. Hotspots colored by BSA burial depth
  6. Docking box visualization
  7. Structural overlay — 7NIO + 6XEZ (two crystal contexts)
  8. Full composite ranking table

Output: notebooks/NSP13-Helicase_3D_7.ipynb
"""

import os, json
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

# ── Paths ──────────────────────────────────────────────────────────────────
BASE       = os.path.expanduser("~/projects/rtc-pan-coronavirus")
NB_DIR     = f"{BASE}/notebooks"
VAL_DIR    = f"{BASE}/02-validation/NSP13-Helicase"
os.makedirs(NB_DIR, exist_ok=True)

NB_PATH    = f"{NB_DIR}/NSP13-Helicase_3D_7.ipynb"

# ── Load ranking data ──────────────────────────────────────────────────────
with open(f"{VAL_DIR}/bsa_alascan_NSP13-Helicase_7.json") as f:
    ranking_data = json.load(f)

top3 = ranking_data["top3_pharmacophores"]

# ── Notebook cells ─────────────────────────────────────────────────────────
cells = []

# ── Cell 0: Title ─────────────────────────────────────────────────────────
cells.append(new_markdown_cell("""# NSP13-Helicase — 3D Interface Visualization
## Pan-Coronavirus RTC Inhibitor Discovery | University of Liège GIGA-VIN Lab

**Complex:** NSP13 homodimer (SARS-CoV-1/2 selective target)
**Primary structures:** 7NIO (2.80 Å), 6XEZ (3.50 Å)
**Key pharmacophore:** ILE480 hydrophobic anchor + LYS414 dual SB donor
**Selectivity:** SARS-CoV-1/2 — ASP580 D→A in MERS, D→T in 229E/NL63

### 8 Interactive Views
1. Full homodimer — all hotspots
2. Dual salt bridge zoomed (LYS414→ASP580 + ASP583) ★
3. ILE480 primary hydrophobic pharmacophore ★
4. Hotspots colored by composite score
5. Hotspots colored by BSA burial depth
6. Docking box
7. Structural overlay (7NIO + 6XEZ)
8. Full composite ranking table
"""))

# ── Cell 1: Setup ─────────────────────────────────────────────────────────
cells.append(new_code_cell("""
import nglview as nv
import numpy as np
import pandas as pd
from IPython.display import display, HTML
import os, json

BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = f"{BASE}/00-reference/pdb_structures"
AF3_DIR = f"{BASE}/01-alphafold3/NSP13-Helicase"
VAL_DIR = f"{BASE}/02-validation/NSP13-Helicase"

PDB_7NIO = f"{PDB_DIR}/7NIO.pdb"
PDB_6XEZ = f"{PDB_DIR}/6XEZ.pdb"
AF3_PDB  = f"{AF3_DIR}/NSP13_Helicase_best_model.pdb"

# Hotspot definitions
ALL_HOTSPOTS     = [116,413,414,415,477,478,479,480,481,482,
                    551,552,553,579,580,583,584]
PRIMARY_HOTSPOTS = [414,477,480,482,579,580,583,584]
SB_RESIDUES      = [414,580,583]
HYDROPHOBIC_ANCHORS = [480,479,482]  # ILE480 top, VAL479, HIS482

DOCKING_BOX = {
    "center": [-30.151, 14.648, -9.240],
    "size"  : [37.237, 31.585, 33.862]
}

# Load ranking
with open(f"{VAL_DIR}/bsa_alascan_NSP13-Helicase_7.json") as f:
    rank_data = json.load(f)

scores = {r["res_num"]: r["composite_score"]
          for r in rank_data["composite_ranking"]}
bsa    = {r["res_num"]: r["bsa_A2"]
          for r in rank_data["composite_ranking"]}

max_score = max(scores.values()) if scores else 1.0
max_bsa   = max(bsa.values())    if bsa    else 1.0

def sel(chain, residues):
    return " or ".join(f"({chain} and {r})" for r in residues)

def score_to_color(score, max_val):
    t = score / max_val if max_val > 0 else 0
    r = int(255 * t)
    b = int(255 * (1 - t))
    return f"#{r:02x}00{b:02x}"

print("Setup complete. nglview ready.")
print(f"Hotspots: {len(ALL_HOTSPOTS)} total, {len(PRIMARY_HOTSPOTS)} primary")
"""))

# ── Cell 2: View 1 — Full homodimer ───────────────────────────────────────
cells.append(new_markdown_cell("## View 1 — Full NSP13 Homodimer | All Hotspots"))
cells.append(new_code_cell("""
view1 = nv.show_file(PDB_7NIO)
view1.clear_representations()

# Chain A — cartoon blue
view1.add_representation("cartoon", selection=":A", color="skyblue",   opacity=0.7)
# Chain E — cartoon teal
view1.add_representation("cartoon", selection=":E", color="lightcoral",opacity=0.7)

# All hotspots — licorice orange
for rn in ALL_HOTSPOTS:
    color = "#d62728" if rn in SB_RESIDUES else ("#ff7f0e" if rn in PRIMARY_HOTSPOTS else "#1f77b4")
    view1.add_representation("licorice", selection=f":A and {rn}",
                              color=color, radius=0.25)

view1.camera = "orthographic"
view1.layout.width  = "900px"
view1.layout.height = "500px"
display(view1)
print("Blue=chain A, Coral=chain E")
print("Red=SB anchors, Orange=primary hotspots, Blue=secondary hotspots")
"""))

# ── Cell 3: View 2 — Dual salt bridge ─────────────────────────────────────
cells.append(new_markdown_cell("""## View 2 — Dual Salt Bridge ★
**LYS414(A) → ASP580(E) [4.49 Å] + ASP583(E) [4.47 Å]**
Unique in this project — single donor bridges two acceptors simultaneously"""))
cells.append(new_code_cell("""
view2 = nv.show_file(PDB_7NIO)
view2.clear_representations()

# Cartoon context
view2.add_representation("cartoon", selection=":A", color="lightgrey", opacity=0.3)
view2.add_representation("cartoon", selection=":E", color="lightgrey", opacity=0.3)

# SB residues — ball+stick prominent
view2.add_representation("ball+stick", selection=":A and 414",
                         color="#e31a1c", radius=0.4)   # LYS414 red
view2.add_representation("ball+stick", selection=":E and 580",
                         color="#1f78b4", radius=0.4)   # ASP580 blue
view2.add_representation("ball+stick", selection=":E and 583",
                         color="#33a02c", radius=0.4)   # ASP583 green

# Neighbourhood context
for rn in [413,415,477,478,479,480,482,579,584]:
    view2.add_representation("licorice", selection=f":A and {rn}",
                              color="#ff7f0e", radius=0.2, opacity=0.7)

view2.camera = "orthographic"
view2.layout.width  = "900px"
view2.layout.height = "500px"
display(view2)
print("RED   = LYS414 (dual SB donor, chain A)")
print("BLUE  = ASP580 (SB acceptor, chain E) — D->A in MERS")
print("GREEN = ASP583 (SB acceptor, chain E)")
"""))

# ── Cell 4: View 3 — ILE480 hydrophobic pharmacophore ─────────────────────
cells.append(new_markdown_cell("""## View 3 — ILE480 Primary Hydrophobic Pharmacophore ★
**ILE480: rank #1, score=0.4823, BSA=59.2 Å², 50 hydrophobic contacts lost**
Adjacent hydrophobic cluster: VAL479, HIS482"""))
cells.append(new_code_cell("""
view3 = nv.show_file(PDB_7NIO)
view3.clear_representations()

view3.add_representation("cartoon", selection=":A", color="lightgrey", opacity=0.25)
view3.add_representation("cartoon", selection=":E", color="lightgrey", opacity=0.25)

# ILE480 — gold prominent
view3.add_representation("ball+stick", selection=":A and 480",
                         color="gold", radius=0.5)
# VAL479 + HIS482 — orange
view3.add_representation("ball+stick", selection=":A and 479",
                         color="#ff7f0e", radius=0.35)
view3.add_representation("ball+stick", selection=":A and 482",
                         color="#ff7f0e", radius=0.35)

# LYS414 for context
view3.add_representation("ball+stick", selection=":A and 414",
                         color="#d62728", radius=0.35)

# Chain E interface residues in region
for rn in [477,478,479,480,481,482,413,414,415]:
    view3.add_representation("licorice", selection=f":E and {rn}",
                              color="#1f77b4", radius=0.2, opacity=0.6)

view3.camera = "orthographic"
view3.layout.width  = "900px"
view3.layout.height = "500px"
display(view3)
print("GOLD  = ILE480 (rank #1, hydrophobic anchor, 50 contacts)")
print("ORANGE = VAL479 + HIS482 (hydrophobic neighbourhood)")
print("RED   = LYS414 (rank #2, dual SB donor)")
print("BLUE  = chain E partner residues")
"""))

# ── Cell 5: View 4 — Hotspots by composite score ──────────────────────────
cells.append(new_markdown_cell("## View 4 — Hotspots Colored by Composite Score\nRed = highest score → Blue = lowest score"))
cells.append(new_code_cell("""
view4 = nv.show_file(PDB_7NIO)
view4.clear_representations()

view4.add_representation("cartoon", selection=":A", color="lightgrey", opacity=0.25)
view4.add_representation("cartoon", selection=":E", color="lightgrey", opacity=0.25)

for rn in ALL_HOTSPOTS:
    sc    = scores.get(rn, 0.0)
    color = score_to_color(sc, max_score)
    view4.add_representation("licorice", selection=f":A and {rn}",
                              color=color, radius=0.3)

view4.camera = "orthographic"
view4.layout.width  = "900px"
view4.layout.height = "500px"
display(view4)
print("Color scale: Red=highest composite score -> Blue=lowest")
for rn in sorted(scores, key=lambda x: -scores[x])[:5]:
    print(f"  {rn}: score={scores[rn]:.4f}")
"""))

# ── Cell 6: View 5 — Hotspots by BSA ──────────────────────────────────────
cells.append(new_markdown_cell("## View 5 — Hotspots Colored by BSA Burial Depth\nRed = most buried → Blue = least buried"))
cells.append(new_code_cell("""
view5 = nv.show_file(PDB_7NIO)
view5.clear_representations()

view5.add_representation("cartoon", selection=":A", color="lightgrey", opacity=0.25)
view5.add_representation("cartoon", selection=":E", color="lightgrey", opacity=0.25)

for rn in ALL_HOTSPOTS:
    b     = bsa.get(rn, 0.0)
    color = score_to_color(b, max_bsa)
    view5.add_representation("licorice", selection=f":A and {rn}",
                              color=color, radius=0.3)

view5.camera = "orthographic"
view5.layout.width  = "900px"
view5.layout.height = "500px"
display(view5)
print("Color scale: Red=most buried -> Blue=least buried")
for rn in sorted(bsa, key=lambda x: -bsa[x])[:5]:
    print(f"  {rn}: BSA={bsa[rn]:.1f} A2")
"""))

# ── Cell 7: View 6 — Docking box ──────────────────────────────────────────
cells.append(new_markdown_cell("""## View 6 — Docking Box
Center: (-30.151, 14.648, -9.240) | Size: 37.2 × 31.6 × 33.9 Å | Volume: 39,826 Å³"""))
cells.append(new_code_cell("""
view6 = nv.show_file(PDB_7NIO)
view6.clear_representations()

view6.add_representation("cartoon", selection=":A", color="skyblue",    opacity=0.5)
view6.add_representation("cartoon", selection=":E", color="lightcoral", opacity=0.5)

for rn in PRIMARY_HOTSPOTS:
    color = "#d62728" if rn in SB_RESIDUES else "#ff7f0e"
    view6.add_representation("ball+stick", selection=f":A and {rn}",
                              color=color, radius=0.3)

# Docking box as shape
cx, cy, cz = DOCKING_BOX["center"]
sx, sy, sz = DOCKING_BOX["size"]
shape = view6.shape
shape.add_box([cx, cy, cz], [0,0,0,1], [sx/2, sy/2, sz/2],
              [0.2, 0.8, 0.2], "DockingBox")

view6.camera = "orthographic"
view6.layout.width  = "900px"
view6.layout.height = "500px"
display(view6)
print(f"Docking box: center={DOCKING_BOX['center']}")
print(f"             size  ={DOCKING_BOX['size']}")
print(f"             volume=39,826 A3")
"""))

# ── Cell 8: View 7 — Structural overlay ───────────────────────────────────
cells.append(new_markdown_cell("""## View 7 — Structural Overlay: 7NIO + 6XEZ
Two crystal contexts: pure NSP13 dimer (7NIO) vs NSP13 in RdRp complex (6XEZ)"""))
cells.append(new_code_cell("""
view7 = nv.show_file(PDB_7NIO)
view7.add_component(PDB_6XEZ)
view7.clear_representations()

# 7NIO chains
view7.add_representation("cartoon", selection=":A", color="skyblue",
                         opacity=0.7, component=0)
view7.add_representation("cartoon", selection=":E", color="steelblue",
                         opacity=0.7, component=0)

# 6XEZ NSP13 chains
view7.add_representation("cartoon", selection=":E", color="salmon",
                         opacity=0.6, component=1)
view7.add_representation("cartoon", selection=":F", color="coral",
                         opacity=0.6, component=1)

# Hotspots on 7NIO
for rn in PRIMARY_HOTSPOTS:
    view7.add_representation("licorice", selection=f":A and {rn}",
                              color="#ff7f0e", radius=0.25, component=0)

view7.camera = "orthographic"
view7.layout.width  = "900px"
view7.layout.height = "500px"
display(view7)
print("Blue shades = 7NIO (pure NSP13 dimer, 2.80 A)")
print("Red shades  = 6XEZ (NSP13 in RdRp complex, 3.50 A)")
print("Orange sticks = primary hotspots on 7NIO chain A")
"""))

# ── Cell 9: View 8 — Ranking table ────────────────────────────────────────
cells.append(new_markdown_cell("## View 8 — Full Composite Ranking Table"))
cells.append(new_code_cell("""
df = pd.DataFrame(rank_data["composite_ranking"])
df = df.sort_values("composite_score", ascending=False)

def highlight_row(row):
    if row["is_sb_residue"]:
        return ["background-color: #ffe0e0"] * len(row)
    elif row["is_primary"]:
        return ["background-color: #fff3e0"] * len(row)
    return [""] * len(row)

cols = ["rank","res_name","res_num","res_aa","conservation",
        "contact_count","bsa_A2","total_loss","composite_score",
        "is_primary","is_sb_residue"]

styled = (df[cols].style
          .apply(highlight_row, axis=1)
          .format({
              "conservation"   : "{:.3f}",
              "bsa_A2"         : "{:.1f}",
              "composite_score": "{:.4f}",
          })
          .set_caption("NSP13-Helicase Composite Hotspot Ranking")
          .set_table_styles([
              {"selector":"caption",
               "props":[("font-size","14px"),("font-weight","bold"),
                        ("text-align","left")]},
              {"selector":"th",
               "props":[("background-color","#2c3e50"),
                        ("color","white"),("font-size","11px")]}
          ]))

display(styled)

print("\\nScientific summary:")
print("  #1  ILE480 — hydrophobic anchor (50 contacts, BSA=59.2 A2)")
print("  #2  LYS414 — dual SB donor (23 SB contacts) ★ PRIMARY")
print("  #3  ASP580 — SB acceptor (22 contacts) — D->A in MERS ★")
print("  #4  HIS482 — burial anchor (BSA=90.7 A2)")
print("  #5  ASP583 — secondary SB acceptor ★")
print("\\nDrug design: aromatic/hydrophobic core (ILE480)")
print("  + charged group engaging LYS414-ASP580/583 cluster")
print("  SARS-CoV-1/2 selective compound strategy")
"""))

# ── Write notebook ─────────────────────────────────────────────────────────
nb = new_notebook(cells=cells)
nb.metadata["kernelspec"] = {
    "display_name": "Python 3 (rtc-discovery)",
    "language"    : "python",
    "name"        : "rtc-discovery"
}

with open(NB_PATH, "w") as f:
    nbformat.write(nb, f)

print(f"Notebook written: {NB_PATH}")
print(f"Cells: {len(cells)}")
print(f"\nTo open:")
print(f"  conda activate rtc-discovery")
print(f"  cd ~/projects/rtc-pan-coronavirus")
print(f"  jupyter notebook --port=8889 --no-browser")
print(f"  Then open: notebooks/NSP13-Helicase_3D_7.ipynb")
