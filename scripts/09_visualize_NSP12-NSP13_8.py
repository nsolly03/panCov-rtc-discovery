#!/usr/bin/env python3
"""
Script 09_visualize_NSP12-NSP13_8.py

Publication figures for NSP12-NSP13 interface analysis.

Fig1: Conservation bar chart — NSP12 + NSP13 hotspots
Fig2: Conservation heatmap — residues x species
Fig3: Contact type distribution — 3 crystal structures

Output: results/
  Fig1_NSP12-NSP13_conservation_bars_8.png
  Fig2_NSP12-NSP13_conservation_heatmap_8.png
  Fig3_NSP12-NSP13_contact_types_8.png
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

# ── Paths ──────────────────────────────────────────────────────────────────
BASE        = os.path.expanduser("~/projects/rtc-pan-coronavirus")
VAL_DIR     = f"{BASE}/02-validation/NSP12-NSP13"
RESULTS_DIR = f"{BASE}/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────
with open(f"{VAL_DIR}/conservation_summary_8.json") as f:
    cons_data = json.load(f)
with open(f"{VAL_DIR}/interface_analysis_8.json") as f:
    iface_data = json.load(f)

nsp12_cons = cons_data["nsp12_conservation"]
nsp13_cons = cons_data["nsp13_conservation"]
coronaviruses = cons_data["coronaviruses"]

# ── Color scheme (project-standard) ───────────────────────────────────────
COL_SB      = "#d62728"   # red    — salt bridge / primary
COL_HY      = "#ff7f0e"   # orange — hydrophobic primary
COL_CONS    = "#2ca02c"   # green  — conserved
COL_VAR     = "#9467bd"   # purple — variable
COL_NSP12   = "#1f77b4"   # blue
COL_NSP13   = "#e377c2"   # pink

PRIMARY_NSP12 = {901, 902}
PRIMARY_NSP13 = {94}
SB_NSP12      = {901}
SB_NSP13      = {94}

# ══════════════════════════════════════════════════════════════════════════
# FIG 1 — Conservation bars (dual panel: NSP12 top, NSP13 bottom)
# ══════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                                gridspec_kw={"hspace": 0.55})
fig.patch.set_facecolor("white")

def bar_color(res_num, primary_set, sb_set):
    if res_num in sb_set:     return COL_SB
    if res_num in primary_set: return COL_HY
    return COL_NSP12

# NSP12
labels12 = [f"{d['res_name_3']}{d['pdb_res_num']}" for d in nsp12_cons]
scores12  = [d["conservation"] for d in nsp12_cons]
colors12  = [bar_color(d["pdb_res_num"], PRIMARY_NSP12, SB_NSP12)
             for d in nsp12_cons]

bars = ax1.bar(labels12, scores12, color=colors12, edgecolor="black",
               linewidth=0.8, zorder=3)
ax1.axhline(0.8, color="gray", linestyle="--", linewidth=1.2,
            alpha=0.7, label="Conservation threshold (0.8)", zorder=2)
ax1.set_ylim(0, 1.15)
ax1.set_ylabel("Conservation score", fontsize=11)
ax1.set_title("NSP12 C-terminal tail — Interface hotspot conservation",
              fontsize=12, fontweight="bold", pad=8)
ax1.set_facecolor("#f8f8f8")
ax1.grid(axis="y", alpha=0.4, zorder=1)
for bar, score in zip(bars, scores12):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f"{score:.3f}", ha="center", va="bottom", fontsize=9,
             fontweight="bold")

# NSP13
labels13 = [f"{d['res_name_3']}{d['pdb_res_num']}" for d in nsp13_cons]
scores13  = [d["conservation"] for d in nsp13_cons]
colors13  = [COL_SB if d["pdb_res_num"] in SB_NSP13
             else COL_CONS if d["conservation"] >= 0.8
             else COL_VAR
             for d in nsp13_cons]

bars2 = ax2.bar(labels13, scores13, color=colors13, edgecolor="black",
                linewidth=0.8, zorder=3)
ax2.axhline(0.8, color="gray", linestyle="--", linewidth=1.2,
            alpha=0.7, zorder=2)
ax2.set_ylim(0, 1.15)
ax2.set_ylabel("Conservation score", fontsize=11)
ax2.set_title("NSP13 N-terminal region — Interface hotspot conservation",
              fontsize=12, fontweight="bold", pad=8)
ax2.set_facecolor("#f8f8f8")
ax2.grid(axis="y", alpha=0.4, zorder=1)
for bar, score in zip(bars2, scores13):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f"{score:.3f}", ha="center", va="bottom", fontsize=9,
             fontweight="bold")

# Legend
legend_patches = [
    mpatches.Patch(color=COL_SB,   label="Salt bridge anchor (primary)"),
    mpatches.Patch(color=COL_HY,   label="Hydrophobic anchor (primary)"),
    mpatches.Patch(color=COL_CONS, label="Conserved (>=0.8)"),
    mpatches.Patch(color=COL_VAR,  label="Variable (<0.8)"),
]
fig.legend(handles=legend_patches, loc="lower center", ncol=4,
           fontsize=9, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.02))

fig.suptitle("NSP12–NSP13 Interface: Conservation Across 5 Coronaviruses\n"
             "SB: ASP901–LYS94  |  Hydrophobic: MET902  |  "
             "SARS-CoV-1/2 selective (MET902) + pan-cov SB (LYS94)",
             fontsize=11, y=1.02, style="italic")

out1 = f"{RESULTS_DIR}/Fig1_NSP12-NSP13_conservation_bars_8.png"
fig.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {out1}")

# ══════════════════════════════════════════════════════════════════════════
# FIG 2 — Conservation heatmap
# ══════════════════════════════════════════════════════════════════════════
all_cons = nsp12_cons + nsp13_cons
res_labels = [f"{'NSP12:' if i<len(nsp12_cons) else 'NSP13:'}"
              f"{d['res_name_3']}{d['pdb_res_num']}"
              for i, d in enumerate(all_cons)]

# Build matrix: residues x species
species_keys = [f"aa_{cov.replace('-','_').replace('/','_')}"
                for cov in coronaviruses]
species_display = [cov.replace("HCoV-","") for cov in coronaviruses]

matrix = []
for d in all_cons:
    row = []
    for sk in species_keys:
        aa = d.get(sk, "-")
        row.append(aa)
    matrix.append(row)

# Color by conservation: same as SARS-CoV-2 = green, different = red, gap = gray
sars2_key = f"aa_{coronaviruses[0].replace('-','_').replace('/','_')}"
ref_aas   = [d.get(sars2_key, "-") for d in all_cons]

color_matrix = []
for i, row in enumerate(matrix):
    ref = ref_aas[i]
    crow = []
    for aa in row:
        if aa == "-":
            crow.append(0.5)
        elif aa == ref:
            crow.append(1.0)
        else:
            crow.append(0.0)
    color_matrix.append(crow)

cmap = LinearSegmentedColormap.from_list(
    "cons_map", ["#d62728", "#f5f5f5", "#2ca02c"], N=256
)

fig, ax = plt.subplots(figsize=(max(8, len(coronaviruses)*1.5),
                                 max(6, len(all_cons)*0.45)))
fig.patch.set_facecolor("white")

im = ax.imshow(color_matrix, cmap=cmap, vmin=0, vmax=1,
               aspect="auto", interpolation="nearest")

# Annotate with amino acid letters
for i in range(len(all_cons)):
    for j in range(len(coronaviruses)):
        aa = matrix[i][j]
        cv = color_matrix[i][j]
        text_color = "white" if cv < 0.2 or cv > 0.8 else "black"
        ax.text(j, i, aa, ha="center", va="center",
                fontsize=10, fontweight="bold", color=text_color)

ax.set_xticks(range(len(coronaviruses)))
ax.set_xticklabels(species_display, rotation=30, ha="right", fontsize=10)
ax.set_yticks(range(len(res_labels)))
ax.set_yticklabels(res_labels, fontsize=9)

# Divider between NSP12 and NSP13
ax.axhline(len(nsp12_cons) - 0.5, color="black", linewidth=2.5)
ax.text(len(coronaviruses) - 0.4, len(nsp12_cons)/2 - 0.5,
        "NSP12", va="center", fontsize=9, fontweight="bold",
        color=COL_NSP12, rotation=90)
ax.text(len(coronaviruses) - 0.4,
        len(nsp12_cons) + len(nsp13_cons)/2 - 0.5,
        "NSP13", va="center", fontsize=9, fontweight="bold",
        color=COL_NSP13, rotation=90)

# Mark primary residues
primary_rows = [i for i, d in enumerate(all_cons)
                if (i < len(nsp12_cons) and d["pdb_res_num"] in PRIMARY_NSP12) or
                   (i >= len(nsp12_cons) and d["pdb_res_num"] in PRIMARY_NSP13)]
for pr in primary_rows:
    ax.get_yticklabels()[pr].set_color(COL_SB)
    ax.get_yticklabels()[pr].set_fontweight("bold")

plt.colorbar(im, ax=ax, label="Conservation (green=identical, red=different)",
             shrink=0.6, pad=0.12)
ax.set_title("NSP12–NSP13 Interface Residue Conservation Across Coronaviruses\n"
             "★ Primary pharmacophore residues in red labels",
             fontsize=11, fontweight="bold", pad=10)

out2 = f"{RESULTS_DIR}/Fig2_NSP12-NSP13_conservation_heatmap_8.png"
fig.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {out2}")

# ══════════════════════════════════════════════════════════════════════════
# FIG 3 — Contact type distribution across 3 crystal structures
# ══════════════════════════════════════════════════════════════════════════
struct_labels = ["6XEZ\n(3.50 Å)", "7CXM\n(3.20 Å)", "7RDY\n(3.10 Å)"]
pdb_ids       = ["6XEZ", "7CXM", "7RDY"]

sb_counts = []
hb_counts = []
hy_counts = []
vdw_counts = []

for pdb_id in pdb_ids:
    r = iface_data["interface_results"][pdb_id]
    sb_counts.append(r["n_salt_bridges"])
    hb_counts.append(r["n_h_bonds"])
    hy_counts.append(r["n_hydrophobic"])
    vdw_counts.append(r["n_total"] - r["n_salt_bridges"]
                      - r["n_h_bonds"] - r["n_hydrophobic"])

x     = np.arange(len(struct_labels))
width = 0.20

fig, ax = plt.subplots(figsize=(9, 6))
fig.patch.set_facecolor("white")
ax.set_facecolor("#f8f8f8")

b1 = ax.bar(x - 1.5*width, sb_counts,  width, label="Salt bridge", color=COL_SB,   edgecolor="black", lw=0.7)
b2 = ax.bar(x - 0.5*width, hb_counts,  width, label="H-bond",      color="#1f77b4", edgecolor="black", lw=0.7)
b3 = ax.bar(x + 0.5*width, hy_counts,  width, label="Hydrophobic",  color=COL_HY,   edgecolor="black", lw=0.7)
b4 = ax.bar(x + 1.5*width, vdw_counts, width, label="VdW",          color="#bcbd22", edgecolor="black", lw=0.7)

for bars in [b1, b2, b3, b4]:
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.5,
                    str(int(h)), ha="center", va="bottom", fontsize=8)

ax.set_xticks(x)
ax.set_xticklabels(struct_labels, fontsize=11)
ax.set_ylabel("Contact count", fontsize=11)
ax.set_title("NSP12–NSP13 Interface Contact Types Across Crystal Structures\n"
             "Hydrophobic-dominated interface | MET902 primary anchor",
             fontsize=11, fontweight="bold")
ax.legend(fontsize=10, framealpha=0.9)
ax.grid(axis="y", alpha=0.4)

# Annotation: MET902 dominance
ax.annotate("MET902: 38–51\ncontacts/structure",
            xy=(x[0] + 0.5*width, hy_counts[0]),
            xytext=(x[1], max(hy_counts) * 0.75),
            fontsize=9, color=COL_HY, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=COL_HY, lw=1.5))

out3 = f"{RESULTS_DIR}/Fig3_NSP12-NSP13_contact_types_8.png"
fig.savefig(out3, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {out3}")

print("\nAll figures saved:")
print(f"  {out1}")
print(f"  {out2}")
print(f"  {out3}")
