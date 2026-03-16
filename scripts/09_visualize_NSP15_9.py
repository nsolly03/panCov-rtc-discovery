"""
Script 09_9 — Publication Figures NSP15 Homodimer
Complex  : NSP15-NSP15 (NendoU dimer interface)
Suffix   : _9
Figures  : Fig1 conservation bars, Fig2 heatmap, Fig3 contact types
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

OUT_DIR  = "results"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────
with open("02-validation/NSP15/conservation_summary_9.json") as f:
    cons = json.load(f)

with open("02-validation/NSP15/interface_analysis_9.json") as f:
    iface = json.load(f)

STRUCTURES = ["6VWW", "6W01", "6WLC", "6WXC"]
ORGANISMS  = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV","HCoV-229E","HCoV-NL63"]

# Primary SB residues (UniProt numbering for conservation, PDB for structure)
PRIMARY_UNI = {39:"ASP39", 61:"ARG61", 90:"ARG90", 266:"GLU266"}
PRIMARY_PDB = {40:"ASP40", 62:"ARG62", 91:"ARG91", 267:"GLU267"}

# ── Figure 1 — Conservation bars ──────────────────────────────
hotspots = cons["hotspot_details"]
hotspots_sorted = sorted(hotspots, key=lambda x: x["score"], reverse=True)

positions = [h["position"] for h in hotspots_sorted]
scores    = [h["score"]    for h in hotspots_sorted]
labels    = ["{}{}".format(h["ref_aa"], h["position"])
             for h in hotspots_sorted]

# Color coding
colors = []
for h in hotspots_sorted:
    pos = h["position"]
    if pos in PRIMARY_UNI:
        colors.append("#C0392B")   # red — primary SB
    elif h["score"] >= 0.8:
        colors.append("#2980B9")   # blue — conserved
    else:
        colors.append("#95A5A6")   # gray — variable

fig, ax = plt.subplots(figsize=(14, 7))
bars = ax.barh(range(len(scores)), scores,
               color=colors, edgecolor='white', linewidth=0.5)

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel("Conservation score (5 coronaviruses)", fontsize=11)
ax.set_title("NSP15 Homodimer Interface — Conservation of Hotspot Residues",
             fontsize=12, fontweight='bold')
ax.axvline(x=0.8, color='navy', linestyle='--', alpha=0.6,
           label='Conservation threshold (0.8)')
ax.axvline(x=1.0, color='black', linestyle='-', alpha=0.3)
ax.set_xlim(0, 1.1)

# Star primary SB residues
for i, h in enumerate(hotspots_sorted):
    if h["position"] in PRIMARY_UNI:
        ax.text(h["score"] + 0.02, i, "★", va='center',
                fontsize=9, color='#C0392B')

legend_elements = [
    mpatches.Patch(color='#C0392B', label='Primary SB anchor ★'),
    mpatches.Patch(color='#2980B9', label='Conserved (≥0.8)'),
    mpatches.Patch(color='#95A5A6', label='Variable (<0.8)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

# Annotation box
ax.text(0.02, 0.02,
    "Interface mean: {:.3f}\nSurface mean: {:.3f}\nConserved: {}/{}".format(
        cons["mean_iface_conservation"],
        cons["mean_surface_conservation"],
        cons["n_conserved_hotspots"],
        cons["n_hotspots"]),
    transform=ax.transAxes, fontsize=9,
    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
path1 = f"{OUT_DIR}/Fig1_NSP15_conservation_bars_9.png"
plt.savefig(path1, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {path1}")

# ── Figure 2 — Conservation heatmap ───────────────────────────
# Build matrix: hotspots x organisms
matrix = []
row_labels = []
for h in hotspots_sorted:
    row = []
    for org in ORGANISMS:
        aa_here = h["per_org"].get(org, "-")
        ref_aa  = h["ref_aa"]
        row.append(1.0 if aa_here == ref_aa
                   else (0.5 if aa_here not in ["-","X"] else 0.0))
    matrix.append(row)
    row_labels.append("{}{}".format(h["ref_aa"], h["position"]))

matrix = np.array(matrix)

fig, ax = plt.subplots(figsize=(8, max(8, len(row_labels)*0.22)))
im = ax.imshow(matrix, cmap='RdYlGn', aspect='auto',
               vmin=0, vmax=1)

ax.set_xticks(range(len(ORGANISMS)))
ax.set_xticklabels([o.replace("HCoV-","") for o in ORGANISMS],
                   rotation=30, ha='right', fontsize=9)
ax.set_yticks(range(len(row_labels)))
ax.set_yticklabels(row_labels, fontsize=7)

# Add AA identity in each cell
for i in range(len(row_labels)):
    for j, org in enumerate(ORGANISMS):
        h = hotspots_sorted[i]
        aa = h["per_org"].get(org, "-")
        ax.text(j, i, aa, ha='center', va='center',
                fontsize=6, color='black')

# Mark primary SB rows with red border
for i, h in enumerate(hotspots_sorted):
    if h["position"] in PRIMARY_UNI:
        for spine_side in ['top','bottom','left','right']:
            pass  # mark with text instead
        ax.text(len(ORGANISMS)+0.1, i, "★", va='center',
                fontsize=8, color='#C0392B')

ax.set_title("NSP15 Interface Conservation Heatmap\n(green=identical, red=different)",
             fontsize=11, fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.4, label='Identity score')
plt.tight_layout()
path2 = f"{OUT_DIR}/Fig2_NSP15_conservation_heatmap_9.png"
plt.savefig(path2, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {path2}")

# ── Figure 3 — Contact types ───────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))

sb_vals = [iface["per_structure"][s]["sb"] for s in STRUCTURES]
hb_vals = [iface["per_structure"][s]["hb"] for s in STRUCTURES]
hy_vals = [iface["per_structure"][s]["hy"] for s in STRUCTURES]

x   = np.arange(len(STRUCTURES))
w   = 0.25

ax.bar(x - w,   sb_vals, w, label='Salt bridges', color='#E74C3C')
ax.bar(x,       hb_vals, w, label='H-bonds',       color='#3498DB')
ax.bar(x + w,   hy_vals, w, label='Hydrophobic',   color='#2ECC71')

ax.set_xticks(x)
ax.set_xticklabels(STRUCTURES, fontsize=11)
ax.set_ylabel("Contact count", fontsize=11)
ax.set_title("NSP15 Homodimer Interface — Contact Types per Structure",
             fontsize=12, fontweight='bold')
ax.legend(fontsize=10)

# Annotate SB pairs
sb_text = "Primary SBs (all 4 structures):\nASP40--ARG91 (2.65-3.02 Å) ★\nARG62--GLU267 (3.49-3.73 Å) ★"
ax.text(0.98, 0.95, sb_text,
        transform=ax.transAxes, fontsize=8, va='top', ha='right',
        bbox=dict(boxstyle='round', facecolor='#FADBD8', alpha=0.8))

# Annotate interface character
ax.text(0.02, 0.95,
        "Interface: HY-dominated\n(~87% hydrophobic contacts)",
        transform=ax.transAxes, fontsize=8, va='top',
        bbox=dict(boxstyle='round', facecolor='#D5F5E3', alpha=0.8))

plt.tight_layout()
path3 = f"{OUT_DIR}/Fig3_NSP15_contact_types_9.png"
plt.savefig(path3, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {path3}")

print("\nScript 09_9 complete — 3 figures generated")
print("Next: Script 10_9 (BSA + AlaScan + composite ranking)")
