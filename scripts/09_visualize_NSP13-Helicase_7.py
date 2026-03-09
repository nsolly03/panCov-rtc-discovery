#!/usr/bin/env python3
"""
Script 09_visualize_NSP13-Helicase_7.py

Publication figures for NSP13-Helicase homodimer interface.

Fig1: Conservation bar chart — all 17 hotspot residues
Fig2: Conservation heatmap — 5 coronaviruses x hotspot residues
Fig3: Contact types bar chart — 7NIO + 6XEZ + salt bridge inventory

Output: results/
  Fig1_NSP13-Helicase_conservation_bars_7.png
  Fig2_NSP13-Helicase_conservation_heatmap_7.png
  Fig3_NSP13-Helicase_contact_types_7.png
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ── Paths ──────────────────────────────────────────────────────────────────
BASE        = os.path.expanduser("~/projects/rtc-pan-coronavirus")
VAL_DIR     = f"{BASE}/02-validation/NSP13-Helicase"
RESULTS_DIR = f"{BASE}/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────
with open(f"{VAL_DIR}/conservation_summary_7.json") as f:
    cons_data = json.load(f)

with open(f"{VAL_DIR}/interface_analysis_7.json") as f:
    iface_data = json.load(f)

cons_df = pd.read_csv(f"{VAL_DIR}/conservation_NSP13.csv")

# ── Data preparation ───────────────────────────────────────────────────────
ORGS = ["SARS-CoV-2", "SARS-CoV-1", "MERS-CoV", "HCoV-229E", "HCoV-NL63"]

ORG_COLS = {
    "SARS-CoV-2" : "aa_SARS_CoV_2",
    "SARS-CoV-1" : "aa_SARS_CoV_1",
    "MERS-CoV"   : "aa_MERS_CoV",
    "HCoV-229E"  : "aa_HCoV_229E",
    "HCoV-NL63"  : "aa_HCoV_NL63",
}

AA_FULL = {
    "ALA":"Alanine",   "ARG":"Arginine",  "ASN":"Asparagine",
    "ASP":"Aspartate", "CYS":"Cysteine",  "GLN":"Glutamine",
    "GLU":"Glutamate", "GLY":"Glycine",   "HIS":"Histidine",
    "ILE":"Isoleucine","LEU":"Leucine",   "LYS":"Lysine",
    "MET":"Methionine","PHE":"Phenylalanine","PRO":"Proline",
    "SER":"Serine",    "THR":"Threonine", "TRP":"Tryptophan",
    "TYR":"Tyrosine",  "VAL":"Valine",
}

PRIMARY_HOTSPOTS = {414, 477, 480, 482, 579, 580, 583, 584}
SB_RESIDUES      = {414, 580, 583}

# Sort hotspots by position
cons_df = cons_df.sort_values("local_pos").reset_index(drop=True)

# Residue labels — full AA name + position
cons_df["label"] = cons_df.apply(
    lambda r: f"{AA_FULL.get(r['res_name_3'], r['res_name_3'])}\n{r['local_pos']}", axis=1)

# ── FIGURE 1: Conservation bar chart ──────────────────────────────────────
print("Generating Fig1: conservation bar chart ...")

fig, ax = plt.subplots(figsize=(16, 6))

colors = []
for _, row in cons_df.iterrows():
    rn = row["local_pos"]
    if rn in SB_RESIDUES:
        colors.append("#d62728")   # red — salt bridge anchor
    elif rn in PRIMARY_HOTSPOTS:
        colors.append("#ff7f0e")   # orange — primary hotspot
    else:
        colors.append("#1f77b4")   # blue — secondary hotspot

bars = ax.bar(range(len(cons_df)), cons_df["conservation"],
              color=colors, edgecolor="white", linewidth=0.5, zorder=3)

# Threshold line
ax.axhline(0.8, color="black", linestyle="--", linewidth=1.2,
           label="Conservation threshold (0.8)", zorder=4)

# Annotations — conservation score on bar
for i, (_, row) in enumerate(cons_df.iterrows()):
    score = row["conservation"]
    va    = "bottom" if score < 0.85 else "top"
    ypos  = score + 0.02 if score < 0.85 else score - 0.04
    ax.text(i, ypos, f"{score:.2f}", ha="center", va=va,
            fontsize=7.5, fontweight="bold", color="black")

# Star markers for salt bridge residues
for i, (_, row) in enumerate(cons_df.iterrows()):
    if row["local_pos"] in SB_RESIDUES:
        ax.text(i, row["conservation"] + 0.07, "★",
                ha="center", fontsize=11, color="#d62728")

ax.set_xticks(range(len(cons_df)))
ax.set_xticklabels(cons_df["label"], fontsize=8, rotation=0)
ax.set_ylabel("Conservation score\n(1 − normalized Shannon entropy)", fontsize=11)
ax.set_xlabel("Hotspot residue (7NIO chain A numbering)", fontsize=11)
ax.set_title("NSP13 Homodimer Interface — Conservation Across 5 Coronaviruses",
             fontsize=13, fontweight="bold", pad=12)
ax.set_ylim(0, 1.25)
ax.set_xlim(-0.7, len(cons_df) - 0.3)
ax.grid(axis="y", alpha=0.3, zorder=0)
ax.spines[["top","right"]].set_visible(False)

# Legend
legend_patches = [
    mpatches.Patch(color="#d62728", label="Salt bridge anchor (★)"),
    mpatches.Patch(color="#ff7f0e", label="Primary hotspot"),
    mpatches.Patch(color="#1f77b4", label="Secondary hotspot"),
    plt.Line2D([0],[0], color="black", linestyle="--", label="Threshold (0.8)"),
]
ax.legend(handles=legend_patches, loc="upper right", fontsize=9,
          framealpha=0.9, edgecolor="gray")

# Scientific note
ax.text(0.01, 0.97,
        "SARS-CoV-1 and SARS-CoV-2 identical at all 17 positions\n"
        "ASP580: D→A (MERS), D→T (229E/NL63) — charge loss",
        transform=ax.transAxes, fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow",
                  ec="goldenrod", alpha=0.9))

plt.tight_layout()
fig1_path = f"{RESULTS_DIR}/Fig1_NSP13-Helicase_conservation_bars_7.png"
plt.savefig(fig1_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig1_path}")

# ── FIGURE 2: Conservation heatmap ────────────────────────────────────────
print("Generating Fig2: conservation heatmap ...")

# Build matrix: organisms × hotspot residues
heatmap_data = []
aa_annotations = []

for org in ORGS:
    col     = ORG_COLS[org]
    row_dat = []
    row_aa  = []
    for _, r in cons_df.iterrows():
        aa = r.get(col, "-")
        row_aa.append(str(aa) if pd.notna(aa) else "-")
        # Color by conservation score (same for all organisms at that column)
        row_dat.append(r["conservation"])
    heatmap_data.append(row_dat)
    aa_annotations.append(row_aa)

heatmap_arr = np.array(heatmap_data)
annot_arr   = np.array(aa_annotations)

fig, ax = plt.subplots(figsize=(18, 5))

cmap = sns.color_palette("RdYlGn", as_cmap=True)
sns.heatmap(heatmap_arr,
            ax=ax,
            cmap=cmap,
            vmin=0, vmax=1,
            annot=annot_arr,
            fmt="",
            annot_kws={"size": 9, "weight": "bold"},
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "Conservation score", "shrink": 0.8})

ax.set_xticks(np.arange(len(cons_df)) + 0.5)
ax.set_xticklabels(
    [f"{r['res_name_3']}{r['local_pos']}" for _, r in cons_df.iterrows()],
    rotation=45, ha="right", fontsize=9)
ax.set_yticklabels(ORGS, rotation=0, fontsize=10)
ax.set_title("NSP13 Homodimer Interface — Residue Identity Across 5 Coronaviruses",
             fontsize=13, fontweight="bold", pad=10)

# Mark primary hotspot columns
for i, (_, row) in enumerate(cons_df.iterrows()):
    if row["local_pos"] in SB_RESIDUES:
        ax.add_patch(plt.Rectangle((i, 0), 1, len(ORGS),
                                    fill=False, edgecolor="#d62728",
                                    linewidth=2.5, zorder=5))
    elif row["local_pos"] in PRIMARY_HOTSPOTS:
        ax.add_patch(plt.Rectangle((i, 0), 1, len(ORGS),
                                    fill=False, edgecolor="#ff7f0e",
                                    linewidth=1.5, zorder=5))

# AA code key
aa_key = ("One-letter codes: K=Lys, D=Asp, G=Gly, I=Ile, H=His, "
          "R=Arg, T=Thr, E=Glu, V=Val, A=Ala, N=Asn, S=Ser, Q=Gln")
fig.text(0.5, -0.04, aa_key, ha="center", fontsize=8,
         style="italic", color="dimgray")

plt.tight_layout()
fig2_path = f"{RESULTS_DIR}/Fig2_NSP13-Helicase_conservation_heatmap_7.png"
plt.savefig(fig2_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig2_path}")

# ── FIGURE 3: Contact types bar chart ─────────────────────────────────────
print("Generating Fig3: contact types bar chart ...")

structures = {
    "7NIO\n(A vs E)\n2.80 Å" : iface_data["interface_7NIO"],
    "6XEZ\n(E vs F)\n3.50 Å" : iface_data["interface_6XEZ"],
}

labels    = list(structures.keys())
sb_counts = [s["n_salt_bridges"] for s in structures.values()]
hb_counts = [s["n_h_bonds"]      for s in structures.values()]
hy_counts = [s["n_hydrophobic"]  for s in structures.values()]
tot_counts= [s["n_total"]        for s in structures.values()]

x     = np.arange(len(labels))
width = 0.2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                gridspec_kw={"width_ratios": [2, 1]})

# Left panel: grouped bars
b1 = ax1.bar(x - width,     sb_counts, width, label="Salt bridges",  color="#d62728", alpha=0.85)
b2 = ax1.bar(x,             hb_counts, width, label="H-bonds",       color="#1f77b4", alpha=0.85)
b3 = ax1.bar(x + width,     hy_counts, width, label="Hydrophobic",   color="#2ca02c", alpha=0.85)

# Total labels on top
for i, tot in enumerate(tot_counts):
    ax1.text(i, max(sb_counts[i], hb_counts[i], hy_counts[i]) + 0.5,
             f"Total={tot}", ha="center", fontsize=9,
             fontweight="bold", color="black")

ax1.set_xticks(x)
ax1.set_xticklabels(labels, fontsize=11)
ax1.set_ylabel("Number of contacts", fontsize=11)
ax1.set_title("Contact Types per Structure", fontsize=12, fontweight="bold")
ax1.legend(fontsize=10, framealpha=0.9)
ax1.grid(axis="y", alpha=0.3)
ax1.spines[["top","right"]].set_visible(False)

# Right panel: salt bridge inventory
ax2.axis("off")
sb_text = "Salt Bridge Inventory\n" + "─" * 28 + "\n\n"

all_sbs = iface_data.get("salt_bridges_all", [])
if all_sbs:
    for sb in all_sbs:
        sb_text += (f"[{sb['structure']}]\n"
                    f"  {sb['res1']} ── {sb['res2']}\n"
                    f"  {sb['dist']:.2f} Å\n\n")
else:
    sb_text += "No salt bridges detected\n"

sb_text += "─" * 28 + "\n"
sb_text += "\nNote: LYS414 forms dual\n"
sb_text += "simultaneous SBs with\n"
sb_text += "ASP580 AND ASP583 —\n"
sb_text += "unique in this project"

ax2.text(0.05, 0.95, sb_text, transform=ax2.transAxes,
         fontsize=9, va="top", fontfamily="monospace",
         bbox=dict(boxstyle="round,pad=0.6", fc="#fff8f0",
                   ec="#d62728", linewidth=1.5))

fig.suptitle("NSP13 Homodimer — Interface Contact Analysis",
             fontsize=13, fontweight="bold", y=1.02)

plt.tight_layout()
fig3_path = f"{RESULTS_DIR}/Fig3_NSP13-Helicase_contact_types_7.png"
plt.savefig(fig3_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig3_path}")

# ── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("FIGURES COMPLETE")
print("=" * 60)
print(f"  Fig1: {fig1_path}")
print(f"  Fig2: {fig2_path}")
print(f"  Fig3: {fig3_path}")
print(f"\n  Key scientific message:")
print(f"  SARS-CoV-1/2 selective interface — LYS414 dual SB")
print(f"  ASP580 charge loss in MERS (D->A) and HCoV (D->T)")
print(f"  4 backbone residues (GLY415/478, THR552, ALA553) conserved")
print(f"  pan-coronavirus but not druggable anchors")
