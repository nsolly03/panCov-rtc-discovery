#!/usr/bin/env python3
"""
Script 09_visualize_NSP12-NSP13_8.py  (IMPROVED)

Publication figures for NSP12-NSV13 interface analysis.
Style aligned with NSP10-NSP14 Script 09_2 (project standard).

Fig1 — Conservation bar charts (NSP12 + NSP13, horizontal)
Fig2 — Conservation heatmap with AA code legend
Fig3 — Contact type distribution across 3 crystal structures

conda activate rtc-discovery
cd ~/projects/rtc-pan-coronavirus
python scripts/09_visualize_NSP12-NSP13_8.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from Bio import AlignIO
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
VAL_DIR = PROJECT / "02-validation" / "NSP12-NSP13"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "results"
RES_DIR.mkdir(exist_ok=True)

# ── Global style (project standard) ───────────────────────────────────────
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.titlesize":    13,
    "axes.titleweight":  "bold",
    "axes.labelsize":    11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "xtick.major.size":  4,
    "ytick.major.size":  4,
    "legend.frameon":    False,
    "legend.fontsize":   9,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

BLUE   = "#2C5F8A"
GREEN  = "#2E8B57"
RED    = "#C0392B"
ORANGE = "#E67E22"
PURPLE = "#7D3C98"
GREY   = "#95A5A6"
LGREY  = "#ECF0F1"
PINK   = "#C9507A"

AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS",
    "Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE",
    "L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
    "S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV","HCoV-229E","HCoV-NL63"]

# Hotspot positions (6XEZ/7RDY PDB numbering)
NSP12_HOTSPOTS = [900, 901, 902, 903, 904]
NSP13_HOTSPOTS = [90, 91, 92, 93, 94, 95, 96]

PRIMARY_NSP12  = {901, 902}   # ASP901(SB) + MET902(hydrophobic)
PRIMARY_NSP13  = {93, 94}     # TYR93(#1) + LYS94(SB pan-cov)
SB_NSP12       = {901}
SB_NSP13       = {94}

# Contact counts from Script 05_8
CONTACT_COUNTS = {
    "6XEZ\n(3.50 Å)" : {"Salt Bridge": 0, "H-Bond": 1, "Hydrophobic": 42, "VdW": 21},
    "7CXM\n(3.20 Å)" : {"Salt Bridge": 1, "H-Bond": 1, "Hydrophobic": 55, "VdW": 38},
    "7RDY\n(3.10 Å)" : {"Salt Bridge": 1, "H-Bond": 2, "Hydrophobic": 53, "VdW": 52},
}

SB_ANNOTATIONS = {
    "6XEZ\n(3.50 Å)" : ("ASP901–LYS94", "absent",  False),
    "7CXM\n(3.20 Å)" : ("ASP901–LYS94", "4.12 Å",  True),
    "7RDY\n(3.10 Å)" : ("ASP901–LYS94", "3.95 Å",  True),
}


# ── Load data ──────────────────────────────────────────────────────────────
def load_data():
    cons12 = pd.read_csv(VAL_DIR / "conservation_NSP12.csv")
    cons13 = pd.read_csv(VAL_DIR / "conservation_NSP13.csv")
    with open(VAL_DIR / "interface_analysis_8.json") as f:
        iface = json.load(f)
    return cons12, cons13, iface


# ── Figure 1 — Conservation bar charts ────────────────────────────────────
def fig_conservation_bars(cons12, cons13):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6),
                             constrained_layout=True)
    fig.suptitle(
        "NSP12–NSP13 Interface: Hotspot Conservation Across 5 Coronaviruses",
        fontsize=14, fontweight="bold", y=1.02)

    datasets = [
        (axes[0], cons12, "NSP12", NSP12_HOTSPOTS, PRIMARY_NSP12, SB_NSP12,
         "C-terminal tail (residues 900–904)"),
        (axes[1], cons13, "NSP13", NSP13_HOTSPOTS, PRIMARY_NSP13, SB_NSP13,
         "N-terminal region (residues 90–96)"),
    ]

    for ax, df, nsp, hotspots, primary_set, sb_set, subtitle in datasets:
        # Filter to hotspot positions only
        hdf = df[df["pdb_res_num"].isin(hotspots)].copy()
        hdf = hdf.sort_values("conservation", ascending=True)

        colors = []
        for _, row in hdf.iterrows():
            rn = row["pdb_res_num"]
            if rn in sb_set:
                colors.append(RED)
            elif rn in primary_set:
                colors.append(ORANGE)
            elif row["conservation"] >= 0.8:
                colors.append(BLUE)
            else:
                colors.append(GREY)

        y    = np.arange(len(hdf))
        ax.barh(y, hdf["conservation"], color=colors,
                edgecolor="white", linewidth=0.6, height=0.72)
        ax.axvline(0.8, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6, label="Threshold (0.8)")

        for i, (_, row) in enumerate(hdf.iterrows()):
            ax.text(row["conservation"] + 0.01, i,
                    f"{row['conservation']:.3f}",
                    va="center", fontsize=8.5)

        ylabels = []
        for _, row in hdf.iterrows():
            rn      = row["pdb_res_num"]
            resname = row["res_name_3"]
            lbl     = f"{resname}{rn}"
            if rn in sb_set:
                lbl = f"★  {lbl}  (SB)"
            elif rn in primary_set:
                lbl = f"★  {lbl}"
            ylabels.append(lbl)
        ax.set_yticks(y)
        ax.set_yticklabels(ylabels, fontsize=9.5)

        ax.set_xlim(0, 1.18)
        ax.set_xlabel("Conservation Score  (1.0 = identical in all 5 species)")
        n_cons = len(hdf[hdf["conservation"] >= 0.8])
        ax.set_title(
            f"{nsp}  —  {subtitle}\n({n_cons}/{len(hdf)} residues conserved ≥ 0.8)",
            pad=8)
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)

        patches = [
            mpatches.Patch(color=RED,    label="Salt bridge anchor ★"),
            mpatches.Patch(color=ORANGE, label="Primary pharmacophore ★"),
            mpatches.Patch(color=BLUE,   label="Conserved ≥ 0.8"),
            mpatches.Patch(color=GREY,   label="Variable < 0.8"),
        ]
        ax.legend(handles=patches, loc="lower right", fontsize=8.5)

    fig.text(
        0.5, -0.03,
        "★  ASP901 / MET902 (NSP12) primary pharmacophores  |  "
        "TYR93 (#1 ranked, AlaLoss=29) / LYS94 (SB, pan-cov, cons=1.000) (NSP13)\n"
        "MET902: SARS-CoV-1/2 selective (M→S in MERS/229E/NL63)  |  "
        "LYS94: K conserved in all 5 species",
        ha="center", fontsize=8.5, style="italic", color="#555555")

    out = RES_DIR / "Fig1_NSP12-NSP13_conservation_bars_8.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Figure 2 — Conservation heatmap ───────────────────────────────────────
def fig_heatmap(cons12, cons13):
    cmap = LinearSegmentedColormap.from_list(
        "cons", ["#FDFEFE","#AED6F1","#2980B9","#1A5276"], N=256)

    fig = plt.figure(figsize=(20, 6), constrained_layout=True)
    gs  = fig.add_gridspec(1, 3, width_ratios=[8, 10, 1.6])
    ax12  = fig.add_subplot(gs[0])
    ax13  = fig.add_subplot(gs[1])
    ax_leg = fig.add_subplot(gs[2])

    fig.suptitle(
        "Amino Acid Identity Across 5 Coronaviruses\n"
        "NSP12–NSP13 Interface Hotspot Residues",
        fontsize=13, fontweight="bold")

    im_last = None

    for ax, nsp, df, hotspots, primary_set, sb_set in [
        (ax12,  "NSP12", cons12, NSP12_HOTSPOTS, PRIMARY_NSP12, SB_NSP12),
        (ax13,  "NSP13", cons13, NSP13_HOTSPOTS, PRIMARY_NSP13, SB_NSP13),
    ]:
        # Build species key mapping
        species_keys = {}
        for cov in CORONAVIRUSES:
            key = f"aa_{cov.replace('-','_').replace('/','_').replace(' ','_')}"
            species_keys[cov] = key

        score_map = dict(zip(df["pdb_res_num"], df["conservation"]))

        # Build matrix
        matrix = []
        scores = []
        for cov in CORONAVIRUSES:
            key = species_keys[cov]
            row_aa  = []
            row_sc  = []
            for p in hotspots:
                sub = df[df["pdb_res_num"] == p]
                if sub.empty or key not in sub.columns:
                    row_aa.append("-")
                    row_sc.append(0.0)
                else:
                    aa = sub[key].values[0]
                    row_aa.append(aa if pd.notna(aa) else "-")
                    row_sc.append(score_map.get(p, 0.0))
            matrix.append(row_aa)
            scores.append(row_sc)

        scores_arr = np.array(scores)
        im = ax.imshow(scores_arr, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        im_last = im

        # Annotate cells with AA letters
        for i in range(len(CORONAVIRUSES)):
            for j in range(len(hotspots)):
                aa = matrix[i][j]
                tc = "white" if scores_arr[i, j] > 0.65 else "black"
                ax.text(j, i, aa, ha="center", va="center",
                        fontsize=12, fontweight="bold", color=tc)

        # Red borders on primary residues
        for j, p in enumerate(hotspots):
            if p in sb_set:
                color = RED
            elif p in primary_set:
                color = ORANGE
            else:
                continue
            ax.add_patch(plt.Rectangle(
                (j-0.5, -0.5), 1, len(CORONAVIRUSES),
                fill=False, edgecolor=color, linewidth=2.5, zorder=3))
            marker = "★" if p in sb_set else "◆"
            ax.text(j, -0.9, marker, ha="center", va="center",
                    fontsize=11, color=color)

        # X-axis labels: 3-letter AA + position (from SARS-CoV-2 row)
        sars2_key = species_keys["SARS-CoV-2"]
        col_labels = []
        for p in hotspots:
            sub = df[df["pdb_res_num"] == p]
            if not sub.empty and sars2_key in sub.columns:
                aa1   = sub[sars2_key].values[0]
                three = AA3.get(str(aa1), str(aa1)) if pd.notna(aa1) else "?"
                col_labels.append(f"{three}{p}")
            else:
                rname = sub["res_name_3"].values[0] if not sub.empty else "?"
                col_labels.append(f"{rname}{p}")

        ax.set_xticks(range(len(hotspots)))
        ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(CORONAVIRUSES)))
        ax.set_yticklabels(CORONAVIRUSES if ax == ax12 else [],
                           fontsize=10)
        ax.set_title(nsp, pad=20)
        ax.tick_params(length=0)

        for x in np.arange(-0.5, len(hotspots), 1):
            ax.axvline(x, color="white", linewidth=1.0)
        for y_ in np.arange(-0.5, len(CORONAVIRUSES), 1):
            ax.axhline(y_, color="white", linewidth=1.0)

    # Colorbar
    cbar = fig.colorbar(im_last, ax=[ax12, ax13], shrink=0.75,
                        pad=0.015, label="Conservation Score")
    cbar.ax.tick_params(labelsize=9)

    # AA code legend panel
    ax_leg.axis("off")
    entries = [
        ("A","ALA — Alanine"),    ("D","ASP — Aspartate"),
        ("E","GLU — Glutamate"),  ("F","PHE — Phenylalanine"),
        ("G","GLY — Glycine"),    ("K","LYS — Lysine"),
        ("L","LEU — Leucine"),    ("M","MET — Methionine"),
        ("N","ASN — Asparagine"), ("P","PRO — Proline"),
        ("S","SER — Serine"),     ("T","THR — Threonine"),
        ("V","VAL — Valine"),     ("Y","TYR — Tyrosine"),
    ]
    ax_leg.text(0.05, 0.99, "AA Code Key",
                transform=ax_leg.transAxes,
                fontsize=9, fontweight="bold", va="top")
    for i, (code, name) in enumerate(entries):
        ax_leg.text(0.05, 0.93 - i*0.063,
                    f"{code}  =  {name}",
                    transform=ax_leg.transAxes,
                    fontsize=7.5, va="top", color="#2C3E50",
                    fontfamily="monospace")

    fig.text(
        0.45, -0.04,
        "★  Salt bridge anchors  |  ◆  Primary pharmacophores  |  "
        "Red border = SB residue  |  Orange border = primary",
        ha="center", fontsize=8.5, style="italic", color="#555555")

    out = RES_DIR / "Fig2_NSP12-NSP13_conservation_heatmap_8.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Figure 3 — Contact type distribution ──────────────────────────────────
def fig_contacts():
    structs = list(CONTACT_COUNTS.keys())
    ctypes  = ["Salt Bridge", "H-Bond", "Hydrophobic", "VdW"]
    colors  = [RED, BLUE, ORANGE, LGREY]

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    fig.suptitle(
        "NSP12–NSP13 Interface: Contact Types Across Three Crystal Structures",
        fontsize=13, fontweight="bold")

    x       = np.arange(len(structs))
    width   = 0.19
    offsets = np.array([-0.285, -0.095, 0.095, 0.285])

    for ct, col, off in zip(ctypes, colors, offsets):
        vals = [CONTACT_COUNTS[s][ct] for s in structs]
        bars = ax.bar(x + off, vals, width, label=ct,
                      color=col, edgecolor="white", linewidth=0.6, zorder=3)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 0.8,
                        str(v), ha="center", va="bottom",
                        fontsize=9, fontweight="bold")

    # Salt bridge annotation boxes
    for i, s in enumerate(structs):
        name, dist, present = SB_ANNOTATIONS[s]
        sb_y  = CONTACT_COUNTS[s]["Salt Bridge"]
        mark  = "✔" if present else "✘"
        mc    = GREEN if present else GREY
        label = f"{mark} {name}\n{dist}"
        ax.text(x[i] + offsets[0],
                sb_y + 2.5,
                label,
                ha="center", va="bottom",
                fontsize=8, color=mc, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25",
                          facecolor="white", edgecolor=mc,
                          linewidth=1.2, alpha=0.9))

    # MET902 dominance annotation
    ax.annotate(
        "MET902:\n38–51 contacts\n(hydrophobic)",
        xy=(x[1] + offsets[2], CONTACT_COUNTS["7CXM\n(3.20 Å)"]["Hydrophobic"]),
        xytext=(x[2] - 0.1, 65),
        fontsize=8.5, color=ORANGE, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.5))

    ax.set_xticks(x)
    ax.set_xticklabels(structs, fontsize=10)
    ax.set_ylabel("Number of contacts")
    ax.set_ylim(0, 80)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.grid(axis="y", alpha=0.25, linewidth=0.7, zorder=0)
    ax.legend(title="Contact type", loc="upper left",
              fontsize=9, title_fontsize=9)

    # Total contacts
    totals = {"6XEZ\n(3.50 Å)": 64,
              "7CXM\n(3.20 Å)": 95,
              "7RDY\n(3.10 Å)": 108}
    for i, s in enumerate(structs):
        ax.text(x[i], 76, f"Total: {totals[s]}",
                ha="center", fontsize=9,
                color="#2C3E50", fontweight="bold")

    fig.text(
        0.5, -0.03,
        "Hydrophobic-dominated interface  |  "
        "MET902(NSP12) primary anchor  |  "
        "ASP901–LYS94 SB absent in 6XEZ (3.50 Å resolution artifact)",
        ha="center", fontsize=8.5, style="italic", color="#555555")

    out = RES_DIR / "Fig3_NSP12-NSP13_contact_types_8.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("\n" + "="*55)
    print("  Script 09_8 (IMPROVED): Figures — NSP12-NSP13")
    print("="*55 + "\n")

    cons12, cons13, iface = load_data()

    print("  Generating Figure 1 — Conservation bar charts ...")
    fig_conservation_bars(cons12, cons13)

    print("  Generating Figure 2 — Conservation heatmap ...")
    fig_heatmap(cons12, cons13)

    print("  Generating Figure 3 — Contact types ...")
    fig_contacts()

    print(f"\n  All 3 figures saved to: results/")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
