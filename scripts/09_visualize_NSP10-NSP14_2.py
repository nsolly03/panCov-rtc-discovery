"""
Script 09_2: Publication Figures — NSP10-NSP14
===============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Fig 1 — Hotspot conservation bar chart (full 3-letter AA names)
Fig 2 — Conservation heatmap with one-letter code legend
Fig 3 — Top20 hotspot scores across 3 structures

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/09_visualize_NSP10-NSP14_2.py
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
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
OUT_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
OUT_DIR.mkdir(exist_ok=True)

# ── Global style ───────────────────────────────────────────
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
GREEN  = "#27AE60"
GREEN  = "#2E8B57"
RED    = "#C0392B"
ORANGE = "#E67E22"
PURPLE = "#7D3C98"
GREY   = "#95A5A6"
LGREY  = "#ECF0F1"

# One-letter → three-letter amino acid map
AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS",
    "Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE",
    "L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
    "S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}

CONS_NSP10 = {5,19,21,40,42,44,45,80,93}
CONS_NSP14 = {4,7,8,9,10,20,25,27,127,201}
CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

# Known contact counts from Script 05 output (WORKLOG entry 010)
CONTACT_COUNTS = {
    "7DIY\nSARS-CoV-2\n(2.69Å)": {"H-Bond":17, "Hydrophobic":41,
                                    "Salt Bridge":1, "VdW":115},
    "5C8T\nSARS-CoV-1\n(3.20Å)": {"H-Bond":17, "Hydrophobic":39,
                                    "Salt Bridge":2, "VdW":108},
    "AF3 model\n(iptm=0.89)":     {"H-Bond":16, "Hydrophobic":42,
                                    "Salt Bridge":1, "VdW":113},
}


def load_data():
    cons10 = pd.read_csv(RES_DIR / "conservation_NSP10.csv")
    cons14 = pd.read_csv(RES_DIR / "conservation_NSP14.csv")
    with open(RES_DIR / "interface_analysis.json") as f:
        interface = json.load(f)
    return cons10, cons14, interface


# ── Figure 1 — Conservation bar chart ─────────────────────
def fig_conservation_bars(cons10, cons14):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                             constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP14 Interface: Hotspot Conservation Across 5 Coronaviruses",
        fontsize=14, fontweight="bold", y=1.02)

    for ax, df, nsp, cons_set, primary_set in [
        (axes[0], cons10, "NSP10", CONS_NSP10, {80, 93}),
        (axes[1], cons14, "NSP14", CONS_NSP14, {126, 127}),
    ]:
        hdf = df[df["is_hotspot"]].copy()
        hdf = hdf.sort_values("conservation", ascending=True)

        colors = []
        for _, row in hdf.iterrows():
            if row["position"] in primary_set:
                colors.append(RED)
            elif row["conservation"] >= 0.8:
                colors.append(BLUE)
            else:
                colors.append(GREY)

        y    = np.arange(len(hdf))
        bars = ax.barh(y, hdf["conservation"], color=colors,
                       edgecolor="white", linewidth=0.6, height=0.72)

        ax.axvline(0.8, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6)

        # Score annotations
        for i, (_, row) in enumerate(hdf.iterrows()):
            ax.text(row["conservation"] + 0.01, i,
                    f"{row['conservation']:.2f}",
                    va="center", fontsize=8.5)

        # Y-axis: full 3-letter labels e.g. HIS80, LYS93
        ax.set_yticks(y)
        ylabels = []
        for _, row in hdf.iterrows():
            three = AA3.get(row["aa_SARS2"], row["aa_SARS2"])
            pos   = int(row["position"])
            lbl   = f"{three}{pos}"
            if row["position"] in primary_set:
                lbl = f"★  {lbl}"
            ylabels.append(lbl)
        ax.set_yticklabels(ylabels, fontsize=9.5,
                           fontfamily="DejaVu Sans")

        ax.set_xlim(0, 1.13)
        ax.set_xlabel("Conservation Score  (1.0 = identical in all 5 coronaviruses)")
        n_cons = len(hdf[hdf["conservation"] >= 0.8])
        ax.set_title(f"{nsp}  ({n_cons} / {len(hdf)} residues conserved ≥ 0.8)",
                     pad=8)
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)

        patches = [
            mpatches.Patch(color=RED,  label="Primary salt bridge residue (★)"),
            mpatches.Patch(color=BLUE, label="Conserved ≥ 0.8"),
            mpatches.Patch(color=GREY, label="Variable < 0.8"),
        ]
        ax.legend(handles=patches, loc="lower right", fontsize=8.5)

    fig.text(0.5, -0.02,
             "★  HIS80 / LYS93 (NSP10)  and  ASP126 / THR127 (NSP14)"
             " — primary salt bridge anchor residues",
             ha="center", fontsize=9, style="italic", color="#555555")

    out = OUT_DIR / "Fig1_NSP10-NSP14_conservation_bars_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Figure 2 — Conservation heatmap ───────────────────────
def fig_heatmap(cons10, cons14):
    HOTSPOTS = {
        "NSP10": [5,19,21,40,42,44,45,80,93],
        "NSP14": [4,7,8,9,10,20,25,27,127,201],
    }
    cmap = LinearSegmentedColormap.from_list(
        "cons", ["#FDFEFE","#AED6F1","#2980B9","#1A5276"], N=256)

    fig = plt.figure(figsize=(18, 7), constrained_layout=True)
    # Main axes + legend strip
    gs   = fig.add_gridspec(1, 3, width_ratios=[9, 9, 1.4])
    axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    ax_leg = fig.add_subplot(gs[2])

    fig.suptitle(
        "Amino Acid Identity Across 5 Coronaviruses\n"
        "NSP10–NSP14 Conserved Hotspot Residues",
        fontsize=13, fontweight="bold")

    im_last = None
    for ax, nsp, df, hotspots in [
        (axes[0], "NSP10", cons10, HOTSPOTS["NSP10"]),
        (axes[1], "NSP14", cons14, HOTSPOTS["NSP14"]),
    ]:
        aln     = AlignIO.read(SEQ_DIR / f"{nsp}_aligned.fasta","fasta")
        ref_seq = next(r for r in aln if "SARS-CoV-2" in r.id)

        pos_map, ref_pos = {}, 0
        for i, aa in enumerate(str(ref_seq.seq)):
            if aa != "-":
                ref_pos += 1
                pos_map[ref_pos] = i

        matrix, row_labels = [], []
        for cov in CORONAVIRUSES:
            rec = next(
                (r for r in aln
                 if cov.replace("-","_") in r.id
                 or r.id.startswith(cov.replace(" ","_")[:6])),
                aln[CORONAVIRUSES.index(cov)])
            matrix.append([
                str(rec.seq)[pos_map[p]] if p in pos_map else "-"
                for p in hotspots])
            row_labels.append(cov)

        score_map = dict(zip(df["position"], df["conservation"]))
        scores = np.array([
            [score_map.get(p, 0) for p in hotspots]
            for _ in CORONAVIRUSES])

        im = ax.imshow(scores, cmap=cmap, vmin=0, vmax=1,
                       aspect="auto")
        im_last = im

        # Cell AA letters + 3-letter tooltip in title
        for i in range(len(CORONAVIRUSES)):
            for j, p in enumerate(hotspots):
                aa1 = matrix[i][j]
                tc  = "white" if scores[i][j] > 0.65 else "black"
                ax.text(j, i, aa1, ha="center", va="center",
                        fontsize=12, fontweight="bold", color=tc)

        # Primary target column red borders
        primary_pos = {80,93} if nsp == "NSP10" else {126,127}
        for j, p in enumerate(hotspots):
            if p in primary_pos:
                ax.add_patch(plt.Rectangle(
                    (j-0.5, -0.5), 1, len(CORONAVIRUSES),
                    fill=False, edgecolor=RED, linewidth=2.5,
                    zorder=3))
                ax.text(j, -0.9, "★", ha="center", va="center",
                        fontsize=12, color=RED)

        # X-axis: 3-letter AA + position
        sars2_seq = str(next(
            r for r in aln if "SARS-CoV-2" in r.id).seq)
        col_labels = []
        for p in hotspots:
            if p in pos_map:
                aa1   = sars2_seq[pos_map[p]]
                three = AA3.get(aa1, aa1)
                col_labels.append(f"{three}{p}")
            else:
                col_labels.append(str(p))

        ax.set_xticks(range(len(col_labels)))
        ax.set_xticklabels(col_labels, rotation=45,
                           ha="right", fontsize=9)
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels if ax == axes[0] else [],
                           fontsize=10)
        ax.set_title(nsp, pad=18)
        ax.tick_params(length=0)

        for x in np.arange(-0.5, len(hotspots), 1):
            ax.axvline(x, color="white", linewidth=1.0)
        for y in np.arange(-0.5, len(CORONAVIRUSES), 1):
            ax.axhline(y, color="white", linewidth=1.0)

    # Colourbar
    cbar = fig.colorbar(im_last, ax=axes, shrink=0.75,
                        pad=0.015, label="Conservation Score")
    cbar.ax.tick_params(labelsize=9)

    # One-letter code legend panel
    ax_leg.axis("off")
    legend_text = "One-letter\nAmino Acid\nCode\n\n"
    entries = [
        ("A", "ALA — Alanine"),
        ("D", "ASP — Aspartate"),
        ("F", "PHE — Phenylalanine"),
        ("H", "HIS — Histidine"),
        ("I", "ILE — Isoleucine"),
        ("K", "LYS — Lysine"),
        ("L", "LEU — Leucine"),
        ("M", "MET — Methionine"),
        ("N", "ASN — Asparagine"),
        ("P", "PRO — Proline"),
        ("S", "SER — Serine"),
        ("T", "THR — Threonine"),
        ("V", "VAL — Valine"),
        ("W", "TRP — Tryptophan"),
        ("Y", "TYR — Tyrosine"),
    ]
    ax_leg.text(0.05, 0.98, "AA Code Key",
                transform=ax_leg.transAxes,
                fontsize=9, fontweight="bold", va="top")
    for i, (code, name) in enumerate(entries):
        ax_leg.text(0.05, 0.92 - i * 0.057,
                    f"{code}  =  {name}",
                    transform=ax_leg.transAxes,
                    fontsize=7.5, va="top", color="#2C3E50",
                    fontfamily="monospace")

    fig.text(0.42, -0.03,
             "★  Primary drug target residues  |  "
             "Red borders = salt bridge anchor columns",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig2_NSP10-NSP14_conservation_heatmap_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Figure 3 — Contact types bar chart ────────────────────
def fig_contacts():
    structs  = list(CONTACT_COUNTS.keys())
    ctypes   = ["H-Bond", "Hydrophobic", "Salt Bridge", "VdW"]
    colors   = [BLUE, ORANGE, RED, LGREY]

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP14 Interface: Contact Types Across Three Structures",
        fontsize=13, fontweight="bold")

    x       = np.arange(len(structs))
    width   = 0.19
    offsets = np.array([-0.285, -0.095, 0.095, 0.285])

    for ct, col, off in zip(ctypes, colors, offsets):
        vals = [CONTACT_COUNTS[s][ct] for s in structs]
        bars = ax.bar(x + off, vals, width, label=ct,
                      color=col, edgecolor="white",
                      linewidth=0.6, zorder=3)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 1.5,
                        str(v), ha="center", va="bottom",
                        fontsize=9, fontweight="bold")

    # Salt bridge per-structure annotation
    # From Script 05: 7DIY=3.65A, 5C8T=2.59A, AF3=2.91A
    sb_data = {
        "7DIY\nSARS-CoV-2\n(2.69Å)": ("HIS80–ASP126", "3.65 Å", True),
        "5C8T\nSARS-CoV-1\n(3.20Å)": ("HIS80–ASP126", "2.59 Å", True),
        "AF3 model\n(iptm=0.89)":      ("HIS80–ASP126", "2.91 Å", True),
    }
    for i, s in enumerate(structs):
        name, dist, present = sb_data[s]
        sb_y = CONTACT_COUNTS[s]["Salt Bridge"]
        marker = "✔" if present else "✘"
        mc     = GREEN if present else RED
        # Label above the salt bridge bar
        ax.text(x[i] + offsets[2],
                sb_y + 3.5,
                f"{marker} {name}\n{dist}",
                ha="center", va="bottom",
                fontsize=8, color=mc,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25",
                          facecolor="white",
                          edgecolor=mc,
                          linewidth=1.2,
                          alpha=0.9))

    ax.set_xticks(x)
    ax.set_xticklabels(structs, fontsize=10)
    ax.set_ylabel("Number of contacts")
    ax.set_ylim(0, 148)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
    ax.grid(axis="y", alpha=0.25, linewidth=0.7, zorder=0)
    ax.legend(title="Contact type", loc="upper right",
              fontsize=9, title_fontsize=9)

    # Total contacts annotation
    totals = {"7DIY\nSARS-CoV-2\n(2.69Å)": 174,
              "5C8T\nSARS-CoV-1\n(3.20Å)": 166,
              "AF3 model\n(iptm=0.89)":     172}
    for i, s in enumerate(structs):
        ax.text(x[i], 142, f"Total: {totals[s]}",
                ha="center", fontsize=9,
                color="#2C3E50", fontweight="bold")

    out = OUT_DIR / "Fig3_NSP10-NSP14_contact_types_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────
def main():
    print("\n" + "="*55)
    print("  Script 09_2: Publication Figures — NSP10-NSP14")
    print("="*55 + "\n")

    cons10, cons14, interface = load_data()

    print("  Generating Figure 1 — Conservation bar charts...")
    fig_conservation_bars(cons10, cons14)

    print("  Generating Figure 2 — Conservation heatmap...")
    fig_heatmap(cons10, cons14)

    print("  Generating Figure 3 — Contact types...")
    fig_contacts()

    print(f"\n  All figures saved to: results/")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
