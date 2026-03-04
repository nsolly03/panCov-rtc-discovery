"""
Script 09_2: Publication Figures — NSP10-NSP16
===============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Fig 1 — Conservation bar charts (NSP10 + NSP16)
Fig 2 — Conservation heatmap (5 coronaviruses × hotspots)
Fig 3 — Contact types + 3 salt bridges + Zn1 annotation

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/09_visualize_NSP10-NSP16_2.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from Bio import AlignIO
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
OUT_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
OUT_DIR.mkdir(exist_ok=True)

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
    "legend.frameon":    False,
    "legend.fontsize":   9,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

BLUE   = "#2C5F8A"
GREEN  = "#27AE60"
RED    = "#C0392B"
ORANGE = "#E67E22"
PURPLE = "#7D3C98"
GOLD   = "#F39C12"
GREY   = "#95A5A6"

AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS",
    "Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE",
    "L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
    "S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

# Primary salt bridge residues (AF3 local)
PRIMARY_NSP10 = {80, 93, 95}
PRIMARY_NSP16 = {102, 106}

# Zn1 coordination neighbours
ZN1_NEIGHBOURS = {71, 74, 76, 77, 78, 80, 83, 90, 93}

# Contact data from Script 05_2
CONTACT_COUNTS = {
    "6W4H\nSARS-CoV-2\n(1.80Å)": {
        "H-Bond":49, "Hydrophobic":98,
        "Salt Bridge":2, "VdW":341},
    "6WVN\nSARS-CoV-1\n(1.95Å)": {
        "H-Bond":52, "Hydrophobic":103,
        "Salt Bridge":3, "VdW":359},
    "6WKQ\nSARS-CoV-2\n(2.05Å)": {
        "H-Bond":55, "Hydrophobic":108,
        "Salt Bridge":3, "VdW":371},
    "AF3 model\n(iptm=0.79)": {
        "H-Bond":57, "Hydrophobic":115,
        "Salt Bridge":4, "VdW":479},
}

SALT_BRIDGE_DATA = {
    "6W4H\nSARS-CoV-2\n(1.80Å)": [
        ("ASP106–LYS93", "3.90Å", True),
        ("ASP106–LYS95", "3.90Å", True),
    ],
    "6WVN\nSARS-CoV-1\n(1.95Å)": [
        ("ASP102–HIS80",  "3.78Å", True),
        ("ASP106–LYS93", "3.89Å", True),
        ("ASP106–LYS95", "3.95Å", True),
    ],
    "6WKQ\nSARS-CoV-2\n(2.05Å)": [
        ("ASP102–HIS80",  "3.85Å", True),
        ("ASP106–LYS93", "3.87Å", True),
        ("ASP106–LYS95", "3.86Å", True),
    ],
    "AF3 model\n(iptm=0.79)": [
        ("GLU6–LYS76",   "3.53Å", False),
        ("ASP102–HIS80",  "3.33Å", True),
        ("ASP106–LYS93", "3.90Å", True),
        ("ASP106–LYS95", "3.89Å", True),
    ],
}


def load_data():
    cons10 = pd.read_csv(RES_DIR / "conservation_NSP10.csv")
    cons16 = pd.read_csv(RES_DIR / "conservation_NSP16.csv")
    with open(RES_DIR / "interface_analysis.json") as f:
        interface = json.load(f)
    return cons10, cons16, interface


# ── Fig 1 — Conservation bars ──────────────────────────────
def fig_conservation_bars(cons10, cons16):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                             constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP16 Interface: Hotspot Conservation "
        "Across 5 Coronaviruses",
        fontsize=14, fontweight="bold", y=1.02)

    for ax, df, nsp, primary_set, zn_set in [
        (axes[0], cons10, "NSP10",
         PRIMARY_NSP10, ZN1_NEIGHBOURS),
        (axes[1], cons16, "NSP16",
         PRIMARY_NSP16, set()),
    ]:
        hdf = df[df["is_hotspot"]].copy()
        hdf = hdf.sort_values("conservation", ascending=True)

        colors = []
        for _, row in hdf.iterrows():
            if row["position"] in primary_set:
                colors.append(RED)
            elif int(row["position"]) in zn_set:
                colors.append(GOLD)
            elif row["conservation"] >= 0.8:
                colors.append(BLUE)
            else:
                colors.append(GREY)

        y = np.arange(len(hdf))
        ax.barh(y, hdf["conservation"], color=colors,
                edgecolor="white", linewidth=0.6, height=0.72)
        ax.axvline(0.8, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6)

        for i, (_, row) in enumerate(hdf.iterrows()):
            ax.text(row["conservation"] + 0.01, i,
                    f"{row['conservation']:.3f}",
                    va="center", fontsize=8.5)

        ax.set_yticks(y)
        ylabels = []
        for _, row in hdf.iterrows():
            three = AA3.get(row["aa_SARS2"], row["aa_SARS2"])
            pos   = int(row["position"])
            lbl   = f"{three}{pos}"
            if pos in primary_set:
                lbl = f"★  {lbl}"
            elif pos in zn_set and nsp == "NSP10":
                lbl = f"⬡  {lbl}"
            ylabels.append(lbl)
        ax.set_yticklabels(ylabels, fontsize=9.5)
        ax.set_xlim(0, 1.13)
        ax.set_xlabel(
            "Conservation Score  "
            "(1.0 = identical in all 5 coronaviruses)")
        n_cons = len(hdf[hdf["conservation"] >= 0.8])
        ax.set_title(
            f"{nsp}  ({n_cons}/{len(hdf)} residues ≥ 0.8)",
            pad=8)
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)

        patches = [
            mpatches.Patch(color=RED,
                label="Primary salt bridge (★)"),
        ]
        if nsp == "NSP10":
            patches.append(mpatches.Patch(
                color=GOLD,
                label="Near Zn1 finger (⬡)"))
        patches += [
            mpatches.Patch(color=BLUE,
                label="Conserved ≥ 0.8"),
            mpatches.Patch(color=GREY,
                label="Variable < 0.8"),
        ]
        ax.legend(handles=patches, loc="lower right",
                  fontsize=8.5)

    fig.text(0.5, -0.02,
             "★ = salt bridge anchors  |  "
             "⬡ = residues flanking Zn1 zinc finger  |  "
             "HIS80 sits inside Zn1 loop (CYS77–HIS83)",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig1_NSP10-NSP16_conservation_bars_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Fig 2 — Conservation heatmap ───────────────────────────
def fig_heatmap(cons10, cons16):
    HOTSPOTS = {
        "NSP10": sorted([5,40,42,43,44,45,71,76,78,
                         80,93,94,95,96]),
        "NSP16": sorted([40,41,44,76,83,87,102,104,
                         106,244,247]),
    }
    cmap = LinearSegmentedColormap.from_list(
        "cons",
        ["#FDFEFE","#AED6F1","#2980B9","#1A5276"], N=256)

    fig  = plt.figure(figsize=(20, 7),
                      constrained_layout=True)
    gs   = fig.add_gridspec(1, 3, width_ratios=[10, 8, 1.4])
    axes = [fig.add_subplot(gs[0]),
            fig.add_subplot(gs[1])]
    ax_leg = fig.add_subplot(gs[2])

    fig.suptitle(
        "Amino Acid Identity Across 5 Coronaviruses\n"
        "NSP10–NSP16 Conserved Hotspot Residues",
        fontsize=13, fontweight="bold")

    im_last = None
    for ax, nsp, df, hotspots in [
        (axes[0], "NSP10", cons10, HOTSPOTS["NSP10"]),
        (axes[1], "NSP16", cons16, HOTSPOTS["NSP16"]),
    ]:
        try:
            aln     = AlignIO.read(
                SEQ_DIR / f"{nsp}_aligned.fasta","fasta")
            ref_seq = next(r for r in aln
                           if "SARS-CoV-2" in r.id or
                           "SARS_CoV_2" in r.id)
        except Exception:
            print(f"  Warning: could not load {nsp} alignment")
            continue

        pos_map, ref_pos = {}, 0
        for i, aa in enumerate(str(ref_seq.seq)):
            if aa != "-":
                ref_pos += 1
                pos_map[ref_pos] = i

        matrix, row_labels = [], []
        for cov in CORONAVIRUSES:
            rec = next(
                (r for r in aln
                 if cov in r.id.replace("_"," ")),
                aln[0])
            matrix.append([
                str(rec.seq)[pos_map[p]]
                if p in pos_map else "-"
                for p in hotspots])
            row_labels.append(cov)

        score_map = dict(zip(df["position"],
                             df["conservation"]))
        scores = np.array([
            [score_map.get(p, 0) for p in hotspots]
            for _ in CORONAVIRUSES])

        im = ax.imshow(scores, cmap=cmap,
                       vmin=0, vmax=1, aspect="auto")
        im_last = im

        for i in range(len(CORONAVIRUSES)):
            for j in range(len(hotspots)):
                aa1 = matrix[i][j]
                tc  = ("white" if scores[i][j] > 0.65
                        else "black")
                ax.text(j, i, aa1, ha="center",
                        va="center", fontsize=11,
                        fontweight="bold", color=tc)

        # Primary target borders — red
        primary = (PRIMARY_NSP10 if nsp == "NSP10"
                   else PRIMARY_NSP16)
        for j, p in enumerate(hotspots):
            if p in primary:
                ax.add_patch(plt.Rectangle(
                    (j-0.5, -0.5), 1,
                    len(CORONAVIRUSES),
                    fill=False, edgecolor=RED,
                    linewidth=2.5, zorder=3))
                ax.text(j, -0.9, "★", ha="center",
                        va="center", fontsize=11,
                        color=RED)
            # Gold border for Zn1 neighbours (NSP10 only)
            if nsp == "NSP10" and p in ZN1_NEIGHBOURS:
                ax.add_patch(plt.Rectangle(
                    (j-0.5, -0.5), 1,
                    len(CORONAVIRUSES),
                    fill=False, edgecolor=GOLD,
                    linewidth=1.2, zorder=2,
                    linestyle="--"))

        # X-axis labels
        col_labels = []
        for p in hotspots:
            if p in pos_map:
                aa1   = str(ref_seq.seq)[pos_map[p]]
                three = AA3.get(aa1, aa1)
                col_labels.append(f"{three}{p}")
            else:
                col_labels.append(f"?{p}")
        ax.set_xticks(range(len(col_labels)))
        ax.set_xticklabels(col_labels, rotation=45,
                           ha="right", fontsize=8.5)
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(
            row_labels if ax == axes[0] else [],
            fontsize=10)
        ax.set_title(nsp, pad=18)
        ax.tick_params(length=0)

        for x in np.arange(-0.5, len(hotspots), 1):
            ax.axvline(x, color="white", linewidth=1.0)
        for y in np.arange(-0.5, len(CORONAVIRUSES), 1):
            ax.axhline(y, color="white", linewidth=1.0)

    if im_last is not None:
        cbar = fig.colorbar(im_last, ax=axes,
                            shrink=0.75, pad=0.015,
                            label="Conservation Score")
        cbar.ax.tick_params(labelsize=9)

    # Legend panel
    ax_leg.axis("off")
    ax_leg.text(0.05, 0.98, "AA Code Key",
                transform=ax_leg.transAxes,
                fontsize=9, fontweight="bold", va="top")
    entries = [
        ("A","ALA"),("D","ASP"),("F","PHE"),("G","GLY"),
        ("H","HIS"),("I","ILE"),("K","LYS"),("L","LEU"),
        ("M","MET"),("N","ASN"),("P","PRO"),("R","ARG"),
        ("S","SER"),("T","THR"),("V","VAL"),("Y","TYR"),
    ]
    for i, (code, name) in enumerate(entries):
        ax_leg.text(0.05, 0.92 - i*0.054,
                    f"{code}  =  {name}",
                    transform=ax_leg.transAxes,
                    fontsize=7.5, va="top",
                    color="#2C3E50",
                    fontfamily="monospace")

    fig.text(0.42, -0.03,
             "★ Red = primary salt bridge  |  "
             "⬡ Gold dashes = Zn1 finger region (NSP10)",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig2_NSP10-NSP16_conservation_heatmap_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Fig 3 — Contact types + salt bridges + Zn annotation ──
def fig_contacts():
    structs = list(CONTACT_COUNTS.keys())
    ctypes  = ["H-Bond","Hydrophobic","Salt Bridge","VdW"]
    colors  = [BLUE, ORANGE, RED, "#BDC3C7"]

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(16, 6),
        gridspec_kw={"width_ratios":[2.5,1]},
        constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP16 Interface: Contact Analysis",
        fontsize=13, fontweight="bold")

    # ── Left: grouped bar chart ────────────────────────────
    x       = np.arange(len(structs))
    width   = 0.19
    offsets = np.array([-0.285,-0.095,0.095,0.285])

    for ct, col, off in zip(ctypes, colors, offsets):
        vals = [CONTACT_COUNTS[s][ct] for s in structs]
        bars = ax.bar(x+off, vals, width, label=ct,
                      color=col, edgecolor="white",
                      linewidth=0.6, zorder=3)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x()+bar.get_width()/2,
                        bar.get_height()+3,
                        str(v), ha="center", va="bottom",
                        fontsize=8, fontweight="bold")

    totals = {s: sum(CONTACT_COUNTS[s].values())
              for s in structs}
    for i, s in enumerate(structs):
        ax.text(x[i], max(totals.values())*0.97,
                f"Total: {totals[s]}",
                ha="center", fontsize=9,
                color="#2C3E50", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(structs, fontsize=9.5)
    ax.set_ylabel("Number of contacts")
    ax.set_ylim(0, max(totals.values())*1.08)
    ax.grid(axis="y", alpha=0.25, linewidth=0.7, zorder=0)
    ax.legend(title="Contact type", loc="upper left",
              fontsize=9, title_fontsize=9)

    # Zn1 annotation on ax
    ax.annotate(
        "HIS80 inside\nZn1 finger loop",
        xy=(x[0]+offsets[2], CONTACT_COUNTS[structs[0]]
            ["Salt Bridge"]+2),
        xytext=(x[0]+offsets[2]-0.6,
                CONTACT_COUNTS[structs[0]]
                ["Salt Bridge"]+55),
        fontsize=8.5, color=GOLD, ha="center",
        arrowprops=dict(arrowstyle="->", color=GOLD,
                        lw=1.2))

    # ── Right: salt bridge presence table ─────────────────
    ax2.axis("off")
    ax2.set_title("Salt Bridge Inventory", pad=8,
                  fontsize=11, fontweight="bold")

    col_labels = ["Structure","Salt Bridge","Dist","✓"]
    all_rows   = []
    for s, sbs in SALT_BRIDGE_DATA.items():
        s_short = s.split("\n")[0]
        for name, dist, confirmed in sbs:
            all_rows.append([s_short, name, dist,
                             "✅" if confirmed else "⚠️"])

    y_pos = 0.95
    header_cols = ["Struct", "Salt Bridge", "Dist", ""]
    col_x = [0.0, 0.18, 0.72, 0.92]
    for i, h in enumerate(header_cols):
        ax2.text(col_x[i], y_pos, h,
                 transform=ax2.transAxes,
                 fontsize=8.5, fontweight="bold",
                 color="#2C3E50")
    y_pos -= 0.05
    ax2.axhline(y=y_pos+0.02,
                xmin=0, xmax=1,
                color="#BDC3C7", linewidth=0.8,
)

    prev_struct = ""
    for row in all_rows:
        struct, sb, dist, ok = row
        color = RED if ok == "✅" else ORANGE
        s_lbl = struct if struct != prev_struct else ""
        prev_struct = struct
        ax2.text(col_x[0], y_pos, s_lbl,
                 transform=ax2.transAxes,
                 fontsize=7.5, color="#555555")
        ax2.text(col_x[1], y_pos, sb,
                 transform=ax2.transAxes,
                 fontsize=7.5, color=color,
                 fontweight="bold")
        ax2.text(col_x[2], y_pos, dist,
                 transform=ax2.transAxes,
                 fontsize=7.5, color="#555555")
        ax2.text(col_x[3], y_pos, ok,
                 transform=ax2.transAxes,
                 fontsize=8)
        y_pos -= 0.072

    ax2.text(0.0, y_pos-0.04,
             "⚠️ = AF3-predicted only\n"
             "✅ = crystallographically confirmed",
             transform=ax2.transAxes,
             fontsize=7, color="#777777",
             style="italic")

    out = OUT_DIR / "Fig3_NSP10-NSP16_contact_types_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────
def main():
    print("\n" + "="*55)
    print("  Script 09_2: Publication Figures — NSP10-NSP16")
    print("="*55 + "\n")

    cons10, cons16, interface = load_data()

    print("  Fig 1 — Conservation bar charts...")
    fig_conservation_bars(cons10, cons16)

    print("  Fig 2 — Conservation heatmap...")
    fig_heatmap(cons10, cons16)

    print("  Fig 3 — Contact types + salt bridges...")
    fig_contacts()

    print(f"\n  All figures saved to: results/")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
