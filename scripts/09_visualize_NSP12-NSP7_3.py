"""
Script 09_3: Publication Figures — NSP12-NSP7
=============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Fig 1 — Conservation bar charts (NSP12 + NSP7)
Fig 2 — Conservation heatmap (5 coronaviruses × hotspots)
Fig 3 — Contact types + salt bridge inventory (all 4 structures)

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/09_visualize_NSP12-NSP7_3.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from Bio import AlignIO
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
RES_DIR = PROJECT / "02-validation" / "NSP12-NSP7"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
OUT_DIR = PROJECT / "02-validation" / "NSP12-NSP7"
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

CORONAVIRUSES   = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                   "HCoV-229E","HCoV-NL63"]
PRIMARY_NSP12   = {431}
PRIMARY_NSP7    = {1}
HYDROPHOBIC_CORE= {412,413,415,420,440,442,843}

HOTSPOTS_NSP12  = [440,412,442,443,420,843,409,
                    40,33,41,37,413,415,14,23]
HOTSPOTS_NSP7   = [40,14,33,41,37,11,23,5,15,29,12,4,1]

# Contact data from Script 05_3 (4 structures)
CONTACT_COUNTS = {
    "7BV2\nSARS-CoV-2\n(2.90Å)": {
        "H-Bond":6,  "Hydrophobic":8,  "Salt Bridge":1},
    "6NUR\nSARS-CoV-2\n(3.10Å)": {
        "H-Bond":4,  "Hydrophobic":10, "Salt Bridge":0},
    "7C2K\nSARS-CoV-2\n(2.50Å)": {
        "H-Bond":6,  "Hydrophobic":11, "Salt Bridge":0},
    "AF3 model\n(iptm=0.81)": {
        "H-Bond":6,  "Hydrophobic":10, "Salt Bridge":1},
}

SALT_BRIDGE_DATA = {
    "7BV2\nSARS-CoV-2\n(2.90Å)": [
        ("GLU431–LYS1(NSP7)", "4.24Å", True),
    ],
    "6NUR\nSARS-CoV-2\n(3.10Å)": [],
    "7C2K\nSARS-CoV-2\n(2.50Å)": [],
    "AF3 model\n(iptm=0.81)": [
        ("GLU431–LYS1(NSP7)", "4.26Å", False),
    ],
}


def load_data():
    cons12 = pd.read_csv(RES_DIR / "conservation_NSP12.csv")
    cons7  = pd.read_csv(RES_DIR / "conservation_NSP7.csv")
    with open(RES_DIR / "interface_analysis_3.json") as f:
        interface = json.load(f)
    return cons12, cons7, interface


# ── Fig 1 — Conservation bars ──────────────────────────────
def fig_conservation_bars(cons12, cons7):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                             constrained_layout=True)
    fig.suptitle(
        "NSP12–NSP7 Interface: Hotspot Conservation "
        "Across 5 Coronaviruses",
        fontsize=14, fontweight="bold", y=1.02)

    for ax, df, nsp, primary_set, hcore in [
        (axes[0], cons12, "NSP12", PRIMARY_NSP12,
         HYDROPHOBIC_CORE),
        (axes[1], cons7,  "NSP7",  PRIMARY_NSP7,  set()),
    ]:
        hdf = df[df["is_hotspot"]].copy()
        hdf = hdf.sort_values("conservation", ascending=True)

        colors = []
        for _, row in hdf.iterrows():
            if row["position"] in primary_set:
                colors.append(RED)
            elif int(row["position"]) in hcore:
                colors.append(ORANGE)
            elif row["conservation"] >= 0.8:
                colors.append(BLUE if nsp=="NSP12"
                               else GREEN)
            else:
                colors.append(GREY)

        y = np.arange(len(hdf))
        ax.barh(y, hdf["conservation"], color=colors,
                edgecolor="white", linewidth=0.6,
                height=0.72)
        ax.axvline(0.8, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6)

        for i, (_, row) in enumerate(hdf.iterrows()):
            ax.text(row["conservation"]+0.01, i,
                    f"{row['conservation']:.3f}",
                    va="center", fontsize=8.5)

        ax.set_yticks(y)
        ylabels = []
        for _, row in hdf.iterrows():
            three = AA3.get(row["aa_SARS2"],
                            row["aa_SARS2"])
            pos   = int(row["position"])
            lbl   = f"{three}{pos}"
            if pos in primary_set:
                lbl = f"★  {lbl}"
            elif pos in hcore and nsp=="NSP12":
                lbl = f"◆  {lbl}"
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
        if nsp == "NSP12":
            patches.append(mpatches.Patch(
                color=ORANGE,
                label="Hydrophobic core (◆)"))
        patches += [
            mpatches.Patch(
                color=BLUE if nsp=="NSP12" else GREEN,
                label="Conserved ≥ 0.8"),
            mpatches.Patch(color=GREY,
                label="Variable < 0.8"),
        ]
        ax.legend(handles=patches, loc="lower right",
                  fontsize=8.5)

    fig.text(
        0.5, -0.02,
        "★ = GLU431(NSP12)–LYS1(NSP7) salt bridge  |  "
        "◆ = hydrophobic core residues (PHE412/415/420/"
        "440/442, PHE843)  |  "
        "Note: LYS1(NSP7) SARS-CoV-1/2 only",
        ha="center", fontsize=9,
        style="italic", color="#555555")

    out = OUT_DIR / "Fig1_NSP12-NSP7_conservation_bars_3.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Fig 2 — Conservation heatmap ───────────────────────────
def fig_heatmap(cons12, cons7):
    HOTSPOTS = {
        "NSP12": sorted(HOTSPOTS_NSP12),
        "NSP7":  sorted(HOTSPOTS_NSP7),
    }
    cmap = LinearSegmentedColormap.from_list(
        "cons",
        ["#FDFEFE","#AED6F1","#2980B9","#1A5276"], N=256)

    fig  = plt.figure(figsize=(22, 7),
                      constrained_layout=True)
    gs   = fig.add_gridspec(1, 3,
                             width_ratios=[12, 6, 1.4])
    axes = [fig.add_subplot(gs[0]),
            fig.add_subplot(gs[1])]
    ax_leg = fig.add_subplot(gs[2])

    fig.suptitle(
        "Amino Acid Identity Across 5 Coronaviruses\n"
        "NSP12–NSP7 Conserved Hotspot Residues",
        fontsize=13, fontweight="bold")

    im_last = None
    for ax, nsp, df, hotspots in [
        (axes[0], "NSP12", cons12, HOTSPOTS["NSP12"]),
        (axes[1], "NSP7",  cons7,  HOTSPOTS["NSP7"]),
    ]:
        try:
            aln     = AlignIO.read(
                SEQ_DIR / f"{nsp}_aligned.fasta","fasta")
            ref_seq = next(
                r for r in aln
                if "SARS-CoV-2" in r.id.replace("_"," "))
        except Exception as e:
            print(f"  Warning: {nsp} alignment: {e}")
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
                tc  = ("white" if scores[i,j] > 0.65
                        else "black")
                ax.text(j, i, aa1, ha="center",
                        va="center", fontsize=10,
                        fontweight="bold", color=tc)

        primary = (PRIMARY_NSP12 if nsp=="NSP12"
                   else PRIMARY_NSP7)
        hcore   = (HYDROPHOBIC_CORE if nsp=="NSP12"
                   else set())
        for j, p in enumerate(hotspots):
            if p in primary:
                ax.add_patch(plt.Rectangle(
                    (j-0.5,-0.5), 1,
                    len(CORONAVIRUSES),
                    fill=False, edgecolor=RED,
                    linewidth=2.5, zorder=3))
                ax.text(j, -0.9, "★", ha="center",
                        va="center", fontsize=11,
                        color=RED)
            if p in hcore:
                ax.add_patch(plt.Rectangle(
                    (j-0.5,-0.5), 1,
                    len(CORONAVIRUSES),
                    fill=False, edgecolor=ORANGE,
                    linewidth=1.2, zorder=2,
                    linestyle="--"))

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
            row_labels if ax==axes[0] else [],
            fontsize=10)
        ax.set_title(nsp, pad=18)
        ax.tick_params(length=0)
        for x in np.arange(-0.5,len(hotspots),1):
            ax.axvline(x,color="white",linewidth=1.0)
        for y in np.arange(-0.5,len(CORONAVIRUSES),1):
            ax.axhline(y,color="white",linewidth=1.0)

    if im_last is not None:
        cbar = fig.colorbar(im_last, ax=axes,
                            shrink=0.75, pad=0.015,
                            label="Conservation Score")
        cbar.ax.tick_params(labelsize=9)

    ax_leg.axis("off")
    ax_leg.text(0.05, 0.98, "AA Code Key",
                transform=ax_leg.transAxes,
                fontsize=9, fontweight="bold", va="top")
    entries = [
        ("A","ALA"),("C","CYS"),("D","ASP"),("E","GLU"),
        ("F","PHE"),("G","GLY"),("H","HIS"),("I","ILE"),
        ("K","LYS"),("L","LEU"),("M","MET"),("N","ASN"),
        ("P","PRO"),("R","ARG"),("S","SER"),("T","THR"),
        ("V","VAL"),("Y","TYR"),
    ]
    for i,(code,name) in enumerate(entries):
        ax_leg.text(0.05, 0.92-i*0.048,
                    f"{code}  =  {name}",
                    transform=ax_leg.transAxes,
                    fontsize=7.5, va="top",
                    color="#2C3E50",
                    fontfamily="monospace")

    fig.text(0.42, -0.03,
             "★ Red = GLU431–LYS1 salt bridge  |  "
             "◆ Orange dashes = hydrophobic core  |  "
             "Note: LYS1(NSP7) absent in MERS/HCoV",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig2_NSP12-NSP7_conservation_heatmap_3.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Fig 3 — Contact types + salt bridge inventory ─────────
def fig_contacts():
    structs = list(CONTACT_COUNTS.keys())
    ctypes  = ["H-Bond","Hydrophobic","Salt Bridge"]
    colors  = [BLUE, ORANGE, RED]

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(16, 6),
        gridspec_kw={"width_ratios":[2.5,1]},
        constrained_layout=True)
    fig.suptitle(
        "NSP12–NSP7 Interface: Contact Analysis "
        "(All 4 Structures)",
        fontsize=13, fontweight="bold")

    x       = np.arange(len(structs))
    width   = 0.22
    offsets = np.array([-0.22, 0.0, 0.22])

    for ct, col, off in zip(ctypes, colors, offsets):
        vals = [CONTACT_COUNTS[s][ct] for s in structs]
        bars = ax.bar(x+off, vals, width, label=ct,
                      color=col, edgecolor="white",
                      linewidth=0.6, zorder=3)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(
                    bar.get_x()+bar.get_width()/2,
                    bar.get_height()+0.15,
                    str(v), ha="center", va="bottom",
                    fontsize=9, fontweight="bold")

    totals = {s: sum(CONTACT_COUNTS[s].values())
              for s in structs}
    for i, s in enumerate(structs):
        ax.text(x[i],
                max(totals.values())*0.97,
                f"Total: {totals[s]}",
                ha="center", fontsize=9,
                color="#2C3E50", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(structs, fontsize=9.5)
    ax.set_ylabel("Number of contacts")
    ax.set_ylim(0, max(totals.values())*1.1)
    ax.grid(axis="y", alpha=0.25,
            linewidth=0.7, zorder=0)
    ax.legend(title="Contact type", loc="upper left",
              fontsize=9, title_fontsize=9)

    ax.annotate(
        "Hydrophobic-dominated\ninterface",
        xy=(x[0]+offsets[1],
            CONTACT_COUNTS[structs[0]]["Hydrophobic"]+0.2),
        xytext=(x[0]+offsets[1]-0.5,
                CONTACT_COUNTS[structs[0]]
                ["Hydrophobic"]+5),
        fontsize=8.5, color=ORANGE, ha="center",
        arrowprops=dict(arrowstyle="->",
                        color=ORANGE, lw=1.2))

    # Salt bridge table
    ax2.axis("off")
    ax2.set_title("Salt Bridge Inventory", pad=8,
                  fontsize=11, fontweight="bold")

    col_x = [0.0, 0.20, 0.72, 0.92]
    y_pos = 0.95
    for i, h in enumerate(["Struct","Salt Bridge",
                            "Dist",""]):
        ax2.text(col_x[i], y_pos, h,
                 transform=ax2.transAxes,
                 fontsize=8.5, fontweight="bold",
                 color="#2C3E50")
    y_pos -= 0.05
    ax2.axhline(y=y_pos+0.02, xmin=0, xmax=1,
                color="#BDC3C7", linewidth=0.8)

    prev_struct = ""
    has_any = False
    for s, sbs in SALT_BRIDGE_DATA.items():
        s_short = s.split("\n")[0]
        if not sbs:
            continue
        for name, dist, confirmed in sbs:
            has_any = True
            color   = RED if confirmed else ORANGE
            s_lbl   = (s_short
                       if s_short != prev_struct else "")
            prev_struct = s_short
            ax2.text(col_x[0], y_pos, s_lbl,
                     transform=ax2.transAxes,
                     fontsize=7.5, color="#555555")
            ax2.text(col_x[1], y_pos, name,
                     transform=ax2.transAxes,
                     fontsize=7.5, color=color,
                     fontweight="bold")
            ax2.text(col_x[2], y_pos, dist,
                     transform=ax2.transAxes,
                     fontsize=7.5, color="#555555")
            ax2.text(col_x[3], y_pos,
                     "✅" if confirmed else "⚠️",
                     transform=ax2.transAxes,
                     fontsize=8)
            y_pos -= 0.09

    if not has_any:
        ax2.text(0.1, 0.7, "No salt bridges\ndetected",
                 transform=ax2.transAxes,
                 fontsize=9, color=GREY,
                 style="italic")

    ax2.text(0.0, y_pos-0.06,
             "⚠️ = AF3-predicted only\n"
             "✅ = crystallographically confirmed\n"
             "Note: 6NUR and 7C2K: no SB detected\n"
             "(resolution/conformation dependent)",
             transform=ax2.transAxes,
             fontsize=7, color="#777777",
             style="italic")

    out = OUT_DIR / "Fig3_NSP12-NSP7_contact_types_3.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────
def main():
    print("\n" + "="*55)
    print("  Script 09_3: Publication Figures — NSP12-NSP7")
    print("="*55 + "\n")

    cons12, cons7, interface = load_data()

    print("  Fig 1 — Conservation bar charts...")
    fig_conservation_bars(cons12, cons7)

    print("  Fig 2 — Conservation heatmap...")
    fig_heatmap(cons12, cons7)

    print("  Fig 3 — Contact types + salt bridges...")
    fig_contacts()

    print(f"\n  All figures → 02-validation/NSP12-NSP7/")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
