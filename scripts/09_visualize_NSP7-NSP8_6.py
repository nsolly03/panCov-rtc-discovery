"""
Script 09_6: Publication Figures — NSP7-NSP8
=============================================
Fig 1 — Conservation bars (NSP8 + NSP7, dual mode)
Fig 2 — Conservation heatmap (5 coronaviruses)
Fig 3 — Contact types + dual mode comparison
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
RES_DIR = PROJECT / "02-validation" / "NSP7-NSP8"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"

plt.rcParams.update({
    "font.family":"DejaVu Sans","font.size":11,
    "axes.spines.top":False,"axes.spines.right":False,
    "axes.linewidth":0.8,"figure.dpi":150,
    "savefig.dpi":300,"savefig.bbox":"tight",
})
BLUE="#2C5F8A"; GREEN="#27AE60"; RED="#C0392B"
ORANGE="#E67E22"; GOLD="#F39C12"; GREY="#95A5A6"
PURPLE="#7D3C98"; TEAL="#148F77"

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                  "HCoV-229E","HCoV-NL63"]
AA3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS",
       "Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE",
       "L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
       "S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}

CONTACT_COUNTS = {
    "7BV2\ncrystal\nMode A": {"H-Bond":2,"Hydrophobic":0,"Salt Bridge":0},
    "6NUR\ncrystal\nMode A": {"H-Bond":4,"Hydrophobic":0,"Salt Bridge":0},
    "AF3\nMode B\n(iptm=0.87)":{"H-Bond":8,"Hydrophobic":34,"Salt Bridge":1},
}
SB_DATA = {
    "7BV2\ncrystal\nMode A": [],
    "6NUR\ncrystal\nMode A": [],
    "AF3\nMode B\n(iptm=0.87)": [("ARG190–GLU50","4.36A",False)],
}

def fig1_conservation_bars():
    df8 = pd.read_csv(RES_DIR/"conservation_NSP8.csv")
    df7 = pd.read_csv(RES_DIR/"conservation_NSP7.csv")
    fig,axes = plt.subplots(1,2,figsize=(18,8),
                             constrained_layout=True)
    fig.suptitle(
        "NSP7–NSP8 Interface: Conservation (Dual Binding Mode)\n"
        "Mode A = Crystal | Mode B = AF3 physiological",
        fontsize=13, fontweight="bold")

    for ax, df, nsp, primary in [
        (axes[0], df8, "NSP8", {190}),
        (axes[1], df7, "NSP7", {50}),
    ]:
        df = df[df["is_hotspot"]].copy()
        df = df.sort_values("conservation", ascending=True)
        colors = []
        for _,row in df.iterrows():
            if row["position"] in primary:
                colors.append(RED)
            elif row.get("mode","B") == "A":
                colors.append(PURPLE)
            elif row["conservation"] >= 0.8:
                colors.append(BLUE if nsp=="NSP8"
                               else GREEN)
            else:
                colors.append(GREY)

        y = np.arange(len(df))
        ax.barh(y, df["conservation"], color=colors,
                edgecolor="white", linewidth=0.5,
                height=0.72)
        ax.axvline(0.8, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6)
        for i,(_,row) in enumerate(df.iterrows()):
            ax.text(row["conservation"]+0.01, i,
                    f"{row['conservation']:.3f}",
                    va="center", fontsize=8)

        ax.set_yticks(y)
        ylabels = []
        for _,row in df.iterrows():
            three = AA3.get(str(row["aa_SARS2"]),
                             str(row["aa_SARS2"]))
            pos   = int(row["position"])
            mode  = row.get("mode","B")
            sb    = row["position"] in primary
            lbl   = f"[{mode}] {three}{pos}"
            if sb: lbl = f"★  {lbl}"
            ylabels.append(lbl)
        ax.set_yticklabels(ylabels, fontsize=8.5)
        ax.set_xlim(0,1.15)
        ax.set_xlabel("Conservation Score")
        n_c = len(df[df["conservation"]>=0.8])
        ax.set_title(f"{nsp}  ({n_c}/{len(df)} ≥ 0.8)",
                     pad=8)
        ax.grid(axis="x", alpha=0.2)
        patches = [
            mpatches.Patch(color=RED,    label="★ AF3 SB anchor"),
            mpatches.Patch(color=PURPLE, label="Mode A (crystal)"),
            mpatches.Patch(color=BLUE if nsp=="NSP8" else GREEN,
                           label="Mode B cons ≥ 0.8"),
            mpatches.Patch(color=GREY,   label="Variable < 0.8"),
        ]
        ax.legend(handles=patches, fontsize=8.5,
                  loc="lower right")

    fig.text(0.5,-0.01,
             "★ = ARG190(NSP8)–GLU50(NSP7) SB (AF3 Mode B) | "
             "Mode A crystal interface undruggable (fpocket=0.276) | "
             "Mode B AF3 primary target (fpocket=0.531)",
             ha="center", fontsize=8.5,
             style="italic", color="#555555")
    out = RES_DIR/"Fig1_NSP7-NSP8_conservation_bars_6.png"
    fig.savefig(out); plt.close()
    print(f"  Saved: {out.name}")


def fig2_heatmap():
    cmap = LinearSegmentedColormap.from_list(
        "cons",["#FDFEFE","#AED6F1","#2980B9","#1A5276"],N=256)
    fig  = plt.figure(figsize=(24,7),constrained_layout=True)
    gs   = fig.add_gridspec(1,3,width_ratios=[13,9,1.4])
    axes = [fig.add_subplot(gs[0]),
            fig.add_subplot(gs[1])]
    ax_leg = fig.add_subplot(gs[2])
    fig.suptitle(
        "NSP7–NSP8: Amino Acid Identity Across 5 Coronaviruses\n"
        "(Mode A = Crystal | Mode B = AF3)",
        fontsize=13, fontweight="bold")

    for ax, nsp, csv_f, hotspots, primary in [
        (axes[0], "NSP8", "conservation_NSP8.csv",
         [87,91,92,94,98,110,111,116,120,150,163,178,179,180,190],
         {190}),
        (axes[1], "NSP7", "conservation_NSP7.csv",
         [2,6,16,24,26,27,49,50,56],
         {50}),
    ]:
        df = pd.read_csv(RES_DIR/csv_f)
        try:
            aln     = AlignIO.read(
                SEQ_DIR/f"{nsp}_aligned.fasta","fasta")
            ref_seq = next(r for r in aln
                if "SARS-CoV-2" in r.id.replace("_"," "))
        except Exception as e:
            print(f"  Warning {nsp}: {e}"); continue

        pos_map, ref_pos = {}, 0
        for i,aa in enumerate(str(ref_seq.seq)):
            if aa != "-":
                ref_pos += 1
                pos_map[ref_pos] = i

        hotspots = [p for p in hotspots if p in pos_map]
        score_map = dict(zip(df["position"],
                              df["conservation"]))
        matrix, row_labels = [], []
        for cov in CORONAVIRUSES:
            rec = next((r for r in aln
                if cov in r.id.replace("_"," ")), aln[0])
            matrix.append([str(rec.seq)[pos_map[p]]
                            if p in pos_map else "-"
                            for p in hotspots])
            row_labels.append(cov)

        scores = np.array([[score_map.get(p,0)
                             for p in hotspots]
                            for _ in CORONAVIRUSES])
        im = ax.imshow(scores, cmap=cmap,
                       vmin=0, vmax=1, aspect="auto")
        for i in range(len(CORONAVIRUSES)):
            for j in range(len(hotspots)):
                tc = "white" if scores[i,j]>0.65 else "black"
                ax.text(j,i,matrix[i][j],ha="center",
                        va="center",fontsize=10,
                        fontweight="bold",color=tc)
        for j,p in enumerate(hotspots):
            if p in primary:
                ax.add_patch(plt.Rectangle(
                    (j-0.5,-0.5),1,len(CORONAVIRUSES),
                    fill=False,edgecolor=RED,
                    linewidth=2.5,zorder=3))
        mode_row = df.set_index("position")["mode"].to_dict() \
            if "mode" in df.columns else {}
        col_labels = []
        for p in hotspots:
            aa1   = str(ref_seq.seq)[pos_map[p]] \
                    if p in pos_map else "?"
            three = AA3.get(aa1,aa1)
            mode  = mode_row.get(p,"B")
            col_labels.append(f"[{mode}]{three}{p}")
        ax.set_xticks(range(len(col_labels)))
        ax.set_xticklabels(col_labels,rotation=45,
                           ha="right",fontsize=8)
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(
            row_labels if ax==axes[0] else [],fontsize=10)
        ax.set_title(nsp,pad=18)
        ax.tick_params(length=0)
        for x in np.arange(-0.5,len(hotspots),1):
            ax.axvline(x,color="white",linewidth=1.0)
        for y in np.arange(-0.5,len(CORONAVIRUSES),1):
            ax.axhline(y,color="white",linewidth=1.0)
        fig.colorbar(im,ax=ax,shrink=0.6,pad=0.015,
                     label="Conservation Score")

    ax_leg.axis("off")
    ax_leg.text(0.05,0.98,"AA Code",
                transform=ax_leg.transAxes,
                fontsize=9,fontweight="bold",va="top")
    for i,(c,n) in enumerate([
            ("A","ALA"),("D","ASP"),("E","GLU"),
            ("F","PHE"),("I","ILE"),("K","LYS"),
            ("L","LEU"),("M","MET"),("N","ASN"),
            ("P","PRO"),("Q","GLN"),("R","ARG"),
            ("S","SER"),("T","THR"),("V","VAL")]):
        ax_leg.text(0.05,0.92-i*0.056,f"{c} = {n}",
                    transform=ax_leg.transAxes,
                    fontsize=7.5,fontfamily="monospace")
    out = RES_DIR/"Fig2_NSP7-NSP8_conservation_heatmap_6.png"
    fig.savefig(out); plt.close()
    print(f"  Saved: {out.name}")


def fig3_contacts():
    structs = list(CONTACT_COUNTS.keys())
    ctypes  = ["H-Bond","Hydrophobic","Salt Bridge"]
    colors  = [BLUE,ORANGE,RED]
    fig,(ax,ax2) = plt.subplots(
        1,2,figsize=(14,6),
        gridspec_kw={"width_ratios":[2,1]},
        constrained_layout=True)
    fig.suptitle(
        "NSP7–NSP8: Contact Analysis\n"
        "Mode A (crystal) vs Mode B (AF3)",
        fontsize=13,fontweight="bold")

    x = np.arange(len(structs))
    for ct,col,off in zip(ctypes,colors,
                           [-0.22,0.0,0.22]):
        vals = [CONTACT_COUNTS[s][ct] for s in structs]
        bars = ax.bar(x+off,vals,0.22,label=ct,
                      color=col,edgecolor="white",
                      linewidth=0.6,zorder=3)
        for bar,v in zip(bars,vals):
            if v>0:
                ax.text(bar.get_x()+bar.get_width()/2,
                        bar.get_height()+0.3,str(v),
                        ha="center",va="bottom",
                        fontsize=9,fontweight="bold")
    totals = {s:sum(CONTACT_COUNTS[s].values())
              for s in structs}
    for i,s in enumerate(structs):
        ax.text(x[i],max(totals.values())*0.97,
                f"Total:{totals[s]}",ha="center",
                fontsize=9,fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(structs,fontsize=9.5)
    ax.set_ylabel("Number of contacts")
    ax.set_ylim(0,max(totals.values())*1.15)
    ax.grid(axis="y",alpha=0.25,zorder=0)
    ax.legend(fontsize=9,loc="upper left")
    ax.annotate("Mode B:\nphysiological\ninterface",
                xy=(x[2],totals[structs[2]]*0.85),
                xytext=(x[2]+0.5,totals[structs[2]]*0.7),
                fontsize=8.5,color=PURPLE,ha="center",
                arrowprops=dict(arrowstyle="->",
                                color=PURPLE,lw=1.2))

    ax2.axis("off")
    ax2.set_title("Salt Bridge Inventory",
                  pad=8,fontsize=11,fontweight="bold")
    y_pos = 0.92
    ax2.text(0.0,y_pos,"Struct",
             transform=ax2.transAxes,
             fontsize=8.5,fontweight="bold")
    ax2.text(0.22,y_pos,"Salt Bridge",
             transform=ax2.transAxes,
             fontsize=8.5,fontweight="bold")
    ax2.text(0.72,y_pos,"Dist",
             transform=ax2.transAxes,
             fontsize=8.5,fontweight="bold")
    y_pos -= 0.06
    ax2.axhline(y=y_pos+0.02,xmin=0,xmax=1,
                color="#BDC3C7",linewidth=0.8)
    prev = ""
    for s,sbs in SB_DATA.items():
        sl = s.split("\n")[0]
        if not sbs:
            lbl = sl if sl!=prev else ""
            prev = sl
            ax2.text(0.0,y_pos,lbl,
                     transform=ax2.transAxes,
                     fontsize=7,color="#555")
            ax2.text(0.22,y_pos,"None",
                     transform=ax2.transAxes,
                     fontsize=7,color=GREY,
                     style="italic")
            y_pos -= 0.09
            continue
        for name,dist,conf in sbs:
            lbl = sl if sl!=prev else ""
            prev = sl
            ax2.text(0.0,y_pos,lbl,
                     transform=ax2.transAxes,
                     fontsize=7,color="#555")
            ax2.text(0.22,y_pos,name,
                     transform=ax2.transAxes,
                     fontsize=7,color=ORANGE,
                     fontweight="bold")
            ax2.text(0.72,y_pos,dist,
                     transform=ax2.transAxes,
                     fontsize=7,color="#555")
            y_pos -= 0.09
    ax2.text(0.0,y_pos-0.04,
             "~ = AF3 Mode B only\n"
             "ARG190 cons=1.000 pan-cov ✅\n"
             "GLU50  cons=0.410 SARS-only ⚠️",
             transform=ax2.transAxes,
             fontsize=7,color="#777",style="italic")
    out = RES_DIR/"Fig3_NSP7-NSP8_contact_types_6.png"
    fig.savefig(out); plt.close()
    print(f"  Saved: {out.name}")


print("\n"+"="*55)
print("  Script 09_6: Publication Figures — NSP7-NSP8")
print("="*55+"\n")
fig1_conservation_bars()
fig2_heatmap()
fig3_contacts()
print(f"\n  All figures -> 02-validation/NSP7-NSP8/")
print("="*55+"\n")
