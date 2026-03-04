"""
Patch AA identity labels across all 5 complexes
- Updates residue_aa column in all ranking CSVs
- Regenerates Fig4 (BSA), Fig5 (AlaScan), Fig6 (Ranking)
- Updates notebook ranking table cells
No BSA/AlaScan recomputation needed.
"""

import json
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from Bio import PDB

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"

plt.rcParams.update({
    "font.family":"DejaVu Sans","font.size":11,
    "axes.spines.top":False,"axes.spines.right":False,
    "figure.dpi":150,"savefig.dpi":300,
    "savefig.bbox":"tight",
})
BLUE="#2C5F8A"; RED="#C0392B"; ORANGE="#E67E22"
GOLD="#F39C12"; GREEN="#27AE60"; GREY="#95A5A6"

AA3_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

# ── Complex definitions ──────────────────────────────────
COMPLEXES = {
    "NSP10-NSP16": {
        "csv":   "02-validation/NSP10-NSP16/composite_ranking_NSP10-NSP16_2.csv",
        "figs":  "02-validation/NSP10-NSP16",
        "nb":    "notebooks/NSP10-NSP16_3D_2.ipynb",
        "pdb":   PDB_DIR / "6W4H.pdb",
        "chains":{
            "A": {"label":"NSP16","offset":6797,
                  "shift":0},
            "B": {"label":"NSP10","offset":4270,
                  "shift":2},   # NSP10_SHIFT=2 from script
        },
        "suffix":"_2",
        "colors":{"NSP16":BLUE,"NSP10":GREEN},
        "primary_col":"primary_sb",
        "pancov_col":None,
    },
    "NSP10-NSP14": {
        "csv":   "02-validation/NSP10-NSP14/composite_ranking_NSP10-NSP14_2.csv",
        "figs":  "02-validation/NSP10-NSP14",
        "nb":    "notebooks/NSP10-NSP14_3D_2.ipynb",
        "pdb":   PDB_DIR / "7DIY.pdb",
        "chains":{
            "A": {"label":"NSP14","offset":5926,
                  "shift":0},
            "B": {"label":"NSP10","offset":4270,
                  "shift":2},
        },
        "suffix":"_2",
        "colors":{"NSP14":BLUE,"NSP10":GREEN},
        "primary_col":"primary_sb",
        "pancov_col":None,
    },
    "NSP12-NSP7": {
        "csv":   "02-validation/NSP12-NSP7/composite_ranking_NSP12-NSP7_3.csv",
        "figs":  "02-validation/NSP12-NSP7",
        "nb":    "notebooks/NSP12-NSP7_3D_3.ipynb",
        "pdb":   PDB_DIR / "7BV2_NSP12-NSP7.pdb",
        "chains":{
            "A": {"label":"NSP12","offset":0,"shift":0},
            "C": {"label":"NSP7", "offset":1,"shift":0},
        },
        "suffix":"_3",
        "colors":{"NSP12":BLUE,"NSP7":GREEN},
        "primary_col":"primary_sb",
        "pancov_col":None,
    },
    "NSP12-NSP8": {
        "csv":   "02-validation/NSP12-NSP8/composite_ranking_NSP12-NSP8_4.csv",
        "figs":  "02-validation/NSP12-NSP8",
        "nb":    "notebooks/NSP12-NSP8_3D_4.ipynb",
        "pdb":   PDB_DIR / "7BV2_NSP12-NSP8.pdb",
        "chains":{
            "A": {"label":"NSP12","offset":0,"shift":0},
            "B": {"label":"NSP8", "offset":77,"shift":0},
        },
        "suffix":"_4",
        "colors":{"NSP12":BLUE,"NSP8":GREEN},
        "primary_col":"primary_sb",
        "pancov_col":None,
    },
    "NSP9-NSP12": {
        "csv":   "02-validation/NSP9-NSP12/composite_ranking_NSP9-NSP12_5.csv",
        "figs":  "02-validation/NSP9-NSP12",
        "nb":    "notebooks/NSP9-NSP12_3D_5.ipynb",
        "pdb":   PDB_DIR / "8SQK_NSP9-NSP12.pdb",
        "chains":{
            "A": {"label":"NSP12","offset":0,"shift":0},
            "G": {"label":"NSP9", "offset":0,"shift":0},
        },
        "suffix":"_5",
        "colors":{"NSP12":BLUE,"NSP9":GREEN},
        "primary_col":"primary_sb",
        "pancov_col":"pan_cov",
    },
}


def build_aa_map(pdb_file, chains):
    """Build {chain_label: {local_pos: AA3}} from PDB."""
    parser = PDB.PDBParser(QUIET=True)
    s      = parser.get_structure("x", pdb_file)
    aa_map = {}
    for chain_id, cfg in chains.items():
        lbl   = cfg["label"]
        off   = cfg["offset"]
        shift = cfg["shift"]
        aa_map[lbl] = {}
        for m in s:
            for c in m:
                if c.id != chain_id:
                    continue
                for r in c:
                    if r.id[0] != " ":
                        continue
                    local = r.id[1] - off + shift
                    aa_map[lbl][local] = \
                        r.resname.strip()
    return aa_map


def get_aa_label(row, aa_map):
    chain = row["chain"]
    pos   = int(row["position"])
    aa3   = aa_map.get(chain,{}).get(pos,"???")
    return f"{chain}-{aa3}{pos}"


def add_aa_to_csv(csv_path, aa_map):
    df = pd.read_csv(csv_path)
    df["residue_aa"] = df.apply(
        lambda r: get_aa_label(r, aa_map), axis=1)
    df.to_csv(csv_path, index=False)
    return df


def get_color(row, cfg):
    primary_col = cfg["primary_col"]
    pancov_col  = cfg.get("pancov_col")
    if row.get(primary_col, False):
        return RED
    if pancov_col and row.get(pancov_col, False):
        return GOLD
    return cfg["colors"].get(row["chain"], GREY)


def fig4_bsa(df, cfg, out_dir, suffix):
    top    = df.nlargest(15,"bsa")
    colors = [get_color(r, cfg)
              for _,r in top.iterrows()]
    fig,ax = plt.subplots(figsize=(11,6))
    y = np.arange(len(top))
    ax.barh(y, top["bsa"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    for i,(_,row) in enumerate(top.iterrows()):
        ax.text(row["bsa"]+0.5,i,
                f"{row['bsa']:.1f} Å²",
                va="center",fontsize=8.5)
    ax.set_yticks(y)
    ylabels = []
    for _,row in top.iterrows():
        lbl = row["residue_aa"]
        if row.get(cfg["primary_col"],False):
            lbl = f"★  {lbl}"
        elif cfg.get("pancov_col") and \
                row.get(cfg["pancov_col"],False):
            lbl = f"◆  {lbl}"
        ylabels.append(lbl)
    ax.set_yticklabels(ylabels, fontsize=9.5)
    ax.set_xlabel("Buried Surface Area (Å²)")
    name = out_dir.name
    ax.set_title(f"{name}: Top 15 Residues by BSA",
                 fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    patches = [
        mpatches.Patch(color=RED,  label="Primary SB (★)"),
        mpatches.Patch(color=GOLD, label="Pan-cov (◆)"),
    ]
    for lbl,col in cfg["colors"].items():
        patches.append(mpatches.Patch(
            color=col, label=lbl))
    ax.legend(handles=patches, fontsize=9,
              loc="lower right")
    plt.tight_layout()
    out = out_dir/f"Fig4_{name}_BSA{suffix}.png"
    plt.savefig(out); plt.close()
    print(f"    Saved: {out.name}")


def fig5_alascan(df, cfg, out_dir, suffix):
    top = df[df["total_loss"]>0].nlargest(
        15,"total_loss")
    if top.empty:
        top = df.nlargest(15,"bsa")
    x   = np.arange(len(top))
    fig,ax = plt.subplots(figsize=(13,5))
    ax.bar(x-0.25, top["sb_loss"]*3, 0.25,
           label="Salt bridge (x3)",
           color=RED, alpha=0.85)
    ax.bar(x,      top["hb_loss"]*2, 0.25,
           label="H-bond (x2)",
           color=BLUE, alpha=0.85)
    ax.bar(x+0.25, top["hy_loss"],   0.25,
           label="Hydrophobic (x1)",
           color=ORANGE, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(top["residue_aa"],
                       rotation=45, ha="right",
                       fontsize=9)
    ax.set_ylabel("Estimated energy loss (AU)")
    name = out_dir.name
    ax.set_title(
        f"{name}: Alanine Scanning",
        fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    out = out_dir/f"Fig5_{name}_AlaScan{suffix}.png"
    plt.savefig(out); plt.close()
    print(f"    Saved: {out.name}")


def fig6_ranking(df, cfg, out_dir, suffix):
    top    = df.head(15)
    colors = [get_color(r,cfg)
              for _,r in top.iterrows()]
    fig,ax = plt.subplots(figsize=(11,6))
    y = np.arange(len(top))
    ax.barh(y, top["composite"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    for i,(_,row) in enumerate(top.iterrows()):
        ax.text(row["composite"]+0.005,i,
                f"{row['composite']:.4f}",
                va="center",fontsize=8.5)
    ax.set_yticks(y)
    ylabels = []
    for _,row in top.iterrows():
        lbl = row["residue_aa"]
        if row.get(cfg["primary_col"],False):
            lbl = f"★  {lbl}"
        elif cfg.get("pancov_col") and \
                row.get(cfg["pancov_col"],False):
            lbl = f"◆  {lbl}"
        ylabels.append(lbl)
    ax.set_yticklabels(ylabels, fontsize=9.5)
    ax.set_xlabel("Composite Drug Target Score")
    name = out_dir.name
    ax.set_title(
        f"{name}: Composite Hotspot Ranking\n"
        "BSA (35%) + Contact loss (35%) + "
        "Conservation (20%) + Primary SB (10%)",
        fontweight="bold")
    ax.set_xlim(0,1.12)
    ax.grid(axis="x",alpha=0.2)
    patches = [
        mpatches.Patch(color=RED,  label="Primary SB (★)"),
        mpatches.Patch(color=GOLD, label="Pan-cov (◆)"),
    ]
    for lbl,col in cfg["colors"].items():
        patches.append(mpatches.Patch(
            color=col,label=lbl))
    ax.legend(handles=patches,fontsize=9,
              loc="lower right")
    plt.tight_layout()
    out = out_dir / \
          f"Fig6_{name}_composite_ranking{suffix}.png"
    plt.savefig(out); plt.close()
    print(f"    Saved: {out.name}")


def patch_notebook(nb_path, df):
    """Replace ranking table cell in notebook."""
    if not nb_path.exists():
        print(f"    Notebook not found: {nb_path.name}")
        return
    nb = json.load(open(nb_path))
    table_lines = [
        "| Rank | Residue | BSA | Loss | Cons "
        "| Score | Flags |\n",
        "|------|---------|-----|------|------"
        "|-------|-------|\n",
    ]
    for rank,(_,row) in enumerate(
            df.head(15).iterrows(),1):
        flags = ""
        if row.get("primary_sb",False): flags+="★SB "
        if row.get("pan_cov",   False): flags+="◆PC "
        table_lines.append(
            f"| {rank} | {row['residue_aa']} "
            f"| {row['bsa']:.1f} "
            f"| {int(row['total_loss'])} "
            f"| {row['conservation']:.3f} "
            f"| {row['composite']:.4f} "
            f"| {flags.strip()} |\n")
    new_src = "".join(table_lines)

    patched = False
    for cell in nb["cells"]:
        src = "".join(cell.get("source",""))
        if ("Rank" in src and "Residue" in src and
                "Score" in src and
                cell["cell_type"]=="markdown"):
            cell["source"] = new_src
            patched = True
            break
    if patched:
        json.dump(nb, open(nb_path,"w"), indent=2)
        print(f"    Notebook patched: {nb_path.name}")
    else:
        print(f"    No ranking table found in: "
              f"{nb_path.name}")


def main():
    print("\n" + "="*60)
    print("  Patch: AA labels for all 5 complexes")
    print("="*60)

    for name, cfg in COMPLEXES.items():
        print(f"\n  ── {name} ──")
        csv_path = PROJECT / cfg["csv"]
        out_dir  = PROJECT / cfg["figs"]
        nb_path  = PROJECT / cfg["nb"]
        suffix   = cfg["suffix"]

        if not csv_path.exists():
            print(f"    CSV not found — skipping")
            continue
        if not cfg["pdb"].exists():
            print(f"    PDB not found — skipping")
            continue

        # Build AA map
        aa_map = build_aa_map(cfg["pdb"],
                               cfg["chains"])
        total  = sum(len(v) for v in aa_map.values())
        print(f"    AA map: {total} residues loaded")

        # Patch CSV
        df = add_aa_to_csv(csv_path, aa_map)
        print(f"    CSV updated: {csv_path.name}")

        # Print top 5
        print(f"    Top 5:")
        for rank,(_,row) in enumerate(
                df.head(5).iterrows(),1):
            flags = ""
            if row.get("primary_sb",False):
                flags+="★SB "
            if row.get("pan_cov",False):
                flags+="◆PC "
            print(f"      {rank}. {row['residue_aa']:<22}"
                  f" score={row['composite']:.4f} "
                  f"{flags}")

        # Regenerate Figs 4-6
        fig4_bsa(df, cfg, out_dir, suffix)
        fig5_alascan(df, cfg, out_dir, suffix)
        fig6_ranking(df, cfg, out_dir, suffix)

        # Patch notebook
        patch_notebook(nb_path, df)

    print("\n" + "="*60)
    print("  All complexes patched ✅")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
