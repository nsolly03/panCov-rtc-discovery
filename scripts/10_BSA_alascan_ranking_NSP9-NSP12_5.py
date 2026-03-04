"""
Script 10_5: BSA + Alanine Scanning + Composite Ranking — NSP9-NSP12
=====================================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Primary PDB: 8SQK (Chain A=NSP12, Chain G=NSP9)
Output: Figs 4-6, composite_ranking CSV, BSA/AlaScan JSON

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/10_BSA_alascan_ranking_NSP9-NSP12_5.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley
from pathlib import Path

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
RES_DIR  = PROJECT / "02-validation" / "NSP9-NSP12"
OUT_DIR  = PROJECT / "02-validation" / "NSP9-NSP12"
OUT_DIR.mkdir(exist_ok=True)

PDB_FILE       = PDB_DIR / "8SQK_NSP9-NSP12.pdb"
CHAIN_NSP12    = "A"
CHAIN_NSP9     = "G"

HOTSPOTS_NSP12 = [38,1,3,4,96,733,202,103,
                   221,233,291,2,223]
HOTSPOTS_NSP9  = [38,1,3,4,96,103,2,97]
PRIMARY_NSP12  = {740,744}
PRIMARY_NSP9   = {36}
PAN_COV_NSP12  = {733,202,221}
PAN_COV_NSP9   = {97,103}

HYDROPHOBIC    = {"ALA","VAL","ILE","LEU","MET",
                  "PHE","TRP","TYR","PRO"}
CHARGED_POS    = {"LYS","ARG","HIS"}
CHARGED_NEG    = {"ASP","GLU"}
HBOND_ATOMS    = {"N","O","NZ","NH1","NH2","NE",
                  "OG","OG1","OE1","OE2","ND1",
                  "ND2","NE2","OH"}

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

BLUE   = "#2C5F8A"
RED    = "#C0392B"
ORANGE = "#E67E22"
GOLD   = "#F39C12"
GREEN  = "#27AE60"
GREY   = "#95A5A6"


def compute_bsa(structure):
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    sasa_complex = {}
    for m in structure:
        for c in m:
            if c.id not in [CHAIN_NSP12, CHAIN_NSP9]:
                continue
            for r in c:
                if r.id[0] != " ":
                    continue
                sasa_complex[(c.id, r.id[1])] = r.sasa

    parser   = PDB.PDBParser(QUIET=True)
    sasa_mono = {}
    import tempfile, os
    for chain_id in [CHAIN_NSP12, CHAIN_NSP9]:
        class ChainSelect(PDB.Select):
            def accept_chain(self, c):
                return c.id == chain_id
        tmp = tempfile.NamedTemporaryFile(
            suffix=".pdb", delete=False)
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(tmp.name, ChainSelect())
        tmp.close()
        mono_s = parser.get_structure("mono", tmp.name)
        sr.compute(mono_s, level="R")
        for m2 in mono_s:
            for c2 in m2:
                for r2 in c2:
                    if r2.id[0] != " ":
                        continue
                    sasa_mono[(chain_id,
                               r2.id[1])] = r2.sasa
        os.unlink(tmp.name)

    bsa = {}
    for key in sasa_complex:
        mono     = sasa_mono.get(key, 0)
        bsa[key] = max(0.0,
                       float(mono -
                              sasa_complex[key]))
    return bsa


def count_contacts(res_a, res_b_dict):
    rname = res_a.resname.strip()
    lost  = {"salt_bridge":0,"hbond":0,
             "hydrophobic":0}
    for rb in res_b_dict.values():
        rnb = rb.resname.strip()
        for a1 in res_a.get_atoms():
            for a2 in rb.get_atoms():
                try:
                    d = a1 - a2
                except Exception:
                    continue
                if d <= 5.0:
                    if ((rname in CHARGED_POS and
                         rnb   in CHARGED_NEG) or
                        (rname in CHARGED_NEG and
                         rnb   in CHARGED_POS)):
                        lost["salt_bridge"] += 1
                if (d <= 3.5 and
                    a1.name in HBOND_ATOMS and
                    a2.name in HBOND_ATOMS):
                    lost["hbond"] += 1
                if (d <= 4.5 and
                    rname in HYDROPHOBIC and
                    rnb   in HYDROPHOBIC and
                    a1.element == "C" and
                    a2.element == "C"):
                    lost["hydrophobic"] += 1
    return lost


def alanine_scan(structure):
    res_12 = {r.id[1]: r
              for m in structure for c in m
              if c.id == CHAIN_NSP12
              for r in c if r.id[0] == " "}
    res_9  = {r.id[1]: r
              for m in structure for c in m
              if c.id == CHAIN_NSP9
              for r in c if r.id[0] == " "}
    results = {}
    all12 = set(HOTSPOTS_NSP12) | PRIMARY_NSP12
    all9  = set(HOTSPOTS_NSP9)  | PRIMARY_NSP9
    for rn,res in res_12.items():
        if rn in all12:
            results[(CHAIN_NSP12, rn)] = \
                count_contacts(res, res_9)
    for rn,res in res_9.items():
        if rn in all9:
            results[(CHAIN_NSP9, rn)] = \
                count_contacts(res, res_12)
    return results


def composite_ranking(bsa, alascan, cons12, cons9):
    cons_map = {}
    for _,row in cons12.iterrows():
        cons_map[(CHAIN_NSP12,
                  int(row["position"]))] = \
            float(row["conservation"])
    for _,row in cons9.iterrows():
        cons_map[(CHAIN_NSP9,
                  int(row["position"]))] = \
            float(row["conservation"])

    all_pos = (
        [(CHAIN_NSP12, p)
         for p in set(HOTSPOTS_NSP12)|PRIMARY_NSP12]+
        [(CHAIN_NSP9, p)
         for p in set(HOTSPOTS_NSP9)|PRIMARY_NSP9])

    rows = []
    for (chain, pos) in set(all_pos):
        bsa_val    = float(bsa.get((chain,pos),0))
        cons       = float(cons_map.get((chain,pos),0))
        scan       = alascan.get((chain,pos),{})
        sb_loss    = scan.get("salt_bridge",0)
        hb_loss    = scan.get("hbond",0)
        hy_loss    = scan.get("hydrophobic",0)
        total_loss = sb_loss*3+hb_loss*2+hy_loss
        primary    = (pos in PRIMARY_NSP12
                      if chain==CHAIN_NSP12
                      else pos in PRIMARY_NSP9)
        pancov     = (pos in PAN_COV_NSP12
                      if chain==CHAIN_NSP12
                      else pos in PAN_COV_NSP9)
        # Pan-coronavirus bonus
        pc_bonus   = 1.20 if pancov else 1.0

        rows.append({
            "residue":      (("NSP12-"
                               if chain==CHAIN_NSP12
                               else "NSP9-")+str(pos)),
            "chain":        ("NSP12"
                              if chain==CHAIN_NSP12
                              else "NSP9"),
            "position":     pos,
            "bsa":          round(bsa_val,1),
            "contact_score":float(total_loss),
            "conservation": cons,
            "total_loss":   total_loss,
            "sb_loss":      sb_loss,
            "hb_loss":      hb_loss,
            "hy_loss":      hy_loss,
            "primary_sb":   primary,
            "pan_cov":      pancov,
            "pc_bonus":     pc_bonus,
            "bsa_score":    bsa_val,
            "cons_score":   cons,
        })

    df      = pd.DataFrame(rows)
    max_bsa = df["bsa_score"].max()
    max_ct  = df["contact_score"].max()
    df["bsa_norm"] = (df["bsa_score"]/max_bsa
                      if max_bsa>0 else 0.0)
    df["ct_norm"]  = (df["contact_score"]/max_ct
                      if max_ct>0  else 0.0)
    df["composite"] = (
        0.35*df["bsa_norm"]  +
        0.35*df["ct_norm"]   +
        0.20*df["cons_score"]+
        0.10*df["primary_sb"].astype(float)
    ) * df["pc_bonus"]
    max_c = df["composite"].max()
    if max_c > 0:
        df["composite"] = df["composite"]/max_c
    return df.sort_values("composite",
                           ascending=False)


def fig_bsa(df):
    top = df.nlargest(15,"bsa")
    colors = []
    for _,row in top.iterrows():
        if row["primary_sb"]:
            colors.append(RED)
        elif row["pan_cov"]:
            colors.append(GOLD)
        elif row["chain"]=="NSP12":
            colors.append(BLUE)
        else:
            colors.append(GREEN)

    fig,ax = plt.subplots(figsize=(10,6))
    y = np.arange(len(top))
    ax.barh(y, top["bsa"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    for i,(_,row) in enumerate(top.iterrows()):
        ax.text(row["bsa"]+0.5, i,
                f"{row['bsa']:.1f} A2",
                va="center", fontsize=8.5)
    ax.set_yticks(y)
    ylabels = []
    for _,row in top.iterrows():
        lbl = row["residue"]
        if row["primary_sb"]:
            lbl = f"★  {lbl}"
        elif row["pan_cov"]:
            lbl = f"◆  {lbl}"
        ylabels.append(lbl)
    ax.set_yticklabels(ylabels, fontsize=9.5)
    ax.set_xlabel("Buried Surface Area (A2)")
    ax.set_title(
        "NSP9-NSP12: Top 15 Residues by BSA\n"
        "8SQK (2.90 A)", fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    patches = [
        mpatches.Patch(color=RED,
            label="AF3 SB anchor (★)"),
        mpatches.Patch(color=GOLD,
            label="Pan-coronavirus (◆)"),
        mpatches.Patch(color=BLUE,label="NSP12"),
        mpatches.Patch(color=GREEN,label="NSP9"),
    ]
    ax.legend(handles=patches, fontsize=9,
              loc="lower right")
    plt.tight_layout()
    out = OUT_DIR/"Fig4_NSP9-NSP12_BSA_5.png"
    plt.savefig(out); plt.close()
    print(f"  Saved: {out.name}")


def fig_alascan(df):
    top = df[df["total_loss"]>0].nlargest(
        15,"total_loss")
    if top.empty:
        top = df.nlargest(15,"bsa")
    x   = np.arange(len(top))
    fig,ax = plt.subplots(figsize=(12,5))
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
    ax.set_xticklabels(top["residue"],
                       rotation=45, ha="right",
                       fontsize=9)
    ax.set_ylabel("Estimated energy loss (AU)")
    ax.set_title(
        "NSP9-NSP12: Alanine Scanning\n"
        "Estimated contact loss per residue",
        fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    out = OUT_DIR/"Fig5_NSP9-NSP12_AlaScan_5.png"
    plt.savefig(out); plt.close()
    print(f"  Saved: {out.name}")


def fig_ranking(df):
    top = df.head(15)
    colors = []
    for _,row in top.iterrows():
        if row["primary_sb"]:
            colors.append(RED)
        elif row["pan_cov"]:
            colors.append(GOLD)
        elif row["chain"]=="NSP12":
            colors.append(BLUE)
        else:
            colors.append(GREEN)

    fig,ax = plt.subplots(figsize=(10,6))
    y = np.arange(len(top))
    ax.barh(y, top["composite"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    for i,(_,row) in enumerate(top.iterrows()):
        ax.text(row["composite"]+0.005, i,
                f"{row['composite']:.4f}",
                va="center", fontsize=8.5)
    ax.set_yticks(y)
    ylabels = []
    for _,row in top.iterrows():
        lbl = row["residue"]
        if row["primary_sb"]:
            lbl = f"★  {lbl}"
        elif row["pan_cov"]:
            lbl = f"◆  {lbl}"
        ylabels.append(lbl)
    ax.set_yticklabels(ylabels, fontsize=9.5)
    ax.set_xlabel("Composite Drug Target Score")
    ax.set_title(
        "NSP9-NSP12: Composite Hotspot Ranking\n"
        "BSA (35%) + Contact loss (35%) + "
        "Conservation (20%) + Primary SB (10%)\n"
        "Pan-coronavirus bonus x1.20",
        fontweight="bold")
    ax.set_xlim(0,1.12)
    ax.grid(axis="x", alpha=0.2)
    patches = [
        mpatches.Patch(color=RED,
            label="AF3 SB anchor (★)"),
        mpatches.Patch(color=GOLD,
            label="Pan-coronavirus (◆) x1.20"),
        mpatches.Patch(color=BLUE,label="NSP12"),
        mpatches.Patch(color=GREEN,label="NSP9"),
    ]
    ax.legend(handles=patches, fontsize=9,
              loc="lower right")
    plt.tight_layout()
    out = OUT_DIR / \
          "Fig6_NSP9-NSP12_composite_ranking_5.png"
    plt.savefig(out); plt.close()
    print(f"  Saved: {out.name}")


def main():
    print("\n" + "="*57)
    print("  Script 10_5: BSA+AlaScan+Ranking — NSP9-NSP12")
    print("  PDB: 8SQK (Chain G=NSP9)")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("8SQK", PDB_FILE)
    cons12    = pd.read_csv(
        RES_DIR/"conservation_NSP12.csv")
    cons9     = pd.read_csv(
        RES_DIR/"conservation_NSP9.csv")

    # AA lookup for residue labels
    res_aa_12 = {r.id[1]: r.resname.strip()
                 for m in structure for c in m
                 if c.id==CHAIN_NSP12
                 for r in c if r.id[0]==" "}
    res_aa_9  = {r.id[1]: r.resname.strip()
                 for m in structure for c in m
                 if c.id==CHAIN_NSP9
                 for r in c if r.id[0]==" "}

    print("  1. Computing BSA...")
    bsa = compute_bsa(structure)
    for (ch,pos),val in sorted(
            bsa.items(),key=lambda x:-x[1])[:10]:
        print(f"     {'NSP12' if ch==CHAIN_NSP12 else 'NSP9 '}"
              f"-{pos}: {val:.1f} A2")

    print("\n  2. Alanine scanning...")
    alascan = alanine_scan(structure)
    for (ch,pos),lost in sorted(
            alascan.items(),
            key=lambda x:-(x[1]["salt_bridge"]*3+
                            x[1]["hbond"]*2+
                            x[1]["hydrophobic"]))[:10]:
        total = (lost["salt_bridge"]*3+
                 lost["hbond"]*2+
                 lost["hydrophobic"])
        if total > 0:
            print(f"     {'NSP12' if ch==CHAIN_NSP12 else 'NSP9 '}"
                  f"-{pos}: loss={total} "
                  f"(SB={lost['salt_bridge']} "
                  f"HB={lost['hbond']} "
                  f"HY={lost['hydrophobic']})")

    print("\n  3. Composite ranking...")
    df = composite_ranking(bsa, alascan, cons12, cons9)
    print(f"\n  {'Rank':<5} {'Residue':<15} "
          f"{'BSA':>7} {'Loss':>5} "
          f"{'Cons':>6} {'Score':>7} {'Flags'}")
    print(f"  {'-'*62}")
    for rank,(_,row) in enumerate(
            df.head(10).iterrows(),1):
        flags = ""
        if row["primary_sb"]: flags += "★SB "
        if row["pan_cov"]:    flags += "◆PC "
        pos = int(row['position'])
        aa  = (res_aa_12.get(pos,"???")
               if row['chain']=="NSP12"
               else res_aa_9.get(pos,"???"))
        lbl = f"{row['chain']}-{aa}{pos}"
        print(f"  {rank:<5} {lbl:<20} "
              f"{row['bsa']:>7.1f} "
              f"{int(row['total_loss']):>5} "
              f"{row['conservation']:>6.3f} "
              f"{row['composite']:>7.4f} {flags}")

    print("\n  4. Generating figures...")
    fig_bsa(df)
    fig_alascan(df)
    fig_ranking(df)

    csv_out = OUT_DIR / \
              "composite_ranking_NSP9-NSP12_5.csv"
    df.to_csv(csv_out, index=False)
    print(f"\n  Saved: {csv_out.name}")

    bsa_out = OUT_DIR / \
              "bsa_alascan_NSP9-NSP12_5.json"
    with open(bsa_out,"w") as f:
        json.dump({
            "bsa":     {f"{ch},{pos}":v
                        for (ch,pos),v in bsa.items()},
            "alascan": {f"{ch},{pos}":lost
                        for (ch,pos),lost in
                        alascan.items()},
        }, f, indent=2)
    print(f"  Saved: {bsa_out.name}")

    print(f"\n  Top pharmacophore: "
          f"{df.iloc[0]['residue']} "
          f"(score={df.iloc[0]['composite']:.4f})")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
