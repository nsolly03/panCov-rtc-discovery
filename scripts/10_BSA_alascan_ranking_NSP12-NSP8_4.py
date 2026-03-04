"""
Script 10_4: BSA + Alanine Scanning + Composite Ranking — NSP12-NSP8
=====================================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Primary PDB: 7BV2_NSP12-NSP8 (2.90 Å)
Output: Figs 4-6, composite_ranking CSV, BSA/AlaScan JSON

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/10_BSA_alascan_ranking_NSP12-NSP8_4.py
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
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP8"
OUT_DIR  = PROJECT / "02-validation" / "NSP12-NSP8"
OUT_DIR.mkdir(exist_ok=True)

PDB_FILE       = PDB_DIR / "7BV2_NSP12-NSP8.pdb"
CHAIN_NSP12    = "A"
CHAIN_NSP8     = "B"
NSP12_OFFSET   = 0
NSP8_OFFSET    = 77

HOTSPOTS_NSP12 = [387,129,389,271,330,131,380,523,
                   91,87,332,95,117,517,99]
HOTSPOTS_NSP8  = [117,129,80,115,131,112,91,87,
                   116,95,98,83,128,90,121]
PRIMARY_NSP12  = {523,332,517}
PRIMARY_NSP8   = {80,99,79}

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


def pdb_to_local(gpos, chain):
    if chain == CHAIN_NSP12:
        return gpos - NSP12_OFFSET
    return gpos - NSP8_OFFSET


def compute_bsa(structure):
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    sasa_complex = {}
    for m in structure:
        for c in m:
            if c.id not in [CHAIN_NSP12, CHAIN_NSP8]:
                continue
            for r in c:
                if r.id[0] != " ":
                    continue
                local = pdb_to_local(r.id[1], c.id)
                sasa_complex[(c.id, local)] = r.sasa

    parser = PDB.PDBParser(QUIET=True)
    sasa_mono = {}
    import tempfile, os
    for chain_id in [CHAIN_NSP12, CHAIN_NSP8]:
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
                    local = pdb_to_local(
                        r2.id[1], chain_id)
                    sasa_mono[(chain_id, local)] = \
                        r2.sasa
        os.unlink(tmp.name)

    bsa = {}
    for key in sasa_complex:
        mono = sasa_mono.get(key, 0)
        bsa[key] = max(0.0,
                       float(mono - sasa_complex[key]))
    return bsa


def count_contacts(res_a, res_b_dict):
    rname = res_a.resname.strip()
    lost  = {"salt_bridge":0,"hbond":0,"hydrophobic":0}
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
    res_8  = {r.id[1]: r
              for m in structure for c in m
              if c.id == CHAIN_NSP8
              for r in c if r.id[0] == " "}
    results = {}
    for rn, res in res_12.items():
        local = pdb_to_local(rn, CHAIN_NSP12)
        if local not in set(HOTSPOTS_NSP12) | \
                PRIMARY_NSP12:
            continue
        results[(CHAIN_NSP12, local)] = \
            count_contacts(res, res_8)
    for rn, res in res_8.items():
        local = pdb_to_local(rn, CHAIN_NSP8)
        if local not in set(HOTSPOTS_NSP8) | \
                PRIMARY_NSP8:
            continue
        results[(CHAIN_NSP8, local)] = \
            count_contacts(res, res_12)
    return results


def composite_ranking(bsa, alascan, cons12, cons8):
    cons_map = {}
    for _, row in cons12.iterrows():
        cons_map[(CHAIN_NSP12,
                  int(row["position"]))] = \
            float(row["conservation"])
    for _, row in cons8.iterrows():
        cons_map[(CHAIN_NSP8,
                  int(row["position"]))] = \
            float(row["conservation"])

    all_pos = (
        [(CHAIN_NSP12, p)
         for p in set(HOTSPOTS_NSP12)|PRIMARY_NSP12] +
        [(CHAIN_NSP8,  p)
         for p in set(HOTSPOTS_NSP8) |PRIMARY_NSP8])

    rows = []
    for (chain, pos) in set(all_pos):
        bsa_val    = float(bsa.get((chain,pos), 0))
        cons       = float(cons_map.get((chain,pos), 0))
        scan       = alascan.get((chain,pos), {})
        sb_loss    = scan.get("salt_bridge", 0)
        hb_loss    = scan.get("hbond", 0)
        hy_loss    = scan.get("hydrophobic", 0)
        total_loss = sb_loss*3 + hb_loss*2 + hy_loss
        primary    = (pos in PRIMARY_NSP12
                      if chain==CHAIN_NSP12
                      else pos in PRIMARY_NSP8)

        rows.append({
            "residue":      (("NSP12-" if
                               chain==CHAIN_NSP12
                               else "NSP8-") + str(pos)),
            "chain":        ("NSP12" if
                              chain==CHAIN_NSP12
                              else "NSP8"),
            "position":     pos,
            "bsa":          round(bsa_val, 1),
            "contact_score":float(total_loss),
            "conservation": cons,
            "total_loss":   total_loss,
            "sb_loss":      sb_loss,
            "hb_loss":      hb_loss,
            "hy_loss":      hy_loss,
            "primary_sb":   primary,
            "bsa_score":    bsa_val,
            "cons_score":   cons,
        })

    df      = pd.DataFrame(rows)
    max_bsa = df["bsa_score"].max()
    max_ct  = df["contact_score"].max()
    df["bsa_norm"] = (df["bsa_score"]/max_bsa
                      if max_bsa > 0 else 0.0)
    df["ct_norm"]  = (df["contact_score"]/max_ct
                      if max_ct  > 0 else 0.0)
    df["composite"] = (
        0.35*df["bsa_norm"]  +
        0.35*df["ct_norm"]   +
        0.20*df["cons_score"]+
        0.10*df["primary_sb"].astype(float))
    max_c = df["composite"].max()
    if max_c > 0:
        df["composite"] = df["composite"] / max_c
    return df.sort_values("composite", ascending=False)


def fig_bsa(df):
    top = df.nlargest(15,"bsa")
    colors = [RED if r["primary_sb"]
              else (BLUE if r["chain"]=="NSP12"
                    else GREEN)
              for _,r in top.iterrows()]
    fig, ax = plt.subplots(figsize=(10,6))
    y = np.arange(len(top))
    ax.barh(y, top["bsa"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    for i,(_,row) in enumerate(top.iterrows()):
        ax.text(row["bsa"]+0.5, i,
                f"{row['bsa']:.1f} Å²",
                va="center", fontsize=8.5)
    ax.set_yticks(y)
    ylabels = [f"★  {r['residue']}"
               if r["primary_sb"]
               else r["residue"]
               for _,r in top.iterrows()]
    ax.set_yticklabels(ylabels, fontsize=9.5)
    ax.set_xlabel("Buried Surface Area (Å²)")
    ax.set_title(
        "NSP12–NSP8: Top 15 Residues by BSA\n"
        "7BV2 (2.90 Å)", fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    patches = [
        mpatches.Patch(color=RED,
            label="Primary SB (★)"),
        mpatches.Patch(color=BLUE, label="NSP12"),
        mpatches.Patch(color=GREEN,label="NSP8"),
    ]
    ax.legend(handles=patches, fontsize=9,
              loc="lower right")
    plt.tight_layout()
    out = OUT_DIR / "Fig4_NSP12-NSP8_BSA_4.png"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def fig_alascan(df):
    top = df[df["total_loss"]>0].nlargest(
        15,"total_loss")
    if top.empty:
        top = df.nlargest(15,"contact_score")
    x   = np.arange(len(top))
    w   = 0.25
    fig, ax = plt.subplots(figsize=(12,5))
    ax.bar(x-w, top["sb_loss"]*3, w,
           label="Salt bridge (x3)",
           color=RED, alpha=0.85)
    ax.bar(x,   top["hb_loss"]*2, w,
           label="H-bond (x2)",
           color=BLUE, alpha=0.85)
    ax.bar(x+w, top["hy_loss"],   w,
           label="Hydrophobic (x1)",
           color=ORANGE, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(top["residue"],
                       rotation=45, ha="right",
                       fontsize=9)
    ax.set_ylabel("Estimated energy loss (AU)")
    ax.set_title(
        "NSP12–NSP8: Alanine Scanning\n"
        "Estimated contact loss per residue",
        fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    out = OUT_DIR / "Fig5_NSP12-NSP8_AlaScan_4.png"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def fig_ranking(df):
    top = df.head(15)
    colors = [RED if r["primary_sb"]
              else (BLUE if r["chain"]=="NSP12"
                    else GREEN)
              for _,r in top.iterrows()]
    fig, ax = plt.subplots(figsize=(10,6))
    y = np.arange(len(top))
    ax.barh(y, top["composite"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    for i,(_,row) in enumerate(top.iterrows()):
        ax.text(row["composite"]+0.005, i,
                f"{row['composite']:.4f}",
                va="center", fontsize=8.5)
    ax.set_yticks(y)
    ylabels = [f"★  {r['residue']}"
               if r["primary_sb"]
               else r["residue"]
               for _,r in top.iterrows()]
    ax.set_yticklabels(ylabels, fontsize=9.5)
    ax.set_xlabel("Composite Drug Target Score")
    ax.set_title(
        "NSP12–NSP8: Composite Hotspot Ranking\n"
        "BSA (35%) + Contact loss (35%) + "
        "Conservation (20%) + Primary SB (10%)",
        fontweight="bold")
    ax.set_xlim(0, 1.12)
    ax.grid(axis="x", alpha=0.2)
    patches = [
        mpatches.Patch(color=RED,
            label="Primary SB (ASP523/LYS332/ASP517) (★)"),
        mpatches.Patch(color=BLUE, label="NSP12"),
        mpatches.Patch(color=GREEN,label="NSP8"),
    ]
    ax.legend(handles=patches, fontsize=9,
              loc="lower right")
    plt.tight_layout()
    out = OUT_DIR / \
          "Fig6_NSP12-NSP8_composite_ranking_4.png"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def main():
    print("\n" + "="*57)
    print("  Script 10_4: BSA+AlaScan+Ranking — NSP12-NSP8")
    print("  PDB: 7BV2_NSP12-NSP8")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("7BV2", PDB_FILE)
    cons12    = pd.read_csv(
        RES_DIR / "conservation_NSP12.csv")
    cons8     = pd.read_csv(
        RES_DIR / "conservation_NSP8.csv")

    print("  1. Computing BSA...")
    bsa = compute_bsa(structure)
    for (ch,pos),val in sorted(
            bsa.items(), key=lambda x:-x[1])[:10]:
        print(f"     {'NSP12' if ch==CHAIN_NSP12 else 'NSP8 '}"
              f"-{pos}: {val:.1f} Å²")

    print("\n  2. Alanine scanning...")
    alascan = alanine_scan(structure)
    for (ch,pos),lost in sorted(
            alascan.items(),
            key=lambda x:-(x[1]["salt_bridge"]*3+
                            x[1]["hbond"]*2+
                            x[1]["hydrophobic"]))[:10]:
        total = (lost["salt_bridge"]*3 +
                 lost["hbond"]*2 +
                 lost["hydrophobic"])
        if total > 0:
            print(f"     {'NSP12' if ch==CHAIN_NSP12 else 'NSP8 '}"
                  f"-{pos}: loss={total} "
                  f"(SB={lost['salt_bridge']} "
                  f"HB={lost['hbond']} "
                  f"HY={lost['hydrophobic']})")

    print("\n  3. Composite ranking...")
    df = composite_ranking(bsa, alascan, cons12, cons8)
    print(f"\n  {'Rank':<5} {'Residue':<15} "
          f"{'BSA':>7} {'Loss':>5} "
          f"{'Cons':>6} {'Score':>7} {'Flags'}")
    print(f"  {'-'*60}")
    for rank,(_,row) in enumerate(
            df.head(10).iterrows(),1):
        flags = "★SB" if row["primary_sb"] else ""
        print(f"  {rank:<5} {row['residue']:<15} "
              f"{row['bsa']:>7.1f} "
              f"{int(row['total_loss']):>5} "
              f"{row['conservation']:>6.3f} "
              f"{row['composite']:>7.4f} {flags}")

    print("\n  4. Generating figures...")
    fig_bsa(df)
    fig_alascan(df)
    fig_ranking(df)

    csv_out = OUT_DIR / \
              "composite_ranking_NSP12-NSP8_4.csv"
    df.to_csv(csv_out, index=False)
    print(f"\n  Saved: {csv_out.name}")

    bsa_out = OUT_DIR / "bsa_alascan_NSP12-NSP8_4.json"
    with open(bsa_out,"w") as f:
        json.dump({
            "bsa":     {f"{ch},{pos}": v
                        for (ch,pos),v in bsa.items()},
            "alascan": {f"{ch},{pos}": lost
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
