"""
Script 10_2: BSA + Alanine Scanning + Composite Ranking — NSP10-NSP16
======================================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Analyses:
  1. Buried Surface Area (BSA) per hotspot residue
  2. Computational alanine scanning (energetic impact)
  3. Composite hotspot ranking (BSA + conservation + energy)
  4. Zn1 finger proximity bonus scoring

Primary PDB: 6W4H (1.80 Å) — highest resolution
Output: Figs 4-6, composite_ranking CSV, BSA/AlaScan JSON

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/10_BSA_alascan_ranking_NSP10-NSP16_2.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP10-NSP16"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
OUT_DIR = PROJECT / "results"
OUT_DIR.mkdir(exist_ok=True)

PDB_FILE       = PDB_DIR / "6W4H.pdb"
CHAIN_NSP10    = "B"
CHAIN_NSP16    = "A"
NSP10_OFFSET   = 4270
NSP10_SHIFT    = 17
NSP16_OFFSET   = 6797

HOTSPOTS_NSP10 = [40,42,43,44,45,71,76,78,80,93,94,95,96]
HOTSPOTS_NSP16 = [40,41,44,76,83,87,102,104,106,244,247]
PRIMARY_NSP10  = {80, 93, 95}
PRIMARY_NSP16  = {102, 106}
ZN1_COORD      = {74, 77, 83, 90}
ZN1_NEIGHBOURS = {71,74,76,77,78,80,83,90,93}

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
TEAL   = "#148F77"
GREY   = "#95A5A6"


def genome_to_local(gpos, chain):
    if chain == CHAIN_NSP10:
        return gpos - NSP10_OFFSET + NSP10_SHIFT
    else:
        return gpos - NSP16_OFFSET


# ── 1. BSA ─────────────────────────────────────────────────

def compute_bsa(structure):
    """
    BSA = SASA(isolated chain) - SASA(in complex).
    Returns dict: (chain, local_pos) -> BSA in Å²
    """
    sr = ShrakeRupley()

    # SASA of full complex
    sr.compute(structure, level="R")
    sasa_complex = {}
    for m in structure:
        for c in m:
            for r in c:
                if r.id[0] != " ":
                    continue
                local = genome_to_local(r.id[1], c.id)
                sasa_complex[(c.id, local)] = r.sasa

    # SASA of isolated chains
    bsa = {}
    for target_chain in [CHAIN_NSP10, CHAIN_NSP16]:
        # Build single-chain structure
        class ChainSelect(PDB.Select):
            def accept_chain(self, chain):
                return chain.id == target_chain
            def accept_residue(self, res):
                return res.id[0] == " "

        io  = PDB.PDBIO()
        tmp = Path("/tmp/chain_tmp.pdb")
        io.set_structure(structure)
        io.save(str(tmp), ChainSelect())

        parser2 = PDB.PDBParser(QUIET=True)
        single  = parser2.get_structure("s", tmp)
        sr.compute(single, level="R")

        for m in single:
            for c in m:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    local  = genome_to_local(r.id[1],
                                             target_chain)
                    key    = (target_chain, local)
                    sasa_i = r.sasa
                    sasa_c = sasa_complex.get(key, sasa_i)
                    bsa[key] = max(0.0,
                                   round(sasa_i - sasa_c, 2))
    return bsa


# ── 2. Alanine scanning ────────────────────────────────────

def alanine_scan(structure):
    """
    Estimate energetic loss upon Ala mutation.
    Counts contacts (H-bond + salt bridge + hydrophobic)
    that each hotspot residue makes across the interface.
    """
    res_nsp10 = {
        genome_to_local(r.id[1], CHAIN_NSP10): r
        for m in structure for c in m
        if c.id == CHAIN_NSP10
        for r in c if r.id[0] == " "}

    res_nsp16 = {
        genome_to_local(r.id[1], CHAIN_NSP16): r
        for m in structure for c in m
        if c.id == CHAIN_NSP16
        for r in c if r.id[0] == " "}

    losses = {}

    def count_contacts(rdict_a, rdict_b, hotspots_a):
        for pos_a in hotspots_a:
            if pos_a not in rdict_a:
                continue
            res_a  = rdict_a[pos_a]
            rname  = res_a.resname.strip()
            loss   = 0
            detail = {"hbond":0,"salt":0,"hydrophobic":0}
            for pos_b, res_b in rdict_b.items():
                rname_b = res_b.resname.strip()
                for a1 in res_a.get_atoms():
                    for a2 in res_b.get_atoms():
                        try:
                            d = a1 - a2
                        except Exception:
                            continue
                        # Salt bridge
                        if d <= 4.0:
                            if ((rname in CHARGED_POS and
                                 rname_b in CHARGED_NEG) or
                                (rname in CHARGED_NEG and
                                 rname_b in CHARGED_POS)):
                                loss += 3
                                detail["salt"] += 1
                        # H-bond
                        if (d <= 3.5 and
                            a1.name in HBOND_ATOMS and
                                a2.name in HBOND_ATOMS):
                            loss += 2
                            detail["hbond"] += 1
                        # Hydrophobic
                        if (d <= 4.5 and
                            rname in HYDROPHOBIC and
                                rname_b in HYDROPHOBIC and
                                a1.element == "C" and
                                a2.element == "C"):
                            loss += 1
                            detail["hydrophobic"] += 1
            key = f"{rname}{pos_a}"
            losses[key] = {
                "total_loss":  loss,
                "detail":      detail,
                "resname":     rname,
                "position":    pos_a,
                "is_hotspot":  pos_a in hotspots_a,
            }

    count_contacts(res_nsp10, res_nsp16, HOTSPOTS_NSP10)
    count_contacts(res_nsp16, res_nsp10, HOTSPOTS_NSP16)
    return losses


# ── 3. Composite ranking ───────────────────────────────────

def composite_ranking(bsa, alascan, cons10, cons16):
    """
    Score = (contact_score/10)
          × conservation
          × (0.5 + 0.5 × burial_norm)
          × (0.5 + 0.5 × energy_norm)
          × zn1_bonus

    zn1_bonus = 1.2 for NSP10 residues near Zn1 finger
    """
    rows = []

    # Gather all data
    all_losses = [v["total_loss"]
                  for v in alascan.values()]
    max_loss   = max(all_losses) if all_losses else 1
    all_bsa    = list(bsa.values())
    max_bsa    = max(all_bsa) if all_bsa else 1

    cons10_map = dict(zip(cons10["position"],
                          cons10["conservation"]))
    cons16_map = dict(zip(cons16["position"],
                          cons16["conservation"]))

    with open(RES_DIR / "interface_analysis.json") as f:
        iface = json.load(f)
    top20_10 = dict(
        iface["PDB_6W4H"]["top20_NSP10"])
    top20_16 = dict(
        iface["PDB_6W4H"]["top20_NSP16"])
    max_score = max(
        list(top20_10.values()) +
        list(top20_16.values()) + [1])

    for key, scan in alascan.items():
        pos    = scan["position"]
        rname  = scan["resname"]
        is_10  = pos in HOTSPOTS_NSP10
        chain  = CHAIN_NSP10 if is_10 else CHAIN_NSP16

        bsa_val  = bsa.get((chain, pos), 0)
        cons_val = (cons10_map.get(pos, 0) if is_10
                    else cons16_map.get(pos, 0))
        cont_map = top20_10 if is_10 else top20_16
        cont_val = cont_map.get(key, 0)

        burial_norm = bsa_val / max_bsa
        energy_norm = scan["total_loss"] / max_loss
        cont_norm   = cont_val / max_score

        # Zn1 bonus for NSP10 residues near zinc finger
        zn1_bonus = (1.2 if (is_10 and
                              pos in ZN1_NEIGHBOURS)
                     else 1.0)

        score = (cont_norm
                 * cons_val
                 * (0.5 + 0.5 * burial_norm)
                 * (0.5 + 0.5 * energy_norm)
                 * zn1_bonus)

        primary_sb = (pos in PRIMARY_NSP10 if is_10
                      else pos in PRIMARY_NSP16)
        zn1_flag   = (is_10 and
                      pos in ZN1_NEIGHBOURS)

        rows.append({
            "residue":       key,
            "chain":         "NSP10" if is_10 else "NSP16",
            "position":      pos,
            "resname":       rname,
            "bsa":           round(bsa_val, 2),
            "contact_score": round(cont_val, 2),
            "conservation":  round(cons_val, 3),
            "total_loss":    scan["total_loss"],
            "burial_norm":   round(burial_norm, 4),
            "energy_norm":   round(energy_norm, 4),
            "zn1_bonus":     zn1_bonus,
            "composite":     round(score, 4),
            "primary_sb":    primary_sb,
            "zn1_region":    zn1_flag,
        })

    df = pd.DataFrame(rows).sort_values(
        "composite", ascending=False)
    return df


# ── Figures ────────────────────────────────────────────────

def fig_bsa(bsa, title=""):
    data = []
    for pos in HOTSPOTS_NSP10:
        v = bsa.get((CHAIN_NSP10, pos), 0)
        data.append({"residue": f"NSP10-{pos}",
                     "bsa": v, "chain":"NSP10",
                     "pos": pos})
    for pos in HOTSPOTS_NSP16:
        v = bsa.get((CHAIN_NSP16, pos), 0)
        data.append({"residue": f"NSP16-{pos}",
                     "bsa": v, "chain":"NSP16",
                     "pos": pos})
    df = pd.DataFrame(data).sort_values("bsa",
                                         ascending=True)

    fig, ax = plt.subplots(figsize=(10, 10),
                           constrained_layout=True)
    colors = []
    for _, r in df.iterrows():
        pos = r["pos"]
        if r["chain"] == "NSP10":
            if pos in PRIMARY_NSP10:
                colors.append(RED)
            elif pos in ZN1_NEIGHBOURS:
                colors.append(GOLD)
            else:
                colors.append(BLUE)
        else:
            colors.append(
                RED if pos in PRIMARY_NSP16
                else TEAL)

    y = np.arange(len(df))
    ax.barh(y, df["bsa"], color=colors,
            edgecolor="white", linewidth=0.6,
            height=0.72)
    ax.axvline(20, color="black", linestyle="--",
               linewidth=1.0, alpha=0.6,
               label="BSA threshold (20 Å²)")

    for i, (_, r) in enumerate(df.iterrows()):
        ax.text(r["bsa"]+1.5, i,
                f"{r['bsa']:.1f} Å²",
                va="center", fontsize=8.5)

    ax.set_yticks(y)
    ax.set_yticklabels(df["residue"], fontsize=9.5)
    ax.set_xlabel("Buried Surface Area (Å²)")
    ax.set_title(
        "NSP10–NSP16: Buried Surface Area per Hotspot\n"
        "(6W4H, 1.80 Å resolution)",
        fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.2)

    patches = [
        mpatches.Patch(color=RED,
            label="Primary salt bridge (NSP10)"),
        mpatches.Patch(color=GOLD,
            label="Near Zn1 finger (NSP10)"),
        mpatches.Patch(color=BLUE,
            label="Other hotspot (NSP10)"),
        mpatches.Patch(color=RED,
            label="Primary salt bridge (NSP16)"),
        mpatches.Patch(color=TEAL,
            label="Other hotspot (NSP16)"),
    ]
    ax.legend(handles=patches, loc="lower right",
              fontsize=8.5)

    out = OUT_DIR / "Fig4_NSP10-NSP16_BSA_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


def fig_alascan(alascan):
    data = []
    for key, v in alascan.items():
        pos   = v["position"]
        is_10 = pos in HOTSPOTS_NSP10
        data.append({
            "residue":    key,
            "chain":      "NSP10" if is_10 else "NSP16",
            "pos":        pos,
            "total":      v["total_loss"],
            "hbond":      v["detail"]["hbond"],
            "salt":       v["detail"]["salt"],
            "hydrophobic":v["detail"]["hydrophobic"],
        })
    df = pd.DataFrame(data).sort_values("total",
                                         ascending=True)

    fig, ax = plt.subplots(figsize=(10, 10),
                           constrained_layout=True)
    y     = np.arange(len(df))
    width = 0.65
    hb    = np.array(df["hbond"]) * 2
    sb    = np.array(df["salt"])  * 3
    hy    = np.array(df["hydrophobic"]) * 1

    ax.barh(y, hy, width, color=ORANGE,
            label="Hydrophobic (×1)")
    ax.barh(y, hb, width, left=hy, color=BLUE,
            label="H-bond (×2)")
    ax.barh(y, sb, width, left=hy+hb, color=RED,
            label="Salt bridge (×3)")

    ax.axvline(6, color="black", linestyle="--",
               linewidth=1.0, alpha=0.6,
               label="Energetic hotspot threshold")

    for i, (_, r) in enumerate(df.iterrows()):
        ax.text(r["total"]+0.5, i,
                str(int(r["total"])),
                va="center", fontsize=8.5,
                fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels(df["residue"], fontsize=9.5)
    ax.set_xlabel("Estimated Energetic Contact Loss\n"
                  "(upon Ala mutation)")
    ax.set_title(
        "NSP10–NSP16: Computational Alanine Scanning\n"
        "(6W4H, interface contact loss per residue)",
        fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.2)

    out = OUT_DIR / "Fig5_NSP10-NSP16_AlaScan_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


def fig_composite(df_rank):
    top = df_rank.head(15).iloc[::-1]

    fig, ax = plt.subplots(figsize=(11, 8),
                           constrained_layout=True)
    y      = np.arange(len(top))
    colors = []
    for _, r in top.iterrows():
        if r["primary_sb"]:
            colors.append(RED)
        elif r["zn1_region"]:
            colors.append(GOLD)
        elif r["chain"] == "NSP16":
            colors.append(TEAL)
        else:
            colors.append(BLUE)

    bars = ax.barh(y, top["composite"],
                   color=colors, edgecolor="white",
                   linewidth=0.6, height=0.72)

    for i, (_, r) in enumerate(top.iterrows()):
        ax.text(r["composite"]+0.003, i,
                f"{r['composite']:.3f}",
                va="center", fontsize=8.5)
        info = (f"BSA={r['bsa']:.0f}Å²  "
                f"loss={r['total_loss']}  "
                f"cons={r['conservation']:.2f}")
        ax.text(0.001, i-0.32, info,
                va="center", fontsize=7,
                color="#555555")

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r['chain']}-{r['residue']}"
         for _, r in top.iterrows()],
        fontsize=9.5)
    ax.set_xlabel("Composite Drug Target Score")
    ax.set_title(
        "NSP10–NSP16: Top 15 Drug Target Candidates\n"
        "Composite Score = contact × conservation "
        "× burial × energy × Zn1 bonus",
        fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.2)

    patches = [
        mpatches.Patch(color=RED,
            label="Primary salt bridge ★"),
        mpatches.Patch(color=GOLD,
            label="Near Zn1 finger region"),
        mpatches.Patch(color=TEAL,
            label="NSP16 hotspot"),
        mpatches.Patch(color=BLUE,
            label="NSP10 hotspot"),
    ]
    ax.legend(handles=patches, loc="lower right",
              fontsize=9)

    out = OUT_DIR / "Fig6_NSP10-NSP16_composite_ranking_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────

def main():
    print("\n" + "="*57)
    print("  Script 10_2: BSA + AlaScan + Ranking "
          "— NSP10-NSP16")
    print("="*57 + "\n")

    parser_pdb = PDB.PDBParser(QUIET=True)
    structure  = parser_pdb.get_structure(
        "6W4H", PDB_FILE)

    cons10 = pd.read_csv(RES_DIR / "conservation_NSP10.csv")
    cons16 = pd.read_csv(RES_DIR / "conservation_NSP16.csv")

    # ── BSA ───────────────────────────────────────────────
    print("  1. Computing BSA (ShrakeRupley)...")
    bsa = compute_bsa(structure)
    buried = {k:v for k,v in bsa.items() if v >= 20}
    print(f"     Residues with BSA ≥ 20 Å²: "
          f"{len(buried)}/{len(bsa)}")
    print(f"\n     Top BSA residues:")
    for (ch,pos), val in sorted(
            bsa.items(), key=lambda x:-x[1])[:8]:
        nsp = "NSP10" if ch==CHAIN_NSP10 else "NSP16"
        print(f"       {nsp}-{pos:<5} BSA={val:.1f} Å²")

    # ── AlaScan ───────────────────────────────────────────
    print("\n  2. Computational alanine scanning...")
    alascan = alanine_scan(structure)
    energetic = {k:v for k,v in alascan.items()
                 if v["total_loss"] >= 2}
    print(f"     Energetic hotspots (loss ≥ 2): "
          f"{len(energetic)}/{len(alascan)}")
    print(f"\n     Top energetic contributors:")
    for k,v in sorted(
            alascan.items(),
            key=lambda x:-x[1]["total_loss"])[:8]:
        d = v["detail"]
        print(f"       {k:<10} loss={v['total_loss']:>4}  "
              f"(HB:{d['hbond']} SB:{d['salt']} "
              f"HY:{d['hydrophobic']})")

    # ── Composite ranking ─────────────────────────────────
    print("\n  3. Composite ranking...")
    df_rank = composite_ranking(
        bsa, alascan, cons10, cons16)
    print(f"\n  Top 10 drug target candidates:")
    print(f"  {'Rank':<5} {'Residue':<12} {'Chain':<7} "
          f"{'BSA':>6} {'Loss':>5} {'Cons':>6} "
          f"{'Score':>7}  Flags")
    print(f"  {'-'*5} {'-'*12} {'-'*7} {'-'*6} "
          f"{'-'*5} {'-'*6} {'-'*7}  {'-'*15}")
    for i, (_, r) in enumerate(
            df_rank.head(10).iterrows(), 1):
        flags = []
        if r["primary_sb"]:
            flags.append("★ salt bridge")
        if r["zn1_region"]:
            flags.append("⬡ Zn1")
        print(f"  {i:<5} {r['residue']:<12} "
              f"{r['chain']:<7} "
              f"{r['bsa']:>6.1f} "
              f"{int(r['total_loss']):>5} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>7.4f}  "
              f"{' '.join(flags)}")

    # ── Figures ───────────────────────────────────────────
    print("\n  4. Generating figures...")
    fig_bsa(bsa)
    fig_alascan(alascan)
    fig_composite(df_rank)

    # ── Save data ─────────────────────────────────────────
    df_rank.to_csv(
        RES_DIR / "composite_ranking_NSP10-NSP16_2.csv",
        index=False)

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.floating,
                                np.float32,
                                np.float64)):
                return float(obj)
            if isinstance(obj, np.integer):
                return int(obj)
            return super().default(obj)

    with open(RES_DIR /
              "bsa_alascan_NSP10-NSP16_2.json","w") as f:
        json.dump({
            "bsa":     {f"{c},{p}":v
                        for (c,p),v in bsa.items()},
            "alascan": alascan,
        }, f, indent=2, cls=NumpyEncoder)

    print(f"\n  Saved: composite_ranking_NSP10-NSP16_2.csv")
    print(f"  Saved: bsa_alascan_NSP10-NSP16_2.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
