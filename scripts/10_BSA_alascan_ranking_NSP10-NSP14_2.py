"""
Script 10_2: BSA + Alanine Scanning + Composite Ranking — NSP10-NSP14
======================================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  1. Buried Surface Area (BSA) per hotspot residue
     BSA = SASA(unbound chain) - SASA(complex) per residue
     Threshold: BSA > 20 Å² = significantly buried

  2. Computational alanine scanning
     Estimates energetic impact of mutating each hotspot to ALA
     Counts H-bonds + salt bridges + hydrophobic contacts lost
     Energetic hotspot: contact_loss >= 2

  3. Composite hotspot ranking
     Score = interface_score × conservation × burial_norm × energy_norm
     Produces final drug target priority list

  4. Publication figures (Fig4, Fig5, Fig6)

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/10_BSA_alascan_ranking_NSP10-NSP14_2.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
OUT_DIR = PROJECT / "results"
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
GREY   = "#95A5A6"

AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS",
    "Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE",
    "L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
    "S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}

# ── Hotspot residues ───────────────────────────────────────
HOTSPOTS = {
    "NSP10": [5,19,21,40,42,44,45,80,93],   # conserved ≥0.8
    "NSP14": [4,7,8,9,10,20,25,27,127,201],
}
ALL_HOTSPOTS = (
    [("A", r) for r in HOTSPOTS["NSP10"]] +
    [("B", r) for r in HOTSPOTS["NSP14"]]
)
PRIMARY = {("A",80), ("A",93), ("B",126), ("B",127)}

# ── Contact distance thresholds ────────────────────────────
HBOND_DIST     = 3.5
HYDRO_DIST     = 4.5
SALTBR_DIST    = 4.0
CHARGED_POS    = {"LYS","ARG","HIS"}
CHARGED_NEG    = {"ASP","GLU"}
HYDROPHOBIC_AA = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","TYR","PRO"}


# ══════════════════════════════════════════════════════════
# PART 1 — BURIED SURFACE AREA
# ══════════════════════════════════════════════════════════

def compute_bsa(pdb_file, chain_a="A", chain_b="B"):
    """
    Compute BSA per residue for hotspot residues.
    BSA = SASA(isolated chain) - SASA(complex)
    """
    parser = PDBParser(QUIET=True)
    sr     = ShrakeRupley()

    # SASA of full complex
    struct_complex = parser.get_structure("cpx", pdb_file)
    sr.compute(struct_complex, level="R")
    sasa_complex = {}
    for model in struct_complex:
        for chain in model:
            for res in chain:
                if res.id[0] == " ":
                    sasa_complex[(chain.id, res.id[1])] = \
                        res.sasa

    # SASA of each chain in isolation
    sasa_isolated = {}
    for chain_id in [chain_a, chain_b]:
        struct_iso = parser.get_structure("iso", pdb_file)
        # Remove the other chain
        for model in struct_iso:
            chains_to_remove = [
                c.id for c in model
                if c.id != chain_id
            ]
            for cid in chains_to_remove:
                model.detach_child(cid)
        sr.compute(struct_iso, level="R")
        for model in struct_iso:
            for chain in model:
                for res in chain:
                    if res.id[0] == " ":
                        sasa_isolated[
                            (chain.id, res.id[1])] = \
                            res.sasa

    # BSA = isolated - complex (clamp to 0)
    bsa = {}
    for key in sasa_complex:
        iso = sasa_isolated.get(key, 0)
        cpx = sasa_complex[key]
        bsa[key] = max(0.0, iso - cpx)

    return bsa


# ══════════════════════════════════════════════════════════
# PART 2 — COMPUTATIONAL ALANINE SCANNING
# ══════════════════════════════════════════════════════════

def get_residue_name(structure, chain_id, resnum):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for res in chain:
                    if res.id[1] == resnum and res.id[0] == " ":
                        return res.resname.strip()
    return "UNK"


def count_contacts_for_residue(structure, chain_id,
                                resnum, partner_chain):
    """
    Count H-bonds, salt bridges, hydrophobic contacts
    that a residue makes with the partner chain.
    This approximates the contact loss upon Ala mutation.
    """
    parser_local = PDBParser(QUIET=True)
    res_atoms    = []
    partner_atoms = []

    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                if chain.id == chain_id and \
                   res.id[1] == resnum:
                    res_atoms.extend(list(res.get_atoms()))
                elif chain.id == partner_chain:
                    partner_atoms.extend(
                        list(res.get_atoms()))

    if not res_atoms or not partner_atoms:
        return 0, 0, 0

    res_name    = get_residue_name(
        structure, chain_id, resnum)
    hbond_count = 0
    hydro_count = 0
    salt_count  = 0

    hbond_donors    = {"N","O","NZ","NH1","NH2","NE",
                       "OG","OG1","OE1","OE2","ND1",
                       "ND2","NE2"}
    hbond_acceptors = {"O","OD1","OD2","OE1","OE2",
                       "OG","OG1","OH","N","NZ"}

    for a1 in res_atoms:
        for a2 in partner_atoms:
            try:
                dist = a1 - a2
            except Exception:
                continue

            # H-bond
            if (dist <= HBOND_DIST and
                a1.name in hbond_donors and
                    a2.name in hbond_acceptors):
                hbond_count += 1

            # Salt bridge
            if dist <= SALTBR_DIST:
                a2_res = a2.get_parent().resname.strip()
                if (res_name in CHARGED_POS and
                        a2_res in CHARGED_NEG):
                    salt_count += 1
                elif (res_name in CHARGED_NEG and
                        a2_res in CHARGED_POS):
                    salt_count += 1

            # Hydrophobic
            if (dist <= HYDRO_DIST and
                res_name in HYDROPHOBIC_AA and
                    a1.element == "C" and
                    a2.element == "C"):
                hydro_count += 1

    return hbond_count, salt_count, hydro_count


def alanine_scanning(pdb_file, chain_a="A", chain_b="B"):
    """
    Estimate contact loss for each hotspot upon Ala mutation.
    Returns dict keyed by (chain, resnum).
    """
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_file)
    results   = {}

    for chain_id, resnum in ALL_HOTSPOTS:
        partner = chain_b if chain_id == chain_a else chain_a
        hb, sb, hy = count_contacts_for_residue(
            structure, chain_id, resnum, partner)
        total = hb + sb + hy
        results[(chain_id, resnum)] = {
            "hbond":       hb,
            "salt_bridge": sb,
            "hydrophobic": hy,
            "total_loss":  total,
            "is_hotspot":  total >= 2,
        }
    return results


# ══════════════════════════════════════════════════════════
# PART 3 — COMPOSITE RANKING
# ══════════════════════════════════════════════════════════

def composite_ranking(bsa, ala_scan,
                      cons10, cons14, interface):
    """
    Composite score = interface_score × conservation
                    × burial_norm × energy_norm
    """
    # Interface scores from Script 05 top20
    iface_scores = {}
    for nsp, key in [("NSP10","top20_NSP10"),
                     ("NSP14","top20_NSP14")]:
        chain = "A" if nsp == "NSP10" else "B"
        for entry in interface.get("PDB_7DIY",{}).get(key,[]):
            resname, score = entry[0], entry[1]
            # Extract residue number
            num = int("".join(filter(str.isdigit, resname)))
            iface_scores[(chain, num)] = score

    # Conservation scores
    cons_scores = {}
    cons10_map = dict(zip(cons10["position"],
                          cons10["conservation"]))
    cons14_map = dict(zip(cons14["position"],
                          cons14["conservation"]))
    for r in HOTSPOTS["NSP10"]:
        cons_scores[("A", r)] = cons10_map.get(r, 0)
    for r in HOTSPOTS["NSP14"]:
        cons_scores[("B", r)] = cons14_map.get(r, 0)

    # Normalize BSA and energy loss
    bsa_vals = [bsa.get(k, 0) for k in ALL_HOTSPOTS]
    ala_vals = [ala_scan.get(k, {}).get(
        "total_loss", 0) for k in ALL_HOTSPOTS]
    max_bsa  = max(bsa_vals) if max(bsa_vals) > 0 else 1
    max_ala  = max(ala_vals) if max(ala_vals) > 0 else 1

    rows = []
    for chain, resnum in ALL_HOTSPOTS:
        nsp     = "NSP10" if chain == "A" else "NSP14"
        bsa_val = bsa.get((chain, resnum), 0)
        ala_val = ala_scan.get(
            (chain, resnum), {}).get("total_loss", 0)
        iface   = iface_scores.get((chain, resnum), 1.0)
        cons    = cons_scores.get((chain, resnum), 0)

        burial_norm = bsa_val / max_bsa
        energy_norm = ala_val / max_ala

        composite = (iface / 10.0) * cons * \
                    (0.5 + 0.5 * burial_norm) * \
                    (0.5 + 0.5 * energy_norm)

        # Get AA name
        aa1   = cons10["aa_SARS2"][
            cons10["position"] == resnum].values
        if len(aa1) == 0:
            aa1 = cons14["aa_SARS2"][
                cons14["position"] == resnum].values
        aa3 = AA3.get(aa1[0] if len(aa1) > 0
                      else "?", "???")

        rows.append({
            "chain":        chain,
            "nsp":          nsp,
            "resnum":       resnum,
            "residue":      f"{aa3}{resnum}",
            "iface_score":  round(iface, 2),
            "conservation": round(cons, 3),
            "bsa":          round(bsa_val, 2),
            "burial_norm":  round(burial_norm, 3),
            "ala_loss":     ala_val,
            "energy_norm":  round(energy_norm, 3),
            "composite":    round(composite, 4),
            "is_primary":   (chain, resnum) in PRIMARY,
        })

    df = pd.DataFrame(rows).sort_values(
        "composite", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


# ══════════════════════════════════════════════════════════
# FIGURES
# ══════════════════════════════════════════════════════════

def fig_bsa(bsa, cons10, cons14):
    """Fig 4 — BSA per hotspot residue"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP14: Buried Surface Area per Hotspot Residue\n"
        "Structure: 7DIY (SARS-CoV-2, 2.69 Å)",
        fontsize=13, fontweight="bold", y=1.02)

    for ax, nsp, df, chain, primary_set in [
        (axes[0], "NSP10", cons10, "A", {80,93}),
        (axes[1], "NSP14", cons14, "B", {126,127}),
    ]:
        hdf  = df[df["is_hotspot"]].copy()
        hdf  = hdf.sort_values("position")
        bsas = [bsa.get((chain, int(r)), 0)
                for r in hdf["position"]]

        colors = []
        for r, b in zip(hdf["position"], bsas):
            if r in primary_set:
                colors.append(RED)
            elif b >= 20:
                colors.append(BLUE)
            else:
                colors.append(GREY)

        labels = []
        for _, row in hdf.iterrows():
            three = AA3.get(row["aa_SARS2"], row["aa_SARS2"])
            lbl   = f"{three}{int(row['position'])}"
            if int(row["position"]) in primary_set:
                lbl = f"★  {lbl}"
            labels.append(lbl)

        y    = np.arange(len(hdf))
        ax.barh(y, bsas, color=colors,
                edgecolor="white", linewidth=0.6,
                height=0.72)
        ax.axvline(20, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6,
                   label="BSA threshold (20 Å²)")
        for i, b in enumerate(bsas):
            ax.text(b + 0.3, i, f"{b:.1f}",
                    va="center", fontsize=8.5)

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9.5)
        ax.set_xlabel("Buried Surface Area (Å²)")
        n_bur = sum(1 for b in bsas if b >= 20)
        ax.set_title(
            f"{nsp}  ({n_bur}/{len(hdf)} residues BSA ≥ 20 Å²)",
            pad=8)
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)

        patches = [
            mpatches.Patch(color=RED,  label="Primary target (★)"),
            mpatches.Patch(color=BLUE, label="Significantly buried ≥ 20 Å²"),
            mpatches.Patch(color=GREY, label="Partially buried < 20 Å²"),
        ]
        ax.legend(handles=patches, loc="lower right",
                  fontsize=8.5)

    fig.text(0.5, -0.02,
             "BSA = SASA(isolated chain) − SASA(complex) "
             "per residue  |  ★ = salt bridge anchor",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig4_NSP10-NSP14_BSA_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


def fig_alascan(ala_scan, cons10, cons14):
    """Fig 5 — Alanine scanning contact loss"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP14: Computational Alanine Scanning\n"
        "Estimated Contact Loss Upon Ala Mutation (7DIY)",
        fontsize=13, fontweight="bold", y=1.02)

    for ax, nsp, df, chain, primary_set in [
        (axes[0], "NSP10", cons10, "A", {80,93}),
        (axes[1], "NSP14", cons14, "B", {126,127}),
    ]:
        hdf  = df[df["is_hotspot"]].copy()
        hdf  = hdf.sort_values("position")
        ress = list(hdf["position"])

        hb_vals = [ala_scan.get((chain,r),{}).get("hbond",0)
                   for r in ress]
        sb_vals = [ala_scan.get((chain,r),{}).get(
            "salt_bridge",0) for r in ress]
        hy_vals = [ala_scan.get((chain,r),{}).get(
            "hydrophobic",0) for r in ress]

        labels = []
        for _, row in hdf.iterrows():
            three = AA3.get(row["aa_SARS2"], row["aa_SARS2"])
            lbl   = f"{three}{int(row['position'])}"
            if int(row["position"]) in primary_set:
                lbl = f"★  {lbl}"
            labels.append(lbl)

        y      = np.arange(len(hdf))
        height = 0.25
        ax.barh(y + height, hb_vals, height,
                label="H-bond", color=BLUE,
                edgecolor="white", linewidth=0.5)
        ax.barh(y, sb_vals, height,
                label="Salt bridge", color=RED,
                edgecolor="white", linewidth=0.5)
        ax.barh(y - height, hy_vals, height,
                label="Hydrophobic", color=ORANGE,
                edgecolor="white", linewidth=0.5)

        ax.axvline(2, color="black", linestyle="--",
                   linewidth=1.0, alpha=0.6,
                   label="Hotspot threshold (≥ 2)")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9.5)
        ax.set_xlabel("Contacts lost upon Ala mutation")
        ax.set_title(nsp, pad=8)
        ax.legend(loc="lower right", fontsize=8.5)
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)

    fig.text(0.5, -0.02,
             "Contact loss ≥ 2 = energetic hotspot  |  "
             "★ = primary salt bridge anchor",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig5_NSP10-NSP14_AlaScan_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


def fig_composite(df_rank):
    """Fig 6 — Composite hotspot ranking"""
    fig, ax = plt.subplots(figsize=(12, 8),
                           constrained_layout=True)
    fig.suptitle(
        "NSP10–NSP14: Composite Hotspot Drug Target Ranking\n"
        "Score = Interface × Conservation × Burial × Energy",
        fontsize=13, fontweight="bold")

    top = df_rank.sort_values(
        "composite", ascending=True).tail(20)

    colors = []
    for _, row in top.iterrows():
        if row["is_primary"]:
            colors.append(RED)
        elif row["nsp"] == "NSP10":
            colors.append(BLUE)
        else:
            colors.append(GREEN)

    y     = np.arange(len(top))
    bars  = ax.barh(y, top["composite"], color=colors,
                    edgecolor="white", linewidth=0.6,
                    height=0.72)

    for i, (_, row) in enumerate(top.iterrows()):
        ax.text(row["composite"] + 0.0005, i,
                f"{row['composite']:.3f}  "
                f"(cons={row['conservation']:.2f}, "
                f"BSA={row['bsa']:.0f}Å², "
                f"loss={row['ala_loss']})",
                va="center", fontsize=8)

    ax.set_yticks(y)
    ylabels = []
    for _, row in top.iterrows():
        lbl = f"#{int(row['rank'])}  {row['residue']}  ({row['nsp']})"
        if row["is_primary"]:
            lbl = f"★  {lbl}"
        ylabels.append(lbl)
    ax.set_yticklabels(ylabels, fontsize=9.5)

    ax.set_xlabel(
        "Composite Score  (interface × conservation × burial × energy)")
    ax.grid(axis="x", alpha=0.2, linewidth=0.6)

    patches = [
        mpatches.Patch(color=RED,   label="Primary salt bridge residue (★)"),
        mpatches.Patch(color=BLUE,  label="NSP10 hotspot"),
        mpatches.Patch(color=GREEN, label="NSP14 hotspot"),
    ]
    ax.legend(handles=patches, loc="lower right", fontsize=9)

    fig.text(0.5, -0.02,
             "Top-ranked residues = highest priority for "
             "docking box definition and virtual screening",
             ha="center", fontsize=9,
             style="italic", color="#555555")

    out = OUT_DIR / "Fig6_NSP10-NSP14_composite_ranking_2.png"
    fig.savefig(out)
    print(f"  Saved: {out.name}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════

def main():
    print("\n" + "="*55)
    print("  Script 10_2: BSA + Ala Scan + Ranking")
    print("  Structure: 7DIY (SARS-CoV-2, 2.69 Å)")
    print("="*55 + "\n")

    pdb_file = PDB_DIR / "7DIY.pdb"
    cons10   = pd.read_csv(RES_DIR / "conservation_NSP10.csv")
    cons14   = pd.read_csv(RES_DIR / "conservation_NSP14.csv")
    with open(RES_DIR / "interface_analysis.json") as f:
        interface = json.load(f)

    # ── Part 1: BSA ───────────────────────────────────────
    print("  Part 1: Computing Buried Surface Area...")
    bsa = compute_bsa(pdb_file)
    print(f"  BSA computed for {len(bsa)} residues")
    print("\n  Hotspot BSA values (Å²):")
    for chain, resnum in ALL_HOTSPOTS:
        b = bsa.get((chain, resnum), 0)
        flag = " ✅" if b >= 20 else ""
        print(f"    Chain {chain} Res {resnum:>4} : "
              f"{b:>7.2f}{flag}")

    # ── Part 2: Alanine scanning ──────────────────────────
    print("\n  Part 2: Computational alanine scanning...")
    ala = alanine_scanning(pdb_file)
    print("\n  Contact loss upon Ala mutation:")
    print(f"  {'Residue':<12} {'HBond':>6} "
          f"{'SaltBr':>7} {'Hydro':>7} "
          f"{'Total':>7}  Hotspot?")
    print(f"  {'-'*12} {'-'*6} {'-'*7} "
          f"{'-'*7} {'-'*7}  {'-'*8}")
    for chain, resnum in ALL_HOTSPOTS:
        r    = ala.get((chain, resnum), {})
        flag = " ✅" if r.get("is_hotspot") else ""
        print(f"  {chain}{resnum:<11} "
              f"{r.get('hbond',0):>6} "
              f"{r.get('salt_bridge',0):>7} "
              f"{r.get('hydrophobic',0):>7} "
              f"{r.get('total_loss',0):>7}"
              f"{flag}")

    # ── Part 3: Composite ranking ─────────────────────────
    print("\n  Part 3: Composite ranking...")
    df_rank = composite_ranking(bsa, ala, cons10,
                                cons14, interface)
    print(f"\n  Top 10 drug target candidates:")
    print(f"  {'Rank':<5} {'Residue':<10} {'NSP':<7} "
          f"{'Composite':>10} {'Cons':>6} "
          f"{'BSA':>8} {'Loss':>6}")
    print(f"  {'-'*5} {'-'*10} {'-'*7} "
          f"{'-'*10} {'-'*6} {'-'*8} {'-'*6}")
    for _, row in df_rank.head(10).iterrows():
        star = "★" if row["is_primary"] else " "
        print(f"  {star}{int(row['rank']):<4} "
              f"{row['residue']:<10} "
              f"{row['nsp']:<7} "
              f"{row['composite']:>10.4f} "
              f"{row['conservation']:>6.3f} "
              f"{row['bsa']:>7.1f}Å² "
              f"{int(row['ala_loss']):>5}")

    # ── Save results ──────────────────────────────────────
    df_rank.to_csv(
        RES_DIR / "composite_ranking_NSP10-NSP14_2.csv",
        index=False)
    bsa_out = {f"{c}{r}": round(v,2)
               for (c,r), v in bsa.items()}
    ala_out = {f"{c}{r}": v
               for (c,r), v in ala.items()}
    with open(RES_DIR /
              "bsa_alascan_NSP10-NSP14_2.json","w") as f:
        json.dump({"bsa": bsa_out, "ala_scan": ala_out},
                  f, indent=2)

    # ── Figures ───────────────────────────────────────────
    print("\n  Generating Figure 4 — BSA...")
    fig_bsa(bsa, cons10, cons14)

    print("  Generating Figure 5 — Alanine scanning...")
    fig_alascan(ala, cons10, cons14)

    print("  Generating Figure 6 — Composite ranking...")
    fig_composite(df_rank)

    print(f"\n  Files saved:")
    print(f"  02-validation/NSP10-NSP14/"
          f"composite_ranking_NSP10-NSP14_2.csv")
    print(f"  02-validation/NSP10-NSP14/"
          f"bsa_alascan_NSP10-NSP14_2.json")
    print(f"  results/Fig4_NSP10-NSP14_BSA_2.png")
    print(f"  results/Fig5_NSP10-NSP14_AlaScan_2.png")
    print(f"  results/Fig6_NSP10-NSP14_composite_ranking_2.png")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
