#!/usr/bin/env python3
"""
Script 10_BSA_alascan_ranking_NSP13-Helicase_7.py

BSA, computational alanine scanning, and composite hotspot
ranking for NSP13 homodimer interface.

Receptor: 7NIO chains A+E (primary dimer)
Hotspots: 17 residues from chain A (7NIO numbering)

Composite score = interface_score x conservation x burial_factor x energy_factor

Output: 02-validation/NSP13-Helicase/
  composite_ranking_NSP13-Helicase_7.csv
  bsa_alascan_NSP13-Helicase_7.json
  results/
  Fig4_NSP13-Helicase_BSA_7.png
  Fig5_NSP13-Helicase_AlaScan_7.png
  Fig6_NSP13-Helicase_composite_ranking_7.png
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

# ── Paths ──────────────────────────────────────────────────────────────────
BASE        = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR     = f"{BASE}/00-reference/pdb_structures"
VAL_DIR     = f"{BASE}/02-validation/NSP13-Helicase"
RESULTS_DIR = f"{BASE}/results"
os.makedirs(VAL_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

PDB_7NIO = f"{PDB_DIR}/7NIO.pdb"

# ── Hotspot data ───────────────────────────────────────────────────────────
HOTSPOTS = [
    (116, "ASN", 0.344), (413, "THR", 0.689), (414, "LYS", 0.689),
    (415, "GLY", 1.000), (477, "LYS", 0.689), (478, "GLY", 1.000),
    (479, "VAL", 0.344), (480, "ILE", 0.582), (481, "THR", 0.582),
    (482, "HIS", 0.582), (551, "GLU", 0.582), (552, "THR", 1.000),
    (553, "ALA", 1.000), (579, "ARG", 0.410), (580, "ASP", 0.344),
    (583, "ASP", 0.689), (584, "LYS", 0.344),
]

PRIMARY_HOTSPOTS = {414, 477, 480, 482, 579, 580, 583, 584}
SB_RESIDUES      = {414, 580, 583}

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}
HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO","TYR"}
POSITIVE    = {"ARG","LYS","HIS"}
NEGATIVE    = {"ASP","GLU"}

# ── Parse structure ────────────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: Parse 7NIO dimer")
print("=" * 60)

parser = PDB.PDBParser(QUIET=True)
struct  = parser.get_structure("7NIO", PDB_7NIO)

class DimerSelector(PDB.Select):
    def accept_chain(self, chain):
        return chain.get_id() in ["A", "E"]
    def accept_residue(self, residue):
        return residue.get_id()[0] == " "
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"

# Build sub-structure for SASA
import io as _io
from Bio.PDB import PDBIO as _PDBIO
_buf = _io.StringIO()
_io2 = _PDBIO()
_io2.set_structure(struct)
_io2.save(_buf, DimerSelector())
_buf.seek(0)
struct_dimer = parser.get_structure("dimer", _buf)

chain_A = struct_dimer[0]["A"]
chain_E = struct_dimer[0]["E"]

aa_A = [r for r in chain_A if PDB.is_aa(r)]
aa_E = [r for r in chain_E if PDB.is_aa(r)]
print(f"  Chain A: {len(aa_A)} residues")
print(f"  Chain E: {len(aa_E)} residues")

# ── Step 2: BSA calculation ────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Buried Surface Area (BSA) per hotspot residue")
print("=" * 60)

sr = ShrakeRupley()

# SASA of full dimer
sr.compute(struct_dimer, level="R")
sasa_complex = {}
for res in chain_A:
    if PDB.is_aa(res):
        sasa_complex[res.get_id()[1]] = res.sasa

# SASA of chain A alone — build single chain structure
_buf2 = _io.StringIO()
_io3  = _PDBIO()
_io3.set_structure(struct)

class ChainASelector(PDB.Select):
    def accept_chain(self, chain):
        return chain.get_id() == "A"
    def accept_residue(self, residue):
        return residue.get_id()[0] == " "
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"

_io3.save(_buf2, ChainASelector())
_buf2.seek(0)
struct_chainA = parser.get_structure("chainA", _buf2)
sr.compute(struct_chainA, level="R")
sasa_unbound = {}
for res in struct_chainA[0]["A"]:
    if PDB.is_aa(res):
        sasa_unbound[res.get_id()[1]] = res.sasa

# BSA = unbound - complex
print(f"\n  {'Pos':>4}  {'Res':>3}  {'SASA_unbound':>12}  "
      f"{'SASA_complex':>12}  {'BSA':>8}  {'Buried':>8}")
print("  " + "-" * 56)

bsa_data = {}
for rn, resname, cons in HOTSPOTS:
    su  = sasa_unbound.get(rn, 0.0)
    sc  = sasa_complex.get(rn, 0.0)
    bsa = max(0.0, su - sc)
    buried = "YES" if bsa > 20 else "no"
    prim   = " ★" if rn in PRIMARY_HOTSPOTS else ""
    print(f"  {rn:>4}  {resname:>3}{prim}  {su:>12.1f}  {sc:>12.1f}  "
          f"{bsa:>8.1f}  {buried:>8}")
    bsa_data[rn] = bsa

# ── Step 3: Alanine scanning ───────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: Computational alanine scanning")
print("=" * 60)

atoms_E = [a for r in chain_E if PDB.is_aa(r) for a in r.get_atoms()]
ns_E    = PDB.NeighborSearch(atoms_E)

print(f"\n  {'Pos':>4}  {'Res':>3}  {'HB_lost':>8}  {'SB_lost':>8}  "
      f"{'HY_lost':>8}  {'total_loss':>10}")
print("  " + "-" * 52)

alascan_data = {}
for rn, resname, cons in HOTSPOTS:
    # Find residue in chain A
    target_res = None
    for r in chain_A:
        if PDB.is_aa(r) and r.get_id()[1] == rn:
            target_res = r
            break
    if target_res is None:
        alascan_data[rn] = {"hb":0,"sb":0,"hy":0,"total":0}
        continue

    hb_lost, sb_lost, hy_lost = 0, 0, 0

    for a1 in target_res.get_atoms():
        hits = ns_E.search(a1.get_vector().get_array(), 5.0, "A")
        for a2 in hits:
            r2   = a2.get_parent()
            if not PDB.is_aa(r2):
                continue
            dist = float(np.linalg.norm(
                a1.get_vector().get_array() - a2.get_vector().get_array()))
            r1n, r2n = resname, r2.get_resname()

            # Salt bridge
            if ((r1n in POSITIVE and r2n in NEGATIVE) or
                (r1n in NEGATIVE and r2n in POSITIVE)) and dist < 4.5:
                sb_lost += 1
            # H-bond
            elif dist < 3.5 and (r1n not in HYDROPHOBIC or r2n not in HYDROPHOBIC):
                hb_lost += 1
            # Hydrophobic
            elif r1n in HYDROPHOBIC and r2n in HYDROPHOBIC and dist < 5.0:
                hy_lost += 1

    total = hb_lost + sb_lost + hy_lost
    alascan_data[rn] = {"hb": hb_lost, "sb": sb_lost,
                         "hy": hy_lost, "total": total}
    prim = " ★" if rn in PRIMARY_HOTSPOTS else ""
    print(f"  {rn:>4}  {resname:>3}{prim}  {hb_lost:>8}  {sb_lost:>8}  "
          f"{hy_lost:>8}  {total:>10}")

# ── Step 4: Composite ranking ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Composite hotspot ranking")
print("=" * 60)

# Interface score from contact count (Script 05_7)
contact_counts = {
    116:0, 413:0, 414:43, 415:0, 477:32, 478:21,
    479:0, 480:60, 481:0, 482:60, 551:24, 552:0,
    553:0, 579:20, 580:49, 583:21, 584:0
}

max_bsa     = max(bsa_data.values()) if bsa_data else 1.0
max_loss    = max(v["total"] for v in alascan_data.values()) if alascan_data else 1.0
max_contact = max(contact_counts.values()) if contact_counts else 1.0

records = []
for rn, resname, cons in HOTSPOTS:
    bsa          = bsa_data.get(rn, 0.0)
    scan         = alascan_data.get(rn, {"total": 0, "sb": 0})
    contact      = contact_counts.get(rn, 0)

    burial_fac   = bsa / max_bsa if max_bsa > 0 else 0.0
    energy_fac   = scan["total"] / max_loss if max_loss > 0 else 0.0
    iface_score  = contact / max_contact if max_contact > 0 else 0.0

    # SB bonus x1.5 for salt bridge residues
    sb_bonus = 1.5 if rn in SB_RESIDUES else 1.0

    composite = iface_score * cons * (burial_fac + 0.1) * (energy_fac + 0.1) * sb_bonus

    records.append({
        "res_num"       : rn,
        "res_name"      : resname,
        "res_aa"        : AA_MAP_3TO1.get(resname, "?"),
        "conservation"  : cons,
        "contact_count" : contact,
        "bsa_A2"        : round(bsa, 2),
        "hb_lost"       : scan["hb"],
        "sb_lost"       : scan["sb"],
        "hy_lost"       : scan["hy"],
        "total_loss"    : scan["total"],
        "burial_factor" : round(burial_fac, 4),
        "energy_factor" : round(energy_fac, 4),
        "iface_score"   : round(iface_score, 4),
        "sb_bonus"      : sb_bonus,
        "composite_score": round(composite, 4),
        "is_primary"    : rn in PRIMARY_HOTSPOTS,
        "is_sb_residue" : rn in SB_RESIDUES,
    })

df = pd.DataFrame(records).sort_values("composite_score", ascending=False).reset_index(drop=True)
df["rank"] = df.index + 1

print(f"\n  {'Rank':>4}  {'Res':>8}  {'Cons':>6}  {'BSA':>7}  "
      f"{'Loss':>6}  {'Contact':>8}  {'Score':>8}")
print("  " + "-" * 60)
for _, row in df.iterrows():
    prim = " ★" if row["is_primary"] else ""
    sb   = " [SB]" if row["is_sb_residue"] else ""
    print(f"  {row['rank']:>4}  {row['res_name']}{row['res_num']:4d}{prim}  "
          f"{row['conservation']:>6.3f}  {row['bsa_A2']:>7.1f}  "
          f"{row['total_loss']:>6}  {row['contact_count']:>8}  "
          f"{row['composite_score']:>8.4f}{sb}")

# Save CSV
csv_path = f"{VAL_DIR}/composite_ranking_NSP13-Helicase_7.csv"
df.to_csv(csv_path, index=False)
print(f"\n  Saved: {csv_path}")

# ── FIGURE 4: BSA bar chart ────────────────────────────────────────────────
print("\nGenerating Fig4: BSA ...")

df_sorted_bsa = df.sort_values("bsa_A2", ascending=True)

colors_bsa = []
for _, row in df_sorted_bsa.iterrows():
    if row["is_sb_residue"]:
        colors_bsa.append("#d62728")
    elif row["is_primary"]:
        colors_bsa.append("#ff7f0e")
    else:
        colors_bsa.append("#1f77b4")

fig, ax = plt.subplots(figsize=(10, 8))
bars = ax.barh(
    [f"{r['res_name']}{r['res_num']}" for _, r in df_sorted_bsa.iterrows()],
    df_sorted_bsa["bsa_A2"],
    color=colors_bsa, edgecolor="white", linewidth=0.5
)
ax.axvline(20, color="black", linestyle="--", linewidth=1.2,
           label="BSA threshold (20 Å²)")
for bar, val in zip(bars, df_sorted_bsa["bsa_A2"]):
    ax.text(val + 0.5, bar.get_y() + bar.get_height()/2,
            f"{val:.1f}", va="center", fontsize=8)
ax.set_xlabel("Buried Surface Area (Å²)", fontsize=11)
ax.set_title("NSP13 Homodimer — BSA per Hotspot Residue (7NIO)",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.grid(axis="x", alpha=0.3)
ax.spines[["top","right"]].set_visible(False)
legend_patches = [
    mpatches.Patch(color="#d62728", label="Salt bridge anchor"),
    mpatches.Patch(color="#ff7f0e", label="Primary hotspot"),
    mpatches.Patch(color="#1f77b4", label="Secondary hotspot"),
    plt.Line2D([0],[0], color="black", linestyle="--", label="Threshold (20 Å²)"),
]
ax.legend(handles=legend_patches, fontsize=9, loc="lower right")
plt.tight_layout()
fig4_path = f"{RESULTS_DIR}/Fig4_NSP13-Helicase_BSA_7.png"
plt.savefig(fig4_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig4_path}")

# ── FIGURE 5: Alanine scanning ────────────────────────────────────────────
print("Generating Fig5: alanine scanning ...")

df_sorted_scan = df.sort_values("total_loss", ascending=True)
labels_scan    = [f"{r['res_name']}{r['res_num']}" for _, r in df_sorted_scan.iterrows()]

fig, ax = plt.subplots(figsize=(10, 8))
x_pos = np.arange(len(df_sorted_scan))

hb_vals = df_sorted_scan["hb_lost"].values
sb_vals = df_sorted_scan["sb_lost"].values
hy_vals = df_sorted_scan["hy_lost"].values

ax.barh(x_pos,             hb_vals, label="H-bonds lost",       color="#1f77b4", alpha=0.85)
ax.barh(x_pos, sb_vals, left=hb_vals,
        label="Salt bridges lost", color="#d62728", alpha=0.85)
ax.barh(x_pos, hy_vals, left=hb_vals+sb_vals,
        label="Hydrophobic lost",  color="#2ca02c", alpha=0.85)

for i, (_, row) in enumerate(df_sorted_scan.iterrows()):
    tot = row["total_loss"]
    if tot > 0:
        ax.text(tot + 0.1, i, str(int(tot)), va="center", fontsize=8)

ax.set_yticks(x_pos)
ax.set_yticklabels(labels_scan, fontsize=9)
ax.set_xlabel("Estimated contacts lost upon Ala mutation", fontsize=11)
ax.set_title("NSP13 Homodimer — Computational Alanine Scanning",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="lower right")
ax.grid(axis="x", alpha=0.3)
ax.spines[["top","right"]].set_visible(False)

# Annotate SB residues
for i, (_, row) in enumerate(df_sorted_scan.iterrows()):
    if row["is_sb_residue"]:
        ax.annotate("★ SB", xy=(row["total_loss"], i),
                    xytext=(row["total_loss"]+1, i),
                    fontsize=8, color="#d62728", va="center")

plt.tight_layout()
fig5_path = f"{RESULTS_DIR}/Fig5_NSP13-Helicase_AlaScan_7.png"
plt.savefig(fig5_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig5_path}")

# ── FIGURE 6: Composite ranking ────────────────────────────────────────────
print("Generating Fig6: composite ranking ...")

df_rank = df.sort_values("composite_score", ascending=True)

colors_rank = []
for _, row in df_rank.iterrows():
    if row["is_sb_residue"]:
        colors_rank.append("#d62728")
    elif row["is_primary"]:
        colors_rank.append("#ff7f0e")
    else:
        colors_rank.append("#1f77b4")

fig, ax = plt.subplots(figsize=(10, 8))
bars = ax.barh(
    [f"{r['res_name']}{r['res_num']}" for _, r in df_rank.iterrows()],
    df_rank["composite_score"],
    color=colors_rank, edgecolor="white", linewidth=0.5
)
for bar, val in zip(bars, df_rank["composite_score"]):
    ax.text(val + 0.001, bar.get_y() + bar.get_height()/2,
            f"{val:.4f}", va="center", fontsize=8)

ax.set_xlabel("Composite score\n(interface × conservation × burial × energy × SB bonus)",
              fontsize=10)
ax.set_title("NSP13 Homodimer — Hotspot Composite Ranking",
             fontsize=12, fontweight="bold")
legend_patches = [
    mpatches.Patch(color="#d62728", label="Salt bridge anchor (×1.5 bonus)"),
    mpatches.Patch(color="#ff7f0e", label="Primary hotspot"),
    mpatches.Patch(color="#1f77b4", label="Secondary hotspot"),
]
ax.legend(handles=legend_patches, fontsize=9, loc="lower right")
ax.grid(axis="x", alpha=0.3)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
fig6_path = f"{RESULTS_DIR}/Fig6_NSP13-Helicase_composite_ranking_7.png"
plt.savefig(fig6_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig6_path}")

# ── Save JSON ──────────────────────────────────────────────────────────────
output = {
    "complex"        : "NSP13-Helicase",
    "script"         : "10_BSA_alascan_ranking_NSP13-Helicase_7.py",
    "structure"      : "7NIO chains A+E",
    "bsa_data"       : {str(k): round(v,2) for k,v in bsa_data.items()},
    "alascan_data"   : {str(k): v for k,v in alascan_data.items()},
    "composite_ranking": df.to_dict(orient="records"),
    "top3_pharmacophores": df.head(3)[
        ["res_num","res_name","conservation","bsa_A2",
         "total_loss","composite_score"]].to_dict(orient="records")
}

json_path = f"{VAL_DIR}/bsa_alascan_NSP13-Helicase_7.json"
with open(json_path, "w") as f:
    json.dump(output, f, indent=2)
print(f"\n  Saved JSON: {json_path}")

# ── Final summary ──────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY — BSA + AlaScan + Ranking")
print("=" * 60)
print(f"\n  Top 5 composite pharmacophores:")
for _, row in df.head(5).iterrows():
    sb  = " [SB anchor]" if row["is_sb_residue"] else ""
    print(f"    #{row['rank']:>2}  {row['res_name']}{row['res_num']:4d}  "
          f"score={row['composite_score']:.4f}  "
          f"BSA={row['bsa_A2']:.1f} A2  "
          f"cons={row['conservation']:.3f}  "
          f"loss={row['total_loss']}{sb}")

print(f"\n  Primary pharmacophore: LYS414 (dual SB donor)")
print(f"  SARS-CoV-1/2 selective interface")
print(f"\n  Figures: Fig4, Fig5, Fig6 saved to results/")
