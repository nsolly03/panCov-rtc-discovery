#!/usr/bin/env python3
"""
Script 10_BSA_alascan_ranking_NSP12-NSP13_8.py

BSA calculation, alanine scanning, and composite hotspot
ranking for NSP12-NSP13 interface.

Primary receptor: 7RDY (best resolution 3.10 A)
Supporting: 6XEZ, 7CXM for cross-validation

Hotspots (6XEZ PDB numbering):
  NSP12: LEU900, ASP901★(SB), MET902★(hydrophobic), TYR903, SER904
  NSP13: PHE90, GLY91, LEU92, TYR93, LYS94★(SB pan-cov), ASN95, THR96

Output: 02-validation/NSP12-NSP13/
  composite_ranking_NSP12-NSP13_8.csv
  bsa_alascan_NSP12-NSP13_8.json
  results/Fig4_NSP12-NSP13_BSA_8.png
  results/Fig5_NSP12-NSP13_AlaScan_8.png
  results/Fig6_NSP12-NSP13_composite_ranking_8.png
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
VAL_DIR     = f"{BASE}/02-validation/NSP12-NSP13"
RESULTS_DIR = f"{BASE}/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

PDB_7RDY = f"{PDB_DIR}/7RDY.pdb"
PDB_6XEZ = f"{PDB_DIR}/6XEZ.pdb"
PDB_7CXM = f"{PDB_DIR}/7CXM.pdb"

NSP12_HOTSPOTS = [
    (900, "LEU", 1.000, False),
    (901, "ASP", 0.582, True),   # SB primary
    (902, "MET", 0.582, True),   # hydrophobic primary
    (903, "TYR", 0.582, False),
    (904, "SER", 1.000, False),
]
NSP13_HOTSPOTS = [
    (90,  "PHE", 1.000, False),
    (91,  "GLY", 1.000, False),
    (92,  "LEU", 1.000, False),
    (93,  "TYR", 1.000, False),
    (94,  "LYS", 1.000, True),   # SB primary pan-cov
    (95,  "ASN", 0.689, False),
    (96,  "THR", 0.344, False),
]

SB_NSP12  = {901}
SB_NSP13  = {94}
PRIM_NSP12 = {901, 902}
PRIM_NSP13 = {94}

COL_SB  = "#d62728"
COL_HY  = "#ff7f0e"
COL_DEF = "#1f77b4"
COL_NSP13 = "#e377c2"

# ── Helper: SASA on complex vs isolated chain ──────────────────────────────
def compute_bsa(pdb_file, chain1_id, chain2_id, res_list_chain1):
    """
    BSA = SASA(isolated chain1) - SASA(chain1 in complex)
    Returns dict {res_num: bsa_value}
    """
    parser = PDB.PDBParser(QUIET=True)
    sr     = ShrakeRupley()

    # Complex
    struct_complex = parser.get_structure("cx", pdb_file)
    class TwoChain(PDB.Select):
        def accept_chain(self, c):
            return c.get_id() in (chain1_id, chain2_id)
        def accept_residue(self, r):
            return r.get_id()[0] == " "
    io = PDB.PDBIO()
    io.set_structure(struct_complex)
    tmp_cx = "/tmp/bsa_complex.pdb"
    io.save(tmp_cx, TwoChain())
    s_cx = parser.get_structure("cx2", tmp_cx)
    sr.compute(s_cx, level="R")
    sasa_complex = {r.get_id()[1]: r.sasa
                    for r in s_cx[0][chain1_id] if PDB.is_aa(r)}

    # Isolated chain1
    class OneChain(PDB.Select):
        def accept_chain(self, c):
            return c.get_id() == chain1_id
        def accept_residue(self, r):
            return r.get_id()[0] == " "
    tmp_iso = "/tmp/bsa_isolated.pdb"
    io.set_structure(struct_complex)
    io.save(tmp_iso, OneChain())
    s_iso = parser.get_structure("iso", tmp_iso)
    sr.compute(s_iso, level="R")
    sasa_isolated = {r.get_id()[1]: r.sasa
                     for r in s_iso[0][chain1_id] if PDB.is_aa(r)}

    bsa = {}
    for rn in res_list_chain1:
        iso = sasa_isolated.get(rn, 0.0)
        cx  = sasa_complex.get(rn, 0.0)
        bsa[rn] = max(0.0, round(iso - cx, 2))
    return bsa

# ── Step 1: BSA calculation ────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: BSA calculation (7RDY primary, 6XEZ + 7CXM validation)")
print("=" * 60)

nsp12_resnums = [r[0] for r in NSP12_HOTSPOTS]
nsp13_resnums = [r[0] for r in NSP13_HOTSPOTS]

print("\n  Computing NSP12 BSA (chain A in complex with E) ...")
bsa_nsp12_7rdy = compute_bsa(PDB_7RDY, "A", "E", nsp12_resnums)
bsa_nsp12_6xez = compute_bsa(PDB_6XEZ, "A", "E", nsp12_resnums)
bsa_nsp12_7cxm = compute_bsa(PDB_7CXM, "A", "E", nsp12_resnums)

print("\n  Computing NSP13 BSA (chain E in complex with A) ...")
bsa_nsp13_7rdy = compute_bsa(PDB_7RDY, "E", "A", nsp13_resnums)
bsa_nsp13_6xez = compute_bsa(PDB_6XEZ, "E", "A", nsp13_resnums)
bsa_nsp13_7cxm = compute_bsa(PDB_7CXM, "E", "A", nsp13_resnums)

# Average BSA across 3 structures
def avg_bsa(bsa_dicts, rn):
    vals = [d.get(rn, 0.0) for d in bsa_dicts]
    return round(np.mean(vals), 2)

print(f"\n  NSP12 BSA (avg 3 structures):")
print(f"  {'Residue':>10}  {'7RDY':>8}  {'6XEZ':>8}  {'7CXM':>8}  {'Mean':>8}")
nsp12_bsa_avg = {}
for rn, resname, _, _ in NSP12_HOTSPOTS:
    v7rdy = bsa_nsp12_7rdy.get(rn, 0.0)
    v6xez = bsa_nsp12_6xez.get(rn, 0.0)
    v7cxm = bsa_nsp12_7cxm.get(rn, 0.0)
    avg   = round(np.mean([v7rdy, v6xez, v7cxm]), 2)
    nsp12_bsa_avg[rn] = avg
    print(f"  {resname+str(rn):>10}  {v7rdy:>8.1f}  {v6xez:>8.1f}  {v7cxm:>8.1f}  {avg:>8.1f}")

print(f"\n  NSP13 BSA (avg 3 structures):")
print(f"  {'Residue':>10}  {'7RDY':>8}  {'6XEZ':>8}  {'7CXM':>8}  {'Mean':>8}")
nsp13_bsa_avg = {}
for rn, resname, _, _ in NSP13_HOTSPOTS:
    v7rdy = bsa_nsp13_7rdy.get(rn, 0.0)
    v6xez = bsa_nsp13_6xez.get(rn, 0.0)
    v7cxm = bsa_nsp13_7cxm.get(rn, 0.0)
    avg   = round(np.mean([v7rdy, v6xez, v7cxm]), 2)
    nsp13_bsa_avg[rn] = avg
    print(f"  {resname+str(rn):>10}  {v7rdy:>8.1f}  {v6xez:>8.1f}  {v7cxm:>8.1f}  {avg:>8.1f}")

# ── Step 2: Alanine scanning (contact loss) ────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Alanine scanning (7RDY primary)")
print("=" * 60)

def alanine_scan_contacts(pdb_file, chain_ref, chain_partner,
                           hotspot_resnums, cutoff=5.0):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_file)
    chain_r = struct[0][chain_ref]
    chain_p = struct[0][chain_partner]

    atoms_p = [a for r in chain_p if PDB.is_aa(r) for a in r.get_atoms()]
    ns      = PDB.NeighborSearch(atoms_p)

    results = {}
    for rn in hotspot_resnums:
        res = None
        for r in chain_r:
            if PDB.is_aa(r) and r.get_id()[1] == rn:
                res = r
                break
        if res is None:
            results[rn] = 0
            continue
        # Contacts made by this residue (excluding backbone CB/CA/N/C/O for ALA)
        sc_atoms = [a for a in res.get_atoms()
                    if a.get_name() not in ("N","CA","C","O")]
        lost = 0
        for a in sc_atoms:
            hits = ns.search(a.get_vector().get_array(), cutoff, "A")
            lost += len(hits)
        results[rn] = lost
    return results

print("\n  NSP12 alanine scan ...")
alascan_nsp12 = alanine_scan_contacts(PDB_7RDY, "A", "E", nsp12_resnums)
print("\n  NSP13 alanine scan ...")
alascan_nsp13 = alanine_scan_contacts(PDB_7RDY, "E", "A", nsp13_resnums)

print(f"\n  NSP12 contacts lost on Ala mutation:")
for rn, resname, _, _ in NSP12_HOTSPOTS:
    loss = alascan_nsp12.get(rn, 0)
    prim = " ★" if rn in PRIM_NSP12 else ""
    print(f"    {resname}{rn:4d}  loss={loss:3d}{prim}")

print(f"\n  NSP13 contacts lost on Ala mutation:")
for rn, resname, _, _ in NSP13_HOTSPOTS:
    loss = alascan_nsp13.get(rn, 0)
    prim = " ★" if rn in PRIM_NSP13 else ""
    print(f"    {resname}{rn:4d}  loss={loss:3d}{prim}")

# ── Step 3: Composite ranking ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: Composite ranking")
print("=" * 60)

def composite_score(bsa, ala_loss, conservation, bsa_max, ala_max):
    bsa_n  = bsa      / bsa_max if bsa_max > 0 else 0
    ala_n  = ala_loss / ala_max if ala_max > 0 else 0
    return round(0.4*bsa_n + 0.4*ala_n + 0.2*conservation, 4)

all_bsa  = list(nsp12_bsa_avg.values()) + list(nsp13_bsa_avg.values())
all_ala  = list(alascan_nsp12.values()) + list(alascan_nsp13.values())
bsa_max  = max(all_bsa) if all_bsa else 1
ala_max  = max(all_ala) if all_ala else 1

rows = []
for rn, resname, cons, is_prim in NSP12_HOTSPOTS:
    bsa  = nsp12_bsa_avg.get(rn, 0)
    ala  = alascan_nsp12.get(rn, 0)
    sc   = composite_score(bsa, ala, cons, bsa_max, ala_max)
    rows.append({"nsp": "NSP12", "res_num": rn, "res_name": resname,
                 "bsa": bsa, "ala_loss": ala, "conservation": cons,
                 "is_primary": is_prim, "composite": sc})

for rn, resname, cons, is_prim in NSP13_HOTSPOTS:
    bsa  = nsp13_bsa_avg.get(rn, 0)
    ala  = alascan_nsp13.get(rn, 0)
    sc   = composite_score(bsa, ala, cons, bsa_max, ala_max)
    rows.append({"nsp": "NSP13", "res_num": rn, "res_name": resname,
                 "bsa": bsa, "ala_loss": ala, "conservation": cons,
                 "is_primary": is_prim, "composite": sc})

df = pd.DataFrame(rows).sort_values("composite", ascending=False).reset_index(drop=True)
df["rank"] = df.index + 1

print(f"\n  {'Rank':>4}  {'NSP':>5}  {'Residue':>10}  {'BSA':>7}  "
      f"{'AlaLoss':>8}  {'Cons':>6}  {'Composite':>10}  {'Primary':>8}")
print("  " + "-"*65)
for _, row in df.iterrows():
    prim = "★" if row["is_primary"] else ""
    print(f"  {int(row['rank']):>4}  {row['nsp']:>5}  "
          f"{row['res_name']+str(row['res_num']):>10}  "
          f"{row['bsa']:>7.1f}  {int(row['ala_loss']):>8}  "
          f"{row['conservation']:>6.3f}  {row['composite']:>10.4f}  {prim:>8}")

# ── Save CSV + JSON ────────────────────────────────────────────────────────
csv_path = f"{VAL_DIR}/composite_ranking_NSP12-NSP13_8.csv"
df.to_csv(csv_path, index=False)
print(f"\n  Saved: {csv_path}")

json_out = {
    "complex"        : "NSP12-NSP13",
    "script"         : "10_BSA_alascan_ranking_NSP12-NSP13_8.py",
    "bsa_nsp12_avg"  : nsp12_bsa_avg,
    "bsa_nsp13_avg"  : nsp13_bsa_avg,
    "alascan_nsp12"  : alascan_nsp12,
    "alascan_nsp13"  : alascan_nsp13,
    "composite_ranking": df.to_dict(orient="records"),
    "top3"           : df.head(3)[["nsp","res_name","res_num","composite"]].to_dict(orient="records"),
}
json_path = f"{VAL_DIR}/bsa_alascan_NSP12-NSP13_8.json"
with open(json_path, "w") as f:
    json.dump(json_out, f, indent=2)
print(f"  Saved: {json_path}")

# ── Figures 4-6 ────────────────────────────────────────────────────────────
def get_color(nsp, rn):
    if (nsp == "NSP12" and rn in SB_NSP12) or (nsp == "NSP13" and rn in SB_NSP13):
        return COL_SB
    if nsp == "NSP12" and rn in PRIM_NSP12:
        return COL_HY
    if nsp == "NSP13":
        return COL_NSP13
    return COL_DEF

# Fig 4 — BSA
fig, ax = plt.subplots(figsize=(12, 5))
fig.patch.set_facecolor("white")
ax.set_facecolor("#f8f8f8")
labels = [f"{r['nsp'].replace('NSP','')[-2:]}:{r['res_name']}{r['res_num']}"
          for _, r in df.iterrows()]
bsa_vals = [r["bsa"] for _, r in df.iterrows()]
colors   = [get_color(r["nsp"], r["res_num"]) for _, r in df.iterrows()]
bars = ax.bar(labels, bsa_vals, color=colors, edgecolor="black", lw=0.7)
for bar, val in zip(bars, bsa_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            f"{val:.1f}", ha="center", va="bottom", fontsize=8)
ax.set_ylabel("BSA (Å²)", fontsize=11)
ax.set_title("NSP12–NSP13: Buried Surface Area per Hotspot Residue\n"
             "(averaged across 6XEZ, 7CXM, 7RDY)", fontsize=11, fontweight="bold")
ax.tick_params(axis="x", rotation=30)
legend_patches = [
    mpatches.Patch(color=COL_SB,  label="Salt bridge anchor (primary)"),
    mpatches.Patch(color=COL_HY,  label="Hydrophobic anchor (primary)"),
    mpatches.Patch(color=COL_NSP13, label="NSP13 residue"),
    mpatches.Patch(color=COL_DEF,  label="NSP12 residue"),
]
ax.legend(handles=legend_patches, fontsize=9, framealpha=0.9)
ax.grid(axis="y", alpha=0.4)
out4 = f"{RESULTS_DIR}/Fig4_NSP12-NSP13_BSA_8.png"
fig.savefig(out4, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  Saved: {out4}")

# Fig 5 — AlaScan
fig, ax = plt.subplots(figsize=(12, 5))
fig.patch.set_facecolor("white")
ax.set_facecolor("#f8f8f8")
ala_vals = [r["ala_loss"] for _, r in df.iterrows()]
bars = ax.bar(labels, ala_vals, color=colors, edgecolor="black", lw=0.7)
for bar, val in zip(bars, ala_vals):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.3,
                str(int(val)), ha="center", va="bottom", fontsize=8)
ax.set_ylabel("Contacts lost (Ala mutation)", fontsize=11)
ax.set_title("NSP12–NSP13: Alanine Scanning — Sidechain Contacts Lost\n"
             "(7RDY, 3.10 Å resolution)", fontsize=11, fontweight="bold")
ax.tick_params(axis="x", rotation=30)
ax.legend(handles=legend_patches, fontsize=9, framealpha=0.9)
ax.grid(axis="y", alpha=0.4)
out5 = f"{RESULTS_DIR}/Fig5_NSP12-NSP13_AlaScan_8.png"
fig.savefig(out5, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  Saved: {out5}")

# Fig 6 — Composite ranking
fig, ax = plt.subplots(figsize=(12, 5))
fig.patch.set_facecolor("white")
ax.set_facecolor("#f8f8f8")
comp_vals = [r["composite"] for _, r in df.iterrows()]
bars = ax.bar(labels, comp_vals, color=colors, edgecolor="black", lw=0.7)
for bar, val in zip(bars, comp_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
            f"{val:.4f}", ha="center", va="bottom", fontsize=7.5, rotation=45)
ax.set_ylabel("Composite score (BSA×0.4 + AlaScan×0.4 + Cons×0.2)", fontsize=10)
ax.set_title("NSP12–NSP13: Composite Hotspot Ranking\n"
             "Primary pharmacophores: MET902 (hydrophobic) + ASP901–LYS94 (SB)",
             fontsize=11, fontweight="bold")
ax.tick_params(axis="x", rotation=30)
ax.legend(handles=legend_patches, fontsize=9, framealpha=0.9)
ax.grid(axis="y", alpha=0.4)
out6 = f"{RESULTS_DIR}/Fig6_NSP12-NSP13_composite_ranking_8.png"
fig.savefig(out6, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  Saved: {out6}")

print("\n  Top 3 hotspots:")
for _, row in df.head(3).iterrows():
    print(f"    #{int(row['rank'])} {row['nsp']}:{row['res_name']}{row['res_num']}  "
          f"composite={row['composite']:.4f}  BSA={row['bsa']:.1f}  "
          f"AlaLoss={int(row['ala_loss'])}  cons={row['conservation']:.3f}")
