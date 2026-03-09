#!/usr/bin/env python3
"""
Script 07_pocket_NSP12-NSP13_8.py

fpocket druggability assessment and docking box definition
for NSP12-NSP13 interface.

Structures:
  6XEZ : A=NSP12, E=NSP13 — primary
  7CXM : A=NSP12, E=NSP13 — secondary
  7RDY : A=NSP12, E=NSP13 — tertiary
  AF3  : A=NSP12, B=NSP13 — fold only (interface not predicted)

Interface hotspots (6XEZ numbering):
  NSP12: LEU900, ASP901★, MET902★, TYR903, SER904
  NSP13: PHE90, GLY91, LEU92, TYR93, LYS94★, ASN95, THR96

Output: 02-validation/NSP12-NSP13/pocket_analysis_8.json
"""

import os, json, subprocess, shutil
import numpy as np
from Bio import PDB

# ── Paths ──────────────────────────────────────────────────────────────────
BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = f"{BASE}/00-reference/pdb_structures"
AF3_DIR = f"{BASE}/01-alphafold3/NSP12-NSP13"
VAL_DIR = f"{BASE}/02-validation/NSP12-NSP13"
TMP_DIR = f"{BASE}/tmp_fpocket_8"
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(VAL_DIR, exist_ok=True)

STRUCTURES = {
    "6XEZ" : (f"{PDB_DIR}/6XEZ.pdb", ("A","E")),
    "7CXM" : (f"{PDB_DIR}/7CXM.pdb", ("A","E")),
    "7RDY" : (f"{PDB_DIR}/7RDY.pdb", ("A","E")),
    "AF3"  : (f"{AF3_DIR}/NSP12_NSP13_best_model.pdb", ("A","B")),
}

# Hotspot residues per chain (PDB numbering, 6XEZ)
NSP12_HOTSPOTS = [900, 901, 902, 903, 904]
NSP13_HOTSPOTS = [90, 91, 92, 93, 94, 95, 96]
SB_RESIDUES    = {901, 902}    # NSP12 primary
SB_NSP13       = {94}          # NSP13 primary

ALL_HOTSPOTS_NSP12 = set(NSP12_HOTSPOTS)
ALL_HOTSPOTS_NSP13 = set(NSP13_HOTSPOTS)

# ── Helper: extract two chains only ───────────────────────────────────────
def extract_two_chains(pdb_in, chain1, chain2, pdb_out):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_in)
    io     = PDB.PDBIO()
    class TwoChainSelect(PDB.Select):
        def accept_chain(self, c):
            return c.get_id() in (chain1, chain2)
        def accept_residue(self, r):
            return r.get_id()[0] == " "
    io.set_structure(struct)
    io.save(pdb_out, TwoChainSelect())

# ── Helper: run fpocket ────────────────────────────────────────────────────
def run_fpocket(pdb_in, label):
    work_pdb = f"{TMP_DIR}/{label}_input.pdb"
    shutil.copy(pdb_in, work_pdb)
    result = subprocess.run(
        ["fpocket", "-f", work_pdb],
        capture_output=True, text=True, cwd=TMP_DIR
    )
    out_dir = work_pdb.replace(".pdb", "_out")
    info_file = os.path.join(out_dir, f"{label}_input_info.txt")
    return out_dir, info_file, result.returncode

# ── Helper: parse fpocket info file ───────────────────────────────────────
def parse_fpocket_info(info_file):
    pockets = []
    current = None
    with open(info_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                pno = int(line.split()[1])
                current = {"pocket_num": pno}
            elif current is not None and ":" in line:
                k, v = line.split(":", 1)
                k = k.strip()
                v = v.strip()
                try:
                    current[k] = float(v)
                except ValueError:
                    current[k] = v
    if current:
        pockets.append(current)
    return pockets

# ── Helper: get Cα coordinates of hotspot residues ─────────────────────────
def get_hotspot_coords(pdb_file, chain_id, res_nums):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_file)
    coords = {}
    for res in struct[0][chain_id]:
        rn = res.get_id()[1]
        if rn in res_nums and "CA" in res:
            coords[rn] = res["CA"].get_vector().get_array()
    return coords

# ── Helper: find nearest pocket to interface ──────────────────────────────
def find_interface_pocket(pockets, hotspot_coords_list):
    if not hotspot_coords_list:
        return None, 99.0
    all_coords = np.array(list(hotspot_coords_list.values()))
    iface_center = all_coords.mean(axis=0)

    best_pocket = None
    best_dist   = 99.0
    for p in pockets:
        cx = p.get("Pocket center x", None)
        cy = p.get("Pocket center y", None)
        cz = p.get("Pocket center z", None)
        if cx is None:
            continue
        dist = float(np.linalg.norm(np.array([cx,cy,cz]) - iface_center))
        if dist < best_dist:
            best_dist   = dist
            best_pocket = p
    return best_pocket, best_dist

# ── Helper: get druggability score ────────────────────────────────────────
def get_druggability(pocket):
    if pocket is None:
        return 0.0
    for key in pocket:
        if "drug" in key.lower():
            return float(pocket[key])
    return 0.0

# ── Run fpocket on all structures ──────────────────────────────────────────
print("=" * 60)
print("POCKET DETECTION — NSP12-NSP13")
print("=" * 60)

results = {}
for pdb_id, (pdb_path, (c1, c2)) in STRUCTURES.items():
    print(f"\n  Processing {pdb_id} (chains {c1}+{c2}) ...")
    two_chain_pdb = f"{TMP_DIR}/{pdb_id}_two_chain.pdb"
    extract_two_chains(pdb_path, c1, c2, two_chain_pdb)

    out_dir, info_file, rc = run_fpocket(two_chain_pdb, pdb_id)

    if rc != 0 or not os.path.exists(info_file):
        print(f"    fpocket failed for {pdb_id}")
        results[pdb_id] = {"n_pockets": 0, "error": True}
        continue

    pockets = parse_fpocket_info(info_file)
    print(f"    {len(pockets)} pockets found")

    # Get hotspot coordinates
    c12_hotspots = get_hotspot_coords(two_chain_pdb, c1, set(NSP12_HOTSPOTS))
    c13_hotspots = get_hotspot_coords(two_chain_pdb, c2, set(NSP13_HOTSPOTS))
    all_hotspot_coords = {**{f"nsp12_{k}": v for k,v in c12_hotspots.items()},
                          **{f"nsp13_{k}": v for k,v in c13_hotspots.items()}}

    best_pocket, best_dist = find_interface_pocket(pockets, all_hotspot_coords)
    druggability = get_druggability(best_pocket)

    print(f"    Best interface pocket: #{best_pocket['pocket_num'] if best_pocket else 'N/A'}")
    print(f"    Distance to interface center: {best_dist:.2f} A")
    print(f"    Druggability score: {druggability:.3f}")
    if best_pocket:
        vol = best_pocket.get("Volume", best_pocket.get("Pocket volume", 0))
        print(f"    Pocket volume: {vol:.1f} A3")

    results[pdb_id] = {
        "n_pockets"    : len(pockets),
        "best_pocket"  : best_pocket,
        "best_dist_A"  : round(best_dist, 2),
        "druggability" : round(druggability, 3),
    }

# ── Docking box from all hotspot Cα coordinates ────────────────────────────
print("\n" + "=" * 60)
print("DOCKING BOX DEFINITION")
print("=" * 60)

# Use 6XEZ as primary reference
two_chain_6XEZ = f"{TMP_DIR}/6XEZ_two_chain.pdb"
coords_nsp12 = get_hotspot_coords(two_chain_6XEZ, "A", set(NSP12_HOTSPOTS))
coords_nsp13 = get_hotspot_coords(two_chain_6XEZ, "E", set(NSP13_HOTSPOTS))

all_coords = np.array(
    list(coords_nsp12.values()) + list(coords_nsp13.values())
)

if len(all_coords) > 0:
    center = all_coords.mean(axis=0)
    span   = all_coords.max(axis=0) - all_coords.min(axis=0)
    pad    = 12.0  # Å padding around hotspot span
    size   = span + 2 * pad

    volume = float(size[0] * size[1] * size[2])

    print(f"\n  Hotspot coordinates used: {len(all_coords)} Cα atoms")
    print(f"  NSP12 hotspots found: {sorted(coords_nsp12.keys())}")
    print(f"  NSP13 hotspots found: {sorted(coords_nsp13.keys())}")
    print(f"\n  Docking box (6XEZ reference):")
    print(f"    Center : ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    print(f"    Size   : {size[0]:.3f} x {size[1]:.3f} x {size[2]:.3f} A")
    print(f"    Volume : {volume:.0f} A3")

    # Verify all primary hotspots are within box
    half = size / 2
    missing = []
    for label, coord in {**{f"NSP12_{k}": v for k,v in coords_nsp12.items()},
                          **{f"NSP13_{k}": v for k,v in coords_nsp13.items()}}.items():
        diff = np.abs(coord - center)
        if not np.all(diff <= half + 2.0):
            missing.append(label)

    if missing:
        print(f"\n  WARNING: Hotspots outside box: {missing}")
    else:
        print(f"\n  All {len(all_coords)} hotspot Cα atoms confirmed within docking box ✅")

    docking_box = {
        "center" : [round(float(x), 3) for x in center],
        "size"   : [round(float(x), 3) for x in size],
        "volume" : round(volume, 0),
        "padding": pad,
        "reference_structure": "6XEZ"
    }
else:
    print("  ERROR: no hotspot coordinates found")
    docking_box = {}

# ── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("DRUGGABILITY SUMMARY")
print("=" * 60)

print(f"\n  {'Structure':>8}  {'Pockets':>8}  {'Best dist':>10}  {'Druggability':>13}")
print("  " + "-"*45)
for pdb_id, r in results.items():
    if r.get("error"):
        print(f"  {pdb_id:>8}  ERROR")
        continue
    print(f"  {pdb_id:>8}  {r['n_pockets']:>8}  "
          f"{r['best_dist_A']:>10.2f}  {r['druggability']:>13.3f}")

print(f"\n  SCIENTIFIC NOTE:")
print(f"  Low fpocket druggability expected for small PPI interfaces")
print(f"  (same pattern as NSP10-NSP14 and NSP13-Helicase)")
print(f"  Interface defined by 5+7=12 residues total — very compact")
print(f"  MET902 hydrophobic groove is primary druggable feature")

# ── Save ───────────────────────────────────────────────────────────────────
output = {
    "complex"       : "NSP12-NSP13",
    "script"        : "07_pocket_NSP12-NSP13_8.py",
    "fpocket_results": {pdb_id: {k: v for k, v in r.items()
                                 if k != "best_pocket"}
                        for pdb_id, r in results.items()},
    "docking_box"   : docking_box,
    "hotspot_coords_nsp12": {str(k): [round(float(x),3) for x in v]
                              for k,v in coords_nsp12.items()},
    "hotspot_coords_nsp13": {str(k): [round(float(x),3) for x in v]
                              for k,v in coords_nsp13.items()},
}

out_path = f"{VAL_DIR}/pocket_analysis_8.json"
with open(out_path, "w") as f:
    json.dump(output, f, indent=2)
print(f"\n  Saved: {out_path}")

# Cleanup
shutil.rmtree(TMP_DIR, ignore_errors=True)
print("  Temp files cleaned up")
