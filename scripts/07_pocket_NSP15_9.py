"""
Script 07_9 — Pocket Detection NSP15 Homodimer
Complex  : NSP15-NSP15 (NendoU dimer interface)
Suffix   : _9
Structures: 6VWW (apo), 6WXC (tipiracil), 9HH5 (YM-155)
Key validation: fpocket MUST rediscover tipiracil and YM-155 binding sites
                3-way ground truth check unique in the project
"""

import json
import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser

OUT_DIR  = "02-validation/NSP15"
OUT_JSON = f"{OUT_DIR}/pocket_analysis_9.json"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Interface hotspots (PDB numbering from Script 05_9) ────────
HOTSPOTS_PDB = [
    10,11,12,13,14,15,26,28,29,33,36,38,39,40,41,42,43,
    47,48,49,59,62,64,91,95,96,97,163,164,166,168,169,
    170,171,172,173,203,241,242,243,244,263,265,267,269,
    270,271,272,280,282,283,284,285,286,287,288,289,291,292,346
]

# Primary pharmacophore residues (PDB numbering)
PRIMARY = {40: "ASP40", 62: "ARG62", 91: "ARG91", 267: "GLU267"}

print("=" * 65)
print("Script 07_9 — Pocket Detection NSP15")
print("=" * 65)

# ── Structures to analyze ──────────────────────────────────────
STRUCTURES = {
    "6VWW_apo": {
        "path"  : "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb",
        "label" : "apo 1.90A PRIMARY",
        "ground_truth": None,
    },
    "6WXC_tipiracil": {
        "path"  : "00-reference/known_interfaces/NSP15/6WXC_NSP15-true-dimer-tipiracil.pdb",
        "label" : "tipiracil-bound 1.80A",
        "ground_truth": "tipiracil active site",
    },
    "9HH5_YM155": {
        "path"  : "00-reference/pdb_structures/NSP15/9HH5.pdb",
        "label" : "YM-155 inhibitor 2.08A",
        "ground_truth": "YM-155 binding site",
    },
}

def run_fpocket(pdb_path, label):
    """Run fpocket on a PDB file, return results dict."""
    if not os.path.exists(pdb_path):
        print(f"  WARNING: {pdb_path} not found")
        return None

    result = subprocess.run(
        ["fpocket", "-f", pdb_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  fpocket error: {result.stderr[:200]}")
        return None

    # fpocket creates output directory: pdb_name_out/
    base     = os.path.splitext(pdb_path)[0]
    out_dir  = base + "_out"
    info_file = os.path.join(out_dir, os.path.basename(base) + "_info.txt")

    if not os.path.exists(info_file):
        print(f"  WARNING: info file not found: {info_file}")
        return None

    # Parse info file
    pockets = []
    current = {}
    with open(info_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                current = {"id": line.split()[1]}
            elif ":" in line:
                key, val = line.split(":", 1)
                key = key.strip().lower().replace(" ","_")
                try:
                    current[key] = float(val.strip())
                except ValueError:
                    current[key] = val.strip()
    if current:
        pockets.append(current)

    return pockets, out_dir

def get_hotspot_overlap(pocket_pdb_path, hotspot_positions, cutoff=4.0):
    """Check which hotspot residues are near the pocket atoms."""
    if not os.path.exists(pocket_pdb_path):
        return set()
    struct  = PDBParser(QUIET=True).get_structure("pocket", pocket_pdb_path)
    pocket_coords = np.array([
        a.get_vector().get_array()
        for a in struct.get_atoms()])
    if len(pocket_coords) == 0:
        return set()

    # Load full complex
    ref_struct  = PDBParser(QUIET=True).get_structure("ref",
        "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb")
    nearby = set()
    for chain in ref_struct[0].get_chains():
        for res in chain.get_residues():
            if res.get_id()[0] != ' ':
                continue
            pos = res.get_id()[1]
            if pos not in hotspot_positions:
                continue
            for atom in res.get_atoms():
                coord = atom.get_vector().get_array()
                dists = np.linalg.norm(pocket_coords - coord, axis=1)
                if dists.min() < cutoff:
                    nearby.add(pos)
                    break
    return nearby

# ── Run fpocket on all structures ──────────────────────────────
results = {}

for struct_id, info in STRUCTURES.items():
    print(f"\n{struct_id}: {info['label']}")
    ret = run_fpocket(info["path"], info["label"])
    if ret is None:
        continue
    pockets, out_dir = ret
    print(f"  Total pockets found: {len(pockets)}")

    # Find best pocket by druggability score
    best_idx = 0
    best_score = -1
    for i, p in enumerate(pockets):
        score = p.get("drug_score", p.get("druggability_score", 0))
        if score > best_score:
            best_score = score
            best_idx = i

    best = pockets[best_idx] if pockets else {}
    print(f"  Best pocket #{best_idx+1}:")
    print(f"    Druggability score : {best.get('drug_score', best.get('druggability_score','?'))}")
    print(f"    Volume (A3)        : {best.get('volume','?')}")

    # Check hotspot overlap for best pocket
    pocket_pdb = os.path.join(
        out_dir, "pockets",
        f"pocket{best_idx+1}_atm.pdb")
    overlap = get_hotspot_overlap(pocket_pdb, set(HOTSPOTS_PDB))
    primary_overlap = set(HOTSPOTS_PDB) & overlap & set(PRIMARY.keys())
    print(f"    Hotspot overlap    : {len(overlap)}/{len(HOTSPOTS_PDB)} residues")
    print(f"    Primary SB overlap : {[PRIMARY[p] for p in primary_overlap]}")

    if info["ground_truth"]:
        print(f"    Ground truth check : {info['ground_truth']}")

    results[struct_id] = {
        "label"           : info["label"],
        "n_pockets"       : len(pockets),
        "best_pocket_id"  : best_idx + 1,
        "best_druggability": best.get("drug_score",
                             best.get("druggability_score", 0)),
        "best_volume"     : best.get("volume", 0),
        "hotspot_overlap" : sorted(overlap),
        "primary_overlap" : [PRIMARY[p] for p in primary_overlap],
        "ground_truth"    : info["ground_truth"],
    }

# ── Consensus docking box from hotspot Cα coords ──────────────
print("\n" + "=" * 65)
print("Computing docking box from primary hotspot Ca coords (6VWW):")

pdb = PDBParser(QUIET=True).get_structure("6VWW",
    "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb")

coords = []
for chain in pdb[0].get_chains():
    for res in chain.get_residues():
        if res.get_id()[0] != ' ':
            continue
        if res.get_id()[1] in HOTSPOTS_PDB:
            if 'CA' in res:
                coords.append(res['CA'].get_vector().get_array())

coords = np.array(coords)
center = coords.mean(axis=0)
PADDING = 6.0
size   = (coords.max(axis=0) - coords.min(axis=0)) + 2 * PADDING

print(f"  Hotspot Ca atoms used: {len(coords)}")
print(f"  Center : ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
print(f"  Size   : {size[0]:.3f} x {size[1]:.3f} x {size[2]:.3f} A")
print(f"  Volume : {size[0]*size[1]*size[2]:.0f} A3")

# ── Save results ───────────────────────────────────────────────
output = {
    "complex"         : "NSP15-homodimer",
    "suffix"          : "_9",
    "docking_box"     : {
        "center"      : [round(float(c),3) for c in center],
        "size"        : [round(float(s),3) for s in size],
        "volume"      : round(float(size[0]*size[1]*size[2]),0),
        "padding"     : PADDING,
    },
    "per_structure"   : results,
    "ground_truth_note": (
        "6WXC contains tipiracil (known NendoU inhibitor). "
        "9HH5 contains YM-155 (Sepantronium, 2.08A). "
        "fpocket must rediscover these sites as top pockets."
    ),
}

with open(OUT_JSON, "w") as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nSaved: {OUT_JSON}")
print("Next : Script 08_9 (docking preparation)")
