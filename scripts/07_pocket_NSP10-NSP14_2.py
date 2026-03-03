"""
Script 07_2: Pocket Detection — NSP10-NSP14
============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Runs fpocket on 7DIY.pdb and NSP10_NSP14_best_model.pdb.
  Parses all pocket scores and identifies which pocket overlaps
  with our conserved hotspot residues (especially HIS80, LYS93).
  Computes docking box center and dimensions from hotspot coordinates.
  Saves pocket analysis and docking box for Script 08_2.

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/07_pocket_NSP10-NSP14_2.py
"""

import json
import subprocess
import numpy as np
from pathlib import Path
from Bio import PDB

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP10-NSP14"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
RES_DIR.mkdir(parents=True, exist_ok=True)

# ── Structures to analyze ──────────────────────────────────
STRUCTURES = [
    {
        "label":   "PDB_7DIY",
        "file":    PDB_DIR / "7DIY.pdb",
        "chain_a": "A",
        "chain_b": "B",
        "source":  "SARS-CoV-2 crystal 2.69A",
    },
    {
        "label":   "PDB_5C8T",
        "file":    PDB_DIR / "5C8T.pdb",
        "chain_a": "A",
        "chain_b": "B",
        "source":  "SARS-CoV-1 crystal 3.20A",
    },
    {
        "label":   "AF3_model",
        "file":    AF3_DIR / "NSP10_NSP14_best_model.pdb",
        "chain_a": "A",
        "chain_b": "B",
        "source":  "AlphaFold3 iptm=0.89",
    },
]

# ── Conserved hotspot residues from Scripts 05+06 ─────────
HOTSPOTS_NSP10 = {5, 19, 21, 40, 42, 44, 45, 80, 93}   # conserved ≥0.8
HOTSPOTS_NSP14 = {4, 7, 8, 9, 10, 20, 25, 27, 127, 201}

# Primary target residues
PRIMARY_NSP10 = {80, 93}   # HIS80 salt bridge + LYS93
PRIMARY_NSP14 = {126, 127} # ASP126 salt bridge + THR127


def run_fpocket(pdb_file):
    """Run fpocket on a PDB file"""
    out_dir = pdb_file.parent / f"{pdb_file.stem}_out"
    if out_dir.exists():
        print(f"  SKIP  fpocket (output exists): {out_dir.name}")
        return out_dir
    cmd = ["fpocket", "-f", str(pdb_file)]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd=str(pdb_file.parent))
    if result.returncode == 0:
        print(f"  OK    fpocket → {out_dir.name}")
        return out_dir
    else:
        print(f"  FAIL  fpocket: {result.stderr[:200]}")
        return None


def parse_fpocket_info(info_file):
    """Parse fpocket _info.txt into list of pocket dicts"""
    pockets = []
    current = {}
    with open(info_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                current = {"pocket_id": int(line.split()[1])}
            elif ":" in line:
                key, val = line.split(":", 1)
                key = key.strip().lower().replace(" ", "_")
                try:
                    current[key] = float(val.strip())
                except ValueError:
                    current[key] = val.strip()
    if current:
        pockets.append(current)
    return pockets


def get_pocket_residues(pocket_atm_file):
    """Get set of (chain, resnum) from pocket atom file"""
    residues = set()
    if not pocket_atm_file.exists():
        return residues
    with open(pocket_atm_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain  = line[21]
                resnum = int(line[22:26].strip())
                residues.add((chain, resnum))
    return residues


def score_pocket_overlap(pocket_residues, chain_a, chain_b):
    """Score how much a pocket overlaps with our hotspots"""
    res_a = {r for c, r in pocket_residues if c == chain_a}
    res_b = {r for c, r in pocket_residues if c == chain_b}

    primary_a   = len(res_a & PRIMARY_NSP10)
    primary_b   = len(res_b & PRIMARY_NSP14)
    hotspot_a   = len(res_a & HOTSPOTS_NSP10)
    hotspot_b   = len(res_b & HOTSPOTS_NSP14)
    both_chains = len({c for c, r in pocket_residues}) == 2

    overlap_score = (primary_a * 3 + primary_b * 3 +
                     hotspot_a + hotspot_b +
                     (5 if both_chains else 0))
    return {
        "overlap_score":  overlap_score,
        "primary_NSP10":  primary_a,
        "primary_NSP14":  primary_b,
        "hotspot_NSP10":  hotspot_a,
        "hotspot_NSP14":  hotspot_b,
        "spans_interface": both_chains,
        "NSP10_residues": sorted(res_a),
        "NSP14_residues": sorted(res_b),
    }


def get_hotspot_coordinates(pdb_file, chain_a, chain_b):
    """Get 3D coordinates of all conserved hotspot residues"""
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_file)
    coords    = []

    all_hotspots = (
        [(chain_a, r) for r in HOTSPOTS_NSP10] +
        [(chain_b, r) for r in HOTSPOTS_NSP14]
    )

    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                key = (chain.id, res.id[1])
                if key in all_hotspots:
                    for atom in res.get_atoms():
                        if atom.name == "CA":
                            coords.append(atom.coord)

    return np.array(coords) if coords else None


def compute_docking_box(coords, padding=6.0):
    """Compute docking box center and dimensions from coordinates"""
    center = coords.mean(axis=0)
    mins   = coords.min(axis=0)
    maxs   = coords.max(axis=0)
    size   = (maxs - mins) + padding * 2
    return {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x":   round(float(size[0]), 3),
        "size_y":   round(float(size[1]), 3),
        "size_z":   round(float(size[2]), 3),
    }


def analyze_structure(s):
    print(f"\n{'='*55}")
    print(f"  {s['label']}  ({s['file'].name})")
    print(f"{'='*55}")

    if not s["file"].exists():
        print(f"  SKIP — file not found")
        return None

    # ── Run fpocket ───────────────────────────────────────
    out_dir = run_fpocket(s["file"])
    if out_dir is None:
        return None

    # ── Parse pocket info ─────────────────────────────────
    stem      = s["file"].stem
    info_file = out_dir / f"{stem}_info.txt"
    if not info_file.exists():
        # Try alternate naming
        info_files = list(out_dir.glob("*_info.txt"))
        if info_files:
            info_file = info_files[0]
        else:
            print(f"  FAIL — info file not found in {out_dir}")
            return None

    pockets = parse_fpocket_info(info_file)
    print(f"  Found {len(pockets)} pockets")

    # ── Score pocket overlap with hotspots ────────────────
    pockets_dir = out_dir / "pockets"
    best_pocket = None
    best_score  = -1

    for p in pockets:
        pid  = p["pocket_id"]
        atm  = pockets_dir / f"pocket{pid}_atm.pdb"
        res  = get_pocket_residues(atm)
        ovlp = score_pocket_overlap(res, s["chain_a"], s["chain_b"])
        p.update(ovlp)

        if ovlp["overlap_score"] > best_score:
            best_score  = ovlp["overlap_score"]
            best_pocket = p

    # ── Print top 5 pockets ───────────────────────────────
    top5 = sorted(pockets,
                  key=lambda x: x.get("overlap_score", 0),
                  reverse=True)[:5]
    print(f"\n  Top 5 pockets by hotspot overlap:")
    print(f"  {'ID':<5} {'Drug.':>7} {'Score':>7} "
          f"{'Overlap':>8} {'Interface':>10}")
    print(f"  {'-'*5} {'-'*7} {'-'*7} {'-'*8} {'-'*10}")
    for p in top5:
        flag = " ⭐" if p.get("spans_interface") else ""
        print(f"  {int(p['pocket_id']):<5} "
              f"{p.get('druggability_score', 0):>7.3f} "
              f"{p.get('score', 0):>7.3f} "
              f"{p.get('overlap_score', 0):>8.1f} "
              f"{'YES' if p.get('spans_interface') else 'NO':>10}"
              f"{flag}")

    # ── Best pocket details ───────────────────────────────
    print(f"\n  Best pocket: Pocket {int(best_pocket['pocket_id'])}")
    print(f"  Druggability score : {best_pocket.get('druggability_score',0):.3f}")
    print(f"  fpocket score      : {best_pocket.get('score',0):.3f}")
    print(f"  Spans interface    : {best_pocket.get('spans_interface')}")
    print(f"  Primary NSP10 hits : {best_pocket.get('primary_NSP10')} "
          f"(HIS80, LYS93)")
    print(f"  Primary NSP14 hits : {best_pocket.get('primary_NSP14')} "
          f"(ASP126, THR127)")
    print(f"  NSP10 residues     : {best_pocket.get('NSP10_residues')}")
    print(f"  NSP14 residues     : {best_pocket.get('NSP14_residues')}")

    # ── Docking box ───────────────────────────────────────
    coords = get_hotspot_coordinates(
        s["file"], s["chain_a"], s["chain_b"])
    if coords is not None:
        box = compute_docking_box(coords, padding=6.0)
        print(f"\n  Docking box (centered on conserved hotspots):")
        print(f"  Center : ({box['center_x']}, "
              f"{box['center_y']}, {box['center_z']})")
        print(f"  Size   : ({box['size_x']} x "
              f"{box['size_y']} x {box['size_z']}) Å")
    else:
        box = None

    return {
        "label":        s["label"],
        "n_pockets":    len(pockets),
        "best_pocket":  best_pocket,
        "docking_box":  box,
        "all_pockets":  pockets,
    }


def main():
    print("\n" + "="*55)
    print("  Script 07_2: Pocket Detection — NSP10-NSP14")
    print("  Structures: 7DIY (SARS-CoV-2) + AF3 model")
    print("="*55)

    results = {}
    for s in STRUCTURES:
        r = analyze_structure(s)
        if r:
            results[s["label"]] = r

    # ── Save results ──────────────────────────────────────
    # Convert numpy types for JSON
    def clean(obj):
        if isinstance(obj, (np.floating, float)):
            return round(float(obj), 4)
        if isinstance(obj, (np.integer, int)):
            return int(obj)
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [clean(v) for v in obj]
        return obj

    out = RES_DIR / "pocket_analysis_2.json"
    with open(out, "w") as f:
        json.dump(clean(results), f, indent=2)

    print(f"\n  Saved: 02-validation/NSP10-NSP14/pocket_analysis_2.json")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
