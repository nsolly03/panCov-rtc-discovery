"""
Script 07_4: Pocket Detection — NSP12-NSP8
==========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2, 6NUR, 7C2K + AF3
Chain A = NSP12, Chain B = NSP8

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/07_pocket_NSP12-NSP8_4.py
"""

import json
import subprocess
import shutil
import numpy as np
from pathlib import Path
from Bio import PDB

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP8"
TMP_DIR  = PROJECT / "tmp"
TMP_DIR.mkdir(exist_ok=True)

PDB_FILES = {
    "7BV2": PDB_DIR / "7BV2_NSP12-NSP8.pdb",
    "6NUR": PDB_DIR / "6NUR_NSP12-NSP8.pdb",
    "7C2K": PDB_DIR / "7C2K_NSP12-NSP8.pdb",
    "AF3":  AF3_DIR / "NSP12_NSP8_best_model.pdb",
}
PDB_FILE     = PDB_DIR / "7BV2_NSP12-NSP8.pdb"
CHAIN_NSP12  = "A"
CHAIN_NSP8   = "B"

# Hotspots from Script 05_4
HOTSPOTS_NSP12 = [387,129,389,271,330,131,380,523,
                   91,87,332,95,117,517,99]
HOTSPOTS_NSP8  = [117,129,80,115,131,112,91,87,
                   116,95,98,83,128,90,121]
PRIMARY_NSP12  = {523, 332, 517}
PRIMARY_NSP8   = {80, 99, 79}

# NSP12 offset: PDB 31 = local 31 (offset=0 for 7BV2)
NSP12_OFFSET   = 0
NSP8_OFFSET    = 77  # 7BV2 NSP8 starts at PDB 78


def pdb_to_local(gpos, chain):
    if chain == CHAIN_NSP12:
        return gpos - NSP12_OFFSET
    return gpos - NSP8_OFFSET


def run_fpocket(pdb_file):
    stem    = pdb_file.stem
    out_dir = TMP_DIR / f"{stem}_out"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    subprocess.run(
        ["fpocket", "-f", str(pdb_file)],
        capture_output=True, text=True,
        cwd=str(TMP_DIR))
    alt = pdb_file.parent / f"{stem}_out"
    if alt.exists():
        return alt
    return out_dir


def parse_fpocket_info(pdb_file):
    stem      = pdb_file.stem
    info_file = TMP_DIR / f"{stem}_out" / \
                f"{stem}_info.txt"
    if not info_file.exists():
        info_file = pdb_file.parent / \
                    f"{stem}_out" / f"{stem}_info.txt"
    if not info_file.exists():
        return []
    pockets, current = [], {}
    with open(info_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                current = {"id": int(line.split()[1])}
            elif ":" in line:
                k, _, v = line.partition(":")
                try:
                    current[k.strip()] = float(v.strip())
                except ValueError:
                    current[k.strip()] = v.strip()
    if current:
        pockets.append(current)
    return pockets


def find_water_bridges(structure):
    waters = [
        (c.id, r.id[1], list(r.get_atoms()))
        for m in structure for c in m
        for r in c
        if r.resname.strip() in ["HOH","WAT","SOL"]]

    res_12 = {r.id[1]: r
              for m in structure for c in m
              if c.id == CHAIN_NSP12
              for r in c if r.id[0] == " "}
    res_8  = {r.id[1]: r
              for m in structure for c in m
              if c.id == CHAIN_NSP8
              for r in c if r.id[0] == " "}

    bridges = []
    for wch, wrn, wat_atoms in waters:
        c12, c8 = [], []
        for wa in wat_atoms:
            for rn, res in res_12.items():
                for a in res.get_atoms():
                    try:
                        d = wa - a
                        if d <= 3.5:
                            c12.append((
                                pdb_to_local(
                                    rn, CHAIN_NSP12),
                                res.resname.strip(),
                                round(float(d),2)))
                    except Exception:
                        pass
            for rn, res in res_8.items():
                for a in res.get_atoms():
                    try:
                        d = wa - a
                        if d <= 3.5:
                            c8.append((
                                pdb_to_local(
                                    rn, CHAIN_NSP8),
                                res.resname.strip(),
                                round(float(d),2)))
                    except Exception:
                        pass
        if c12 and c8:
            bridges.append({
                "water":        f"HOH{wrn}",
                "nsp12_contacts":list(set(c12))[:3],
                "nsp8_contacts": list(set(c8))[:3],
            })
    return bridges


def define_docking_box(structure, padding=4.0):
    """Box centered on primary SB + top hotspots."""
    # Focus on primary SB pairs only for tight box
    target_12 = PRIMARY_NSP12  # {523,332,517}
    target_8  = PRIMARY_NSP8   # {80,99,79}
    coords = []
    for m in structure:
        for c in m:
            if c.id == CHAIN_NSP12:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    if pdb_to_local(
                            r.id[1], CHAIN_NSP12) \
                            in target_12:
                        coords.extend([
                            a.get_vector().get_array()
                            for a in r.get_atoms()])
            elif c.id == CHAIN_NSP8:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    if pdb_to_local(
                            r.id[1], CHAIN_NSP8) \
                            in target_8:
                        coords.extend([
                            a.get_vector().get_array()
                            for a in r.get_atoms()])
    if not coords:
        return None
    coords = np.array(coords)
    center = coords.mean(axis=0)
    size   = coords.max(axis=0) - \
             coords.min(axis=0) + 2*padding
    return {
        "center_x": round(float(center[0]),3),
        "center_y": round(float(center[1]),3),
        "center_z": round(float(center[2]),3),
        "size_x":   round(float(size[0]),3),
        "size_y":   round(float(size[1]),3),
        "size_z":   round(float(size[2]),3),
    }


def main():
    print("\n" + "="*57)
    print("  Script 07_4: Pocket Detection — NSP12-NSP8")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("7BV2", PDB_FILE)

    # fpocket all structures
    print("  1. fpocket analysis (all structures)...")
    all_pockets = {}
    for sid, pdb_f in PDB_FILES.items():
        if not pdb_f.exists():
            print(f"     {sid}: not found — skipping")
            continue
        run_fpocket(pdb_f)
        pockets = parse_fpocket_info(pdb_f)
        all_pockets[sid] = pockets
        if pockets:
            best  = max(pockets,
                        key=lambda p: p.get(
                            "Druggability Score",
                            p.get("Drug Score",0)))
            score = best.get("Druggability Score",
                             best.get("Drug Score",0))
            vol   = best.get("Volume",0)
            print(f"     {sid}: {len(pockets)} pockets | "
                  f"best druggability={score:.3f} "
                  f"vol={vol:.1f} Å³")
        else:
            print(f"     {sid}: no pockets parsed")

    pockets = all_pockets.get("7BV2",[])
    if pockets:
        print(f"\n     7BV2 top 5 pockets:")
        print(f"     {'Pocket':<8} {'DrugScore':>10} "
              f"{'Volume':>10}")
        print(f"     {'-'*30}")
        for p in sorted(
                pockets,
                key=lambda x: x.get(
                    "Druggability Score",
                    x.get("Drug Score",0)),
                reverse=True)[:5]:
            print(f"     {p.get('id','?'):<8} "
                  f"{p.get('Druggability Score', p.get('Drug Score',0)):>10.3f} "
                  f"{p.get('Volume',0):>10.1f}")

    # Water bridges
    print("\n  2. Water-mediated interface contacts...")
    bridges = find_water_bridges(structure)
    hotspot_b = [
        b for b in bridges
        if any(p in HOTSPOTS_NSP12
               for p,_,_ in b["nsp12_contacts"])
        or any(p in HOTSPOTS_NSP8
               for p,_,_ in b["nsp8_contacts"])]
    print(f"     Total water bridges:   {len(bridges)}")
    print(f"     At hotspot residues:   {len(hotspot_b)}")
    if hotspot_b:
        print(f"\n     Key water bridges:")
        for b in hotspot_b[:6]:
            n12 = [f"{rn}{p}" for p,rn,d
                   in b["nsp12_contacts"]]
            n8  = [f"{rn}{p}" for p,rn,d
                   in b["nsp8_contacts"]]
            print(f"       {b['water']}: "
                  f"NSP12={n12} ↔ NSP8={n8}")

    # Docking box
    print("\n  3. Docking box definition...")
    box = define_docking_box(structure)
    if box:
        vol = (box["size_x"]*box["size_y"]*
               box["size_z"])
        print(f"     Center: ({box['center_x']}, "
              f"{box['center_y']}, {box['center_z']})")
        print(f"     Size:   {box['size_x']} × "
              f"{box['size_y']} × {box['size_z']} Å")
        print(f"     Volume: {vol:.0f} Å³")

    # Save
    out = RES_DIR / "pocket_analysis_4.json"
    with open(out,"w") as f:
        json.dump({
            "complex":        "NSP12-NSP8",
            "pdb":            "7BV2",
            "fpocket_all":    {k: v[:5]
                               for k,v in
                               all_pockets.items()},
            "fpocket_7BV2":   pockets[:10],
            "water_bridges":  bridges,
            "hotspot_water":  hotspot_b,
            "docking_box":    box,
        }, f, indent=2)
    print(f"\n  Saved: 02-validation/NSP12-NSP8/"
          f"pocket_analysis_4.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
