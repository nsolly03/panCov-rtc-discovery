"""
Script 07_3: Pocket Detection — NSP12-NSP7
==========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2 (primary, 2.90 Å)
Chain A = NSP12, Chain C = NSP7

Analyses:
  1. fpocket — druggable pocket detection
  2. Water-mediated contacts at interface
  3. Disulfide bridge check
  4. Docking box definition

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/07_pocket_NSP12-NSP7_3.py
"""

import json
import subprocess
import shutil
import numpy as np
from pathlib import Path
from Bio import PDB

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP7"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP7"
TMP_DIR  = PROJECT / "tmp"
TMP_DIR.mkdir(exist_ok=True)

# All structures for pocket analysis
PDB_FILES = {
    "7BV2": PDB_DIR / "7BV2_NSP12-NSP7.pdb",
    "6NUR": PDB_DIR / "6NUR_NSP12-NSP7.pdb",
    "7C2K": PDB_DIR / "7C2K_NSP12-NSP7.pdb",
    "AF3":  AF3_DIR / "NSP12_NSP7_best_model.pdb",
}
PDB_FILE        = PDB_DIR / "7BV2_NSP12-NSP7.pdb"  # primary
CHAIN_NSP12     = "A"
CHAIN_NSP7      = "C"
NSP12_OFFSET    = 30
NSP7_OFFSET     = 1

HOTSPOTS_NSP12  = [440,412,442,443,420,843,409,
                    40,33,41,37,413,415,14,23]
HOTSPOTS_NSP7   = [40,14,33,41,37,11,23,5,15,29,12,4,1]
PRIMARY_NSP12   = {431}
PRIMARY_NSP7    = {2}
HYDROPHOBIC_CORE= {412,413,415,420,440,442,843}


def genome_to_local(gpos, chain):
    if chain == CHAIN_NSP12:
        return gpos - NSP12_OFFSET
    return gpos - NSP7_OFFSET


# ── 1. fpocket ─────────────────────────────────────────────

def run_fpocket(pdb_file):
    stem    = pdb_file.stem
    out_dir = TMP_DIR / f"{stem}_out"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    subprocess.run(
        ["fpocket", "-f", str(pdb_file)],
        capture_output=True, text=True,
        cwd=str(TMP_DIR))
    # fpocket outputs next to input file
    alt = pdb_file.parent / f"{stem}_out"
    if alt.exists():
        return alt
    return out_dir


def parse_fpocket_info(pdb_file):
    stem      = pdb_file.stem
    info_file = TMP_DIR / f"{stem}_out" / f"{stem}_info.txt"
    if not info_file.exists():
        # Try next to pdb file
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


# ── 2. Water bridges ───────────────────────────────────────

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
    res_7  = {r.id[1]: r
              for m in structure for c in m
              if c.id == CHAIN_NSP7
              for r in c if r.id[0] == " "}

    bridges = []
    for wch, wrn, wat_atoms in waters:
        c12, c7 = [], []
        for wa in wat_atoms:
            for rn, res in res_12.items():
                for a in res.get_atoms():
                    try:
                        d = wa - a
                        if d <= 3.5:
                            local = genome_to_local(
                                rn, CHAIN_NSP12)
                            c12.append(
                                (local,
                                 res.resname.strip(),
                                 round(float(d),2)))
                    except Exception:
                        pass
            for rn, res in res_7.items():
                for a in res.get_atoms():
                    try:
                        d = wa - a
                        if d <= 3.5:
                            local = genome_to_local(
                                rn, CHAIN_NSP7)
                            c7.append(
                                (local,
                                 res.resname.strip(),
                                 round(float(d),2)))
                    except Exception:
                        pass
        if c12 and c7:
            bridges.append({
                "water":         f"HOH{wrn}",
                "nsp12_contacts":list(set(c12))[:3],
                "nsp7_contacts": list(set(c7))[:3],
            })
    return bridges


# ── 3. Disulfide check ─────────────────────────────────────

def find_disulfides(structure):
    cys_all = [
        (c.id, r.id[1], r)
        for m in structure for c in m
        if c.id in [CHAIN_NSP12, CHAIN_NSP7]
        for r in c
        if r.resname.strip() == "CYS"
        and r.id[0] == " "]

    disulfides = []
    for i, (ch1, rn1, r1) in enumerate(cys_all):
        for ch2, rn2, r2 in cys_all[i+1:]:
            try:
                d = r1["SG"] - r2["SG"]
                if d <= 2.1:
                    l1 = genome_to_local(rn1, ch1)
                    l2 = genome_to_local(rn2, ch2)
                    disulfides.append({
                        "res1": f"CYS{l1}({ch1})",
                        "res2": f"CYS{l2}({ch2})",
                        "distance": round(float(d),2),
                    })
            except Exception:
                pass
    return disulfides


# ── 4. Docking box ─────────────────────────────────────────

def define_docking_box(structure, padding=5.0):
    """Box centered on primary + hydrophobic core hotspots."""
    target_12 = PRIMARY_NSP12 | HYDROPHOBIC_CORE
    target_7  = PRIMARY_NSP7

    coords = []
    for m in structure:
        for c in m:
            if c.id == CHAIN_NSP12:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    local = genome_to_local(
                        r.id[1], CHAIN_NSP12)
                    if local in target_12:
                        coords.extend([
                            a.get_vector().get_array()
                            for a in r.get_atoms()])
            elif c.id == CHAIN_NSP7:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    local = genome_to_local(
                        r.id[1], CHAIN_NSP7)
                    if local in target_7:
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


# ── Main ───────────────────────────────────────────────────

def main():
    print("\n" + "="*57)
    print("  Script 07_3: Pocket Detection — NSP12-NSP7")
    print(f"  PDB: 7BV2_NSP12-NSP7")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("7BV2", PDB_FILE)

    # ── fpocket ───────────────────────────────────────────
    print("  1. fpocket analysis (all structures)...")
    all_pockets = {}
    for sid, pdb_f in PDB_FILES.items():
        if not pdb_f.exists():
            print(f"     {sid}: file not found — skipping")
            continue
        run_fpocket(pdb_f)
        pockets = parse_fpocket_info(pdb_f)
        all_pockets[sid] = pockets
        if pockets:
            best = max(pockets,
                       key=lambda p: p.get(
                           "Druggability Score",
                           p.get("Drug Score", 0)))
            score = best.get("Druggability Score",
                             best.get("Drug Score", 0))
            vol   = best.get("Volume", 0)
            print(f"     {sid}: {len(pockets)} pockets | "
                  f"best druggability={score:.3f} "
                  f"vol={vol:.1f} Å³")
        else:
            print(f"     {sid}: no pockets parsed")

    # Use 7BV2 pockets as primary
    pockets = all_pockets.get("7BV2", [])
    if pockets:
        print(f"\n     7BV2 top pockets:")
        print(f"     {'Pocket':<8} {'DrugScore':>10} "
              f"{'Volume':>10}")
        print(f"     {'-'*8} {'-'*10} {'-'*10}")
        for p in sorted(pockets,
                        key=lambda x: x.get(
                            "Druggability Score",
                            x.get("Drug Score",0)),
                        reverse=True)[:5]:
            pid   = p.get("id","?")
            score = p.get("Druggability Score",
                    p.get("Drug Score", 0))
            vol   = p.get("Volume", 0)
            print(f"     {pid:<8} {score:>10.3f} "
                  f"{vol:>10.1f}")

    # ── Water bridges ─────────────────────────────────────
    print("\n  2. Water-mediated interface contacts...")
    bridges = find_water_bridges(structure)
    hotspot_b = [
        b for b in bridges
        if any(p in HOTSPOTS_NSP12
               for p,_,_ in b["nsp12_contacts"])
        or any(p in HOTSPOTS_NSP7
               for p,_,_ in b["nsp7_contacts"])]
    print(f"     Total water bridges:   {len(bridges)}")
    print(f"     At hotspot residues:   {len(hotspot_b)}")
    if hotspot_b:
        print(f"\n     Key water bridges:")
        for b in hotspot_b[:6]:
            n12 = [(f"{rn}{p}",d)
                   for p,rn,d in b["nsp12_contacts"]]
            n7  = [(f"{rn}{p}",d)
                   for p,rn,d in b["nsp7_contacts"]]
            print(f"       {b['water']}: "
                  f"NSP12={n12} ↔ NSP7={n7}")

    # ── Disulfides ────────────────────────────────────────
    print("\n  3. Disulfide bridge check...")
    disulfides = find_disulfides(structure)
    if disulfides:
        for ds in disulfides:
            print(f"     {ds['res1']} — {ds['res2']} "
                  f"({ds['distance']:.2f}Å)")
    else:
        print("     No disulfide bonds found")

    # ── Docking box ───────────────────────────────────────
    print("\n  4. Docking box definition...")
    box = define_docking_box(structure)
    if box:
        print(f"     Center: ({box['center_x']}, "
              f"{box['center_y']}, {box['center_z']})")
        print(f"     Size:   {box['size_x']} × "
              f"{box['size_y']} × {box['size_z']} Å")
        vol = (box['size_x'] * box['size_y'] *
               box['size_z'])
        print(f"     Volume: {vol:.0f} Å³")

    # ── Save ──────────────────────────────────────────────
    out = RES_DIR / "pocket_analysis_3.json"
    with open(out, "w") as f:
        json.dump({
            "complex":               "NSP12-NSP7",
            "pdb":                   "7BV2",
            "fpocket_pockets":       pockets[:10],
            "fpocket_all_structures": {k: v[:5] for k,v in all_pockets.items()},
            "water_bridges":         bridges,
            "hotspot_water_bridges": hotspot_b,
            "disulfides":            disulfides,
            "docking_box":           box,
        }, f, indent=2)

    print(f"\n  Saved: 02-validation/NSP12-NSP7/"
          f"pocket_analysis_3.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
