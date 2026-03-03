"""
Script 07_2: Pocket Detection — NSP10-NSP16
===========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 6W4H (primary), 6WVN, 6WKQ + AF3 model
Analyses:
  1. fpocket — druggable pocket detection
  2. Zinc coordination — hotspots near Zn finger
  3. Water-mediated contacts — bridging waters at interface
  4. Disulfide bridges — structural rigidity check
  5. Docking box definition from best pocket

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/07_pocket_NSP10-NSP16_2.py
"""

import json
import subprocess
import numpy as np
from pathlib import Path
from Bio import PDB

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP10-NSP16"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
TMP_DIR = PROJECT / "tmp"
TMP_DIR.mkdir(exist_ok=True)

# Hotspots in AF3 local numbering
HOTSPOTS_NSP10 = [5,40,42,43,44,45,71,76,78,80,93,94,95,96]
HOTSPOTS_NSP16 = [40,41,44,76,83,87,102,104,106,244,247]
PRIMARY        = {("NSP10",80),("NSP10",93),("NSP10",95),
                  ("NSP16",102),("NSP16",106)}

# Zn1 coordination residues (AF3 local)
ZN1_COORD = {74,77,83,90}   # CYS74,CYS77,HIS83,CYS90
ZN2_COORD = {117,120,128,130}

# PDB genome offset for NSP10
NSP10_OFFSET = 4270
NSP10_SHIFT  = 17
NSP16_OFFSET = 6797


def genome_to_local(gpos, offset, shift=0):
    return gpos - offset + shift


# ── 1. fpocket analysis ────────────────────────────────────

def run_fpocket(pdb_file, label):
    out_dir = TMP_DIR / f"fpocket_{label}"
    if out_dir.exists():
        import shutil
        shutil.rmtree(out_dir)
    result = subprocess.run(
        ["fpocket", "-f", str(pdb_file),
         "-o", str(out_dir)],
        capture_output=True, text=True)
    if result.returncode != 0:
        # Try without -o flag (older fpocket)
        result = subprocess.run(
            ["fpocket", "-f", str(pdb_file)],
            capture_output=True, text=True,
            cwd=str(TMP_DIR))
    return out_dir


def parse_fpocket(pdb_file, label):
    """Run fpocket and parse pocket info file."""
    stem      = pdb_file.stem
    info_file = TMP_DIR / f"{stem}_out" / f"{stem}_info.txt"

    if not info_file.exists():
        run_fpocket(pdb_file, label)
        # Try alternate output location
        alt = pdb_file.parent / f"{stem}_out" / f"{stem}_info.txt"
        if alt.exists():
            info_file = alt

    pockets = []
    if not info_file.exists():
        print(f"    fpocket output not found for {label}")
        return pockets

    current = {}
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


def get_pocket_residues(pdb_file, pocket_id):
    """Get residues in a specific fpocket pocket."""
    stem     = pdb_file.stem
    pock_pdb = (TMP_DIR / f"{stem}_out" / "pockets" /
                f"pocket{pocket_id}_atm.pdb")
    if not pock_pdb.exists():
        return []
    parser  = PDB.PDBParser(QUIET=True)
    struct  = parser.get_structure("p", pock_pdb)
    resnums = [(c.id, r.id[1])
               for m in struct for c in m
               for r in c if r.id[0] == " "]
    return list(set(resnums))


# ── 2. Zinc coordination ───────────────────────────────────

def find_zinc_coordination(structure, chain_nsp10,
                            nsp10_offset, nsp10_shift):
    """Find zinc atoms and coordinating residues."""
    zn_sites = []
    zn_atoms = [
        (c.id, r.id[1], list(r.get_atoms())[0])
        for m in structure for c in m
        for r in c
        if r.resname.strip() in ["ZN","ZIN","ZNC"]]

    all_atoms = [
        (c.id, r, a)
        for m in structure for c in m
        for r in c if r.id[0] == " "
        for a in r]

    for zch, zrn, zn_at in zn_atoms:
        coords = []
        for ch, res, at in all_atoms:
            try:
                d = zn_at - at
                if d <= 2.6:
                    gpos  = res.id[1]
                    local = genome_to_local(
                        gpos, nsp10_offset, nsp10_shift)
                    is_hotspot = local in HOTSPOTS_NSP10
                    is_zn_coord = local in (ZN1_COORD | ZN2_COORD)
                    coords.append({
                        "chain":      ch,
                        "genome_pos": gpos,
                        "local_pos":  local,
                        "resname":    res.resname.strip(),
                        "atom":       at.name,
                        "distance":   round(d, 2),
                        "is_hotspot": is_hotspot,
                        "is_zn_coord": is_zn_coord,
                    })
            except Exception:
                pass
        zn_sites.append({
            "zn_chain":  zch,
            "zn_resnum": zrn,
            "coordinators": sorted(coords,
                                   key=lambda x: x["distance"])
        })
    return zn_sites


# ── 3. Water-mediated contacts ─────────────────────────────

def find_water_bridges(structure, chain_nsp10, chain_nsp16,
                        nsp10_offset, nsp10_shift,
                        nsp16_offset):
    """
    Find water molecules bridging NSP10-NSP16 interface.
    Water bridge: HOH within 3.5Å of both NSP10 and NSP16 residues.
    """
    waters = [
        (c.id, r.id[1], list(r.get_atoms()))
        for m in structure for c in m
        for r in c
        if r.resname.strip() in ["HOH","WAT","SOL"]]

    res_nsp10 = {
        r.id[1]: r for m in structure for c in m
        if c.id == chain_nsp10
        for r in c if r.id[0] == " "}
    res_nsp16 = {
        r.id[1]: r for m in structure for c in m
        if c.id == chain_nsp16
        for r in c if r.id[0] == " "}

    bridges = []
    for wch, wrn, wat_atoms in waters:
        contacts10, contacts16 = [], []
        for wa in wat_atoms:
            for rn10, res10 in res_nsp10.items():
                for a10 in res10.get_atoms():
                    try:
                        d = wa - a10
                        if d <= 3.5:
                            local = genome_to_local(
                                rn10, nsp10_offset,
                                nsp10_shift)
                            contacts10.append(
                                (local, res10.resname.strip(),
                                 round(d,2)))
                    except Exception:
                        pass
            for rn16, res16 in res_nsp16.items():
                for a16 in res16.get_atoms():
                    try:
                        d = wa - a16
                        if d <= 3.5:
                            local = genome_to_local(
                                rn16, nsp16_offset, 0)
                            contacts16.append(
                                (local, res16.resname.strip(),
                                 round(d,2)))
                    except Exception:
                        pass
        if contacts10 and contacts16:
            bridges.append({
                "water": f"HOH{wrn}",
                "nsp10_contacts": list(set(contacts10))[:3],
                "nsp16_contacts": list(set(contacts16))[:3],
            })
    return bridges


# ── 4. Disulfide bridges ───────────────────────────────────

def find_disulfides(structure, chain_nsp10,
                    nsp10_offset, nsp10_shift):
    """Find disulfide bonds in NSP10 (SG-SG < 2.1 Å)."""
    cys_res = [
        r for m in structure for c in m
        if c.id == chain_nsp10
        for r in c
        if r.resname.strip() == "CYS"
        and r.id[0] == " "]

    disulfides = []
    for i, r1 in enumerate(cys_res):
        for r2 in cys_res[i+1:]:
            try:
                sg1 = r1["SG"]
                sg2 = r2["SG"]
                d   = sg1 - sg2
                if d <= 2.1:
                    l1 = genome_to_local(
                        r1.id[1], nsp10_offset, nsp10_shift)
                    l2 = genome_to_local(
                        r2.id[1], nsp10_offset, nsp10_shift)
                    disulfides.append({
                        "res1": f"CYS{l1}",
                        "res2": f"CYS{l2}",
                        "distance": round(d, 2),
                        "near_zn1": (l1 in ZN1_COORD or
                                     l2 in ZN1_COORD),
                    })
            except Exception:
                pass
    return disulfides


# ── 5. Docking box ─────────────────────────────────────────

def define_docking_box(structure, chain_nsp10, chain_nsp16,
                        nsp10_offset, nsp10_shift,
                        nsp16_offset, padding=8.0):
    """
    Define docking box centered on primary hotspot residues.
    Includes Zn1 coordination site for NSP10.
    """
    primary_nsp10 = {80, 93, 95}
    primary_nsp16 = {102, 106}

    coords = []
    for m in structure:
        for c in m:
            if c.id == chain_nsp10:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    local = genome_to_local(
                        r.id[1], nsp10_offset, nsp10_shift)
                    if local in primary_nsp10:
                        coords.extend(
                            [a.get_vector().get_array()
                             for a in r.get_atoms()])
            elif c.id == chain_nsp16:
                for r in c:
                    if r.id[0] != " ":
                        continue
                    local = genome_to_local(
                        r.id[1], nsp16_offset, 0)
                    if local in primary_nsp16:
                        coords.extend(
                            [a.get_vector().get_array()
                             for a in r.get_atoms()])

    if not coords:
        return None
    coords = np.array(coords)
    center = coords.mean(axis=0)
    size   = coords.max(axis=0) - coords.min(axis=0) + 2*padding
    return {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x":   round(float(size[0]),   3),
        "size_y":   round(float(size[1]),   3),
        "size_z":   round(float(size[2]),   3),
    }


# ── Main ───────────────────────────────────────────────────

def main():
    print("\n" + "="*57)
    print("  Script 07_2: Pocket Detection — NSP10-NSP16")
    print("="*57 + "\n")

    parser_pdb = PDB.PDBParser(QUIET=True)
    pdb_file   = PDB_DIR / "6W4H.pdb"
    structure  = parser_pdb.get_structure("6W4H", pdb_file)

    results = {}

    # ── fpocket ───────────────────────────────────────────
    print("  1. fpocket analysis on 6W4H...")
    pockets = parse_fpocket(pdb_file, "6W4H")
    if pockets:
        print(f"     Found {len(pockets)} pockets")
        print(f"\n     {'Pocket':<8} {'Score':>8} "
              f"{'DrugScore':>10} {'Volume':>10}")
        print(f"     {'-'*8} {'-'*8} {'-'*10} {'-'*10}")
        for p in pockets[:10]:
            pid   = p.get("id","?")
            score = p.get("Druggability Score",
                    p.get("Drug Score", 0))
            vol   = p.get("Volume", 0)
            dscore= p.get("Drug Score", score)
            print(f"     {pid:<8} {score:>8.3f} "
                  f"{dscore:>10.3f} {vol:>10.1f}")
        results["fpocket_pockets"] = pockets[:10]
    else:
        print("     No pockets found or fpocket not available")
        results["fpocket_pockets"] = []

    # ── Zinc coordination ─────────────────────────────────
    print("\n  2. Zinc coordination analysis (6W4H)...")
    zn_sites = find_zinc_coordination(
        structure, "B", NSP10_OFFSET, NSP10_SHIFT)
    print(f"     Found {len(zn_sites)} Zn sites")
    for site in zn_sites:
        print(f"\n     Zn at {site['zn_chain']}"
              f"{site['zn_resnum']}:")
        for c in site["coordinators"]:
            hs = " ← HOTSPOT" if c["is_hotspot"] else ""
            zc = " ← Zn coord" if c["is_zn_coord"] else ""
            print(f"       {c['resname']}{c['local_pos']:<5} "
                  f"({c['atom']}) {c['distance']:.2f}Å"
                  f"{hs}{zc}")
    results["zinc_sites"] = zn_sites

    # ── Water bridges ─────────────────────────────────────
    print("\n  3. Water-mediated interface contacts (6W4H)...")
    bridges = find_water_bridges(
        structure, "B", "A",
        NSP10_OFFSET, NSP10_SHIFT, NSP16_OFFSET)
    hotspot_bridges = [
        b for b in bridges
        if any(pos in HOTSPOTS_NSP10
               for pos, _, _ in b["nsp10_contacts"])
        or any(pos in HOTSPOTS_NSP16
               for pos, _, _ in b["nsp16_contacts"])]
    print(f"     Total water bridges: {len(bridges)}")
    print(f"     Bridges at hotspot residues: "
          f"{len(hotspot_bridges)}")
    if hotspot_bridges:
        print(f"\n     Key water bridges:")
        for b in hotspot_bridges[:8]:
            n10 = [(f"{rn}{p}", d)
                   for p,rn,d in b["nsp10_contacts"]]
            n16 = [(f"{rn}{p}", d)
                   for p,rn,d in b["nsp16_contacts"]]
            print(f"       {b['water']}: "
                  f"NSP10={n10} ↔ NSP16={n16}")
    results["water_bridges"]        = bridges
    results["hotspot_water_bridges"] = hotspot_bridges

    # ── Disulfides ────────────────────────────────────────
    print("\n  4. Disulfide bridge check (6W4H)...")
    disulfides = find_disulfides(
        structure, "B", NSP10_OFFSET, NSP10_SHIFT)
    if disulfides:
        for ds in disulfides:
            zn = " ← near Zn1" if ds["near_zn1"] else ""
            print(f"     {ds['res1']} — {ds['res2']} "
                  f"({ds['distance']:.2f}Å){zn}")
    else:
        print("     No disulfide bonds found "
              "(Zn finger = reduced cysteines, expected)")
    results["disulfides"] = disulfides

    # ── Docking box ───────────────────────────────────────
    print("\n  5. Docking box definition...")
    box = define_docking_box(
        structure, "B", "A",
        NSP10_OFFSET, NSP10_SHIFT, NSP16_OFFSET)
    if box:
        print(f"     Center: ({box['center_x']}, "
              f"{box['center_y']}, {box['center_z']})")
        print(f"     Size:   {box['size_x']} × "
              f"{box['size_y']} × {box['size_z']} Å")
        vol = (box['size_x'] * box['size_y'] *
               box['size_z'])
        print(f"     Volume: {vol:.0f} Å³")
    results["docking_box"] = box

    # ── Save ──────────────────────────────────────────────
    out = RES_DIR / "pocket_analysis_2.json"

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.floating, np.float32,
                                np.float64)):
                return float(obj)
            if isinstance(obj, np.integer):
                return int(obj)
            return super().default(obj)

    with open(out, "w") as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)
    print(f"\n  Saved: 02-validation/NSP10-NSP16/"
          f"pocket_analysis_2.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
