"""
Script 08_3: Docking Preparation — NSP12-NSP7
=============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2_NSP12-NSP7 (primary)
Chain A = NSP12, Chain C = NSP7

Steps:
  1. Extract receptor (NSP12 + NSP7 chains only)
  2. Remove waters, ligands, RNA chains
  3. Verify all hotspot residues present
  4. Write VirtualFlow + Vina config files

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/08_docking_prep_NSP12-NSP7_3.py
"""

import json
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, Select

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP7"
DOCK_DIR = PROJECT / "03-virtual-screening" / "NSP12-NSP7_3"
DOCK_DIR.mkdir(parents=True, exist_ok=True)

PDB_FILE     = PDB_DIR / "7BV2_NSP12-NSP7.pdb"
KEEP_CHAINS  = {"A", "C"}  # NSP12, NSP7

HOTSPOTS_NSP12  = [440,412,442,443,420,843,409,
                    40,33,41,37,413,415,14,23]
HOTSPOTS_NSP7   = [40,14,33,41,37,11,23,5,
                    15,29,12,4,1]
PRIMARY_NSP12   = {431}
PRIMARY_NSP7    = {1}  # LYS at PDB res 2 = local 1
NSP12_OFFSET    = 0
NSP7_OFFSET     = 1


class ReceptorSelect(Select):
    def accept_chain(self, chain):
        return chain.id in KEEP_CHAINS

    def accept_residue(self, residue):
        return (residue.id[0] == " " and
                residue.resname.strip() != "HOH")

    def accept_atom(self, atom):
        # Keep only first alternate location
        return (atom.altloc == " " or
                atom.altloc == "A")


def genome_to_local(gpos, chain):
    if chain == "A":
        return gpos - NSP12_OFFSET
    return gpos - NSP7_OFFSET


def main():
    print("\n" + "="*57)
    print("  Script 08_3: Docking Prep — NSP12-NSP7")
    print("  PDB: 7BV2_NSP12-NSP7")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("7BV2", PDB_FILE)

    # ── Receptor stats ────────────────────────────────────
    print("  Receptor composition:")
    total_res, total_atoms = 0, 0
    for m in structure:
        for c in m:
            if c.id not in KEEP_CHAINS:
                continue
            res = [r for r in c if r.id[0] == " "]
            ats = sum(len(list(r.get_atoms()))
                      for r in res)
            label = ("NSP12" if c.id == "A" else "NSP7")
            print(f"    Chain {c.id} ({label}): "
                  f"{len(res)} residues, {ats} atoms")
            total_res  += len(res)
            total_atoms += ats
    print(f"    Total: {total_res} residues, "
          f"{total_atoms} atoms")

    # ── Write receptor PDB ────────────────────────────────
    io  = PDBIO()
    io.set_structure(structure)
    rec = DOCK_DIR / "receptor_NSP12-NSP7_3.pdb"
    io.save(str(rec), ReceptorSelect())
    print(f"\n  Receptor saved: {rec.name}")

    # ── Verify hotspots ───────────────────────────────────
    print("\n  Hotspot residue verification:")
    clean = parser.get_structure("rec", rec)
    res_12 = {genome_to_local(r.id[1], "A"): r.resname
              for m in clean for c in m
              if c.id == "A"
              for r in c if r.id[0] == " "}
    res_7  = {genome_to_local(r.id[1], "C"): r.resname
              for m in clean for c in m
              if c.id == "C"
              for r in c if r.id[0] == " "}

    missing12, missing7 = [], []
    for pos in HOTSPOTS_NSP12 + list(PRIMARY_NSP12):
        if pos not in res_12:
            missing12.append(pos)
    for pos in HOTSPOTS_NSP7 + list(PRIMARY_NSP7):
        if pos not in res_7:
            missing7.append(pos)

    if missing12:
        print(f"    ⚠️  NSP12 missing: {missing12}")
    else:
        print(f"    NSP12: all hotspots present ✅")
    if missing7:
        print(f"    ⚠️  NSP7  missing: {missing7}")
    else:
        print(f"    NSP7:  all hotspots present ✅")

    # Print primary residues
    for pos in sorted(PRIMARY_NSP12):
        aa = res_12.get(pos, "???")
        print(f"    GLU{pos} (NSP12): {aa} ✅"
              if pos in res_12 else
              f"    GLU{pos} (NSP12): MISSING ⚠️")
    for pos in sorted(PRIMARY_NSP7):
        aa = res_7.get(pos, "???")
        print(f"    LYS{pos} (NSP7):  {aa} ✅"
              if pos in res_7 else
              f"    LYS{pos} (NSP7):  MISSING ⚠️")

    # ── Load docking box ──────────────────────────────────
    pocket = json.load(
        open(RES_DIR / "pocket_analysis_3.json"))
    box    = pocket["docking_box"]

    # ── Vina config ───────────────────────────────────────
    vina_cfg = DOCK_DIR / "vina_config_NSP12-NSP7_3.txt"
    with open(vina_cfg, "w") as f:
        f.write(f"receptor = {rec}\n\n")
        f.write(f"center_x = {box['center_x']}\n")
        f.write(f"center_y = {box['center_y']}\n")
        f.write(f"center_z = {box['center_z']}\n\n")
        f.write(f"size_x = {box['size_x']}\n")
        f.write(f"size_y = {box['size_y']}\n")
        f.write(f"size_z = {box['size_z']}\n\n")
        f.write(f"exhaustiveness = 16\n")
        f.write(f"num_modes = 9\n")
        f.write(f"energy_range = 3\n")
    print(f"\n  Vina config: {vina_cfg.name}")

    # ── VirtualFlow config ────────────────────────────────
    vf_cfg = DOCK_DIR / "virtualflow_config_NSP12-NSP7_3.json"
    with open(vf_cfg, "w") as f:
        json.dump({
            "target":    "NSP12-NSP7",
            "receptor":  str(rec),
            "center":    [box["center_x"],
                          box["center_y"],
                          box["center_z"]],
            "size":      [box["size_x"],
                          box["size_y"],
                          box["size_z"]],
            "exhaustiveness": 16,
            "druggability":   0.874,
            "primary_sb":     "GLU431(NSP12)-LYS2(NSP7)",
            "hotspots_NSP12": HOTSPOTS_NSP12,
            "hotspots_NSP7":  HOTSPOTS_NSP7,
        }, f, indent=2)
    print(f"  VF config:   {vf_cfg.name}")

    print(f"\n  Docking box:")
    print(f"    Center: ({box['center_x']}, "
          f"{box['center_y']}, {box['center_z']})")
    print(f"    Size:   {box['size_x']} × "
          f"{box['size_y']} × {box['size_z']} Å")

    print(f"\n  Output directory: "
          f"03-virtual-screening/NSP12-NSP7_3/")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
