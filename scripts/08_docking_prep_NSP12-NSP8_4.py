"""
Script 08_4: Docking Preparation — NSP12-NSP8
=============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2_NSP12-NSP8 (primary, 2.90 Å)
Chain A = NSP12, Chain B = NSP8

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/08_docking_prep_NSP12-NSP8_4.py
"""

import json
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, Select

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP8"
DOCK_DIR = PROJECT / "03-virtual-screening" / "NSP12-NSP8_4"
DOCK_DIR.mkdir(parents=True, exist_ok=True)

PDB_FILE      = PDB_DIR / "7BV2_NSP12-NSP8.pdb"
KEEP_CHAINS   = {"A","B"}
CHAIN_NSP12   = "A"
CHAIN_NSP8    = "B"
NSP12_OFFSET  = 0
NSP8_OFFSET   = 77

HOTSPOTS_NSP12 = [387,129,389,271,330,131,380,523,
                   91,87,332,95,117,517,99]
HOTSPOTS_NSP8  = [117,129,80,115,131,112,91,87,
                   116,95,98,83,128,90,121]
PRIMARY_NSP12  = {523,332,517}
PRIMARY_NSP8   = {80,99,79}


class ReceptorSelect(Select):
    def accept_chain(self, chain):
        return chain.id in KEEP_CHAINS
    def accept_residue(self, residue):
        return (residue.id[0] == " " and
                residue.resname.strip() != "HOH")
    def accept_atom(self, atom):
        return (atom.altloc == " " or
                atom.altloc == "A")


def pdb_to_local(gpos, chain):
    if chain == CHAIN_NSP12:
        return gpos - NSP12_OFFSET
    return gpos - NSP8_OFFSET


def main():
    print("\n" + "="*57)
    print("  Script 08_4: Docking Prep — NSP12-NSP8")
    print("  PDB: 7BV2_NSP12-NSP8")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("7BV2", PDB_FILE)

    # Receptor stats
    print("  Receptor composition:")
    total_res, total_atoms = 0, 0
    for m in structure:
        for c in m:
            if c.id not in KEEP_CHAINS:
                continue
            res = [r for r in c if r.id[0] == " "]
            ats = sum(len(list(r.get_atoms()))
                      for r in res)
            label = ("NSP12" if c.id == "A"
                     else "NSP8")
            print(f"    Chain {c.id} ({label}): "
                  f"{len(res)} residues, {ats} atoms")
            total_res  += len(res)
            total_atoms += ats
    print(f"    Total: {total_res} residues, "
          f"{total_atoms} atoms")

    # Write receptor
    io  = PDBIO()
    io.set_structure(structure)
    rec = DOCK_DIR / "receptor_NSP12-NSP8_4.pdb"
    io.save(str(rec), ReceptorSelect())
    print(f"\n  Receptor saved: {rec.name}")

    # Verify hotspots
    print("\n  Hotspot residue verification:")
    clean  = parser.get_structure("rec", rec)
    res_12 = {pdb_to_local(r.id[1], "A"): r.resname
              for m in clean for c in m
              if c.id == "A"
              for r in c if r.id[0] == " "}
    res_8  = {pdb_to_local(r.id[1], "B"): r.resname
              for m in clean for c in m
              if c.id == "B"
              for r in c if r.id[0] == " "}

    missing12 = [p for p in HOTSPOTS_NSP12
                 if p not in res_12]
    missing8  = [p for p in HOTSPOTS_NSP8
                 if p not in res_8]

    if missing12:
        print(f"    ⚠️  NSP12 missing: {missing12}")
    else:
        print(f"    NSP12: all hotspots present ✅")
    if missing8:
        print(f"    ⚠️  NSP8  missing: {missing8}")
    else:
        print(f"    NSP8:  all hotspots present ✅")

    # Primary SB residues
    print(f"\n  Primary salt bridge residues:")
    for pos in sorted(PRIMARY_NSP12):
        aa = res_12.get(pos,"???")
        flag = "✅" if pos in res_12 else "⚠️  MISSING"
        print(f"    {aa}{pos} (NSP12) {flag}")
    for pos in sorted(PRIMARY_NSP8):
        aa = res_8.get(pos,"???")
        flag = "✅" if pos in res_8 else "⚠️  MISSING"
        print(f"    {aa}{pos} (NSP8)  {flag}")

    # Load docking box
    pocket = json.load(
        open(RES_DIR / "pocket_analysis_4.json"))
    box    = pocket["docking_box"]

    # Vina config
    vina_cfg = DOCK_DIR / "vina_config_NSP12-NSP8_4.txt"
    with open(vina_cfg,"w") as f:
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

    # VirtualFlow config
    vf_cfg = DOCK_DIR / \
             "virtualflow_config_NSP12-NSP8_4.json"
    with open(vf_cfg,"w") as f:
        json.dump({
            "target":       "NSP12-NSP8",
            "receptor":     str(rec),
            "center":       [box["center_x"],
                             box["center_y"],
                             box["center_z"]],
            "size":         [box["size_x"],
                             box["size_y"],
                             box["size_z"]],
            "exhaustiveness": 16,
            "druggability":   0.874,
            "primary_sb":   ["ASP523-ARG80",
                             "LYS332-ASP99",
                             "ASP517-LYS79"],
            "hotspots_NSP12": HOTSPOTS_NSP12,
            "hotspots_NSP8":  HOTSPOTS_NSP8,
        }, f, indent=2)
    print(f"  VF config:   {vf_cfg.name}")

    print(f"\n  Docking box:")
    print(f"    Center: ({box['center_x']}, "
          f"{box['center_y']}, {box['center_z']})")
    print(f"    Size:   {box['size_x']} × "
          f"{box['size_y']} × {box['size_z']} Å")

    print(f"\n  Output: 03-virtual-screening/"
          f"NSP12-NSP8_4/")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
