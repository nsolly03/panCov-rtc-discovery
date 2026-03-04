"""
Script 08_5: Docking Preparation — NSP9-NSP12
=============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 8SQK (Chain A=NSP12, Chain G=NSP9)
Note: Chain G unusual — must specify explicitly

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/08_docking_prep_NSP9-NSP12_5.py
"""

import json
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, Select

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
RES_DIR  = PROJECT / "02-validation" / "NSP9-NSP12"
DOCK_DIR = PROJECT / "03-virtual-screening" / "NSP9-NSP12_5"
DOCK_DIR.mkdir(parents=True, exist_ok=True)

PDB_FILE      = PDB_DIR / "8SQK_NSP9-NSP12.pdb"
KEEP_CHAINS   = {"A","G"}
CHAIN_NSP12   = "A"
CHAIN_NSP9    = "G"

HOTSPOTS_NSP12 = [38,1,3,4,96,733,202,103,
                   221,233,291,2,223]
HOTSPOTS_NSP9  = [38,1,3,4,96,103,2,97]
PRIMARY_NSP12  = {740,744}
PRIMARY_NSP9   = {36}


class ReceptorSelect(Select):
    def accept_chain(self, chain):
        return chain.id in KEEP_CHAINS
    def accept_residue(self, residue):
        return (residue.id[0] == " " and
                residue.resname.strip() != "HOH")
    def accept_atom(self, atom):
        return (atom.altloc == " " or
                atom.altloc == "A")


def main():
    print("\n" + "="*57)
    print("  Script 08_5: Docking Prep — NSP9-NSP12")
    print("  PDB: 8SQK (Chain A=NSP12, Chain G=NSP9)")
    print("="*57 + "\n")

    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("8SQK", PDB_FILE)

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
                     else "NSP9")
            print(f"    Chain {c.id} ({label}): "
                  f"{len(res)} residues, {ats} atoms")
            total_res  += len(res)
            total_atoms += ats
    print(f"    Total: {total_res} residues, "
          f"{total_atoms} atoms")

    # Write receptor
    io  = PDBIO()
    io.set_structure(structure)
    rec = DOCK_DIR / "receptor_NSP9-NSP12_5.pdb"
    io.save(str(rec), ReceptorSelect())
    print(f"\n  Receptor saved: {rec.name}")

    # Verify hotspots
    clean  = parser.get_structure("rec", rec)
    res_12 = {r.id[1]: r.resname
              for m in clean for c in m
              if c.id == "A"
              for r in c if r.id[0] == " "}
    res_9  = {r.id[1]: r.resname
              for m in clean for c in m
              if c.id == "G"
              for r in c if r.id[0] == " "}

    print("\n  Hotspot residue verification:")
    missing12 = [p for p in HOTSPOTS_NSP12
                 if p not in res_12]
    missing9  = [p for p in HOTSPOTS_NSP9
                 if p not in res_9]
    if missing12:
        print(f"    NSP12 missing: {missing12}")
    else:
        print(f"    NSP12: all hotspots present ✅")
    if missing9:
        print(f"    NSP9  missing: {missing9}")
    else:
        print(f"    NSP9:  all hotspots present ✅")

    print(f"\n  Key interface residues:")
    for pos in sorted({733,740,744} &
                       set(res_12.keys())):
        print(f"    {res_12[pos]}{pos} (NSP12) ✅")
    for pos in sorted({36,96,97,103} &
                       set(res_9.keys())):
        print(f"    {res_9[pos]}{pos} (NSP9)  ✅")

    # Load docking box
    pocket = json.load(
        open(RES_DIR / "pocket_analysis_5.json"))
    box    = pocket["docking_box"]

    # Vina config
    vina_cfg = DOCK_DIR / \
               "vina_config_NSP9-NSP12_5.txt"
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
             "virtualflow_config_NSP9-NSP12_5.json"
    with open(vf_cfg,"w") as f:
        json.dump({
            "target":         "NSP9-NSP12",
            "receptor":       str(rec),
            "center":         [box["center_x"],
                               box["center_y"],
                               box["center_z"]],
            "size":           [box["size_x"],
                               box["size_y"],
                               box["size_z"]],
            "exhaustiveness":  16,
            "druggability":    0.895,
            "note":            "novel NiRAN domain "
                               "interface, single PDB",
            "primary_sb_af3":  ["ASP740-LYS36",
                                 "GLU744-LYS36"],
            "pan_cov_hotspots":["ARG733","VAL202",
                                 "ASP221","LEU97",
                                 "LEU103"],
            "hotspots_NSP12":  HOTSPOTS_NSP12,
            "hotspots_NSP9":   HOTSPOTS_NSP9,
        }, f, indent=2)
    print(f"  VF config:   {vf_cfg.name}")

    print(f"\n  Docking box:")
    print(f"    Center: ({box['center_x']}, "
          f"{box['center_y']}, {box['center_z']})")
    print(f"    Size:   {box['size_x']} × "
          f"{box['size_y']} × {box['size_z']} Å")

    print(f"\n  Output: 03-virtual-screening/"
          f"NSP9-NSP12_5/")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
