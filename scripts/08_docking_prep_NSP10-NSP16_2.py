"""
Script 08_2: Docking Preparation — NSP10-NSP16
===============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 6W4H (SARS-CoV-2, 1.80 Å)
Chain A = NSP16, Chain B = NSP10

Key decisions:
  - KEEP ZN atoms: Zn1 finger constrains HIS80 geometry
  - STRIP waters, ligands, HETATM (except ZN)
  - Verify all hotspot + Zn coordination residues present
  - Output Vina + VirtualFlow configs

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/08_docking_prep_NSP10-NSP16_2.py
"""

import json
from pathlib import Path
from Bio import PDB

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
RES_DIR  = PROJECT / "02-validation" / "NSP10-NSP16"
OUT_DIR  = PROJECT / "03-virtual-screening" / "NSP10-NSP16_2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PDB_FILE       = PDB_DIR / "6W4H.pdb"
CHAIN_NSP10    = "B"
CHAIN_NSP16    = "A"
NSP10_OFFSET   = 4270
NSP10_SHIFT    = 17
NSP16_OFFSET   = 6797

HOTSPOTS_NSP10 = [5,40,42,43,44,45,71,76,78,80,93,94,95,96]
HOTSPOTS_NSP16 = [40,41,44,76,83,87,102,104,106,244,247]
ZN1_COORD      = {74,77,83,90}
PRIMARY_NSP10  = {80,93,95}
PRIMARY_NSP16  = {102,106}


def genome_to_local(gpos, offset, shift=0):
    return gpos - offset + shift


class ZNSelecter(PDB.Select):
    """Keep chains A+B protein residues + ZN atoms only."""
    def accept_residue(self, residue):
        chain = residue.get_parent().id
        rname = residue.resname.strip()
        rid   = residue.id[0]
        # Keep protein residues in chains A and B
        if chain in ["A","B"] and rid == " ":
            return True
        # Keep ZN atoms in chain B (NSP10)
        if chain == "B" and rname in ["ZN","ZIN","ZNC"]:
            return True
        return False


def prepare_receptor(pdb_file):
    """Strip waters/ligands, keep protein + ZN, save receptor."""
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("rec", pdb_file)
    io        = PDB.PDBIO()
    io.set_structure(structure)
    out_pdb   = OUT_DIR / "receptor_NSP10-NSP16_2.pdb"
    io.save(str(out_pdb), ZNSelecter())
    return out_pdb, structure


def verify_residues(structure):
    """Verify all critical residues are present in receptor."""
    present_nsp10 = set()
    present_nsp16 = set()
    zn_present    = []

    for m in structure:
        for c in m:
            for r in c:
                rname = r.resname.strip()
                gpos  = r.id[1]
                if c.id == CHAIN_NSP10:
                    if r.id[0] == " ":
                        local = genome_to_local(
                            gpos, NSP10_OFFSET, NSP10_SHIFT)
                        present_nsp10.add(local)
                    elif rname in ["ZN","ZIN","ZNC"]:
                        zn_present.append(gpos)
                elif c.id == CHAIN_NSP16 and r.id[0] == " ":
                    local = genome_to_local(gpos, NSP16_OFFSET)
                    present_nsp16.add(local)

    print("\n  Hotspot residue verification:")
    print(f"  {'Residue':<12} {'Chain':<7} "
          f"{'Present':>8}  Notes")
    print(f"  {'-'*12} {'-'*7} {'-'*8}  {'-'*20}")

    # NSP10 positions < 18 are in the missing N-terminal region
    # (PDB 6W4H starts at AF3 local position 18 due to crystal truncation)
    EXPECTED_MISSING_NSP10 = {p for p in HOTSPOTS_NSP10 if p < 18}

    all_ok = True
    for pos in sorted(HOTSPOTS_NSP10):
        ok    = pos in present_nsp10
        if not ok and pos in EXPECTED_MISSING_NSP10:
            flag = "⚠️  N-term truncation (expected)"
        else:
            flag  = "✅" if ok else "❌ MISSING"
        notes = ""
        if pos in PRIMARY_NSP10:
            notes = "★ primary salt bridge"
        if pos in ZN1_COORD:
            notes += " Zn1 coord"
        print(f"  NSP10-{pos:<6} {'Chain B':<7} {flag:>8}  {notes}")
        if not ok and pos not in EXPECTED_MISSING_NSP10:
            all_ok = False

    for pos in sorted(HOTSPOTS_NSP16):
        ok    = pos in present_nsp16
        flag  = "✅" if ok else "❌ MISSING"
        notes = "★ primary salt bridge" if pos in PRIMARY_NSP16 else ""
        print(f"  NSP16-{pos:<6} {'Chain A':<7} {flag:>8}  {notes}")
        if not ok:
            all_ok = False

    print(f"\n  Zn1 coordination residues:")
    for pos in sorted(ZN1_COORD):
        ok   = pos in present_nsp10
        flag = "✅" if ok else "❌ MISSING"
        print(f"  NSP10-{pos:<6} {'Chain B':<7} {flag:>8}"
              f"  Zn1 coordinator")

    print(f"\n  ZN atoms retained: {len(zn_present)} "
          f"(genome pos: {zn_present})")
    return all_ok


def write_configs(box):
    """Write Vina and VirtualFlow config files."""
    # Vina config
    vina = OUT_DIR / "vina_config_NSP10-NSP16_2.txt"
    with open(vina, "w") as f:
        f.write(f"# AutoDock Vina config — NSP10-NSP16\n")
        f.write(f"# Primary salt bridges: "
                f"ASP106-LYS93, ASP106-LYS95, ASP102-HIS80\n")
        f.write(f"# Zn1 finger: HIS80 between CYS77 and HIS83\n")
        f.write(f"# Water bridges: HIS80 and LYS93 displaceble\n\n")
        f.write(f"receptor = receptor_NSP10-NSP16_2.pdb\n\n")
        f.write(f"center_x = {box['center_x']}\n")
        f.write(f"center_y = {box['center_y']}\n")
        f.write(f"center_z = {box['center_z']}\n\n")
        f.write(f"size_x   = {box['size_x']}\n")
        f.write(f"size_y   = {box['size_y']}\n")
        f.write(f"size_z   = {box['size_z']}\n\n")
        f.write(f"exhaustiveness = 16\n")
        f.write(f"num_modes      = 9\n")
        f.write(f"energy_range   = 3\n")

    # VirtualFlow config
    vf = OUT_DIR / "virtualflow_config_NSP10-NSP16_2.json"
    with open(vf, "w") as f:
        json.dump({
            "complex":         "NSP10-NSP16",
            "pdb":             "6W4H",
            "receptor":        "receptor_NSP10-NSP16_2.pdb",
            "center":          [box["center_x"],
                                box["center_y"],
                                box["center_z"]],
            "size":            [box["size_x"],
                                box["size_y"],
                                box["size_z"]],
            "primary_hotspots": {
                "NSP10": ["HIS80","LYS93","LYS95"],
                "NSP16": ["ASP106","ASP102"],
            },
            "zn1_site": {
                "coordinators": ["CYS74","CYS77","HIS83","CYS90"],
                "note": "HIS80 inside Zn1 loop — ZN retained in receptor"
            },
            "water_bridges": {
                "HIS80":  "HOH7319 — displacement pharmacophore",
                "LYS93":  "HOH4523 — displacement pharmacophore",
            },
            "salt_bridges": [
                "ASP106(NSP16)-LYS93(NSP10) all 4 structures",
                "ASP106(NSP16)-LYS95(NSP10) all 4 structures",
                "ASP102(NSP16)-HIS80(NSP10) 3/4 structures",
            ],
            "fpocket_score":   0.546,
            "exhaustiveness":  16,
            "num_modes":       9,
        }, f, indent=2)

    return vina, vf


def main():
    print("\n" + "="*57)
    print("  Script 08_2: Docking Preparation — NSP10-NSP16")
    print(f"  PDB: 6W4H (1.80 Å) | ZN atoms: RETAINED")
    print("="*57 + "\n")

    # Load pocket analysis for box
    pocket_file = RES_DIR / "pocket_analysis_2.json"
    with open(pocket_file) as f:
        pocket_data = json.load(f)
    box = pocket_data["docking_box"]

    print(f"  Docking box from Script 07_2:")
    print(f"    Center: ({box['center_x']}, "
          f"{box['center_y']}, {box['center_z']})")
    print(f"    Size:   {box['size_x']} × "
          f"{box['size_y']} × {box['size_z']} Å")

    # Prepare receptor
    print(f"\n  Preparing receptor (protein + ZN only)...")
    out_pdb, structure = prepare_receptor(PDB_FILE)

    # Count atoms
    parser2   = PDB.PDBParser(QUIET=True)
    rec_struct = parser2.get_structure("rec", out_pdb)
    n_atoms   = sum(1 for m in rec_struct for c in m
                    for r in c for a in r)
    n_res     = sum(1 for m in rec_struct for c in m
                    for r in c if r.id[0] == " ")
    n_zn      = sum(1 for m in rec_struct for c in m
                    for r in c
                    if r.resname.strip() in ["ZN","ZIN","ZNC"])
    print(f"    Residues: {n_res}")
    print(f"    Atoms:    {n_atoms}")
    print(f"    ZN ions:  {n_zn} (retained ✅)")

    # Verify
    all_ok = verify_residues(structure)
    print(f"\n  All critical residues present: "
          f"{'✅' if all_ok else '❌'}")

    # Write configs
    print(f"\n  Writing Vina + VirtualFlow configs...")
    vina, vf = write_configs(box)
    print(f"    {vina.name}")
    print(f"    {vf.name}")

    print(f"\n  Files saved to: 03-virtual-screening/NSP10-NSP16_2/")
    print(f"    receptor_NSP10-NSP16_2.pdb")
    print(f"    vina_config_NSP10-NSP16_2.txt")
    print(f"    virtualflow_config_NSP10-NSP16_2.json")
    print("="*57 + "\n")

    if not all_ok:
        print("  ⚠️  Missing residues — review before docking")


if __name__ == "__main__":
    main()
