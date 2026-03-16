"""
Script 08_9 — Docking Preparation NSP15 Homodimer
Complex  : NSP15-NSP15 (NendoU dimer interface)
Suffix   : _9
Receptor : 6VWW_NSP15-true-dimer.pdb (chains A+B, biological assembly)
Strategy : Large docking box spanning distributed PPI interface
"""

import json
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

OUT_DIR = "03-virtual-screening/NSP15_9"
os.makedirs(OUT_DIR, exist_ok=True)

IN_PDB   = "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb"
OUT_PDB  = f"{OUT_DIR}/receptor_NSP15_9.pdb"
OUT_JSON = "02-validation/NSP15/pocket_analysis_9.json"

print("=" * 60)
print("Script 08_9 — Docking Preparation NSP15")
print("=" * 60)

# ── Strip waters and HETATM ────────────────────────────────────
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.get_id()[0] == ' '

pdb = PDBParser(QUIET=True).get_structure("6VWW", IN_PDB)
io  = PDBIO()
io.set_structure(pdb)
io.save(OUT_PDB, ProteinSelect())

# ── Verify receptor ────────────────────────────────────────────
rec = PDBParser(QUIET=True).get_structure("rec", OUT_PDB)
chains  = list(rec[0].get_chains())
all_res = [r for c in rec[0].get_chains()
             for r in c.get_residues() if r.get_id()[0]==' ']
all_atm = list(rec[0].get_atoms())

print(f"\nReceptor: {OUT_PDB}")
print(f"  Chains    : {[c.id for c in chains]}")
print(f"  Residues  : {len(all_res)}")
print(f"  Atoms     : {len(all_atm)}")

for c in chains:
    res = [r for r in c.get_residues() if r.get_id()[0]==' ']
    print(f"  Chain {c.id} : {len(res)} residues  "
          f"(pos {res[0].get_id()[1]}-{res[-1].get_id()[1]})")

# ── Verify primary hotspots present ───────────────────────────
PRIMARY = {40:"ASP40", 62:"ARG62", 91:"ARG91", 267:"GLU267"}
print("\nPrimary SB hotspot verification:")
res_map = {}
for c in rec[0].get_chains():
    for r in c.get_residues():
        if r.get_id()[0] == ' ':
            pos = r.get_id()[1]
            if pos not in res_map:
                res_map[pos] = r.get_resname()

for pos, name in PRIMARY.items():
    found = pos in res_map
    aa    = res_map.get(pos, "MISSING")
    status = "✓" if found else "✗ MISSING"
    print(f"  {name:<12}: pos={pos}  AA={aa}  {status}")

# ── Load docking box from pocket analysis ─────────────────────
with open(OUT_JSON) as f:
    pocket_data = json.load(f)

box    = pocket_data["docking_box"]
center = box["center"]
size   = box["size"]
volume = box["volume"]

print(f"\nDocking box (from pocket_analysis_9.json):")
print(f"  Center : ({center[0]}, {center[1]}, {center[2]})")
print(f"  Size   : {size[0]} x {size[1]} x {size[2]} A")
print(f"  Volume : {volume} A3")
print(f"  Strategy: {box.get('strategy','large box for distributed PPI interface')}")

# ── Write Vina config ──────────────────────────────────────────
vina_config = f"""receptor = receptor_NSP15_9.pdbqt

center_x = {center[0]}
center_y = {center[1]}
center_z = {center[2]}

size_x = {size[0]}
size_y = {size[1]}
size_z = {size[2]}

exhaustiveness = 16
num_modes = 9
energy_range = 3
"""

vina_path = f"{OUT_DIR}/vina_config_NSP15_9.txt"
with open(vina_path, "w") as f:
    f.write(vina_config)
print(f"\nVina config saved: {vina_path}")

# ── Write VirtualFlow config ───────────────────────────────────
vf_config = {
    "receptor"    : "receptor_NSP15_9.pdbqt",
    "center_x"    : center[0],
    "center_y"    : center[1],
    "center_z"    : center[2],
    "size_x"      : size[0],
    "size_y"      : size[1],
    "size_z"      : size[2],
    "exhaustiveness": 16,
    "complex"     : "NSP15-homodimer",
    "suffix"      : "_9",
    "strategy"    : "Fragment-based PPI interface disruption",
    "primary_pharmacophore": "ASP40-ARG62 salt bridge (pan-cov)",
}

vf_path = f"{OUT_DIR}/virtualflow_config_NSP15_9.json"
with open(vf_path, "w") as f:
    json.dump(vf_config, f, indent=2)
print(f"VF config saved : {vf_path}")

print(f"\nNSP15 docking prep COMPLETE")
print(f"Next: convert receptor to PDBQT with OpenBabel")
print(f"  cd {OUT_DIR}")
print(f"  obabel receptor_NSP15_9.pdb -O receptor_NSP15_9.pdbqt -p 7.4 --partialcharge gasteiger")
print("Next: Script 09_9 (publication figures)")
