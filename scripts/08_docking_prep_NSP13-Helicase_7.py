#!/usr/bin/env python3
"""
Script 08_docking_prep_NSP13-Helicase_7.py

Prepares receptor and docking configuration files for
NSP13 homodimer virtual screening.

Primary receptor: 7NIO chains A+E (full dimer)
Docking box: centered on LYS414-ASP580/ASP583 dual SB region

Output: 03-virtual-screening/NSP13-Helicase_7/
  receptor_NSP13-Helicase_7.pdb
  vina_config_NSP13-Helicase_7.txt
  virtualflow_config_NSP13-Helicase_7.json
"""

import os, json
from Bio import PDB

# ── Paths ──────────────────────────────────────────────────────────────────
BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = f"{BASE}/00-reference/pdb_structures"
VS_DIR  = f"{BASE}/03-virtual-screening/NSP13-Helicase_7"
VAL_DIR = f"{BASE}/02-validation/NSP13-Helicase"
os.makedirs(VS_DIR, exist_ok=True)

PDB_7NIO = f"{PDB_DIR}/7NIO.pdb"

# ── Docking box (from Script 07_7) ─────────────────────────────────────────
DOCKING_BOX = {
    "center_x" : -30.151,
    "center_y" :  14.648,
    "center_z" :  -9.240,
    "size_x"   :  37.237,
    "size_y"   :  31.585,
    "size_z"   :  33.862,
}

# ── Hotspot residues for verification ─────────────────────────────────────
ALL_HOTSPOTS     = [116,413,414,415,477,478,479,480,481,482,
                    551,552,553,579,580,583,584]
PRIMARY_HOTSPOTS = {414, 477, 480, 482, 579, 580, 583, 584}

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# ── Step 1: Extract receptor — chains A+E, ATOM only ──────────────────────
print("=" * 60)
print("STEP 1: Extract receptor from 7NIO (chains A+E)")
print("=" * 60)

parser = PDB.PDBParser(QUIET=True)
struct  = parser.get_structure("7NIO", PDB_7NIO)

class DimerSelector(PDB.Select):
    def accept_chain(self, chain):
        return chain.get_id() in ["A", "E"]
    def accept_residue(self, residue):
        return residue.get_id()[0] == " "  # ATOM only — no waters/ligands
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"

receptor_path = f"{VS_DIR}/receptor_NSP13-Helicase_7.pdb"
io = PDB.PDBIO()
io.set_structure(struct)
io.save(receptor_path, DimerSelector())
print(f"  Saved receptor: {receptor_path}")

# ── Step 2: Count and verify receptor ─────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Receptor verification")
print("=" * 60)

struct_rec = parser.get_structure("receptor", receptor_path)

chain_res_counts = {}
total_atoms = 0
for chain in struct_rec[0]:
    aa_res = [r for r in chain if PDB.is_aa(r)]
    chain_res_counts[chain.get_id()] = len(aa_res)
    for r in chain:
        total_atoms += len(list(r.get_atoms()))

total_res = sum(chain_res_counts.values())
print(f"  Chains retained : {list(chain_res_counts.keys())}")
for cid, count in chain_res_counts.items():
    print(f"    Chain {cid}: {count} residues")
print(f"  Total residues  : {total_res}")
print(f"  Total atoms     : {total_atoms}")

# ── Step 3: Hotspot verification ───────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: Hotspot residue verification (chain A)")
print("=" * 60)

found_all      = []
found_primary  = []
missing        = []

chain_A = struct_rec[0]["A"]
for rn in ALL_HOTSPOTS:
    found = False
    for res in chain_A:
        if PDB.is_aa(res) and res.get_id()[1] == rn:
            resname = res.get_resname()
            found_all.append((rn, resname))
            if rn in PRIMARY_HOTSPOTS:
                found_primary.append((rn, resname))
            found = True
            break
    if not found:
        missing.append(rn)

print(f"  All hotspots    : {len(found_all)}/{len(ALL_HOTSPOTS)} present")
print(f"  Primary hotspots: {len(found_primary)}/{len(PRIMARY_HOTSPOTS)} present")

print(f"\n  Hotspot residues confirmed:")
for rn, resname in found_all:
    prim = " ★" if rn in PRIMARY_HOTSPOTS else ""
    expected = dict(ALL_HOTSPOTS_EXPECTED := {
        116:"ASN",413:"THR",414:"LYS",415:"GLY",
        477:"LYS",478:"GLY",479:"VAL",480:"ILE",
        481:"THR",482:"HIS",551:"GLU",552:"THR",
        553:"ALA",579:"ARG",580:"ASP",583:"ASP",584:"LYS"
    }).get(rn, "?")
    match = "✅" if resname == expected else f"⚠️  expected {expected}"
    print(f"    {resname}{rn:4d}{prim}  {match}")

if missing:
    print(f"\n  Missing hotspot positions: {missing}")
    print(f"  (Check if PDB chain A starts at residue > min hotspot)")
else:
    print(f"\n  No missing hotspots ✅")

# ── Step 4: Salt bridge distance check ────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Salt bridge geometry verification")
print("=" * 60)

import numpy as np

def get_res(chain, res_num):
    for res in chain:
        if PDB.is_aa(res) and res.get_id()[1] == res_num:
            return res
    return None

def min_dist(res1, res2):
    min_d = 999.0
    for a1 in res1.get_atoms():
        for a2 in res2.get_atoms():
            d = float(np.linalg.norm(
                a1.get_vector().get_array() - a2.get_vector().get_array()))
            if d < min_d:
                min_d = d
    return min_d

chain_A = struct_rec[0]["A"]
chain_E = struct_rec[0]["E"]

sb_pairs = [
    (chain_A, 414, chain_E, 580, "LYS414(A)-ASP580(E)"),
    (chain_A, 414, chain_E, 583, "LYS414(A)-ASP583(E)"),
]

print(f"\n  Salt bridge distances (receptor):")
for c1, rn1, c2, rn2, label in sb_pairs:
    r1 = get_res(c1, rn1)
    r2 = get_res(c2, rn2)
    if r1 and r2:
        d = min_dist(r1, r2)
        flag = "✅" if d < 5.0 else "⚠️  check"
        print(f"    {label}: {d:.2f} A  {flag}")
    else:
        print(f"    {label}: residue not found")

# ── Step 5: Write Vina config ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 5: Write AutoDock Vina config")
print("=" * 60)

vina_config = f"""# AutoDock Vina config — NSP13-Helicase_7
# Target: NSP13 homodimer interface (LYS414-ASP580/ASP583 dual salt bridge)
# Conservation: SARS-CoV-1/2 selective (ASP580 D->A in MERS, D->T in HCoV)
# Primary pharmacophore: LYS414 (dual SB donor), ILE480/HIS482 (top contacts)
# Receptor: 7NIO chains A+E, waters/ligands stripped

receptor = receptor_NSP13-Helicase_7.pdbqt

center_x = {DOCKING_BOX['center_x']}
center_y = {DOCKING_BOX['center_y']}
center_z = {DOCKING_BOX['center_z']}

size_x = {DOCKING_BOX['size_x']}
size_y = {DOCKING_BOX['size_y']}
size_z = {DOCKING_BOX['size_z']}

exhaustiveness = 16
num_modes      = 9
energy_range   = 3
"""

vina_path = f"{VS_DIR}/vina_config_NSP13-Helicase_7.txt"
with open(vina_path, "w") as f:
    f.write(vina_config)
print(f"  Saved: {vina_path}")
print(f"  Box center : ({DOCKING_BOX['center_x']}, {DOCKING_BOX['center_y']}, {DOCKING_BOX['center_z']})")
print(f"  Box size   : ({DOCKING_BOX['size_x']} x {DOCKING_BOX['size_y']} x {DOCKING_BOX['size_z']}) A")

# ── Step 6: Write VirtualFlow config ──────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 6: Write VirtualFlow config")
print("=" * 60)

vf_config = {
    "target"       : "NSP13-Helicase_7",
    "complex"      : "NSP13 homodimer (SARS-CoV-1/2 selective)",
    "receptor"     : "receptor_NSP13-Helicase_7.pdbqt",
    "interface"    : "LYS414-ASP580/ASP583 dual salt bridge + ILE480/HIS482 hydrophobic core",
    "conservation" : "SARS-CoV-1/2 selective — ASP580 D->A(MERS) D->T(229E/NL63)",
    "primary_pharmacophore": "LYS414 (dual SB donor, cons=0.689)",
    "docking_box"  : DOCKING_BOX,
    "hotspot_residues" : {
        "all"     : ALL_HOTSPOTS,
        "primary" : sorted(PRIMARY_HOTSPOTS)
    },
    "salt_bridges" : [
        {"pair": "LYS414(A)-ASP580(E)", "dist_7NIO": 4.49, "present_in": "7NIO"},
        {"pair": "LYS414(A)-ASP583(E)", "dist_7NIO": 4.47, "present_in": "7NIO"},
    ],
    "top_contacts" : [
        {"residue": "ILE480", "contacts": 60, "conservation": 0.582},
        {"residue": "HIS482", "contacts": 60, "conservation": 0.582},
        {"residue": "ASP580", "contacts": 49, "conservation": 0.344},
        {"residue": "LYS414", "contacts": 43, "conservation": 0.689},
    ],
    "vina_exhaustiveness": 16,
    "vina_num_modes"     : 9,
}

vf_path = f"{VS_DIR}/virtualflow_config_NSP13-Helicase_7.json"
with open(vf_path, "w") as f:
    json.dump(vf_config, f, indent=2)
print(f"  Saved: {vf_path}")

# ── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY — NSP13-Helicase_7 DOCKING PREP COMPLETE")
print("=" * 60)
print(f"\n  Receptor    : {total_res} residues, {total_atoms} atoms")
print(f"  Chains      : A (NSP13 chain 1) + E (NSP13 chain 2)")
print(f"  Hotspots    : {len(found_all)}/{len(ALL_HOTSPOTS)} present ✅")
print(f"  Primary SB  : LYS414-ASP580 + LYS414-ASP583 (dual)")
print(f"  Selectivity : SARS-CoV-1/2 selective")
print(f"\n  Output files:")
print(f"    {receptor_path}")
print(f"    {vina_path}")
print(f"    {vf_path}")
print(f"\n  NSP13-Helicase pipeline: steps 04-08 complete ✅")
print(f"  Next: Script 09_7 (publication figures)")
