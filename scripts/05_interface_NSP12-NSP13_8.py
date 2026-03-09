#!/usr/bin/env python3
"""
Script 05_interface_NSP12-NSP13_8.py

NSP12-NSP13 interface analysis from 3 crystal structures.
AF3 excluded from interface analysis (iptm=0.20, PAE=24.95 A).

Structures:
  6XEZ : A=NSP12, E=NSP13 (3.50 A) — primary
  7CXM : A=NSP12, E=NSP13 (3.20 A) — secondary
  7RDY : A=NSP12, E=NSP13 (3.10 A) — tertiary

Key finding from Script 04_8:
  NSP12 interface: C-terminal tail residues 900-904
  NSP13 interface: N-terminal region residues 90-96
  Smallest interface in project — 4-5 residues per chain

Output: 02-validation/NSP12-NSP13/interface_analysis_8.json
"""

import os, json
import numpy as np
from Bio import PDB, Align
from collections import Counter

# ── Paths ──────────────────────────────────────────────────────────────────
BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = f"{BASE}/00-reference/pdb_structures"
VAL_DIR = f"{BASE}/02-validation/NSP12-NSP13"
os.makedirs(VAL_DIR, exist_ok=True)

STRUCTURES = {
    "6XEZ" : (f"{PDB_DIR}/6XEZ.pdb", "A", "E", "3.50 A"),
    "7CXM" : (f"{PDB_DIR}/7CXM.pdb", "A", "E", "3.20 A"),
    "7RDY" : (f"{PDB_DIR}/7RDY.pdb", "A", "E", "3.10 A"),
}

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}
HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO","TYR"}
POSITIVE    = {"ARG","LYS","HIS"}
NEGATIVE    = {"ASP","GLU"}

CUTOFF = 5.0

# ── Contact classification ─────────────────────────────────────────────────
def classify_contact(res1, res2, dist):
    r1, r2 = res1.get_resname(), res2.get_resname()
    if ((r1 in POSITIVE and r2 in NEGATIVE) or
        (r1 in NEGATIVE and r2 in POSITIVE)) and dist < 4.5:
        return "salt_bridge"
    if dist < 3.5 and not (r1 in HYDROPHOBIC and r2 in HYDROPHOBIC):
        return "h_bond"
    if r1 in HYDROPHOBIC and r2 in HYDROPHOBIC and dist < 5.0:
        return "hydrophobic"
    return "vdw"

# ── Interface analysis ─────────────────────────────────────────────────────
def analyze_interface(chain_ref, chain_partner, label):
    print(f"\n  Analyzing {label} ...")

    atoms_partner = [a for r in chain_partner
                     if PDB.is_aa(r) for a in r.get_atoms()]
    ns = PDB.NeighborSearch(atoms_partner)

    contacts      = []
    salt_bridges  = []
    seen_sb       = set()
    iface_ref     = set()
    contact_count = Counter()

    for r1 in chain_ref:
        if not PDB.is_aa(r1):
            continue
        for a1 in r1.get_atoms():
            hits = ns.search(a1.get_vector().get_array(), CUTOFF, "A")
            for a2 in hits:
                r2 = a2.get_parent()
                if not PDB.is_aa(r2):
                    continue
                dist = float(np.linalg.norm(
                    a1.get_vector().get_array() -
                    a2.get_vector().get_array()))
                ctype = classify_contact(r1, r2, dist)
                rn1   = r1.get_id()[1]
                rn2   = r2.get_id()[1]
                iface_ref.add(rn1)
                contact_count[rn1] += 1

                contacts.append({
                    "res1_num" : rn1,
                    "res1_name": r1.get_resname(),
                    "res2_num" : rn2,
                    "res2_name": r2.get_resname(),
                    "dist_ang" : round(dist, 2),
                    "type"     : ctype
                })

                pair = (min(rn1,rn2), max(rn1,rn2))
                if ctype == "salt_bridge" and pair not in seen_sb:
                    seen_sb.add(pair)
                    salt_bridges.append({
                        "res1" : f"{r1.get_resname()}{rn1}",
                        "res2" : f"{r2.get_resname()}{rn2}",
                        "dist" : round(dist, 2)
                    })

    n_sb = len(salt_bridges)
    n_hb = len([c for c in contacts if c["type"]=="h_bond"])
    n_hy = len([c for c in contacts if c["type"]=="hydrophobic"])
    n_tot= len(contacts)

    hotspots = [{"res_num": rn, "res_name": next(
                    (r.get_resname() for r in chain_ref
                     if PDB.is_aa(r) and r.get_id()[1]==rn), "UNK"),
                 "contacts": cnt}
                for rn, cnt in sorted(
                    contact_count.items(), key=lambda x: -x[1])]

    print(f"    SB={n_sb}  HB={n_hb}  HY={n_hy}  total={n_tot}")
    print(f"    Interface residues: {sorted(iface_ref)}")
    if salt_bridges:
        print(f"    Salt bridges:")
        for sb in salt_bridges:
            print(f"      {sb['res1']} -- {sb['res2']} : {sb['dist']} A")
    else:
        print(f"    No salt bridges detected")
    print(f"    Hotspots by contact count:")
    for h in hotspots:
        print(f"      {h['res_name']}{h['res_num']:5d}  contacts={h['contacts']}")

    return {
        "label"             : label,
        "n_salt_bridges"    : n_sb,
        "n_h_bonds"         : n_hb,
        "n_hydrophobic"     : n_hy,
        "n_total"           : n_tot,
        "salt_bridges"      : salt_bridges,
        "interface_residues": sorted(iface_ref),
        "hotspots"          : hotspots,
    }

# ── Run analysis ───────────────────────────────────────────────────────────
print("=" * 60)
print("INTERFACE ANALYSIS — NSP12-NSP13")
print("=" * 60)
print("Note: AF3 excluded — iptm=0.20, PAE=24.95 A")
print("      3 crystal structures only")

parser  = PDB.PDBParser(QUIET=True)
results = {}

for pdb_id, (pdb_path, c12, c13, res) in STRUCTURES.items():
    struct = parser.get_structure(pdb_id, pdb_path)
    label  = f"{pdb_id} NSP12({c12}) vs NSP13({c13}) [{res}]"
    res_   = analyze_interface(struct[0][c12], struct[0][c13], label)
    results[pdb_id] = res_

# ── Cross-structure salt bridge summary ───────────────────────────────────
print("\n" + "=" * 60)
print("SALT BRIDGE SUMMARY (all structures)")
print("=" * 60)

all_sbs = []
for pdb_id, r in results.items():
    for sb in r["salt_bridges"]:
        sb["structure"] = pdb_id
        all_sbs.append(sb)

if all_sbs:
    for sb in all_sbs:
        print(f"  [{sb['structure']}] {sb['res1']} -- {sb['res2']} : {sb['dist']} A")
else:
    print("  No salt bridges detected in any crystal structure")

# ── Consensus hotspots ─────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("CONSENSUS HOTSPOTS (NSP12, in >= 2 structures)")
print("=" * 60)

# Union of all NSP12 interface residues across structures
all_iface_nsp12 = {}
for pdb_id, r in results.items():
    for rn in r["interface_residues"]:
        all_iface_nsp12[rn] = all_iface_nsp12.get(rn, 0) + 1

print(f"\n  NSP12 residue appearance across 3 structures:")
for rn, count in sorted(all_iface_nsp12.items()):
    # Get residue name from 6XEZ
    struct_6XEZ = parser.get_structure("6XEZ", f"{PDB_DIR}/6XEZ.pdb")
    resname = "UNK"
    for r in struct_6XEZ[0]["A"]:
        if PDB.is_aa(r) and r.get_id()[1] == rn:
            resname = r.get_resname()
            break
    marker = " ✅ ALL 3" if count == 3 else f" ({count}/3)"
    print(f"    {resname}{rn:5d}  {marker}")

consensus = [(rn, count) for rn, count in all_iface_nsp12.items() if count >= 2]
print(f"\n  Consensus NSP12 hotspots (>=2 structures): {len(consensus)}")

# NSP13 interface
all_iface_nsp13 = {}
for pdb_id, r in results.items():
    # Get NSP13 interface from partner side
    struct = parser.get_structure(pdb_id, STRUCTURES[pdb_id][0])
    c12    = STRUCTURES[pdb_id][1]
    c13    = STRUCTURES[pdb_id][2]
    _, iface_nsp13 = set(), set()
    atoms_nsp12 = [a for r_ in struct[0][c12]
                   if PDB.is_aa(r_) for a in r_.get_atoms()]
    ns = PDB.NeighborSearch(atoms_nsp12)
    for r_ in struct[0][c13]:
        if not PDB.is_aa(r_):
            continue
        for a in r_.get_atoms():
            hits = ns.search(a.get_vector().get_array(), CUTOFF, "A")
            if hits:
                iface_nsp13.add(r_.get_id()[1])
    for rn in iface_nsp13:
        all_iface_nsp13[rn] = all_iface_nsp13.get(rn, 0) + 1

print(f"\n  NSP13 residue appearance across 3 structures:")
struct_6XEZ = parser.get_structure("6XEZ", f"{PDB_DIR}/6XEZ.pdb")
for rn, count in sorted(all_iface_nsp13.items()):
    resname = "UNK"
    for r in struct_6XEZ[0]["E"]:
        if PDB.is_aa(r) and r.get_id()[1] == rn:
            resname = r.get_resname()
            break
    marker = " ✅ ALL 3" if count == 3 else f" ({count}/3)"
    print(f"    {resname}{rn:5d}  {marker}")

# ── Interface character summary ────────────────────────────────────────────
print("\n" + "=" * 60)
print("INTERFACE CHARACTER SUMMARY")
print("=" * 60)

print(f"\n  {'Structure':>8}  {'SB':>4}  {'HB':>4}  {'HY':>4}  {'Total':>6}")
print("  " + "-"*35)
for pdb_id, r in results.items():
    print(f"  {pdb_id:>8}  {r['n_salt_bridges']:>4}  "
          f"{r['n_h_bonds']:>4}  {r['n_hydrophobic']:>4}  {r['n_total']:>6}")

print(f"\n  SCIENTIFIC NOTE:")
print(f"  This is the smallest interface in the project.")
print(f"  NSP12 C-terminal tail (900-904) contacts NSP13 N-terminal (90-96)")
print(f"  Consistent with transient regulatory interaction")

# ── Save ───────────────────────────────────────────────────────────────────
nsp12_hotspots_all = sorted([
    {"res_num": rn, "n_structures": cnt}
    for rn, cnt in all_iface_nsp12.items()
], key=lambda x: -x["n_structures"])

nsp13_hotspots_all = sorted([
    {"res_num": rn, "n_structures": cnt}
    for rn, cnt in all_iface_nsp13.items()
], key=lambda x: -x["n_structures"])

output = {
    "complex"              : "NSP12-NSP13",
    "script"               : "05_interface_NSP12-NSP13_8.py",
    "af3_note"             : "Excluded — iptm=0.20, PAE=24.95 A",
    "structures_analyzed"  : list(STRUCTURES.keys()),
    "interface_results"    : results,
    "salt_bridges_all"     : all_sbs,
    "nsp12_hotspots_all"   : nsp12_hotspots_all,
    "nsp13_hotspots_all"   : nsp13_hotspots_all,
    "consensus_nsp12"      : [x for x in nsp12_hotspots_all if x["n_structures"]>=2],
    "consensus_nsp13"      : [x for x in nsp13_hotspots_all if x["n_structures"]>=2],
}

out_path = f"{VAL_DIR}/interface_analysis_8.json"
with open(out_path, "w") as f:
    json.dump(output, f, indent=2)
print(f"\n  Saved: {out_path}")
