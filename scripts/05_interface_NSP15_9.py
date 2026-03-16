"""
Script 05_9 — Interface Analysis NSP15 Homodimer
Complex  : NSP15-NSP15 (NendoU dimer interface)
Suffix   : _9
Reference: 6VWW_NSP15-true-dimer.pdb (chains A+B = biological assembly A+A-2)
Structures: 6VWW (apo), 6W01 (apo), 6WLC (UMP), 6WXC (tipiracil)
Note     : AF3 is monomer — used for pLDDT mapping only, not interface contacts
"""

import json
import os
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio import pairwise2
from Bio.PDB.Polypeptide import protein_letters_3to1

OUT_DIR  = "02-validation/NSP15"
OUT_JSON = f"{OUT_DIR}/interface_analysis_9.json"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Structure registry ─────────────────────────────────────────
# For NSP15 all structures are homodimers — extract true dimer
# from biological assembly for each PDB
STRUCTURES = {
    "6VWW": {
        "path"  : "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb",
        "label" : "SARS-CoV-2 apo 1.90A PRIMARY",
        "chains": ("A","B"),
    },
    "6W01": {
        "path"  : "00-reference/known_interfaces/NSP15/6W01_NSP15-true-dimer.pdb",
        "label" : "SARS-CoV-2 apo 2.20A",
        "chains": ("A","B"),
    },
    "6WLC": {
        "path"  : "00-reference/known_interfaces/NSP15/6WLC_NSP15-true-dimer-UMP.pdb",
        "label" : "SARS-CoV-2 +UMP 1.70A",
        "chains": ("A","B"),
    },
    "6WXC": {
        "path"  : "00-reference/known_interfaces/NSP15/6WXC_NSP15-true-dimer-tipiracil.pdb",
        "label" : "SARS-CoV-2 +tipiracil 1.80A",
        "chains": ("A","B"),
    },
}

CUTOFF_CONTACT = 5.0
CUTOFF_SB      = 4.5   # salt bridge
CUTOFF_HB      = 3.5   # hydrogen bond
CUTOFF_HY      = 5.0   # hydrophobic

# Charged residues for salt bridge detection
POS_CHARGED = {'ARG','LYS','HIS'}
NEG_CHARGED = {'ASP','GLU'}
HYDROPHOBIC = {'ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO','TYR'}

print("=" * 65)
print("Script 05_9 — Interface Analysis NSP15 Homodimer")
print("=" * 65)

# ── Helper functions ───────────────────────────────────────────
def get_seq(chain):
    res = [r for r in chain.get_residues() if r.get_id()[0] == ' ']
    return ''.join(protein_letters_3to1.get(r.get_resname(),'X') for r in res)

def classify_contact(rA, rB):
    """Classify contact type between two residues."""
    nA = rA.get_resname()
    nB = rB.get_resname()
    # Salt bridge
    if ((nA in POS_CHARGED and nB in NEG_CHARGED) or
        (nA in NEG_CHARGED and nB in POS_CHARGED)):
        for aA in rA.get_atoms():
            for aB in rB.get_atoms():
                if aA - aB < CUTOFF_SB:
                    return 'SB'
    # Hydrogen bond (N or O atoms within 3.5A)
    for aA in rA.get_atoms():
        if aA.element in ('N','O'):
            for aB in rB.get_atoms():
                if aB.element in ('N','O'):
                    if aA - aB < CUTOFF_HB:
                        return 'HB'
    # Hydrophobic
    if nA in HYDROPHOBIC and nB in HYDROPHOBIC:
        return 'HY'
    return 'HY'  # default to hydrophobic for any other contact

def analyze_interface(pdb_id, path, chains, label):
    """Full interface analysis for one structure."""
    print(f"\n  {pdb_id}: {label}")

    if not os.path.exists(path):
        print(f"    WARNING: {path} not found — skipping")
        return None

    struct = PDBParser(QUIET=True).get_structure(pdb_id, path)
    try:
        cA = struct[0][chains[0]]
        cB = struct[0][chains[1]]
    except KeyError as e:
        print(f"    WARNING: chain {e} not found — skipping")
        return None

    resA = [r for r in cA.get_residues() if r.get_id()[0] == ' ']
    resB = [r for r in cB.get_residues() if r.get_id()[0] == ' ']

    contacts  = []
    iface_A   = set()
    iface_B   = set()
    sb_list   = []
    hb_count  = 0
    hy_count  = 0
    sb_count  = 0

    for rA in resA:
        for rB in resB:
            in_contact = False
            for aA in rA.get_atoms():
                for aB in rB.get_atoms():
                    if aA - aB < CUTOFF_CONTACT:
                        in_contact = True
                        break
                if in_contact:
                    break
            if in_contact:
                ctype = classify_contact(rA, rB)
                posA  = rA.get_id()[1]
                posB  = rB.get_id()[1]
                iface_A.add(posA)
                iface_B.add(posB)
                contacts.append((posA, rA.get_resname(),
                                  posB, rB.get_resname(), ctype))
                if ctype == 'SB':
                    sb_count += 1
                    # Measure distance
                    min_d = min(aA - aB
                                for aA in rA.get_atoms()
                                for aB in rB.get_atoms())
                    sb_list.append({
                        "donor"    : f"{rA.get_resname()}{posA}(A)",
                        "acceptor" : f"{rB.get_resname()}{posB}(B)",
                        "distance" : round(min_d, 2)
                    })
                elif ctype == 'HB':
                    hb_count += 1
                else:
                    hy_count += 1

    total = sb_count + hb_count + hy_count
    print(f"    SB={sb_count}  HB={hb_count}  HY={hy_count}  total={total}")
    print(f"    Interface A: {len(iface_A)} residues")
    print(f"    Interface B: {len(iface_B)} residues")

    if sb_list:
        print(f"    Salt bridges:")
        for sb in sb_list:
            print(f"      {sb['donor']} -- {sb['acceptor']}: {sb['distance']} A")

    return {
        "pdb_id"      : pdb_id,
        "label"       : label,
        "sb_count"    : sb_count,
        "hb_count"    : hb_count,
        "hy_count"    : hy_count,
        "total"       : total,
        "n_iface_A"   : len(iface_A),
        "n_iface_B"   : len(iface_B),
        "iface_A"     : sorted(iface_A),
        "iface_B"     : sorted(iface_B),
        "salt_bridges": sb_list,
        "contacts"    : contacts,
    }

# ── Run analysis on all structures ────────────────────────────
print("\nAnalyzing interface across all structures:")
results = {}
for pdb_id, info in STRUCTURES.items():
    r = analyze_interface(pdb_id, info["path"],
                          info["chains"], info["label"])
    if r:
        results[pdb_id] = r

# ── Consensus hotspots (present in >= 2 structures) ───────────
print("\n" + "=" * 65)
print("Consensus hotspot analysis:")

from collections import Counter
all_A = Counter()
all_B = Counter()
for r in results.values():
    for pos in r["iface_A"]:
        all_A[pos] += 1
    for pos in r["iface_B"]:
        all_B[pos] += 1

n_structs   = len(results)
thresh      = max(2, n_structs // 2)
consensus_A = sorted([p for p,c in all_A.items() if c >= thresh])
consensus_B = sorted([p for p,c in all_B.items() if c >= thresh])

print(f"\n  Structures analyzed : {n_structs}")
print(f"  Consensus threshold : >= {thresh} structures")
print(f"  Consensus A ({len(consensus_A)} res): {consensus_A}")
print(f"  Consensus B ({len(consensus_B)} res): {consensus_B}")

# ── Salt bridge summary ────────────────────────────────────────
print("\nSalt bridge inventory:")
all_sbs = {}
for pdb_id, r in results.items():
    for sb in r["salt_bridges"]:
        key = "{}-{}".format(sb["donor"], sb["acceptor"])
        if key not in all_sbs:
            all_sbs[key] = {}
        all_sbs[key][pdb_id] = sb["distance"]

if all_sbs:
    for pair, structs in all_sbs.items():
        presence = list(structs.keys())
        distances = ["{} {:.2f}A".format(p,d) for p,d in structs.items()]
        marker = " ★ PRIMARY" if len(presence) == n_structs else ""
        print(f"  {pair}: {', '.join(distances)}{marker}")
else:
    print("  No salt bridges detected — hydrophobic/H-bond dominated interface")

# ── Contact summary table ──────────────────────────────────────
print("\nContact summary:")
print("{:<8} {:<6} {:<6} {:<6} {:<8} {:<12} {}".format(
    "PDB","SB","HB","HY","Total","Iface_A_n","Iface_B_n"))
print("-" * 55)
for pdb_id, r in results.items():
    print("{:<8} {:<6} {:<6} {:<6} {:<8} {:<12} {}".format(
        pdb_id, r["sb_count"], r["hb_count"], r["hy_count"],
        r["total"], r["n_iface_A"], r["n_iface_B"]))

# ── Interface character assessment ────────────────────────────
total_sb = sum(r["sb_count"] for r in results.values())
total_hb = sum(r["hb_count"] for r in results.values())
total_hy = sum(r["hy_count"] for r in results.values())
dominant = max([("SB",total_sb),("HB",total_hb),("HY",total_hy)],
               key=lambda x: x[1])[0]

print(f"\nInterface character: {dominant}-dominated")
print(f"  Total across all structures: SB={total_sb} HB={total_hb} HY={total_hy}")

# ── Save results ───────────────────────────────────────────────
output = {
    "complex"            : "NSP15-homodimer",
    "suffix"             : "_9",
    "structures_analyzed": list(results.keys()),
    "n_structures"       : n_structs,
    "consensus_threshold": thresh,
    "consensus_A"        : consensus_A,
    "consensus_B"        : consensus_B,
    "salt_bridge_summary": all_sbs,
    "interface_character": dominant,
    "per_structure"      : {
        k: {
            "sb": v["sb_count"],
            "hb": v["hb_count"],
            "hy": v["hy_count"],
            "total": v["total"],
            "n_iface_A": v["n_iface_A"],
            "n_iface_B": v["n_iface_B"],
            "iface_A": v["iface_A"],
            "iface_B": v["iface_B"],
            "salt_bridges": v["salt_bridges"],
        }
        for k, v in results.items()
    },
    "note": (
        "6W01/6WLC/6WXC use asymmetric unit chains A+B. "
        "6VWW uses biological assembly A+A-2 (true dimer, 75 contacts). "
        "AF3 monomer used for pLDDT only — not for interface contacts."
    )
}

with open(OUT_JSON, "w") as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nSaved: {OUT_JSON}")
print("Next : Script 06_9 (conservation analysis)")
