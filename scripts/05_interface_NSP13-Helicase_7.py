#!/usr/bin/env python3
"""
Script 05_interface_NSP13-Helicase_7.py

NSP13 homodimer interface analysis.
AF3 is a monomer — used only for pLDDT mapping at hotspot positions.

Structures:
  7NIO  : chains A (ref) vs E (partner) — primary  (2.80 A)
  6XEZ  : chains E (ref) vs F (partner) — secondary (3.50 A)
  AF3   : chain A monomer — pLDDT only

Output: 02-validation/NSP13-Helicase/interface_analysis_7.json
"""

import os, json, math
import numpy as np
from Bio import PDB, Align

# ── Paths ──────────────────────────────────────────────────────────────────
BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = f"{BASE}/00-reference/pdb_structures"
AF3_DIR = f"{BASE}/01-alphafold3/NSP13-Helicase"
VAL_DIR = f"{BASE}/02-validation/NSP13-Helicase"
os.makedirs(VAL_DIR, exist_ok=True)

PDB_7NIO = f"{PDB_DIR}/7NIO.pdb"
PDB_6XEZ = f"{PDB_DIR}/6XEZ.pdb"
AF3_PDB  = f"{AF3_DIR}/NSP13_Helicase_best_model.pdb"

CUTOFF = 5.0

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO","TYR"}
POSITIVE    = {"ARG","LYS","HIS"}
NEGATIVE    = {"ASP","GLU"}

# ── Parse structures ───────────────────────────────────────────────────────
parser = PDB.PDBParser(QUIET=True)
struct_7NIO = parser.get_structure("7NIO", PDB_7NIO)
struct_6XEZ = parser.get_structure("6XEZ", PDB_6XEZ)
struct_AF3  = parser.get_structure("AF3",  AF3_PDB)

# ── Contact analysis ───────────────────────────────────────────────────────
def classify_contact(res1, res2, dist):
    """Classify a contact as salt bridge, H-bond, hydrophobic, or VDW."""
    r1, r2 = res1.get_resname(), res2.get_resname()
    # Salt bridge
    if (r1 in POSITIVE and r2 in NEGATIVE) or (r1 in NEGATIVE and r2 in POSITIVE):
        if dist < 4.5:
            return "salt_bridge"
    # H-bond donors/acceptors (N, O atoms within 3.5 A)
    donors_acceptors = {"N","O","ND1","ND2","NE","NE1","NE2","NH1","NH2",
                        "NZ","OD1","OD2","OE1","OE2","OG","OG1","OH","SD"}
    # Check if close atoms include N or O
    if dist < 3.5 and (r1 not in HYDROPHOBIC or r2 not in HYDROPHOBIC):
        return "h_bond"
    # Hydrophobic
    if r1 in HYDROPHOBIC and r2 in HYDROPHOBIC and dist < 5.0:
        return "hydrophobic"
    return "vdw"

def analyze_interface(chain_ref, chain_partner, label, cutoff=5.0):
    """Full interface analysis between two chains."""
    print(f"\n  Analyzing {label} ...")

    atoms_partner = [a for r in chain_partner if PDB.is_aa(r) for a in r.get_atoms()]
    ns = PDB.NeighborSearch(atoms_partner)

    contacts      = []
    salt_bridges  = []
    h_bonds       = []
    hydrophobics  = []
    seen_sb_pairs = set()
    iface_ref     = set()

    for r1 in chain_ref:
        if not PDB.is_aa(r1):
            continue
        for a1 in r1.get_atoms():
            hits = ns.search(a1.get_vector().get_array(), cutoff, "A")
            for a2 in hits:
                r2   = a2.get_parent()
                if not PDB.is_aa(r2):
                    continue
                dist = float(np.linalg.norm(
                    a1.get_vector().get_array() - a2.get_vector().get_array()))
                contact_type = classify_contact(r1, r2, dist)
                rn1 = r1.get_id()[1]
                rn2 = r2.get_id()[1]
                iface_ref.add(rn1)

                contacts.append({
                    "res1_num"  : rn1,
                    "res1_name" : r1.get_resname(),
                    "atom1"     : a1.get_name(),
                    "res2_num"  : rn2,
                    "res2_name" : r2.get_resname(),
                    "atom2"     : a2.get_name(),
                    "dist_ang"  : round(dist, 2),
                    "type"      : contact_type
                })

                pair_key = (min(rn1,rn2), max(rn1,rn2))

                if contact_type == "salt_bridge" and pair_key not in seen_sb_pairs:
                    seen_sb_pairs.add(pair_key)
                    salt_bridges.append({
                        "res1": f"{r1.get_resname()}{rn1}",
                        "res2": f"{r2.get_resname()}{rn2}",
                        "dist": round(dist, 2)
                    })
                elif contact_type == "h_bond":
                    h_bonds.append((rn1, rn2))
                elif contact_type == "hydrophobic":
                    hydrophobics.append((rn1, rn2))

    # Hotspot scoring
    from collections import Counter
    res_contact_count = Counter(c["res1_num"] for c in contacts)

    hotspots = []
    for rn, count in sorted(res_contact_count.items(), key=lambda x: -x[1]):
        # find residue name
        for r in chain_ref:
            if PDB.is_aa(r) and r.get_id()[1] == rn:
                resname = r.get_resname()
                break
        else:
            resname = "UNK"
        hotspots.append({
            "res_num"  : rn,
            "res_name" : resname,
            "contacts" : count
        })

    n_sb  = len(salt_bridges)
    n_hb  = len(set(h_bonds))
    n_hy  = len(set(hydrophobics))
    n_tot = len(contacts)

    print(f"    SB={n_sb}  HB={n_hb}  HY={n_hy}  total={n_tot}")
    print(f"    Interface residues: {sorted(iface_ref)}")

    if salt_bridges:
        print(f"    Salt bridges:")
        for sb in salt_bridges:
            print(f"      {sb['res1']} -- {sb['res2']} : {sb['dist']} A")

    print(f"    Top 10 hotspots:")
    for h in hotspots[:10]:
        print(f"      {h['res_name']}{h['res_num']:4d}  contacts={h['contacts']}")

    return {
        "label"            : label,
        "n_salt_bridges"   : n_sb,
        "n_h_bonds"        : n_hb,
        "n_hydrophobic"    : n_hy,
        "n_total"          : n_tot,
        "salt_bridges"     : salt_bridges,
        "interface_residues": sorted(iface_ref),
        "hotspots"         : hotspots
    }

# ── Run analysis on both crystal structures ────────────────────────────────
print("=" * 60)
print("INTERFACE ANALYSIS — NSP13 HOMODIMER")
print("=" * 60)

res_7NIO = analyze_interface(
    struct_7NIO[0]["A"],
    struct_7NIO[0]["E"],
    "7NIO_A-vs-E"
)

res_6XEZ = analyze_interface(
    struct_6XEZ[0]["E"],
    struct_6XEZ[0]["F"],
    "6XEZ_E-vs-F"
)

# ── Consensus hotspots ─────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("CONSENSUS HOTSPOTS (residues in BOTH crystal structures)")
print("=" * 60)

# Build sequence-alignment mapping between 7NIO chain A and 6XEZ chain E
# so we can compare residue numbers across structures

def get_seq_reslist(chain):
    return [(r.get_id()[1], r.get_resname())
            for r in chain if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1]

reslist_7NIO_A = get_seq_reslist(struct_7NIO[0]["A"])
reslist_6XEZ_E = get_seq_reslist(struct_6XEZ[0]["E"])

seq_7NIO_A = "".join(AA_MAP_3TO1[r[1]] for r in reslist_7NIO_A)
seq_6XEZ_E = "".join(AA_MAP_3TO1[r[1]] for r in reslist_6XEZ_E)

aligner = Align.PairwiseAligner()
aligner.mode           = "global"
aligner.match_score    =  2
aligner.mismatch_score = -1
aligner.open_gap_score    = -5
aligner.extend_gap_score  = -0.5

aln = aligner.align(seq_7NIO_A, seq_6XEZ_E)[0]

# Map 7NIO_A seq_idx -> 6XEZ_E seq_idx
map_7NIO_to_6XEZ = {}
for (r_start, r_end), (q_start, q_end) in zip(aln.aligned[0], aln.aligned[1]):
    for r_pos, q_pos in zip(range(r_start, r_end), range(q_start, q_end)):
        map_7NIO_to_6XEZ[r_pos] = q_pos

# Build pdb_num -> seq_idx for each structure
pdb_to_seqidx_7NIO = {pdb_num: i for i, (pdb_num, _) in enumerate(reslist_7NIO_A)}
pdb_to_seqidx_6XEZ = {pdb_num: i for i, (pdb_num, _) in enumerate(reslist_6XEZ_E)}
seqidx_to_pdb_6XEZ = {i: pdb_num for pdb_num, i in pdb_to_seqidx_6XEZ.items()}

iface_7NIO_set = set(res_7NIO["interface_residues"])
iface_6XEZ_set = set(res_6XEZ["interface_residues"])

consensus = []
only_7NIO = []
only_6XEZ = []

for pdb_num_7NIO in sorted(iface_7NIO_set):
    if pdb_num_7NIO not in pdb_to_seqidx_7NIO:
        continue
    seq_idx_7NIO = pdb_to_seqidx_7NIO[pdb_num_7NIO]
    if seq_idx_7NIO in map_7NIO_to_6XEZ:
        seq_idx_6XEZ = map_7NIO_to_6XEZ[seq_idx_7NIO]
        pdb_num_6XEZ = seqidx_to_pdb_6XEZ.get(seq_idx_6XEZ)
        if pdb_num_6XEZ and pdb_num_6XEZ in iface_6XEZ_set:
            # Find residue name
            for r in struct_7NIO[0]["A"]:
                if PDB.is_aa(r) and r.get_id()[1] == pdb_num_7NIO:
                    resname = r.get_resname()
                    break
            else:
                resname = "UNK"
            consensus.append({
                "res_num_7NIO" : pdb_num_7NIO,
                "res_num_6XEZ" : pdb_num_6XEZ,
                "res_name"     : resname
            })
        else:
            only_7NIO.append(pdb_num_7NIO)
    else:
        only_7NIO.append(pdb_num_7NIO)

for pdb_num_6XEZ in sorted(iface_6XEZ_set):
    if pdb_num_6XEZ not in pdb_to_seqidx_6XEZ:
        continue
    seq_idx_6XEZ = pdb_to_seqidx_6XEZ[pdb_num_6XEZ]
    mapped_to_7NIO = False
    for s7, s6 in map_7NIO_to_6XEZ.items():
        if s6 == seq_idx_6XEZ:
            pdb_num_7 = reslist_7NIO_A[s7][0] if s7 < len(reslist_7NIO_A) else None
            if pdb_num_7 and pdb_num_7 in iface_7NIO_set:
                mapped_to_7NIO = True
            break
    if not mapped_to_7NIO:
        only_6XEZ.append(pdb_num_6XEZ)

print(f"\n  Consensus hotspots (both crystal structures): {len(consensus)}")
for c in consensus:
    print(f"    {c['res_name']:3s} 7NIO-A:{c['res_num_7NIO']:4d}  6XEZ-E:{c['res_num_6XEZ']:4d}")

print(f"\n  Only in 7NIO (primary, 2.80 A): {len(only_7NIO)} residues")
print(f"    {sorted(only_7NIO)}")

print(f"\n  Only in 6XEZ (secondary, 3.50 A): {len(only_6XEZ)} residues")
print(f"    {sorted(only_6XEZ)}")

# All hotspots for conservation analysis (7NIO numbering — primary)
all_hotspots_7NIO = sorted(iface_7NIO_set)
print(f"\n  All 7NIO-A hotspot residues (primary, used for conservation):")
for rn in all_hotspots_7NIO:
    for r in struct_7NIO[0]["A"]:
        if PDB.is_aa(r) and r.get_id()[1] == rn:
            print(f"    {r.get_resname()}{rn}")
            break

# ── AF3 pLDDT at hotspot positions ────────────────────────────────────────
print("\n" + "=" * 60)
print("AF3 pLDDT AT HOTSPOT POSITIONS (monomer — structural quality)")
print("=" * 60)

reslist_AF3 = [(r.get_id()[1], r.get_resname())
               for r in struct_AF3[0]["A"]
               if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1]
seq_AF3 = "".join(AA_MAP_3TO1[r[1]] for r in reslist_AF3)

aln_af3 = aligner.align(seq_7NIO_A, seq_AF3)[0]
map_7NIO_to_AF3 = {}
for (r_start, r_end), (q_start, q_end) in zip(aln_af3.aligned[0], aln_af3.aligned[1]):
    for r_pos, q_pos in zip(range(r_start, r_end), range(q_start, q_end)):
        map_7NIO_to_AF3[r_pos] = q_pos

seqidx_to_res_AF3 = {i: r for i, (_, r) in enumerate(reslist_AF3)
                     if True}
# rebuild properly
seqidx_to_resobj_AF3 = {}
for i, r in enumerate(r for r in struct_AF3[0]["A"]
                       if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1):
    seqidx_to_resobj_AF3[i] = r

for rn_7NIO in all_hotspots_7NIO:
    if rn_7NIO not in pdb_to_seqidx_7NIO:
        continue
    si_7NIO = pdb_to_seqidx_7NIO[rn_7NIO]
    if si_7NIO in map_7NIO_to_AF3:
        si_af3 = map_7NIO_to_AF3[si_7NIO]
        if si_af3 in seqidx_to_resobj_AF3:
            af3_res = seqidx_to_resobj_AF3[si_af3]
            plddt   = af3_res["CA"].get_bfactor() if "CA" in af3_res else 0.0
            for r in struct_7NIO[0]["A"]:
                if PDB.is_aa(r) and r.get_id()[1] == rn_7NIO:
                    rname = r.get_resname()
                    break
            else:
                rname = "UNK"
            flag = " <-- LOW" if plddt < 70 else ""
            print(f"  {rname}{rn_7NIO:4d} -> AF3 pLDDT {plddt:.1f}{flag}")

# ── Salt bridge summary ────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SALT BRIDGE SUMMARY")
print("=" * 60)

all_sbs = []
for sb in res_7NIO["salt_bridges"]:
    sb["structure"] = "7NIO"
    all_sbs.append(sb)
for sb in res_6XEZ["salt_bridges"]:
    sb["structure"] = "6XEZ"
    all_sbs.append(sb)

if all_sbs:
    for sb in all_sbs:
        print(f"  [{sb['structure']}] {sb['res1']} -- {sb['res2']} : {sb['dist']} A")
else:
    print("  No salt bridges detected in either crystal structure")

# ── Save ───────────────────────────────────────────────────────────────────
result = {
    "complex"         : "NSP13-Helicase",
    "script"          : "05_interface_NSP13-Helicase_7.py",
    "structures"      : ["7NIO (A vs E)", "6XEZ (E vs F)", "AF3 monomer (pLDDT only)"],
    "interface_7NIO"  : res_7NIO,
    "interface_6XEZ"  : res_6XEZ,
    "consensus_hotspots": consensus,
    "only_7NIO"       : only_7NIO,
    "only_6XEZ"       : only_6XEZ,
    "all_hotspots_7NIO_numbering": all_hotspots_7NIO,
    "salt_bridges_all": all_sbs
}

out_path = f"{VAL_DIR}/interface_analysis_7.json"
with open(out_path, "w") as f:
    json.dump(result, f, indent=2)
print(f"\n  Saved: {out_path}")
