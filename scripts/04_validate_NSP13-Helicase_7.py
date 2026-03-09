#!/usr/bin/env python3
"""
Script 04_validate_NSP13-Helicase_7.py

NSP13 homodimer interface — PPI disruption target
AF3 was submitted as MONOMER — no iptm score available.

Validation approach:
  1. Read AF3 confidence.json (ptm, has_clash, fraction_disordered)
  2. Identify dimer interface residues from 7NIO (A vs E) and 6XEZ (E vs F)
  3. Align AF3 chain A vs 7NIO chain A — compute RMSD
  4. Report AF3 pLDDT at interface-equivalent positions

Gates:
  ptm > 0.50
  has_clash = 0.0
  RMSD < 3.0 A

Structures:
  7NIO  : NSP13 homodimer, chains A+E (2.80 A)
  6XEZ  : NSP13 dimer in RdRp complex, chains E+F (3.50 A)
  AF3   : NSP13 monomer, chain A

Output: 02-validation/NSP13-Helicase/validation_result_7.json
"""

import os, json, math
import numpy as np
from Bio import PDB, Align

# ── Paths ──────────────────────────────────────────────────────────────────
BASE     = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR  = f"{BASE}/00-reference/pdb_structures"
AF3_DIR  = f"{BASE}/01-alphafold3/NSP13-Helicase"
VAL_DIR  = f"{BASE}/02-validation/NSP13-Helicase"
os.makedirs(VAL_DIR, exist_ok=True)

PDB_7NIO  = f"{PDB_DIR}/7NIO.pdb"
PDB_6XEZ  = f"{PDB_DIR}/6XEZ.pdb"
AF3_PDB   = f"{AF3_DIR}/NSP13_Helicase_best_model.pdb"
CONF_JSON = f"{AF3_DIR}/NSP13_Helicase_confidence.json"

CUTOFF = 5.0  # Angstrom interface cutoff

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# ── Step 1: AF3 confidence scores ──────────────────────────────────────────
print("=" * 60)
print("STEP 1: AF3 confidence scores")
print("=" * 60)

with open(CONF_JSON) as f:
    conf = json.load(f)

print("Raw confidence keys:", list(conf.keys()))

ptm                 = float(conf.get("ptm", conf.get("pTM", 0)))
has_clash           = float(conf.get("has_clash", conf.get("hasClash", 0)))
fraction_disordered = float(conf.get("fraction_disordered", conf.get("fractionDisordered", 0)))
ranking_score       = float(conf.get("ranking_score", conf.get("rankingScore", 0)))

print(f"\n  ptm                 : {ptm:.3f}   (gate > 0.50)")
print(f"  has_clash           : {has_clash}     (gate = 0.0)")
print(f"  fraction_disordered : {fraction_disordered:.3f}   (gate < 0.20)")
print(f"  ranking_score       : {ranking_score:.3f}")

ptm_pass      = ptm > 0.50
clash_pass    = has_clash == 0.0
disorder_pass = fraction_disordered < 0.20

print(f"\n  ptm gate       : {'PASS' if ptm_pass else 'FAIL'}")
print(f"  clash gate     : {'PASS' if clash_pass else 'FAIL'}")
print(f"  disorder gate  : {'PASS' if disorder_pass else 'FAIL (check structure)'}")

# ── Step 2: Parse structures ───────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Parse structures")
print("=" * 60)

parser = PDB.PDBParser(QUIET=True)
struct_7NIO = parser.get_structure("7NIO", PDB_7NIO)
struct_6XEZ = parser.get_structure("6XEZ", PDB_6XEZ)
struct_AF3  = parser.get_structure("AF3",  AF3_PDB)

chain_7NIO_A = struct_7NIO[0]["A"]
chain_7NIO_E = struct_7NIO[0]["E"]
chain_6XEZ_E = struct_6XEZ[0]["E"]
chain_6XEZ_F = struct_6XEZ[0]["F"]
chain_AF3_A  = struct_AF3[0]["A"]

def count_aa(chain):
    return len([r for r in chain if PDB.is_aa(r)])

print(f"  7NIO chain A : {count_aa(chain_7NIO_A)} residues")
print(f"  7NIO chain E : {count_aa(chain_7NIO_E)} residues")
print(f"  6XEZ chain E : {count_aa(chain_6XEZ_E)} residues")
print(f"  6XEZ chain F : {count_aa(chain_6XEZ_F)} residues")
print(f"  AF3  chain A : {count_aa(chain_AF3_A)}  residues")

# ── Step 3: Interface residues ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: Dimer interface residues")
print("=" * 60)

def get_interface_residues(chain1, chain2, cutoff=5.0):
    atoms2 = [a for r in chain2 if PDB.is_aa(r) for a in r.get_atoms()]
    ns = PDB.NeighborSearch(atoms2)
    iface1, iface2 = set(), set()
    for r1 in chain1:
        if not PDB.is_aa(r1):
            continue
        for a1 in r1.get_atoms():
            hits = ns.search(a1.get_vector().get_array(), cutoff, "A")
            if hits:
                iface1.add(r1.get_id()[1])
                for a2 in hits:
                    iface2.add(a2.get_parent().get_id()[1])
    return iface1, iface2

iface_7NIO_A, iface_7NIO_E = get_interface_residues(chain_7NIO_A, chain_7NIO_E)
iface_6XEZ_E, iface_6XEZ_F = get_interface_residues(chain_6XEZ_E, chain_6XEZ_F)

print(f"  7NIO A vs E : chain A={len(iface_7NIO_A)} res, chain E={len(iface_7NIO_E)} res")
print(f"  6XEZ E vs F : chain E={len(iface_6XEZ_E)} res, chain F={len(iface_6XEZ_F)} res")
print(f"\n  7NIO chain A interface residues: {sorted(iface_7NIO_A)}")
print(f"  7NIO chain E interface residues: {sorted(iface_7NIO_E)}")
print(f"\n  6XEZ chain E interface residues: {sorted(iface_6XEZ_E)}")
print(f"  6XEZ chain F interface residues: {sorted(iface_6XEZ_F)}")

# ── Step 4: Sequence extraction + pairwise alignment ──────────────────────
print("\n" + "=" * 60)
print("STEP 4: Sequence alignment AF3 vs 7NIO chain A")
print("=" * 60)

def get_seq_ca_residues(chain):
    residues = [r for r in chain if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1]
    seq = "".join(AA_MAP_3TO1[r.get_resname()] for r in residues)
    ca   = {}
    for i, r in enumerate(residues):
        if "CA" in r:
            ca[i] = r["CA"].get_vector().get_array()
    return seq, ca, residues

seq_7NIO_A, ca_7NIO_A, res_7NIO_A = get_seq_ca_residues(chain_7NIO_A)
seq_AF3_A,  ca_AF3_A,  res_AF3_A  = get_seq_ca_residues(chain_AF3_A)

print(f"  7NIO chain A: {len(seq_7NIO_A)} aa")
print(f"  AF3  chain A: {len(seq_AF3_A)} aa")

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score    =  2
aligner.mismatch_score = -1
aligner.open_gap_score    = -5
aligner.extend_gap_score  = -0.5

alignments = aligner.align(seq_7NIO_A, seq_AF3_A)
best_aln   = alignments[0]
print(f"  Alignment score: {best_aln.score:.1f}")

# Build position mapping ref_seq_idx -> qry_seq_idx
ref_to_qry = {}
for (r_start, r_end), (q_start, q_end) in zip(best_aln.aligned[0], best_aln.aligned[1]):
    for r_pos, q_pos in zip(range(r_start, r_end), range(q_start, q_end)):
        ref_to_qry[r_pos] = q_pos

# Collect paired CA coordinates for superposition
paired_ref, paired_qry = [], []
for r_pos, q_pos in ref_to_qry.items():
    if r_pos in ca_7NIO_A and q_pos in ca_AF3_A:
        paired_ref.append(ca_7NIO_A[r_pos])
        paired_qry.append(ca_AF3_A[q_pos])

print(f"  Paired CA atoms for RMSD: {len(paired_ref)}")

# ── Step 5: Superposition RMSD ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 5: RMSD after superposition (SVD)")
print("=" * 60)

ref_arr = np.array(paired_ref)
qry_arr = np.array(paired_qry)

c_ref = ref_arr.mean(axis=0)
c_qry = qry_arr.mean(axis=0)
ref_c = ref_arr - c_ref
qry_c = qry_arr - c_qry

H = qry_c.T @ ref_c
U, S, Vt = np.linalg.svd(H)
d = np.linalg.det(Vt.T @ U.T)
D = np.diag([1.0, 1.0, d])
R = Vt.T @ D @ U.T
qry_rot = qry_c @ R.T

diff  = ref_c - qry_rot
rmsd  = math.sqrt((diff**2).sum() / len(diff))

print(f"  RMSD (superposed, {len(paired_ref)} CA pairs): {rmsd:.3f} A")
rmsd_pass = rmsd < 3.0
print(f"  RMSD gate (< 3.0 A): {'PASS' if rmsd_pass else 'FAIL'}")

# ── Step 6: pLDDT at interface positions ───────────────────────────────────
print("\n" + "=" * 60)
print("STEP 6: AF3 pLDDT at interface-equivalent positions")
print("=" * 60)

# Map 7NIO PDB residue numbers -> seq index
pdb_num_to_seqidx = {r.get_id()[1]: i for i, r in enumerate(res_7NIO_A)}
# Map AF3 seq index -> residue object
af3_seqidx_to_res = {i: r for i, r in enumerate(res_AF3_A)}

plddt_vals = []
mapped_iface = []
for pdb_num in sorted(iface_7NIO_A):
    if pdb_num in pdb_num_to_seqidx:
        seq_idx_ref = pdb_num_to_seqidx[pdb_num]
        if seq_idx_ref in ref_to_qry:
            af3_idx = ref_to_qry[seq_idx_ref]
            if af3_idx in af3_seqidx_to_res:
                af3_res = af3_seqidx_to_res[af3_idx]
                if "CA" in af3_res:
                    plddt = af3_res["CA"].get_bfactor()
                    plddt_vals.append(plddt)
                    mapped_iface.append((pdb_num, af3_res.get_resname(), round(plddt, 1)))

print(f"  Interface residues mapped to AF3: {len(mapped_iface)} / {len(iface_7NIO_A)}")
if plddt_vals:
    mean_plddt = sum(plddt_vals) / len(plddt_vals)
    min_plddt  = min(plddt_vals)
    print(f"  Mean pLDDT at interface: {mean_plddt:.1f}")
    print(f"  Min  pLDDT at interface: {min_plddt:.1f}")
    print(f"  (>=90 very high | >=70 confident | <70 low)")
    plddt_pass = mean_plddt >= 70.0
    print(f"  pLDDT gate: {'PASS' if plddt_pass else 'LOW -- check structure'}")
    print(f"\n  Per-residue pLDDT at interface positions:")
    for (pdb_num, resname, pl) in mapped_iface:
        flag = " <-- LOW" if pl < 70 else ""
        print(f"    7NIO-A res {pdb_num:4d} ({resname}) -> AF3 pLDDT {pl:.1f}{flag}")
else:
    mean_plddt = 0.0
    plddt_pass = False
    print("  WARNING: could not map interface residues to AF3")

# ── Step 7: Overall gate ───────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 7: OVERALL GATE")
print("=" * 60)

overall_pass = ptm_pass and clash_pass and rmsd_pass

print(f"  ptm > 0.50          : {'PASS' if ptm_pass else 'FAIL'}  ({ptm:.3f})")
print(f"  has_clash = 0       : {'PASS' if clash_pass else 'FAIL'}")
print(f"  fraction_disordered : {'PASS' if disorder_pass else 'CHECK'}  ({fraction_disordered:.3f})")
print(f"  RMSD < 3.0 A        : {'PASS' if rmsd_pass else 'FAIL'}  ({rmsd:.3f} A)")
print(f"  mean pLDDT >= 70    : {'PASS' if plddt_pass else 'LOW'}  ({mean_plddt:.1f})")
print(f"\n  OVERALL: {'PASS -- proceed to Script 05_7' if overall_pass else 'FAIL -- do not proceed'}")

# ── Save ───────────────────────────────────────────────────────────────────
result = {
    "complex"  : "NSP13-Helicase",
    "script"   : "04_validate_NSP13-Helicase_7.py",
    "af3_type" : "monomer_single_chain",
    "structures": {
        "primary_dimer"   : "7NIO (chains A+E, 2.80 A)",
        "secondary_dimer" : "6XEZ (chains E+F, 3.50 A)",
        "af3"             : "NSP13_Helicase_best_model.pdb (chain A)"
    },
    "af3_confidence": {
        "ptm"                : ptm,
        "has_clash"          : has_clash,
        "fraction_disordered": fraction_disordered,
        "ranking_score"      : ranking_score
    },
    "structural_validation": {
        "rmsd_af3_vs_7NIO_chainA_angstrom" : rmsd,
        "aligned_ca_pairs"                  : len(paired_ref),
        "mean_plddt_at_interface"           : round(mean_plddt, 2),
        "interface_residues_7NIO_A"         : sorted([int(x) for x in iface_7NIO_A]),
        "interface_residues_7NIO_E"         : sorted([int(x) for x in iface_7NIO_E]),
        "interface_residues_6XEZ_E"         : sorted([int(x) for x in iface_6XEZ_E]),
        "interface_residues_6XEZ_F"         : sorted([int(x) for x in iface_6XEZ_F]),
    },
    "gates": {
        "ptm_pass"      : ptm_pass,
        "clash_pass"    : clash_pass,
        "disorder_pass" : disorder_pass,
        "rmsd_pass"     : rmsd_pass,
        "plddt_pass"    : plddt_pass,
        "overall_pass"  : overall_pass
    }
}

out_path = f"{VAL_DIR}/validation_result_7.json"
with open(out_path, "w") as f:
    json.dump(result, f, indent=2)
print(f"\n  Saved: {out_path}")
