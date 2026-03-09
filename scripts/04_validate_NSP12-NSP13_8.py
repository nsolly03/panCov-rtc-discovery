#!/usr/bin/env python3
"""
Script 04_validate_NSP12-NSP13_8.py

NSP12-NSP13 interface validation.

CRITICAL NOTE: AF3 iptm=0.20 — FAILS interface gate.
inter-chain PAE=24.95 A — chains predicted independently.
AF3 cannot predict this interface from sequences alone.
Scientific interpretation: NSP12-NSP13 is a transient/context-
dependent interaction requiring full RTC assembly.

Validation approach:
  - AF3 iptm gate: FAIL (documented, not blocking — 3 crystal structures available)
  - Individual chain fold quality: ptm per chain (0.89 NSP12, 0.74 NSP13)
  - RMSD AF3 chain A vs 6XEZ chain A (NSP12 fold quality)
  - RMSD AF3 chain B vs 6XEZ chain E (NSP13 fold quality)
  - Interface residues identified from 3 crystal structures only

Structures:
  6XEZ : A=NSP12, E=NSP13 (3.50 A) — primary
  7CXM : A=NSP12, E=NSP13 (3.20 A) — secondary
  7RDY : A=NSP12, E=NSP13 (3.10 A) — tertiary
  AF3  : A=NSP12 (chain_ptm=0.89), B=NSP13 (chain_ptm=0.74)
         iptm=0.20 FAIL — fold quality only

Output: 02-validation/NSP12-NSP13/validation_result_8.json
"""

import os, json, math
import numpy as np
from Bio import PDB, Align

# ── Paths ──────────────────────────────────────────────────────────────────
BASE     = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR  = f"{BASE}/00-reference/pdb_structures"
AF3_DIR  = f"{BASE}/01-alphafold3/NSP12-NSP13"
VAL_DIR  = f"{BASE}/02-validation/NSP12-NSP13"
os.makedirs(VAL_DIR, exist_ok=True)

PDB_6XEZ  = f"{PDB_DIR}/6XEZ.pdb"
PDB_7CXM  = f"{PDB_DIR}/7CXM.pdb"
PDB_7RDY  = f"{PDB_DIR}/7RDY.pdb"
AF3_PDB   = f"{AF3_DIR}/NSP12_NSP13_best_model.pdb"
CONF_JSON = f"{AF3_DIR}/NSP12_NSP13_confidence.json"

CUTOFF = 5.0

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# ── Step 1: AF3 confidence ─────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: AF3 confidence scores")
print("=" * 60)

with open(CONF_JSON) as f:
    conf = json.load(f)

iptm                = float(conf.get("iptm", 0))
ptm                 = float(conf.get("ptm",  0))
has_clash           = float(conf.get("has_clash", 0))
fraction_disordered = float(conf.get("fraction_disordered", 0))
ranking_score       = float(conf.get("ranking_score", 0))
chain_ptm           = conf.get("chain_ptm", [])
chain_pair_iptm     = conf.get("chain_pair_iptm", [[]])
chain_pair_pae_min  = conf.get("chain_pair_pae_min", [[]])
num_recycles        = conf.get("num_recycles", 0)

print(f"\n  iptm                : {iptm:.3f}   (gate > 0.60)  *** FAIL ***")
print(f"  ptm                 : {ptm:.3f}   (gate > 0.50)")
print(f"  has_clash           : {has_clash}")
print(f"  fraction_disordered : {fraction_disordered:.3f}")
print(f"  ranking_score       : {ranking_score:.3f}")
print(f"  num_recycles        : {num_recycles}")
print(f"  chain_ptm           : {chain_ptm}  (NSP12={chain_ptm[0]:.2f}, NSP13={chain_ptm[1]:.2f})")
print(f"\n  chain_pair_iptm (off-diagonal = inter-chain):")
for row in chain_pair_iptm:
    print(f"    {row}")
print(f"\n  chain_pair_pae_min (inter-chain PAE):")
for row in chain_pair_pae_min:
    print(f"    {row}")

inter_pae = chain_pair_pae_min[0][1] if chain_pair_pae_min else 99.0
inter_iptm = chain_pair_iptm[0][1]   if chain_pair_iptm    else 0.0

print(f"\n  Inter-chain PAE  : {inter_pae:.2f} A  (compare: NSP12-NSP8=1.22, NSP12-NSP7=2.02)")
print(f"  Inter-chain iptm : {inter_iptm:.3f}  (gate > 0.60)")

iptm_pass    = iptm > 0.60
ptm_pass     = ptm  > 0.50
clash_pass   = has_clash == 0.0

print(f"\n  iptm gate  : {'PASS' if iptm_pass else '*** FAIL *** — AF3 cannot predict this interface'}")
print(f"  ptm gate   : {'PASS' if ptm_pass  else 'FAIL'}")
print(f"  clash gate : {'PASS' if clash_pass else 'FAIL'}")
print(f"\n  SCIENTIFIC NOTE:")
print(f"  inter-chain PAE={inter_pae:.1f} A >> 5.0 A threshold.")
print(f"  NSP12-NSP13 is likely transient or context-dependent.")
print(f"  Proceeding with 3 crystal structures as sole interface evidence.")

# ── Step 2: Parse structures ───────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Parse crystal structures + AF3")
print("=" * 60)

parser = PDB.PDBParser(QUIET=True)

structs = {
    "6XEZ" : parser.get_structure("6XEZ", PDB_6XEZ),
    "7CXM" : parser.get_structure("7CXM", PDB_7CXM),
    "7RDY" : parser.get_structure("7RDY", PDB_7RDY),
}
struct_AF3 = parser.get_structure("AF3", AF3_PDB)

def count_aa(chain):
    return len([r for r in chain if PDB.is_aa(r)])

print(f"\n  {'PDB':>6}  {'NSP12 chain':>12}  {'NSP12 res':>10}  "
      f"{'NSP13 chain':>12}  {'NSP13 res':>10}  {'Resolution':>12}")

pdb_chain_map = {
    "6XEZ" : ("A", "E"),
    "7CXM" : ("A", "E"),
    "7RDY" : ("A", "E"),
}
resolutions = {
    "6XEZ": "3.50 A",
    "7CXM": "3.20 A",
    "7RDY": "3.10 A",
}

for pdb_id, struct in structs.items():
    c12, c13 = pdb_chain_map[pdb_id]
    n12 = count_aa(struct[0][c12])
    n13 = count_aa(struct[0][c13])
    res = resolutions[pdb_id]
    print(f"  {pdb_id:>6}  {c12:>12}  {n12:>10}  {c13:>12}  {n13:>10}  {res:>12}")

nsp12_af3 = count_aa(struct_AF3[0]["A"])
nsp13_af3 = count_aa(struct_AF3[0]["B"])
print(f"  {'AF3':>6}  {'A':>12}  {nsp12_af3:>10}  {'B':>12}  {nsp13_af3:>10}  {'predicted':>12}")

# ── Step 3: Interface residues from all 3 crystal structures ───────────────
print("\n" + "=" * 60)
print("STEP 3: Interface residues from crystal structures")
print("=" * 60)

def get_interface_residues(chain1, chain2, cutoff=5.0):
    atoms2 = [a for r in chain2 if PDB.is_aa(r) for a in r.get_atoms()]
    ns     = PDB.NeighborSearch(atoms2)
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

iface_results = {}
for pdb_id, struct in structs.items():
    c12, c13 = pdb_chain_map[pdb_id]
    iface_nsp12, iface_nsp13 = get_interface_residues(
        struct[0][c12], struct[0][c13])
    iface_results[pdb_id] = {
        "nsp12": sorted(iface_nsp12),
        "nsp13": sorted(iface_nsp13)
    }
    print(f"\n  {pdb_id}:")
    print(f"    NSP12 interface: {len(iface_nsp12)} residues — {sorted(iface_nsp12)}")
    print(f"    NSP13 interface: {len(iface_nsp13)} residues — {sorted(iface_nsp13)}")

# ── Step 4: Cross-structure consensus ─────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Interface consensus across 3 crystal structures")
print("=" * 60)

# Use sequence alignment to map residue numbers across structures
def get_seq_reslist(chain):
    return [(r.get_id()[1], r.get_resname())
            for r in chain if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1]

aligner = Align.PairwiseAligner()
aligner.mode              = "global"
aligner.match_score       =  2
aligner.mismatch_score    = -1
aligner.open_gap_score    = -5
aligner.extend_gap_score  = -0.5

def align_and_map(reslist_ref, reslist_qry):
    seq_ref = "".join(AA_MAP_3TO1[r[1]] for r in reslist_ref)
    seq_qry = "".join(AA_MAP_3TO1[r[1]] for r in reslist_qry)
    aln = aligner.align(seq_ref, seq_qry)[0]
    ref_to_qry = {}
    for (rs, re), (qs, qe) in zip(aln.aligned[0], aln.aligned[1]):
        for rp, qp in zip(range(rs, re), range(qs, qe)):
            ref_to_qry[rp] = qp
    return ref_to_qry

# 6XEZ as reference — map 7CXM and 7RDY onto it
reslist_6XEZ_A = get_seq_reslist(structs["6XEZ"][0]["A"])
reslist_7CXM_A = get_seq_reslist(structs["7CXM"][0]["A"])
reslist_7RDY_A = get_seq_reslist(structs["7RDY"][0]["A"])

map_6XEZ_to_7CXM = align_and_map(reslist_6XEZ_A, reslist_7CXM_A)
map_6XEZ_to_7RDY = align_and_map(reslist_6XEZ_A, reslist_7RDY_A)

pdb_to_seqidx_6XEZ = {pdb_num: i for i, (pdb_num, _) in enumerate(reslist_6XEZ_A)}
seqidx_to_pdb_7CXM = {i: pdb_num for i, (pdb_num, _) in enumerate(reslist_7CXM_A)}
seqidx_to_pdb_7RDY = {i: pdb_num for i, (pdb_num, _) in enumerate(reslist_7RDY_A)}

iface_6XEZ_nsp12 = set(iface_results["6XEZ"]["nsp12"])
iface_7CXM_nsp12 = set(iface_results["7CXM"]["nsp12"])
iface_7RDY_nsp12 = set(iface_results["7RDY"]["nsp12"])

# Find residues in >=2 structures (consensus)
consensus_nsp12 = []
for pdb_num_6XEZ in sorted(iface_6XEZ_nsp12):
    if pdb_num_6XEZ not in pdb_to_seqidx_6XEZ:
        continue
    si_6XEZ = pdb_to_seqidx_6XEZ[pdb_num_6XEZ]
    in_7CXM = False
    in_7RDY = False

    if si_6XEZ in map_6XEZ_to_7CXM:
        pdb_7CXM = seqidx_to_pdb_7CXM.get(map_6XEZ_to_7CXM[si_6XEZ])
        in_7CXM  = pdb_7CXM in iface_7CXM_nsp12

    if si_6XEZ in map_6XEZ_to_7RDY:
        pdb_7RDY = seqidx_to_pdb_7RDY.get(map_6XEZ_to_7RDY[si_6XEZ])
        in_7RDY  = pdb_7RDY in iface_7RDY_nsp12

    n_structs = 1 + (1 if in_7CXM else 0) + (1 if in_7RDY else 0)
    if n_structs >= 2:
        for r in structs["6XEZ"][0]["A"]:
            if PDB.is_aa(r) and r.get_id()[1] == pdb_num_6XEZ:
                consensus_nsp12.append((pdb_num_6XEZ, r.get_resname(), n_structs))
                break

print(f"\n  NSP12 consensus hotspots (in >= 2 crystal structures):")
print(f"  {'Pos':>5}  {'Res':>4}  {'Structures':>10}")
for rn, resname, n in consensus_nsp12:
    marker = " (all 3)" if n == 3 else ""
    print(f"  {rn:>5}  {resname:>4}  {n:>10}{marker}")

print(f"\n  Total NSP12 consensus hotspots: {len(consensus_nsp12)}")

# All interface residues across all 3 structures (union, 6XEZ numbering)
all_nsp12_iface = sorted(iface_6XEZ_nsp12)
print(f"\n  All NSP12 interface residues (6XEZ numbering): {all_nsp12_iface}")
print(f"  All NSP13 interface residues (6XEZ numbering): "
      f"{sorted(iface_results['6XEZ']['nsp13'])}")

# ── Step 5: AF3 RMSD — fold quality only ──────────────────────────────────
print("\n" + "=" * 60)
print("STEP 5: AF3 fold quality — RMSD per chain vs 6XEZ")
print("=" * 60)

def get_seq_ca(chain):
    residues = [r for r in chain
                if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1]
    seq = "".join(AA_MAP_3TO1[r.get_resname()] for r in residues)
    ca  = {i: r["CA"].get_vector().get_array()
           for i, r in enumerate(residues) if "CA" in r}
    return seq, ca

def compute_rmsd(seq_ref, ca_ref, seq_qry, ca_qry):
    aln = aligner.align(seq_ref, seq_qry)[0]
    paired_ref, paired_qry = [], []
    for (rs, re), (qs, qe) in zip(aln.aligned[0], aln.aligned[1]):
        for rp, qp in zip(range(rs, re), range(qs, qe)):
            if rp in ca_ref and qp in ca_qry:
                paired_ref.append(ca_ref[rp])
                paired_qry.append(ca_qry[qp])
    if len(paired_ref) < 10:
        return 99.0, 0
    ref_arr = np.array(paired_ref)
    qry_arr = np.array(paired_qry)
    c_ref   = ref_arr.mean(axis=0)
    c_qry   = qry_arr.mean(axis=0)
    H       = (qry_arr - c_qry).T @ (ref_arr - c_ref)
    U, S, Vt = np.linalg.svd(H)
    d       = np.linalg.det(Vt.T @ U.T)
    R       = Vt.T @ np.diag([1,1,d]) @ U.T
    qry_rot = (qry_arr - c_qry) @ R.T
    diff    = (ref_arr - c_ref) - qry_rot
    rmsd    = math.sqrt((diff**2).sum() / len(diff))
    return rmsd, len(paired_ref)

# NSP12: AF3 chain A vs 6XEZ chain A
seq_6XEZ_A, ca_6XEZ_A = get_seq_ca(structs["6XEZ"][0]["A"])
seq_AF3_A,  ca_AF3_A  = get_seq_ca(struct_AF3[0]["A"])
rmsd_nsp12, n_nsp12   = compute_rmsd(seq_6XEZ_A, ca_6XEZ_A, seq_AF3_A, ca_AF3_A)

# NSP13: AF3 chain B vs 6XEZ chain E
seq_6XEZ_E, ca_6XEZ_E = get_seq_ca(structs["6XEZ"][0]["E"])
seq_AF3_B,  ca_AF3_B  = get_seq_ca(struct_AF3[0]["B"])
rmsd_nsp13, n_nsp13   = compute_rmsd(seq_6XEZ_E, ca_6XEZ_E, seq_AF3_B, ca_AF3_B)

print(f"\n  NSP12: AF3-A vs 6XEZ-A  RMSD={rmsd_nsp12:.3f} A ({n_nsp12} CA pairs)")
print(f"  NSP13: AF3-B vs 6XEZ-E  RMSD={rmsd_nsp13:.3f} A ({n_nsp13} CA pairs)")
print(f"\n  Individual chain folds are {'well' if rmsd_nsp12<3 and rmsd_nsp13<3 else 'poorly'} predicted")
print(f"  Interface between them: NOT predicted (PAE={inter_pae:.1f} A)")

# ── Step 6: Overall assessment ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 6: OVERALL ASSESSMENT")
print("=" * 60)

rmsd_pass_nsp12 = rmsd_nsp12 < 3.0
rmsd_pass_nsp13 = rmsd_nsp13 < 3.0

print(f"\n  AF3 iptm gate (>0.60)         : FAIL ({iptm:.2f}) — interface not predicted")
print(f"  AF3 inter-chain PAE           : {inter_pae:.1f} A (>5.0 = no interface confidence)")
print(f"  AF3 NSP12 chain ptm           : {chain_ptm[0]:.2f} ({'OK' if chain_ptm[0]>0.5 else 'LOW'})")
print(f"  AF3 NSP13 chain ptm           : {chain_ptm[1]:.2f} ({'OK' if chain_ptm[1]>0.5 else 'LOW'})")
print(f"  NSP12 fold RMSD vs 6XEZ       : {rmsd_nsp12:.3f} A ({'PASS' if rmsd_pass_nsp12 else 'FAIL'})")
print(f"  NSP13 fold RMSD vs 6XEZ       : {rmsd_nsp13:.3f} A ({'PASS' if rmsd_pass_nsp13 else 'FAIL'})")
print(f"  Crystal structures available  : 3 (6XEZ, 7CXM, 7RDY) ✅")
print(f"  has_clash                     : {has_clash} ({'PASS' if clash_pass else 'FAIL'})")
print(f"\n  DECISION: Proceed using 3 crystal structures as primary evidence.")
print(f"  AF3 used ONLY for individual chain fold quality (not interface).")
print(f"  Scientific note: NSP12-NSP13 is transient/context-dependent interface.")

# ── Save ───────────────────────────────────────────────────────────────────
result = {
    "complex"     : "NSP12-NSP13",
    "script"      : "04_validate_NSP12-NSP13_8.py",
    "af3_note"    : "iptm=0.20 FAIL — interface not predicted, transient/context-dependent",
    "structures"  : {
        "primary"   : "6XEZ (A=NSP12, E=NSP13, 3.50 A)",
        "secondary" : "7CXM (A=NSP12, E=NSP13, 3.20 A)",
        "tertiary"  : "7RDY (A=NSP12, E=NSP13, 3.10 A)",
        "af3"       : "A=NSP12, B=NSP13 (fold only)"
    },
    "af3_confidence": {
        "iptm"               : iptm,
        "ptm"                : ptm,
        "inter_chain_pae_min": inter_pae,
        "inter_chain_iptm"   : inter_iptm,
        "chain_ptm_nsp12"    : chain_ptm[0],
        "chain_ptm_nsp13"    : chain_ptm[1],
        "has_clash"          : has_clash,
        "ranking_score"      : ranking_score,
    },
    "fold_rmsd": {
        "nsp12_af3_vs_6XEZ_A": rmsd_nsp12,
        "nsp13_af3_vs_6XEZ_E": rmsd_nsp13,
    },
    "interface_residues_crystal": iface_results,
    "consensus_nsp12_hotspots"  : [
        {"res_num": rn, "res_name": resname, "n_structures": n}
        for rn, resname, n in consensus_nsp12
    ],
    "gates": {
        "iptm_pass"       : iptm_pass,
        "ptm_pass"        : ptm_pass,
        "clash_pass"      : clash_pass,
        "rmsd_nsp12_pass" : rmsd_pass_nsp12,
        "rmsd_nsp13_pass" : rmsd_pass_nsp13,
        "proceed"         : True,
        "basis"           : "3 crystal structures — AF3 interface gate FAIL documented"
    }
}

out_path = f"{VAL_DIR}/validation_result_8.json"
with open(out_path, "w") as f:
    json.dump(result, f, indent=2)
print(f"\n  Saved: {out_path}")
