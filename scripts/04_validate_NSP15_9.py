"""
Script 04_9 — AF3 Validation NSP15 Homodimer
Complex  : NSP15-NSP15 (NendoU dimer interface)
Suffix   : _9
Reference: 6VWW chains A+B (1.90 A, apo, PRIMARY)
Note     : AF3 predicted monomer (single chain A, 346 aa)
           Same situation as NSP13-Helicase — use ptm+RMSD+pLDDT gates
"""

import json
import os
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio import pairwise2
from Bio.PDB.Polypeptide import protein_letters_3to1

AF3_DIR  = "01-alphafold3/NSP15"
PDB_APO  = "00-reference/known_interfaces/NSP15/6VWW_NSP15-dimer.pdb"
OUT_DIR  = "02-validation/NSP15"
OUT_JSON = f"{OUT_DIR}/validation_result_9.json"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Script 04_9 — AF3 Validation NSP15")
print("=" * 60)

# ── Step 1: Load confidence ────────────────────────────────────
with open(f"{AF3_DIR}/NSP15_confidence.json") as f:
    conf = json.load(f)

ptm        = conf.get("ptm")
has_clash  = conf.get("has_clash", 0.0)
ranking    = conf.get("ranking_score")
n_rec      = conf.get("num_recycles")
frac_dis   = conf.get("fraction_disordered")
chain_ptm  = conf.get("chain_ptm", [None])
cp_iptm    = conf.get("chain_pair_iptm")
cp_pae     = conf.get("chain_pair_pae_min")

print("\nAF3 confidence metrics:")
print(f"  ptm              : {ptm}")
print(f"  ranking_score    : {ranking}")
print(f"  has_clash        : {has_clash}")
print(f"  num_recycles     : {n_rec}")
print(f"  fraction_disord  : {frac_dis}")
print(f"  chain_ptm        : {chain_ptm}")
print(f"  chain_pair_iptm  : {cp_iptm}  <-- 1x1 = monomer prediction")
print(f"  chain_pair_pae   : {cp_pae}")
print("\nNOTE: AF3 predicted monomer (chain A only, 346 aa)")
print("      Same as NSP13-Helicase — using ptm+RMSD+pLDDT gates")

# ── Step 2: Load AF3 structure ─────────────────────────────────
cif_path = f"{AF3_DIR}/NSP15_best_model.cif"
af3 = MMCIFParser(QUIET=True).get_structure("AF3", cif_path)
af3_chain = list(af3[0].get_chains())[0]
af3_res = [r for r in af3_chain.get_residues() if r.get_id()[0] == ' ']
print(f"\nAF3 model: chain {af3_chain.id}, {len(af3_res)} residues")

# ── Step 3: Load crystal reference ────────────────────────────
pdb = PDBParser(QUIET=True).get_structure("6VWW", PDB_APO)
pdb_chain_A = list(pdb[0].get_chains())[0]
pdb_res = [r for r in pdb_chain_A.get_residues() if r.get_id()[0] == ' ']
print(f"6VWW chain A: {len(pdb_res)} residues")

# ── Step 4: Sequence extraction ────────────────────────────────
def get_seq_ca(chain):
    res = [r for r in chain.get_residues() if r.get_id()[0] == ' ']
    seq, ca = "", []
    for r in res:
        one = protein_letters_3to1.get(r.get_resname(), 'X')
        seq += one
        if 'CA' in r:
            ca.append((r.get_id()[1], r['CA']))
    return seq, ca

seq_pdb, ca_pdb = get_seq_ca(pdb_chain_A)
seq_af3, ca_af3 = get_seq_ca(af3_chain)
print(f"\nSequences: 6VWW={len(seq_pdb)}aa  AF3={len(seq_af3)}aa")

# ── Step 5: Alignment-based RMSD ──────────────────────────────
alns = pairwise2.align.globalms(
    seq_pdb, seq_af3, 2, -1, -5, -0.5, one_alignment_only=True)
aln_ref, aln_qry = alns[0][0], alns[0][1]

cas_r, cas_q = [], []
i_r = i_q = 0
for r_aa, q_aa in zip(aln_ref, aln_qry):
    if r_aa != '-' and q_aa != '-':
        if i_r < len(ca_pdb) and i_q < len(ca_af3):
            cas_r.append(ca_pdb[i_r][1])
            cas_q.append(ca_af3[i_q][1])
    if r_aa != '-': i_r += 1
    if q_aa != '-': i_q += 1

coords_r = np.array([c.get_vector().get_array() for c in cas_r])
coords_q = np.array([c.get_vector().get_array() for c in cas_q])
rmsd = float(np.sqrt(((coords_r - coords_q)**2).sum(axis=1).mean()))
print(f"Aligned pairs : {len(cas_r)}")
print(f"RMSD vs 6VWW  : {rmsd:.3f} A")

# ── Step 6: Interface residues (6VWW A vs B) ──────────────────
chain_A = list(pdb[0].get_chains())[0]
chain_B = list(pdb[0].get_chains())[1]
iface_A, iface_B = set(), set()
CUTOFF = 5.0

for rA in chain_A.get_residues():
    if rA.get_id()[0] != ' ': continue
    for rB in chain_B.get_residues():
        if rB.get_id()[0] != ' ': continue
        try:
            for aA in rA.get_atoms():
                for aB in rB.get_atoms():
                    if aA - aB < CUTOFF:
                        iface_A.add(rA.get_id()[1])
                        iface_B.add(rB.get_id()[1])
                        raise StopIteration
        except StopIteration:
            pass

print(f"\nInterface residues (6VWW 5.0A cutoff):")
print(f"  Chain A: {len(iface_A)} residues -> {sorted(iface_A)}")
print(f"  Chain B: {len(iface_B)} residues -> {sorted(iface_B)}")

# ── Step 7: pLDDT at interface ─────────────────────────────────
plddt_map = {}
for r in af3_chain.get_residues():
    if r.get_id()[0] == ' ' and 'CA' in r:
        plddt_map[r.get_id()[1]] = r['CA'].get_bfactor()

# Map 6VWW interface positions to AF3 positions via alignment
pdb_res_list = [r.get_id()[1] for r in chain_A.get_residues()
                if r.get_id()[0] == ' ']
af3_res_list = [r.get_id()[1] for r in af3_chain.get_residues()
                if r.get_id()[0] == ' ']

pdb_to_af3 = {}
i_r = i_q = 0
for r_aa, q_aa in zip(aln_ref, aln_qry):
    if r_aa != '-' and q_aa != '-':
        if i_r < len(pdb_res_list) and i_q < len(af3_res_list):
            pdb_to_af3[pdb_res_list[i_r]] = af3_res_list[i_q]
    if r_aa != '-': i_r += 1
    if q_aa != '-': i_q += 1

plddts = []
for pos in sorted(iface_A):
    af3_pos = pdb_to_af3.get(pos)
    if af3_pos and af3_pos in plddt_map:
        plddts.append(plddt_map[af3_pos])

mean_plddt = float(np.mean(plddts)) if plddts else 0.0
min_plddt  = float(np.min(plddts))  if plddts else 0.0
print(f"\npLDDT at interface residues (AF3):")
print(f"  Mean : {mean_plddt:.1f}")
print(f"  Min  : {min_plddt:.1f}")

# ── Step 8: Gate evaluation ────────────────────────────────────
print("\n" + "=" * 60)
print("Gate evaluation (monomer mode — no iptm gate):")
gates = {
    "ptm >= 0.70"      : ptm is not None and float(ptm) >= 0.70,
    "has_clash = 0"    : float(has_clash) == 0.0,
    "RMSD < 3.0 A"     : rmsd < 3.0,
    "mean pLDDT >= 70" : mean_plddt >= 70.0,
}
all_pass = True
for gate, result in gates.items():
    status = "PASS" if result else "FAIL ***"
    print(f"  {gate:<25}: {status}")
    if not result: all_pass = False

overall = "PASS" if all_pass else "FAIL"
print(f"\n  OVERALL: {overall}")

# ── Step 9: Save results ───────────────────────────────────────
result = {
    "complex"             : "NSP15-homodimer",
    "suffix"              : "_9",
    "reference_pdb"       : "6VWW",
    "af3_mode"            : "monomer (single chain A, 346aa)",
    "af3_confidence"      : {
        "ptm"             : float(ptm) if ptm else None,
        "ranking_score"   : float(ranking) if ranking else None,
        "has_clash"       : float(has_clash),
        "num_recycles"    : float(n_rec) if n_rec else None,
        "fraction_disord" : float(frac_dis) if frac_dis else None,
        "chain_pair_iptm" : cp_iptm,
        "chain_pair_pae"  : cp_pae,
    },
    "rmsd_vs_6VWW_chainA" : round(rmsd, 3),
    "aligned_pairs"       : len(cas_r),
    "interface_A"         : sorted(iface_A),
    "interface_B"         : sorted(iface_B),
    "n_interface_A"       : len(iface_A),
    "n_interface_B"       : len(iface_B),
    "mean_plddt_iface"    : round(mean_plddt, 1),
    "min_plddt_iface"     : round(min_plddt, 1),
    "gates"               : {k: bool(v) for k, v in gates.items()},
    "overall_gate"        : overall,
    "note"                : "AF3 monomer prediction — ptm+RMSD+pLDDT gates used (same as NSP13-Helicase)"
}

with open(OUT_JSON, "w") as f:
    json.dump(result, f, indent=2)

print(f"\nSaved: {OUT_JSON}")
print("Next : Script 05_9 (interface analysis)")
