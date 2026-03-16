"""
Standalone RMSD fix — superimpose AF3 onto 6VWW using BioPython Superimposer
then report correct RMSD. Patches the validation result JSON.
"""
import json
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio import pairwise2
from Bio.PDB.Polypeptide import protein_letters_3to1

PDB_APO  = "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb"
CIF_PATH = "01-alphafold3/NSP15/NSP15_best_model.cif"
OUT_JSON = "02-validation/NSP15/validation_result_9.json"

# Load structures
pdb = PDBParser(QUIET=True).get_structure("6VWW", PDB_APO)
af3 = MMCIFParser(QUIET=True).get_structure("AF3", CIF_PATH)

chain_pdb = list(pdb[0].get_chains())[0]   # chain A of true dimer
chain_af3 = list(af3[0].get_chains())[0]   # chain A of AF3 monomer

# Get sequences and CA atoms
def get_seq_ca(chain):
    res = [r for r in chain.get_residues() if r.get_id()[0] == ' ']
    seq, ca = "", []
    for r in res:
        one = protein_letters_3to1.get(r.get_resname(), 'X')
        seq += one
        if 'CA' in r:
            ca.append(r['CA'])
    return seq, ca

seq_pdb, ca_pdb = get_seq_ca(chain_pdb)
seq_af3, ca_af3 = get_seq_ca(chain_af3)

print(f"6VWW chain A : {len(seq_pdb)} aa, {len(ca_pdb)} CA atoms")
print(f"AF3  chain A : {len(seq_af3)} aa, {len(ca_af3)} CA atoms")

# Sequence alignment to get matched CA pairs
alns = pairwise2.align.globalms(
    seq_pdb, seq_af3, 2, -1, -5, -0.5, one_alignment_only=True)
aln_ref, aln_qry = alns[0][0], alns[0][1]

fixed_atoms, moving_atoms = [], []
i_r = i_q = 0
for r_aa, q_aa in zip(aln_ref, aln_qry):
    if r_aa != '-' and q_aa != '-':
        if i_r < len(ca_pdb) and i_q < len(ca_af3):
            fixed_atoms.append(ca_pdb[i_r])
            moving_atoms.append(ca_af3[i_q])
    if r_aa != '-': i_r += 1
    if q_aa != '-': i_q += 1

print(f"Aligned CA pairs: {len(fixed_atoms)}")

# Superimpose AF3 onto 6VWW
sup = Superimposer()
sup.set_atoms(fixed_atoms, moving_atoms)
sup.apply(list(af3[0].get_atoms()))

rmsd = sup.rms
print(f"RMSD after superimposition: {rmsd:.3f} A")
print(f"(was 85.9 A before — coordinate frame difference confirmed fixed)")

# Gate check
gate_rmsd = rmsd < 3.0
print(f"\nRMSD gate (< 3.0 A): {'PASS' if gate_rmsd else 'FAIL'}")

# Update JSON
with open(OUT_JSON) as f:
    result = json.load(f)

result["rmsd_vs_6VWW_chainA"]      = round(rmsd, 3)
result["rmsd_note"]                 = "SVD superimposition applied (Superimposer)"
result["gates"]["RMSD < 3.0 A"]    = gate_rmsd

all_pass = all(result["gates"].values())
result["overall_gate"] = "PASS" if all_pass else "FAIL"

with open(OUT_JSON, "w") as f:
    json.dump(result, f, indent=2)

print(f"\nAll gates:")
for k, v in result["gates"].items():
    print(f"  {k:<25}: {'PASS' if v else 'FAIL ***'}")
print(f"\nOVERALL: {result['overall_gate']}")
print(f"Updated: {OUT_JSON}")
