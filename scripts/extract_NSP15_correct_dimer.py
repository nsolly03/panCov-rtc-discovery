from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain, Select
import os

outdir = "00-reference/known_interfaces/NSP15"
os.makedirs(outdir, exist_ok=True)

# Load full biological assembly
struct = MMCIFParser(QUIET=True).get_structure(
    '6VWW_asm',
    '00-reference/pdb_structures/NSP15/6VWW_assembly1.cif')

model = struct[0]

# Get the two chains we want: A and A-2
chain_A  = model['A']
chain_A2 = None
for c in model.get_chains():
    if c.id == 'A-2':
        chain_A2 = c
        break

if chain_A2 is None:
    raise ValueError("Chain A-2 not found in assembly")

# Build a brand new structure with just these two chains
new_struct = Structure.Structure('NSP15_dimer')
new_model  = Model.Model(0)
new_struct.add(new_model)

# Copy chain A as chain A
new_chain_A = Chain.Chain('A')
new_model.add(new_chain_A)
for residue in chain_A.get_residues():
    if residue.get_id()[0] == ' ':   # ATOM only, skip HETATM
        new_chain_A.add(residue.copy())

# Copy chain A-2 as chain B
new_chain_B = Chain.Chain('B')
new_model.add(new_chain_B)
for residue in chain_A2.get_residues():
    if residue.get_id()[0] == ' ':
        new_chain_B.add(residue.copy())

# Save
outpath = f"{outdir}/6VWW_NSP15-true-dimer.pdb"
io = PDBIO()
io.set_structure(new_struct)
io.save(outpath)

# Verify
from Bio.PDB import PDBParser
s2 = PDBParser(QUIET=True).get_structure("x", outpath)
chains = list(s2[0].get_chains())
print("Chains in saved dimer:", [c.id for c in chains])
for c in chains:
    res = [r for r in c.get_residues() if r.get_id()[0]==' ']
    print(f"  Chain {c.id}: {len(res)} residues")

# Verify contacts
CUTOFF = 5.0
iface_a, iface_b = set(), set()
contacts = 0
resA = [r for r in chains[0].get_residues() if r.get_id()[0]==' ']
resB = [r for r in chains[1].get_residues() if r.get_id()[0]==' ']
for rA in resA:
    for rB in resB:
        try:
            for aA in rA.get_atoms():
                for aB in rB.get_atoms():
                    if aA - aB < CUTOFF:
                        contacts += 1
                        iface_a.add(rA.get_id()[1])
                        iface_b.add(rB.get_id()[1])
                        raise StopIteration
        except StopIteration:
            pass

print(f"\nContacts    : {contacts}  (expected ~75)")
print(f"Interface A : {len(iface_a)} res -> {sorted(iface_a)}")
print(f"Interface B : {len(iface_b)} res -> {sorted(iface_b)}")
print(f"\nSaved: {outpath}")
