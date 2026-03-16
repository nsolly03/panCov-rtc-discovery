"""
Extract true dimer (biological assembly A+A-2) for all NSP15 structures.
6VWW already done — process 6W01, 6WLC, 6WXC.
"""
import requests
from Bio.PDB import MMCIFParser, PDBIO, Structure, Model, Chain
import os

outdir = "00-reference/known_interfaces/NSP15"
os.makedirs(outdir, exist_ok=True)

pdbs = {
    "6W01": "6W01_NSP15-true-dimer.pdb",
    "6WLC": "6WLC_NSP15-true-dimer-UMP.pdb",
    "6WXC": "6WXC_NSP15-true-dimer-tipiracil.pdb",
}

for pdb_id, outname in pdbs.items():
    print(f"\nProcessing {pdb_id}...")

    # Download biological assembly
    asm_path = f"00-reference/pdb_structures/NSP15/{pdb_id}_assembly1.cif"
    if not os.path.exists(asm_path):
        url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif"
        r = requests.get(url, timeout=60)
        if r.status_code == 200:
            with open(asm_path, "w") as f:
                f.write(r.text)
            print(f"  Downloaded: {asm_path}")
        else:
            print(f"  FAILED download: {r.status_code}")
            continue
    else:
        print(f"  Already exists: {asm_path}")

    # Parse and find all chains
    struct = MMCIFParser(QUIET=True).get_structure(pdb_id, asm_path)
    model  = struct[0]
    chains = list(model.get_chains())
    print(f"  Chains in assembly: {[c.id for c in chains]}")

    # Find chain A and chain A-2
    chain_A  = None
    chain_A2 = None
    for c in chains:
        if c.id == 'A':
            chain_A = c
        if c.id == 'A-2':
            chain_A2 = c

    if chain_A is None or chain_A2 is None:
        print(f"  WARNING: A or A-2 not found — chains: {[c.id for c in chains]}")
        # Try alternative naming
        for c in chains:
            print(f"    Chain {c.id}: {len([r for r in c.get_residues() if r.get_id()[0]==' '])} aa")
        continue

    # Build clean new structure
    new_struct = Structure.Structure(pdb_id)
    new_model  = Model.Model(0)
    new_struct.add(new_model)

    new_chain_A = Chain.Chain('A')
    new_model.add(new_chain_A)
    for res in chain_A.get_residues():
        if res.get_id()[0] == ' ':
            new_chain_A.add(res.copy())

    new_chain_B = Chain.Chain('B')
    new_model.add(new_chain_B)
    for res in chain_A2.get_residues():
        if res.get_id()[0] == ' ':
            new_chain_B.add(res.copy())

    # Save
    outpath = f"{outdir}/{outname}"
    io = PDBIO()
    io.set_structure(new_struct)
    io.save(outpath)

    # Verify contacts
    from Bio.PDB import PDBParser
    s2 = PDBParser(QUIET=True).get_structure("x", outpath)
    ch = list(s2[0].get_chains())
    resA = [r for r in ch[0].get_residues() if r.get_id()[0]==' ']
    resB = [r for r in ch[1].get_residues() if r.get_id()[0]==' ']

    contacts = 0
    iface_a, iface_b = set(), set()
    for rA in resA:
        for rB in resB:
            try:
                for aA in rA.get_atoms():
                    for aB in rB.get_atoms():
                        if aA - aB < 5.0:
                            contacts += 1
                            iface_a.add(rA.get_id()[1])
                            iface_b.add(rB.get_id()[1])
                            raise StopIteration
            except StopIteration:
                pass

    print(f"  Contacts: {contacts} | A:{len(iface_a)} res | B:{len(iface_b)} res")
    print(f"  Saved: {outpath}")

print("\nDone. Update STRUCTURES dict in Script 05_9 to use true-dimer paths.")
