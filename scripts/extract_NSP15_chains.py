from Bio.PDB import PDBParser, PDBIO, Select
import os

class ChainSelect(Select):
    def __init__(self, chains):
        self.chains = chains
    def accept_chain(self, chain):
        return chain.get_id() in self.chains
    def accept_residue(self, residue):
        return residue.get_id()[0] == ' '

parser = PDBParser(QUIET=True)
outdir = "00-reference/known_interfaces/NSP15"
os.makedirs(outdir, exist_ok=True)

extractions = {
    "6VWW": (["A","B"], "NSP15-dimer"),
    "6W01": (["A","B"], "NSP15-dimer"),
    "6WLC": (["A","B"], "NSP15-dimer-UMP"),
    "6WXC": (["A","B"], "NSP15-dimer-tipiracil"),
    "2H85": (["A"],    "NSP15-monomer-SARS1"),
}

for pdb_id, (chains, label) in extractions.items():
    inpath  = f"00-reference/pdb_structures/NSP15/{pdb_id}.pdb"
    outpath = f"{outdir}/{pdb_id}_{label}.pdb"
    struct  = parser.get_structure(pdb_id, inpath)
    io = PDBIO()
    io.set_structure(struct)
    io.save(outpath, ChainSelect(chains))
    s2   = PDBParser(QUIET=True).get_structure("x", outpath)
    res  = [r for c in s2.get_chains()
              for r in c.get_residues() if r.get_id()[0]==' ']
    atoms = list(s2.get_atoms())
    print(f"  {pdb_id} {chains}: {len(res)} residues, {len(atoms)} atoms -> {outpath}")

print("\nChain extraction complete.")
