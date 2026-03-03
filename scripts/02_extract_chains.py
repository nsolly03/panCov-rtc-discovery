"""
Script 02: Extract binary chains from multi-chain PDB structures
================================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Takes multi-chain PDB files (trimers, complexes) and extracts
  only the two chains needed for each binary interface.
  Output files are named: PDBID_COMPLEX.pdb
  Example: 7BV2_NSP12-NSP7.pdb contains only chains A and C

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/02_extract_chains.py
"""

from pathlib import Path
from Bio import PDB

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"

# ── What to extract ────────────────────────────────────────
# format: (source_pdb, chains_to_keep, output_name)
EXTRACTIONS = [
    # NSP12-NSP7 (chains A + C)
    ("7BV2", ["A", "C"], "7BV2_NSP12-NSP7"),
    ("6NUR", ["A", "C"], "6NUR_NSP12-NSP7"),
    ("7C2K", ["A", "C"], "7C2K_NSP12-NSP7"),

    # NSP12-NSP8 (chains A + B)
    ("7BV2", ["A", "B"], "7BV2_NSP12-NSP8"),
    ("6NUR", ["A", "B"], "6NUR_NSP12-NSP8"),
    ("7C2K", ["A", "B"], "7C2K_NSP12-NSP8"),

    # NSP7-NSP8 (chains C + B)
    ("7BV2", ["C", "B"], "7BV2_NSP7-NSP8"),
    ("6NUR", ["C", "B"], "6NUR_NSP7-NSP8"),

    # NSP9-NSP12 (chains A + G)
    ("8SQK", ["A", "G"], "8SQK_NSP9-NSP12"),

    # NSP12-NSP13 (chains A + E)
    ("6XEZ", ["A", "E"], "6XEZ_NSP12-NSP13"),
    ("7CXM", ["A", "E"], "7CXM_NSP12-NSP13"),
    ("7RDY", ["A", "E"], "7RDY_NSP12-NSP13"),

    # NSP13 alone (chain E from 6XEZ, chain A from 7NIO)
    ("6XEZ", ["E"],      "6XEZ_NSP13"),
    ("7NIO", ["A"],      "7NIO_NSP13"),
]


class ChainSelector(PDB.Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_chain(self, chain):
        return chain.id in self.chains


def extract(source_pdb, chains, out_name):
    in_path  = PDB_DIR / f"{source_pdb}.pdb"
    out_path = PDB_DIR / f"{out_name}.pdb"

    if out_path.exists():
        print(f"  SKIP  {out_name}.pdb  (already exists)")
        return True

    if not in_path.exists():
        print(f"  FAIL  {source_pdb}.pdb not found")
        return False

    try:
        parser    = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(source_pdb, in_path)
        io        = PDB.PDBIO()
        io.set_structure(structure)
        io.save(str(out_path), ChainSelector(chains))

        # Quick check — count atoms saved
        n_atoms = sum(
            1 for line in open(out_path)
            if line.startswith("ATOM")
        )
        print(f"  OK    {out_name}.pdb  "
              f"(chains {'+'.join(chains)}, {n_atoms} atoms)")
        return True

    except Exception as e:
        print(f"  FAIL  {out_name}: {e}")
        return False


def main():
    print("\n" + "="*55)
    print("  Script 02: Extracting binary chains")
    print("="*55 + "\n")

    ok = fail = 0
    for source, chains, name in EXTRACTIONS:
        if extract(source, chains, name):
            ok += 1
        else:
            fail += 1

    print(f"\n  Extracted: {ok}  |  Failed: {fail}")
    print("\n  Files saved to: 00-reference/pdb_structures/")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
