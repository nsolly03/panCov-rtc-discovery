"""
Script 04: Interface Residue Validation — NSP10-NSP14
======================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Compares interface residues between:
  - PDB crystal structure (7DIY) — ground truth
  - AF3 predicted structure    — our model

  For each structure it finds all residues within 5.0 Angstroms
  of the partner chain. Then it computes:
  - Precision : fraction of AF3 predictions that are correct
  - Recall    : fraction of true interface residues AF3 found
  - F1 score  : overall accuracy (main metric)

  Gate: F1 >= 0.70 required to proceed to next steps

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/04_validate_interface.py
"""

from pathlib import Path
from Bio import PDB
from Bio.PDB import NeighborSearch

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
RES_DIR.mkdir(parents=True, exist_ok=True)

# ── Configuration ──────────────────────────────────────────
COMPLEX     = "NSP10-NSP14"
PDB_FILE    = PDB_DIR / "7DIY.pdb"
AF3_FILE    = AF3_DIR / "NSP10-NSP14" / "NSP10_NSP14_best_model.pdb"
CUTOFF      = 5.0   # Angstroms

# Chain assignments
# PDB 7DIY : A=NSP10, B=NSP14
# AF3 model: A=NSP10, B=NSP14
PDB_CHAIN_A = "A"   # NSP10
PDB_CHAIN_B = "B"   # NSP14
AF3_CHAIN_A = "A"   # NSP10
AF3_CHAIN_B = "B"   # NSP14


def get_interface_residues(structure, chain_a, chain_b, cutoff):
    """
    Find all residues within cutoff Angstroms of the partner chain.
    Returns two sets: residues from chain_a and chain_b at interface.
    """
    atoms_a = []
    atoms_b = []

    for model in structure:
        for chain in model:
            if chain.id == chain_a:
                atoms_a.extend(list(chain.get_atoms()))
            elif chain.id == chain_b:
                atoms_b.extend(list(chain.get_atoms()))

    if not atoms_a or not atoms_b:
        return set(), set()

    # Search for atoms within cutoff
    ns_a = NeighborSearch(atoms_a)
    ns_b = NeighborSearch(atoms_b)

    iface_a = set()
    iface_b = set()

    # Chain A residues close to chain B
    for atom in atoms_a:
        neighbors = ns_b.search(atom.coord, cutoff, "A")
        if neighbors:
            res = atom.get_parent()
            iface_a.add((res.id[1], res.resname))

    # Chain B residues close to chain A
    for atom in atoms_b:
        neighbors = ns_a.search(atom.coord, cutoff, "A")
        if neighbors:
            res = atom.get_parent()
            iface_b.add((res.id[1], res.resname))

    return iface_a, iface_b


def compute_metrics(pdb_set, af3_set, chain_name):
    """Compute Precision, Recall, F1"""
    # Use only residue numbers for comparison
    pdb_res = {r[0] for r in pdb_set}
    af3_res = {r[0] for r in af3_set}

    tp = len(pdb_res & af3_res)   # in both
    fp = len(af3_res - pdb_res)   # in AF3 only
    fn = len(pdb_res - af3_res)   # in PDB only

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1        = (2 * precision * recall / (precision + recall)
                 if (precision + recall) > 0 else 0.0)

    return {
        "chain":     chain_name,
        "pdb_size":  len(pdb_res),
        "af3_size":  len(af3_res),
        "tp": tp, "fp": fp, "fn": fn,
        "precision": round(precision, 3),
        "recall":    round(recall, 3),
        "f1":        round(f1, 3),
        "pdb_res":   sorted(pdb_res),
        "af3_res":   sorted(af3_res),
        "only_pdb":  sorted(pdb_res - af3_res),
        "only_af3":  sorted(af3_res - pdb_res),
    }


def print_results(m, label):
    print(f"\n  [{label}]")
    print(f"  PDB interface residues : {m['pdb_size']}")
    print(f"  AF3 interface residues : {m['af3_size']}")
    print(f"  TP={m['tp']}  FP={m['fp']}  FN={m['fn']}")
    print(f"  Precision : {m['precision']:.3f}")
    print(f"  Recall    : {m['recall']:.3f}")
    print(f"  F1 Score  : {m['f1']:.3f}")


def main():
    print("\n" + "="*55)
    print(f"  Script 04: Interface Validation — {COMPLEX}")
    print(f"  PDB reference : 7DIY")
    print(f"  AF3 model     : NSP10_NSP14_best_model.pdb")
    print(f"  Cutoff        : {CUTOFF} Angstroms")
    print("="*55)

    # Load structures
    parser = PDB.PDBParser(QUIET=True)

    print("\n  Loading structures...")
    pdb_struct = parser.get_structure("pdb", PDB_FILE)
    af3_struct = parser.get_structure("af3", AF3_FILE)
    print("  OK    7DIY loaded")
    print("  OK    AF3 model loaded")

    # Get interface residues
    print(f"\n  Finding interface residues (cutoff={CUTOFF}A)...")
    pdb_a, pdb_b = get_interface_residues(
        pdb_struct, PDB_CHAIN_A, PDB_CHAIN_B, CUTOFF)
    af3_a, af3_b = get_interface_residues(
        af3_struct, AF3_CHAIN_A, AF3_CHAIN_B, CUTOFF)

    # Compute metrics
    m_a = compute_metrics(pdb_a, af3_a, "NSP10 (chain A)")
    m_b = compute_metrics(pdb_b, af3_b, "NSP14 (chain B)")

    # Overall F1 = average of both chains
    overall_f1 = round((m_a["f1"] + m_b["f1"]) / 2, 3)

    # Print results
    print("\n" + "="*55)
    print("  RESULTS")
    print("="*55)
    print_results(m_a, "NSP10 chain A")
    print_results(m_b, "NSP14 chain B")

    print(f"\n  {'='*55}")
    print(f"  OVERALL F1 SCORE : {overall_f1:.3f}")

    if overall_f1 >= 0.70:
        gate = "✅ PASS — proceed to interface analysis"
    elif overall_f1 >= 0.50:
        gate = "⚠️  CAUTION — proceed carefully"
    else:
        gate = "❌ FAIL — do not proceed, check AF3 model"

    print(f"  GATE STATUS      : {gate}")
    print(f"  {'='*55}")

    # Print interface residues
    print(f"\n  NSP10 interface residues (PDB 7DIY):")
    print(f"  {m_a['pdb_res']}")
    print(f"\n  NSP14 interface residues (PDB 7DIY):")
    print(f"  {m_b['pdb_res']}")
    print(f"\n  Residues in PDB but missed by AF3 (NSP10):")
    print(f"  {m_a['only_pdb']}")
    print(f"\n  Residues in PDB but missed by AF3 (NSP14):")
    print(f"  {m_b['only_pdb']}")

    # Save results
    import json
    results = {
        "complex":     COMPLEX,
        "pdb":         "7DIY",
        "cutoff_A":    CUTOFF,
        "iptm":        0.89,
        "ptm":         0.88,
        "NSP10":       m_a,
        "NSP14":       m_b,
        "overall_f1":  overall_f1,
        "gate":        gate,
    }
    out = RES_DIR / "validation_result.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: 02-validation/NSP10-NSP14/validation_result.json")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
