"""
Script 04_2: AF3 Validation — NSP10-NSP16
==========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Validates AF3 predicted interface against PDB crystal structure.
Primary PDB: 6W4H (SARS-CoV-2, 1.80 Å) — Chain B=NSP10, Chain A=NSP16
Secondary:   6WVN, 6WKQ
Gate: F1 >= 0.70 required to proceed.

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/04_validate_NSP10-NSP16_2.py
"""

import json
import numpy as np
from pathlib import Path
from Bio import PDB

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP10-NSP16"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
RES_DIR.mkdir(parents=True, exist_ok=True)

# PDB: Chain A=NSP16 (genome 6798-7096), Chain B=NSP10 (genome 4271-4386)
# AF3: Chain A=NSP10 (local 1-138),     Chain B=NSP16 (local 1-298)
# Sequence alignment shows:
#   PDB NSP10 starts at AF3 position 23 (missing 22 N-terminal residues)
#   PDB NSP16 starts at AF3 position 1  (full length, genome offset only)
# NSP10_PDB_OFFSET: PDB local 1 = AF3 position 23, so add 22
# NSP16_PDB_OFFSET: genome offset only (PDB 6798 -> local 1)
PDB_STRUCTURES = [
    {"file": PDB_DIR / "6W4H.pdb", "label": "6W4H", "res": "1.80Å",
     "chain_nsp10": "B", "chain_nsp16": "A",
     "nsp10_shift": 22, "nsp16_shift": 0},
    {"file": PDB_DIR / "6WVN.pdb", "label": "6WVN", "res": "1.95Å",
     "chain_nsp10": "B", "chain_nsp16": "A",
     "nsp10_shift": 22, "nsp16_shift": 0},
    {"file": PDB_DIR / "6WKQ.pdb", "label": "6WKQ", "res": "2.05Å",
     "chain_nsp10": "B", "chain_nsp16": "A",
     "nsp10_shift": 22, "nsp16_shift": 0},
]
AF3_FILE        = AF3_DIR / "NSP10_NSP16_best_model.pdb"
CHAIN_NSP10_AF3 = "A"   # AF3 Chain A = NSP10
CHAIN_NSP16_AF3 = "B"   # AF3 Chain B = NSP16
CUTOFF          = 5.0
F1_GATE         = 0.70


def get_sequence(structure, chain_id):
    """Extract sequence and ordered residue numbers from a chain."""
    from Bio.PDB.Polypeptide import PPBuilder
    ppb = PPBuilder()
    resnums, seq = [], []
    for m in structure:
        for c in m:
            if c.id != chain_id:
                continue
            for r in sorted(c.get_residues(),
                            key=lambda x: x.id[1]):
                if r.id[0] != " ":
                    continue
                try:
                    aa = PDB.Polypeptide.index_to_one(
                         PDB.Polypeptide.three_to_index(
                         r.resname))
                    resnums.append(r.id[1])
                    seq.append(aa)
                except Exception:
                    pass
    return resnums, "".join(seq)


def build_resnum_to_af3(pdb_resnums, pdb_seq, af3_seq):
    """
    Pairwise align PDB seq to AF3 seq.
    Returns dict: pdb_resnum -> af3_local_position (1-indexed).
    """
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.mode            = "global"
    aligner.match_score     = 2
    aligner.mismatch_score  = -1
    aligner.open_gap_score  = -3
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(pdb_seq, af3_seq)
    aln = next(iter(alignments))

    # Walk alignment columns
    mapping = {}
    pdb_i, af3_i = 0, 0
    for col in zip(*aln):
        pdb_aa, af3_aa = col
        if pdb_aa != "-" and af3_aa != "-":
            if pdb_i < len(pdb_resnums):
                mapping[pdb_resnums[pdb_i]] = af3_i + 1
            pdb_i += 1
            af3_i += 1
        elif pdb_aa == "-":
            af3_i += 1
        else:
            pdb_i += 1
    return mapping


def get_interface_residues(structure, chain_a, chain_b, cutoff,
                            af3_seq_a=None, af3_seq_b=None):
    """
    Return interface residues mapped to AF3 local numbering.
    af3_seq_a/b: AF3 sequences for sequence-alignment-based mapping.
    If None, uses raw residue numbers (AF3 already local).
    """
    # Build residue number mappings
    resnums_a, seq_a = get_sequence(structure, chain_a)
    resnums_b, seq_b = get_sequence(structure, chain_b)

    if af3_seq_a and seq_a:
        map_a = build_resnum_to_af3(resnums_a, seq_a, af3_seq_a)
    else:
        map_a = {r: r for r in resnums_a}

    if af3_seq_b and seq_b:
        map_b = build_resnum_to_af3(resnums_b, seq_b, af3_seq_b)
    else:
        map_b = {r: r for r in resnums_b}

    residues_a, residues_b = set(), set()
    atoms_a = [a for m in structure for c in m
               if c.id == chain_a
               for r in c if r.id[0] == " "
               for a in r]
    atoms_b = [a for m in structure for c in m
               if c.id == chain_b
               for r in c if r.id[0] == " "
               for a in r]
    for a1 in atoms_a:
        for a2 in atoms_b:
            try:
                if (a1 - a2) <= cutoff:
                    r_a = a1.get_parent().id[1]
                    r_b = a2.get_parent().id[1]
                    if r_a in map_a:
                        residues_a.add(map_a[r_a])
                    if r_b in map_b:
                        residues_b.add(map_b[r_b])
            except Exception:
                continue
    return residues_a, residues_b


def compute_f1(predicted, reference):
    if not predicted or not reference:
        return 0.0, 0.0, 0.0
    tp = len(predicted & reference)
    precision = tp / len(predicted) if predicted else 0
    recall    = tp / len(reference) if reference else 0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return round(precision, 4), round(recall, 4), round(f1, 4)


def load_confidence():
    cf = AF3_DIR / "NSP10_NSP16_confidence.json"
    return json.load(open(cf)) if cf.exists() else {}


def main():
    print("\n" + "="*55)
    print("  Script 04_2: AF3 Validation — NSP10-NSP16")
    print(f"  PDB: 6W4H (SARS-CoV-2, 1.80 Å)")
    print(f"  Cutoff: {CUTOFF} Å  |  F1 gate: {F1_GATE}")
    print("="*55 + "\n")

    parser = PDB.PDBParser(QUIET=True)

    # ── AF3 confidence ────────────────────────────────────
    conf = load_confidence()
    print("  AF3 Confidence scores:")
    for k in ["iptm","ptm","ranking_score",
              "has_clash","fraction_disordered"]:
        print(f"    {k:<22}: {conf.get(k, 'N/A')}")

    # ── Load AF3 once ─────────────────────────────────────
    print(f"\n  Loading AF3: {AF3_FILE.name}")
    af3_struct  = parser.get_structure("af3", AF3_FILE)
    # Extract AF3 sequences for alignment mapping
    _, af3_seq_nsp10 = get_sequence(af3_struct, CHAIN_NSP10_AF3)
    _, af3_seq_nsp16 = get_sequence(af3_struct, CHAIN_NSP16_AF3)
    print(f"  AF3 NSP10 seq length: {len(af3_seq_nsp10)}")
    print(f"  AF3 NSP16 seq length: {len(af3_seq_nsp16)}")
    # AF3 interface uses identity mapping (already local)
    af3_nsp10, af3_nsp16 = get_interface_residues(
        af3_struct, CHAIN_NSP10_AF3, CHAIN_NSP16_AF3, CUTOFF)
    print(f"  AF3 interface: NSP10={len(af3_nsp10)} "
          f"NSP16={len(af3_nsp16)} residues")

    # ── Validate against each PDB ─────────────────────────
    all_results = {}
    print(f"\n  {'PDB':<8} {'Res':<7} "
          f"{'NSP10 F1':>9} {'NSP16 F1':>9} "
          f"{'Overall F1':>11}  Gate")
    print(f"  {'-'*8} {'-'*7} {'-'*9} {'-'*9} {'-'*11}  {'-'*6}")

    for pdb_info in PDB_STRUCTURES:
        lbl  = pdb_info["label"]
        res  = pdb_info["res"]
        if not pdb_info["file"].exists():
            print(f"  {lbl:<8} MISSING")
            continue
        pdb_struct = parser.get_structure(lbl, pdb_info["file"])
        pdb_nsp10, pdb_nsp16 = get_interface_residues(
            pdb_struct,
            pdb_info["chain_nsp10"],
            pdb_info["chain_nsp16"], CUTOFF,
            af3_seq_a=af3_seq_nsp10,
            af3_seq_b=af3_seq_nsp16)

        p10,r10,f10 = compute_f1(af3_nsp10, pdb_nsp10)
        p16,r16,f16 = compute_f1(af3_nsp16, pdb_nsp16)
        _,_,f_all   = compute_f1(
            af3_nsp10|af3_nsp16, pdb_nsp10|pdb_nsp16)

        gate = "✅ PASS" if f_all >= F1_GATE else "❌ FAIL"
        print(f"  {lbl:<8} {res:<7} {f10:>9.3f} "
              f"{f16:>9.3f} {f_all:>11.3f}  {gate}")

        all_results[lbl] = {
            "nsp10_f1": f10, "nsp16_f1": f16,
            "overall_f1": f_all,
            "gate_passed": f_all >= F1_GATE,
            "pdb_nsp10": sorted(pdb_nsp10),
            "pdb_nsp16": sorted(pdb_nsp16),
        }

    # Primary = 6W4H (highest resolution)
    primary = all_results.get("6W4H", {})
    f_primary = primary.get("overall_f1", 0)
    gate_primary = "✅ PASS" if f_primary >= F1_GATE else "❌ FAIL"
    print(f"\n  Primary gate (6W4H): F1={f_primary:.3f}  {gate_primary}")

    # ── Save results ──────────────────────────────────────
    result = {
        "complex":          "NSP10-NSP16",
        "primary_pdb":      "6W4H",
        "all_pdbs":         ["6W4H","6WVN","6WKQ"],
        "cutoff_ang":       CUTOFF,
        "f1_gate":          F1_GATE,
        "confidence":       conf,
        "overall_f1":       f_primary,
        "gate_passed":      f_primary >= F1_GATE,
        "af3_nsp10_residues": sorted(af3_nsp10),
        "af3_nsp16_residues": sorted(af3_nsp16),
        "per_pdb":          all_results,
    }

    out = RES_DIR / "validation_result_2.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved: 02-validation/NSP10-NSP16/validation_result_2.json")
    print("="*55 + "\n")

    if not result["gate_passed"]:
        print("  ⚠️  F1 below gate — do NOT proceed to Script 05_2")


if __name__ == "__main__":
    main()
