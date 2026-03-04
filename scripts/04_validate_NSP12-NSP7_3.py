"""
Script 04_3: AF3 Validation — NSP12-NSP7
=========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2 (NSP12/NSP7/NSP8 trimer, SARS-CoV-2, 2.90 Å)
Chain A = NSP12, Chain B = NSP7, Chain C = NSP8 (ignored)

AF3 chains SWAPPED:
  AF3 Chain A = NSP7  (83 res)
  AF3 Chain B = NSP12 (932 res)

Method: BioPython PairwiseAligner sequence-alignment mapping
        (standard approach — handles numbering + indels)

Gate: F1 >= 0.70 at 5.0 Å cutoff

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/04_validate_NSP12-NSP7_3.py
"""

import json
import numpy as np
from pathlib import Path
from Bio import PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP7"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP7"
RES_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF   = 5.0
F1_GATE  = 0.70

# PDB chain assignments
PDB_CHAIN_NSP12 = "A"
PDB_CHAIN_NSP7  = "C"

# AF3 chain assignments (swapped)
AF3_CHAIN_NSP7  = "A"
AF3_CHAIN_NSP12 = "B"


def get_sequence(structure, chain_id):
    ppb = PPBuilder()
    for m in structure:
        for c in m:
            if c.id == chain_id:
                pp = ppb.build_peptides(c)
                seq = "".join(str(p.get_sequence()) for p in pp)
                res = [r for r in c if r.id[0] == " "]
                resnums = [r.id[1] for r in res]
                return seq, resnums
    return "", []


def build_alignment_map(pdb_seq, pdb_resnums,
                         af3_seq, af3_resnums):
    """
    Align PDB sequence to AF3 sequence.
    Returns dict: pdb_resnum -> af3_resnum
    """
    alns = pairwise2.align.globalms(
        pdb_seq, af3_seq,
        2, -1, -2, -0.5,
        one_alignment_only=True)
    if not alns:
        return {}

    aln_pdb, aln_af3 = alns[0].seqA, alns[0].seqB
    mapping = {}
    pdb_i, af3_i = 0, 0

    for p_aa, a_aa in zip(aln_pdb, aln_af3):
        if p_aa != "-" and a_aa != "-":
            if pdb_i < len(pdb_resnums) and \
               af3_i < len(af3_resnums):
                mapping[pdb_resnums[pdb_i]] = \
                    af3_resnums[af3_i]
            pdb_i += 1
            af3_i += 1
        elif p_aa == "-":
            af3_i += 1
        else:
            pdb_i += 1

    return mapping


def get_interface_residues(structure,
                            chain_a, chain_b,
                            cutoff=5.0):
    """Get residues within cutoff of the partner chain."""
    atoms_a = {r.id[1]: list(r.get_atoms())
               for m in structure for c in m
               if c.id == chain_a
               for r in c if r.id[0] == " "}
    atoms_b = {r.id[1]: list(r.get_atoms())
               for m in structure for c in m
               if c.id == chain_b
               for r in c if r.id[0] == " "}

    iface_a, iface_b = set(), set()
    for ra, ats_a in atoms_a.items():
        for rb, ats_b in atoms_b.items():
            found = False
            for a1 in ats_a:
                if found:
                    break
                for a2 in ats_b:
                    try:
                        if a1 - a2 <= cutoff:
                            iface_a.add(ra)
                            iface_b.add(rb)
                            found = True
                            break
                    except Exception:
                        pass
    return iface_a, iface_b


def compute_f1(pdb_iface, af3_iface_mapped,
               pdb_resnums, label=""):
    """
    Compute Precision/Recall/F1.
    pdb_iface: set of PDB resnums at interface
    af3_iface_mapped: set of PDB resnums (AF3 mapped back)
    pdb_resnums: all residues in chain (denominator)
    """
    tp = len(pdb_iface & af3_iface_mapped)
    fp = len(af3_iface_mapped - pdb_iface)
    fn = len(pdb_iface - af3_iface_mapped)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1        = (2 * precision * recall /
                 (precision + recall)
                 if (precision + recall) > 0 else 0)

    print(f"    {label}: TP={tp} FP={fp} FN={fn} "
          f"P={precision:.3f} R={recall:.3f} "
          f"F1={f1:.3f}")
    return precision, recall, f1


def validate_structure(pdb_file, pdb_label,
                        af3_structure,
                        af3_seq_nsp12, af3_resnums_nsp12,
                        af3_seq_nsp7,  af3_resnums_nsp7):
    """Validate one PDB structure against AF3."""
    parser = PDB.PDBParser(QUIET=True)
    pdb_s  = parser.get_structure(pdb_label, pdb_file)

    # Get PDB sequences + resnums
    pdb_seq_nsp12, pdb_rns_nsp12 = get_sequence(
        pdb_s, PDB_CHAIN_NSP12)
    pdb_seq_nsp7,  pdb_rns_nsp7  = get_sequence(
        pdb_s, PDB_CHAIN_NSP7)

    print(f"\n  {pdb_label}:")
    print(f"    PDB NSP12: {len(pdb_seq_nsp12)} aa, "
          f"res {pdb_rns_nsp12[0]}-{pdb_rns_nsp12[-1]}")
    print(f"    PDB NSP7:  {len(pdb_seq_nsp7)} aa, "
          f"res {pdb_rns_nsp7[0]}-{pdb_rns_nsp7[-1]}")

    # Build alignment maps
    map_nsp12 = build_alignment_map(
        pdb_seq_nsp12, pdb_rns_nsp12,
        af3_seq_nsp12, af3_resnums_nsp12)
    map_nsp7  = build_alignment_map(
        pdb_seq_nsp7,  pdb_rns_nsp7,
        af3_seq_nsp7,  af3_resnums_nsp7)

    print(f"    Alignment map: NSP12={len(map_nsp12)} "
          f"NSP7={len(map_nsp7)} residues")

    # PDB interface residues
    pdb_iface_nsp12, pdb_iface_nsp7 = \
        get_interface_residues(
            pdb_s, PDB_CHAIN_NSP12,
            PDB_CHAIN_NSP7, CUTOFF)

    print(f"    PDB interface: NSP12={len(pdb_iface_nsp12)} "
          f"NSP7={len(pdb_iface_nsp7)} residues")

    # AF3 interface residues (in AF3 local numbering)
    af3_iface_nsp12_local, af3_iface_nsp7_local = \
        get_interface_residues(
            af3_structure,
            AF3_CHAIN_NSP12, AF3_CHAIN_NSP7, CUTOFF)

    # Map AF3 interface back to PDB numbering
    rev_map_nsp12 = {v:k for k,v in map_nsp12.items()}
    rev_map_nsp7  = {v:k for k,v in map_nsp7.items()}

    af3_iface_nsp12_pdb = {
        rev_map_nsp12[r]
        for r in af3_iface_nsp12_local
        if r in rev_map_nsp12}
    af3_iface_nsp7_pdb  = {
        rev_map_nsp7[r]
        for r in af3_iface_nsp7_local
        if r in rev_map_nsp7}

    # F1 scores
    _, _, f1_nsp12 = compute_f1(
        pdb_iface_nsp12, af3_iface_nsp12_pdb,
        pdb_rns_nsp12, "NSP12")
    _, _, f1_nsp7  = compute_f1(
        pdb_iface_nsp7,  af3_iface_nsp7_pdb,
        pdb_rns_nsp7,  "NSP7 ")

    n_total = len(pdb_iface_nsp12) + len(pdb_iface_nsp7)
    tp_total = len(
        (pdb_iface_nsp12 & af3_iface_nsp12_pdb) |
        (pdb_iface_nsp7  & af3_iface_nsp7_pdb))
    fp_total = len(
        (af3_iface_nsp12_pdb - pdb_iface_nsp12) |
        (af3_iface_nsp7_pdb  - pdb_iface_nsp7))
    fn_total = len(
        (pdb_iface_nsp12 - af3_iface_nsp12_pdb) |
        (pdb_iface_nsp7  - af3_iface_nsp7_pdb))

    p_ov = tp_total/(tp_total+fp_total) \
        if (tp_total+fp_total) > 0 else 0
    r_ov = tp_total/(tp_total+fn_total) \
        if (tp_total+fn_total) > 0 else 0
    f1_overall = (2*p_ov*r_ov/(p_ov+r_ov)
                  if (p_ov+r_ov) > 0 else 0)

    gate = "✅ PASS" if f1_overall >= F1_GATE else "❌ FAIL"
    print(f"    Overall F1={f1_overall:.3f} {gate}")

    return {
        "pdb":           pdb_label,
        "f1_nsp12":      round(f1_nsp12, 3),
        "f1_nsp7":       round(f1_nsp7,  3),
        "f1_overall":    round(f1_overall, 3),
        "pass":          f1_overall >= F1_GATE,
        "iface_nsp12":   sorted(pdb_iface_nsp12),
        "iface_nsp7":    sorted(pdb_iface_nsp7),
        "af3_iface_nsp12": sorted(af3_iface_nsp12_local),
        "af3_iface_nsp7":  sorted(af3_iface_nsp7_local),
    }


def main():
    print("\n" + "="*57)
    print("  Script 04_3: AF3 Validation — NSP12-NSP7")
    print("  PDB: 7BV2 | Gate: F1 >= 0.70 @ 5.0 Å")
    print("="*57)

    parser = PDB.PDBParser(QUIET=True)

    # Load AF3
    af3_file = AF3_DIR / "NSP12_NSP7_best_model.pdb"
    af3_s    = parser.get_structure("AF3", af3_file)

    # AF3 sequences (chains swapped)
    af3_seq_nsp7,  af3_rns_nsp7  = get_sequence(
        af3_s, AF3_CHAIN_NSP7)
    af3_seq_nsp12, af3_rns_nsp12 = get_sequence(
        af3_s, AF3_CHAIN_NSP12)

    print(f"\n  AF3 model:")
    print(f"    Chain A (NSP7):  {len(af3_seq_nsp7)} aa")
    print(f"    Chain B (NSP12): {len(af3_seq_nsp12)} aa")
    print(f"    iptm=0.81  ptm=0.91  ranking=0.84")

    # Load confidence
    conf = json.load(open(AF3_DIR / "confidence.json"))
    print(f"    Confirmed: iptm={conf.get('iptm')} "
          f"ptm={conf.get('ptm')}")

    # Validate against all 3 PDB structures
    print(f"\n  Validating against PDB structures...")
    results = []

    # Use pre-extracted NSP12-NSP7 chain files
    pdb_files = {
        "7BV2": PDB_DIR / "7BV2_NSP12-NSP7.pdb",
        "6NUR": PDB_DIR / "6NUR_NSP12-NSP7.pdb",
        "7C2K": PDB_DIR / "7C2K_NSP12-NSP7.pdb",
    }
    for pdb_id, pdb_file in pdb_files.items():
        if not pdb_file.exists():
            print(f"\n  ⚠️  {pdb_file.name} not found — skipping")
            continue
        r = validate_structure(
            pdb_file, pdb_id,
            af3_s,
            af3_seq_nsp12, af3_rns_nsp12,
            af3_seq_nsp7,  af3_rns_nsp7)
        results.append(r)

    # Summary
    print(f"\n  {'PDB':<8} {'NSP12 F1':>9} "
          f"{'NSP7 F1':>8} {'Overall':>8}  Gate")
    print(f"  {'-'*8} {'-'*9} {'-'*8} {'-'*8}  {'-'*6}")
    for r in results:
        gate = "✅ PASS" if r["pass"] else "❌ FAIL"
        print(f"  {r['pdb']:<8} {r['f1_nsp12']:>9.3f} "
              f"{r['f1_nsp7']:>8.3f} "
              f"{r['f1_overall']:>8.3f}  {gate}")

    primary = results[0]
    print(f"\n  Primary gate (7BV2): "
          f"F1={primary['f1_overall']:.3f} "
          f"{'✅ PASS' if primary['pass'] else '❌ FAIL'}")

    # Save
    out = RES_DIR / "validation_result_3.json"
    with open(out, "w") as f:
        json.dump({
            "complex":       "NSP12-NSP7",
            "af3_iptm":      conf.get("iptm"),
            "af3_ptm":       conf.get("ptm"),
            "af3_ranking":   conf.get("ranking_score"),
            "cutoff_A":      CUTOFF,
            "f1_gate":       F1_GATE,
            "primary_pdb":   "7BV2",
            "primary_f1":    primary["f1_overall"],
            "primary_pass":  primary["pass"],
            "results":       results,
        }, f, indent=2)
    print(f"\n  Saved: 02-validation/NSP12-NSP7/"
          f"validation_result_3.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
