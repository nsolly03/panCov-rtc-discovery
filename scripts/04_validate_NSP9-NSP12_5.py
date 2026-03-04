"""
Script 04_5: AF3 Validation — NSP9-NSP12
=========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB structures:
  8SQK: Chain A=NSP12 (929 res, 1-929), Chain G=NSP9 (113 res, 1-113)
AF3:  Chain A=NSP12 (932 res), Chain B=NSP9 (113 res)

Note: Only 1 PDB structure available (8SQK)
      Chain G in PDB maps to Chain B in AF3
      PAE=3.52 A — higher uncertainty, novel interface

Validation gate: F1 >= 0.70 at 5.0 A cutoff

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/04_validate_NSP9-NSP12_5.py
"""

import json
import warnings
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2
warnings.filterwarnings("ignore")

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP9-NSP12"
RES_DIR  = PROJECT / "02-validation" / "NSP9-NSP12"
RES_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF           = 5.0
F1_GATE          = 0.70
CHAIN_NSP12_PDB  = "A"
CHAIN_NSP9_PDB   = "G"   # unusual chain in 8SQK
CHAIN_NSP12_AF3  = "A"
CHAIN_NSP9_AF3   = "B"

PDB_FILES = {
    "8SQK": PDB_DIR / "8SQK_NSP9-NSP12.pdb",
}
AF3_FILE = AF3_DIR / "NSP9_NSP12_best_model.pdb"


def get_interface_residues(structure, chain1_id,
                            chain2_id, cutoff):
    res1 = [r for m in structure for c in m
            if c.id == chain1_id
            for r in c if r.id[0] == " "]
    res2 = [r for m in structure for c in m
            if c.id == chain2_id
            for r in c if r.id[0] == " "]
    iface1, iface2 = set(), set()
    for r1 in res1:
        for r2 in res2:
            for a1 in r1.get_atoms():
                for a2 in r2.get_atoms():
                    try:
                        if a1 - a2 <= cutoff:
                            iface1.add(r1.id[1])
                            iface2.add(r2.id[1])
                    except Exception:
                        pass
    return iface1, iface2


def get_seq_resnums(structure, chain_id):
    ppb = PPBuilder()
    for m in structure:
        for c in m:
            if c.id == chain_id:
                pp  = ppb.build_peptides(c)
                seq = "".join(str(p.get_sequence())
                              for p in pp)
                rns = [r.id[1] for r in c
                       if r.id[0] == " "]
                return seq, rns
    return "", []


def align_residue_map(pdb_seq, pdb_rns,
                       af3_seq, af3_rns):
    alns = pairwise2.align.globalms(
        pdb_seq, af3_seq, 2,-1,-2,-0.5,
        one_alignment_only=True)
    if not alns:
        return {}
    m  = {}
    pi, ai = 0, 0
    for p_aa, a_aa in zip(alns[0].seqA,
                           alns[0].seqB):
        if p_aa != "-" and a_aa != "-":
            if pi < len(pdb_rns) and ai < len(af3_rns):
                m[pdb_rns[pi]] = af3_rns[ai]
            pi += 1; ai += 1
        elif p_aa == "-":
            ai += 1
        else:
            pi += 1
    return m


def f1_score(pred, true):
    if not pred or not true:
        return 0.0, 0.0, 0.0
    tp = len(pred & true)
    p  = tp/len(pred) if pred else 0.0
    r  = tp/len(true) if true else 0.0
    f1 = (2*p*r/(p+r)) if (p+r) > 0 else 0.0
    return round(p,3), round(r,3), round(f1,3)


def validate(pdb_id, pdb_file, af3_structure,
              af3_i12, af3_i9):
    parser = PDB.PDBParser(QUIET=True)
    pdb_s  = parser.get_structure(pdb_id, pdb_file)

    pdb_i12, pdb_i9 = get_interface_residues(
        pdb_s, CHAIN_NSP12_PDB, CHAIN_NSP9_PDB,
        CUTOFF)

    # Alignment maps
    af3_seq12, af3_rns12 = get_seq_resnums(
        af3_structure, CHAIN_NSP12_AF3)
    af3_seq9,  af3_rns9  = get_seq_resnums(
        af3_structure, CHAIN_NSP9_AF3)
    pdb_seq12, pdb_rns12 = get_seq_resnums(
        pdb_s, CHAIN_NSP12_PDB)
    pdb_seq9,  pdb_rns9  = get_seq_resnums(
        pdb_s, CHAIN_NSP9_PDB)

    map12 = align_residue_map(
        pdb_seq12, pdb_rns12, af3_seq12, af3_rns12)
    map9  = align_residue_map(
        pdb_seq9,  pdb_rns9,  af3_seq9,  af3_rns9)

    pdb_i12_af3 = {map12[r] for r in pdb_i12
                   if r in map12}
    pdb_i9_af3  = {map9[r]  for r in pdb_i9
                   if r in map9}

    p12,r12,f12 = f1_score(af3_i12, pdb_i12_af3)
    p9, r9, f9  = f1_score(af3_i9,  pdb_i9_af3)

    all_af3 = af3_i12 | af3_i9
    all_pdb = pdb_i12_af3 | pdb_i9_af3
    _,_,f_overall = f1_score(all_af3, all_pdb)

    status = "✅ PASS" if f_overall >= F1_GATE \
             else "❌ FAIL"
    print(f"  {pdb_id}: NSP12 F1={f12:.3f}  "
          f"NSP9 F1={f9:.3f}  "
          f"Overall F1={f_overall:.3f}  {status}")
    print(f"    PDB interface: NSP12={len(pdb_i12)} "
          f"NSP9={len(pdb_i9)} residues")
    print(f"    AF3 interface: NSP12={len(af3_i12)} "
          f"NSP9={len(af3_i9)} residues")

    return {
        "pdb_id":     pdb_id,
        "f1_nsp12":   f12,
        "f1_nsp9":    f9,
        "f1_overall": f_overall,
        "pass":       f_overall >= F1_GATE,
        "n_af3":      len(all_af3),
        "n_pdb":      len(all_pdb),
    }


def main():
    print("\n" + "="*57)
    print("  Script 04_5: AF3 Validation — NSP9-NSP12")
    print(f"  Cutoff={CUTOFF} Å  Gate F1≥{F1_GATE}")
    print("  Note: single PDB (8SQK), Chain G=NSP9")
    print("  PAE=3.52 Å — novel interface, higher uncertainty")
    print("="*57 + "\n")

    parser = PDB.PDBParser(QUIET=True)
    af3_s  = parser.get_structure("AF3", AF3_FILE)

    af3_i12, af3_i9 = get_interface_residues(
        af3_s, CHAIN_NSP12_AF3, CHAIN_NSP9_AF3, CUTOFF)

    results = []
    for pdb_id, pdb_file in PDB_FILES.items():
        results.append(
            validate(pdb_id, pdb_file,
                     af3_s, af3_i12, af3_i9))

    primary = results[0]
    gate    = "✅ PASS" if primary["pass"] \
              else "❌ FAIL"
    print(f"\n  Primary gate (8SQK): "
          f"F1={primary['f1_overall']:.3f}  {gate}")

    out = RES_DIR / "validation_result_5.json"
    with open(out,"w") as f:
        json.dump({
            "complex":         "NSP9-NSP12",
            "primary_pdb":     "8SQK",
            "cutoff_angstrom": CUTOFF,
            "f1_gate":         F1_GATE,
            "primary_f1":      primary["f1_overall"],
            "primary_pass":    primary["pass"],
            "note":            "single PDB, chain G=NSP9, "
                               "PAE=3.52A novel interface",
            "results":         results,
        }, f, indent=2)
    print(f"\n  Saved: 02-validation/NSP9-NSP12/"
          f"validation_result_5.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
