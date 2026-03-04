"""
Script 04_6: AF3 Validation — NSP7-NSP8
========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB structures:
  7BV2: Chain B=NSP8 (114 res, 78-191), Chain C=NSP7 (63 res, 2-64)
  6NUR: Chain B=NSP8 (115 res, 77-191), Chain C=NSP7 (70 res, 2-71)
AF3:  Chain A=NSP8 (198 res, 1-198), Chain B=NSP7 (83 res, 1-83)

Note: AF3 chains swapped vs PDB (A=NSP8, B=NSP7)
      Extra gate: fpocket druggability >= 0.5

Validation gate: F1 >= 0.70 at 5.0 A cutoff

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/04_validate_NSP7-NSP8_6.py
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
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP7-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP7-NSP8"
RES_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF           = 5.0
F1_GATE          = 0.70
FPOCKET_GATE     = 0.50

# PDB chain assignments
CHAIN_NSP8_PDB   = "B"
CHAIN_NSP7_PDB   = "C"
# AF3 chain assignments (swapped)
CHAIN_NSP8_AF3   = "A"
CHAIN_NSP7_AF3   = "B"

PDB_FILES = {
    "7BV2": PDB_DIR / "7BV2_NSP7-NSP8.pdb",
    "6NUR": PDB_DIR / "6NUR_NSP7-NSP8.pdb",
}
AF3_FILE = AF3_DIR / "NSP7_NSP8_best_model.pdb"


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


def align_map(pdb_seq, pdb_rns, af3_seq, af3_rns):
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
            if pi < len(pdb_rns) and \
                    ai < len(af3_rns):
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
    f1 = 2*p*r/(p+r) if (p+r) > 0 else 0.0
    return round(p,3), round(r,3), round(f1,3)


def validate(pdb_id, pdb_file, af3_s,
              af3_i8, af3_i7):
    parser = PDB.PDBParser(QUIET=True)
    pdb_s  = parser.get_structure(pdb_id, pdb_file)

    pdb_i8, pdb_i7 = get_interface_residues(
        pdb_s, CHAIN_NSP8_PDB, CHAIN_NSP7_PDB,
        CUTOFF)

    # Alignment maps: PDB→AF3
    af3_seq8, af3_rns8 = get_seq_resnums(
        af3_s, CHAIN_NSP8_AF3)
    af3_seq7, af3_rns7 = get_seq_resnums(
        af3_s, CHAIN_NSP7_AF3)
    pdb_seq8, pdb_rns8 = get_seq_resnums(
        pdb_s, CHAIN_NSP8_PDB)
    pdb_seq7, pdb_rns7 = get_seq_resnums(
        pdb_s, CHAIN_NSP7_PDB)

    map8 = align_map(pdb_seq8, pdb_rns8,
                     af3_seq8, af3_rns8)
    map7 = align_map(pdb_seq7, pdb_rns7,
                     af3_seq7, af3_rns7)

    pdb_i8_af3 = {map8[r] for r in pdb_i8
                  if r in map8}
    pdb_i7_af3 = {map7[r] for r in pdb_i7
                  if r in map7}

    _,_,f8       = f1_score(af3_i8, pdb_i8_af3)
    _,_,f7       = f1_score(af3_i7, pdb_i7_af3)
    all_af3      = af3_i8 | af3_i7
    all_pdb      = pdb_i8_af3 | pdb_i7_af3
    _,_,f_overall = f1_score(all_af3, all_pdb)

    status = "✅ PASS" if f_overall >= F1_GATE \
             else "❌ FAIL"
    print(f"  {pdb_id}: NSP8 F1={f8:.3f}  "
          f"NSP7 F1={f7:.3f}  "
          f"Overall F1={f_overall:.3f}  {status}")
    print(f"    PDB interface: NSP8={len(pdb_i8)} "
          f"NSP7={len(pdb_i7)} residues")
    print(f"    AF3 interface: NSP8={len(af3_i8)} "
          f"NSP7={len(af3_i7)} residues")

    return {
        "pdb_id":      pdb_id,
        "f1_nsp8":     f8,
        "f1_nsp7":     f7,
        "f1_overall":  f_overall,
        "pass":        f_overall >= F1_GATE,
        "n_pdb_nsp8":  len(pdb_i8),
        "n_pdb_nsp7":  len(pdb_i7),
        "n_af3_nsp8":  len(af3_i8),
        "n_af3_nsp7":  len(af3_i7),
    }


def main():
    print("\n" + "="*57)
    print("  Script 04_6: AF3 Validation — NSP7-NSP8")
    print(f"  Cutoff={CUTOFF} Å  Gate F1≥{F1_GATE}")
    print(f"  Extra gate: fpocket ≥ {FPOCKET_GATE}")
    print("  AF3: Chain A=NSP8, Chain B=NSP7 (swapped)")
    print("="*57 + "\n")

    # AF3 quality
    conf = json.load(
        open(AF3_DIR / "NSP7_NSP8_confidence.json"))
    print(f"  AF3 quality:")
    print(f"    iptm={conf.get('iptm')}  "
          f"ptm={conf.get('ptm')}  "
          f"ranking={conf.get('ranking_score')}")
    print(f"    has_clash={conf.get('has_clash')}  "
          f"disordered={conf.get('fraction_disordered')}")
    print(f"    PAE_min="
          f"{conf.get('chain_pair_pae_min')}\n")

    parser = PDB.PDBParser(QUIET=True)
    af3_s  = parser.get_structure("AF3", AF3_FILE)

    af3_i8, af3_i7 = get_interface_residues(
        af3_s, CHAIN_NSP8_AF3, CHAIN_NSP7_AF3,
        CUTOFF)

    results = []
    for pdb_id, pdb_file in PDB_FILES.items():
        results.append(
            validate(pdb_id, pdb_file,
                     af3_s, af3_i8, af3_i7))

    primary = results[0]
    gate    = "✅ PASS" if primary["pass"] \
              else "❌ FAIL"
    print(f"\n  Primary gate (7BV2): "
          f"F1={primary['f1_overall']:.3f}  {gate}")
    print(f"  fpocket gate: ≥{FPOCKET_GATE} "
          f"(checked in Script 07_6)")

    out = RES_DIR / "validation_result_6.json"
    with open(out,"w") as f:
        json.dump({
            "complex":         "NSP7-NSP8",
            "primary_pdb":     "7BV2",
            "cutoff_angstrom": CUTOFF,
            "f1_gate":         F1_GATE,
            "fpocket_gate":    FPOCKET_GATE,
            "primary_f1":      primary["f1_overall"],
            "primary_pass":    primary["pass"],
            "af3_quality":     {
                "iptm":    conf.get("iptm"),
                "ptm":     conf.get("ptm"),
                "ranking": conf.get("ranking_score"),
                "pae_min": conf.get("chain_pair_pae_min"),
            },
            "note": "AF3 chains swapped: A=NSP8, B=NSP7",
            "results": results,
        }, f, indent=2)
    print(f"\n  Saved: 02-validation/NSP7-NSP8/"
          f"validation_result_6.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
