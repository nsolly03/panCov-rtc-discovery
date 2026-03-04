"""
Script 05_3: Interface Analysis — NSP12-NSP7
============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2_NSP12-NSP7 (primary), 6NUR_NSP12-NSP7, 7C2K_NSP12-NSP7 + AF3
Chain A = NSP12, Chain C = NSP7
AF3: Chain A = NSP7, Chain B = NSP12

Analyses:
  1. Salt bridges (< 4.0 Å, charged residue pairs)
  2. Hydrophobic contacts (< 4.5 Å, C-C)
  3. H-bond contacts (< 3.5 Å, donor-acceptor)
  4. Consensus hotspots (≥ 2 structures)
  5. Top 20 contact residues per chain

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/05_interface_NSP12-NSP7_3.py
"""

import json
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP7"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP7"
RES_DIR.mkdir(parents=True, exist_ok=True)

PDB_CHAIN_NSP12 = "A"
PDB_CHAIN_NSP7  = "C"
AF3_CHAIN_NSP12 = "B"
AF3_CHAIN_NSP7  = "A"

NSP12_OFFSET = 30   # PDB res 31 = local 1
NSP7_OFFSET  = 1    # PDB res 2  = local 1

CHARGED_POS = {"LYS","ARG","HIS"}
CHARGED_NEG = {"ASP","GLU"}
HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET",
               "PHE","TRP","TYR","PRO"}
HBOND_ATOMS = {"N","O","NZ","NH1","NH2","NE",
               "OG","OG1","OE1","OE2","ND1",
               "ND2","NE2","OH","OD1","OD2"}

STRUCTURES = {
    "7BV2": (PDB_DIR/"7BV2_NSP12-NSP7.pdb",
             PDB_CHAIN_NSP12, PDB_CHAIN_NSP7),
    "6NUR": (PDB_DIR/"6NUR_NSP12-NSP7.pdb",
             PDB_CHAIN_NSP12, PDB_CHAIN_NSP7),
    "7C2K": (PDB_DIR/"7C2K_NSP12-NSP7.pdb",
             PDB_CHAIN_NSP12, PDB_CHAIN_NSP7),
    "AF3":  (AF3_DIR/"NSP12_NSP7_best_model.pdb",
             AF3_CHAIN_NSP12, AF3_CHAIN_NSP7),
}


# Alignment maps: pdb_resnum -> af3_local
# Built once at startup, reused per structure
_ALIGN_MAPS = {}

def build_align_maps(struct_id, pdb_file,
                      ch_nsp12, ch_nsp7, af3_s):
    from Bio import pairwise2
    ppb = PPBuilder()
    parser2 = PDB.PDBParser(QUIET=True)
    pdb_s   = parser2.get_structure(struct_id, pdb_file)

    def get_seq_resnums(s, chain_id):
        for m in s:
            for c in m:
                if c.id == chain_id:
                    pp  = ppb.build_peptides(c)
                    seq = "".join(str(p.get_sequence())
                                  for p in pp)
                    rns = [r.id[1] for r in c
                           if r.id[0] == " "]
                    return seq, rns
        return "", []

    def align_map(pdb_seq, pdb_rns,
                   af3_seq, af3_rns):
        alns = pairwise2.align.globalms(
            pdb_seq, af3_seq,
            2,-1,-2,-0.5,
            one_alignment_only=True)
        if not alns:
            return {}
        aln_p, aln_a = alns[0].seqA, alns[0].seqB
        m = {}
        pi, ai = 0, 0
        for p_aa, a_aa in zip(aln_p, aln_a):
            if p_aa != "-" and a_aa != "-":
                if pi < len(pdb_rns) and                    ai < len(af3_rns):
                    m[pdb_rns[pi]] = af3_rns[ai]
                pi += 1; ai += 1
            elif p_aa == "-":
                ai += 1
            else:
                pi += 1
        return m

    # AF3 sequences
    af3_seq_7,  af3_rns_7  = get_seq_resnums(
        af3_s, AF3_CHAIN_NSP7)
    af3_seq_12, af3_rns_12 = get_seq_resnums(
        af3_s, AF3_CHAIN_NSP12)

    pdb_seq_12, pdb_rns_12 = get_seq_resnums(
        pdb_s, ch_nsp12)
    pdb_seq_7,  pdb_rns_7  = get_seq_resnums(
        pdb_s, ch_nsp7)

    m12 = align_map(pdb_seq_12, pdb_rns_12,
                    af3_seq_12, af3_rns_12)
    m7  = align_map(pdb_seq_7,  pdb_rns_7,
                    af3_seq_7,  af3_rns_7)

    _ALIGN_MAPS[struct_id] = {
        "nsp12": m12, "nsp7": m7}


def pdb_to_local(resnum, chain, struct_id):
    """Convert PDB resnum to AF3 local position."""
    if struct_id == "AF3":
        return resnum
    maps = _ALIGN_MAPS.get(struct_id, {})
    if chain in [PDB_CHAIN_NSP12, AF3_CHAIN_NSP12]:
        return maps.get("nsp12", {}).get(resnum, resnum)
    else:
        return maps.get("nsp7",  {}).get(resnum, resnum)


def analyze_interface(struct_id, pdb_file,
                       chain_nsp12, chain_nsp7):
    parser = PDB.PDBParser(QUIET=True)
    s      = parser.get_structure(struct_id, pdb_file)

    res_12 = {r.id[1]: r
              for m in s for c in m
              if c.id == chain_nsp12
              for r in c if r.id[0] == " "}
    res_7  = {r.id[1]: r
              for m in s for c in m
              if c.id == chain_nsp7
              for r in c if r.id[0] == " "}

    salt_bridges  = []
    hydrophobic   = []
    hbonds        = []
    contact_score = {}  # resnum -> score

    for rn12, r12 in res_12.items():
        rname12 = r12.resname.strip()
        for rn7, r7 in res_7.items():
            rname7 = r7.resname.strip()
            for a1 in r12.get_atoms():
                for a2 in r7.get_atoms():
                    try:
                        d = a1 - a2
                    except Exception:
                        continue

                    # Salt bridge
                    if d <= 4.0:
                        if ((rname12 in CHARGED_POS and
                             rname7  in CHARGED_NEG) or
                            (rname12 in CHARGED_NEG and
                             rname7  in CHARGED_POS)):
                            pair = (f"{rname12}"
                                    f"{pdb_to_local(rn12, chain_nsp12, struct_id)}",
                                    f"{rname7}"
                                    f"{pdb_to_local(rn7, chain_nsp7, struct_id)}",
                                    round(d, 2))
                            if pair not in salt_bridges:
                                salt_bridges.append(pair)
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+3
                            contact_score[rn7]  = \
                                contact_score.get(rn7, 0)+3

                    # H-bond
                    if (d <= 3.5 and
                        a1.name in HBOND_ATOMS and
                            a2.name in HBOND_ATOMS):
                        hbonds.append((
                            f"{rname12}"
                            f"{pdb_to_local(rn12, chain_nsp12, struct_id)}",
                            f"{rname7}"
                            f"{pdb_to_local(rn7, chain_nsp7, struct_id)}",
                            round(d,2)))
                        contact_score[rn12] = \
                            contact_score.get(rn12,0)+2
                        contact_score[rn7]  = \
                            contact_score.get(rn7, 0)+2

                    # Hydrophobic
                    if (d <= 4.5 and
                        rname12 in HYDROPHOBIC and
                        rname7  in HYDROPHOBIC and
                        a1.element == "C" and
                            a2.element == "C"):
                        hydrophobic.append((
                            f"{rname12}"
                            f"{pdb_to_local(rn12, chain_nsp12, struct_id)}",
                            f"{rname7}"
                            f"{pdb_to_local(rn7, chain_nsp7, struct_id)}",
                            round(d,2)))
                        contact_score[rn12] = \
                            contact_score.get(rn12,0)+1
                        contact_score[rn7]  = \
                            contact_score.get(rn7, 0)+1

    # Top 20 per chain
    nsp12_scores = {
        pdb_to_local(k, chain_nsp12, struct_id): v
        for k,v in contact_score.items()
        if k in res_12}
    nsp7_scores  = {
        pdb_to_local(k, chain_nsp7, struct_id): v
        for k,v in contact_score.items()
        if k in res_7}

    top20_12 = sorted(nsp12_scores.items(),
                      key=lambda x:-x[1])[:20]
    top20_7  = sorted(nsp7_scores.items(),
                      key=lambda x:-x[1])[:20]

    # Unique salt bridge pairs
    sb_unique = list({(a,b):d
                      for a,b,d in salt_bridges}.items())

    total = (len(salt_bridges)*3 +
             len(hbonds) +
             len(hydrophobic))

    print(f"  {struct_id}: SB={len(sb_unique)} "
          f"HB={len(set((a,b) for a,b,d in hbonds))} "
          f"HY={len(set((a,b) for a,b,d in hydrophobic))} "
          f"total={total}")

    if sb_unique:
        for (a,b),d in sb_unique[:6]:
            print(f"    SB: {a} — {b} ({d} Å)")

    return {
        "struct":       struct_id,
        "salt_bridges": [(a,b,d)
                         for (a,b),d in sb_unique],
        "n_hbonds":     len(set((a,b)
                               for a,b,d in hbonds)),
        "n_hydrophobic":len(set((a,b)
                               for a,b,d in hydrophobic)),
        "total_score":  total,
        "top20_NSP12":  top20_12,
        "top20_NSP7":   top20_7,
    }


def consensus_hotspots(results):
    """Residues appearing in ≥ 2 structures."""
    from collections import Counter
    nsp12_counts = Counter()
    nsp7_counts  = Counter()

    for r in results:
        nsp12_counts.update(
            [pos for pos,_ in r["top20_NSP12"]])
        nsp7_counts.update(
            [pos for pos,_ in r["top20_NSP7"]])

    hot12 = sorted(
        [p for p,c in nsp12_counts.items() if c >= 2],
        key=lambda x: -nsp12_counts[x])
    hot7  = sorted(
        [p for p,c in nsp7_counts.items()  if c >= 2],
        key=lambda x: -nsp7_counts[x])

    return hot12, hot7


def main():
    print("\n" + "="*57)
    print("  Script 05_3: Interface Analysis — NSP12-NSP7")
    print("="*57 + "\n")

    # Load AF3 structure for alignment maps
    parser_af3 = PDB.PDBParser(QUIET=True)
    af3_s = parser_af3.get_structure(
        "AF3", AF3_DIR / "NSP12_NSP7_best_model.pdb")

    results = []
    for sid, (pdb_file, ch12, ch7) in STRUCTURES.items():
        if sid != "AF3":
            build_align_maps(sid, pdb_file,
                             ch12, ch7, af3_s)
        results.append(
            analyze_interface(sid, pdb_file, ch12, ch7))

    hot12, hot7 = consensus_hotspots(results)

    print(f"\n  Consensus hotspots (≥ 2 structures):")
    print(f"  NSP12 ({len(hot12)}): {hot12[:15]}")
    print(f"  NSP7  ({len(hot7)}):  {hot7[:15]}")

    # Salt bridges confirmed across structures
    print(f"\n  Salt bridge summary:")
    all_sbs = {}
    for r in results:
        for a, b, d in r["salt_bridges"]:
            key = f"{a}—{b}"
            all_sbs.setdefault(key, []).append(
                (r["struct"], d))
    for pair, occurrences in sorted(
            all_sbs.items(),
            key=lambda x: -len(x[1])):
        structs = [s for s,_ in occurrences]
        print(f"    {pair:<30} "
              f"[{', '.join(structs)}]")

    # Save
    out = RES_DIR / "interface_analysis_3.json"
    with open(out, "w") as f:
        json.dump({
            "complex":          "NSP12-NSP7",
            "primary_pdb":      "7BV2",
            "structures":       [r["struct"]
                                 for r in results],
            "hotspots_NSP12":   hot12,
            "hotspots_NSP7":    hot7,
            "salt_bridges_all": all_sbs,
            "per_structure":    {
                r["struct"]: {
                    "salt_bridges": r["salt_bridges"],
                    "n_hbonds":     r["n_hbonds"],
                    "n_hydrophobic":r["n_hydrophobic"],
                    "total_score":  r["total_score"],
                    "top20_NSP12":  r["top20_NSP12"],
                    "top20_NSP7":   r["top20_NSP7"],
                } for r in results},
        }, f, indent=2)

    print(f"\n  Saved: 02-validation/NSP12-NSP7/"
          f"interface_analysis_3.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
