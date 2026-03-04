"""
Script 05_5: Interface Analysis — NSP9-NSP12
============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 8SQK (primary, only available)
  Chain A = NSP12, Chain G = NSP9
AF3: Chain A = NSP12, Chain B = NSP9

Note: Single PDB — hotspots from 8SQK + AF3 only

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/05_interface_NSP9-NSP12_5.py
"""

import json
import warnings
import numpy as np
from pathlib import Path
from Bio import PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
warnings.filterwarnings("ignore")

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP9-NSP12"
RES_DIR  = PROJECT / "02-validation" / "NSP9-NSP12"
RES_DIR.mkdir(parents=True, exist_ok=True)

CHAIN_NSP12_PDB = "A"
CHAIN_NSP9_PDB  = "G"
CHAIN_NSP12_AF3 = "A"
CHAIN_NSP9_AF3  = "B"

CHARGED_POS = {"LYS","ARG","HIS"}
CHARGED_NEG = {"ASP","GLU"}
HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET",
               "PHE","TRP","TYR","PRO"}
HBOND_ATOMS = {"N","O","NZ","NH1","NH2","NE",
               "OG","OG1","OE1","OE2","ND1",
               "ND2","NE2","OH","OD1","OD2"}
CHARGED_POS_ATOMS = {
    "LYS":["NZ"],"ARG":["NH1","NH2","NE"],
    "HIS":["ND1","NE2"]}
CHARGED_NEG_ATOMS = {
    "ASP":["OD1","OD2"],"GLU":["OE1","OE2"]}
SB_CUTOFF = 5.0

STRUCTURES = {
    "8SQK": (PDB_DIR/"8SQK_NSP9-NSP12.pdb",
              CHAIN_NSP12_PDB, CHAIN_NSP9_PDB),
    "AF3":  (AF3_DIR/"NSP9_NSP12_best_model.pdb",
              CHAIN_NSP12_AF3, CHAIN_NSP9_AF3),
}

_ALIGN_MAP_8SQK = {}


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


def build_align_map(pdb_s, af3_s):
    af3_seq12, af3_rns12 = get_seq_resnums(
        af3_s, CHAIN_NSP12_AF3)
    af3_seq9,  af3_rns9  = get_seq_resnums(
        af3_s, CHAIN_NSP9_AF3)
    pdb_seq12, pdb_rns12 = get_seq_resnums(
        pdb_s, CHAIN_NSP12_PDB)
    pdb_seq9,  pdb_rns9  = get_seq_resnums(
        pdb_s, CHAIN_NSP9_PDB)

    def amap(pseq, prns, aseq, arns):
        alns = pairwise2.align.globalms(
            pseq, aseq, 2,-1,-2,-0.5,
            one_alignment_only=True)
        if not alns:
            return {}
        m  = {}
        pi, ai = 0, 0
        for p_aa, a_aa in zip(alns[0].seqA,
                               alns[0].seqB):
            if p_aa != "-" and a_aa != "-":
                if pi < len(prns) and ai < len(arns):
                    m[prns[pi]] = arns[ai]
                pi += 1; ai += 1
            elif p_aa == "-":
                ai += 1
            else:
                pi += 1
        return m

    _ALIGN_MAP_8SQK["nsp12"] = amap(
        pdb_seq12, pdb_rns12, af3_seq12, af3_rns12)
    _ALIGN_MAP_8SQK["nsp9"]  = amap(
        pdb_seq9,  pdb_rns9,  af3_seq9,  af3_rns9)


def pdb_to_local(resnum, chain, struct_id):
    if struct_id == "AF3":
        return resnum
    if chain == CHAIN_NSP12_PDB:
        return _ALIGN_MAP_8SQK.get(
            "nsp12",{}).get(resnum, resnum)
    return _ALIGN_MAP_8SQK.get(
        "nsp9", {}).get(resnum, resnum)


def analyze_interface(struct_id, pdb_file,
                       ch_nsp12, ch_nsp9):
    parser = PDB.PDBParser(QUIET=True)
    s      = parser.get_structure(struct_id, pdb_file)

    res_12 = {r.id[1]: r
              for m in s for c in m
              if c.id == ch_nsp12
              for r in c if r.id[0] == " "}
    res_9  = {r.id[1]: r
              for m in s for c in m
              if c.id == ch_nsp9
              for r in c if r.id[0] == " "}

    salt_bridges = []
    hbonds       = []
    hydrophobic  = []
    contact_score = {}
    sb_seen = set()
    hb_seen = set()
    hy_seen = set()

    # Salt bridges
    for rn12, r12 in res_12.items():
        rname12 = r12.resname.strip()
        for rn9, r9 in res_9.items():
            rname9  = r9.resname.strip()
            pos_at12 = CHARGED_POS_ATOMS.get(rname12,[])
            neg_at12 = CHARGED_NEG_ATOMS.get(rname12,[])
            pos_at9  = CHARGED_POS_ATOMS.get(rname9, [])
            neg_at9  = CHARGED_NEG_ATOMS.get(rname9, [])
            pairs = []
            if pos_at12 and neg_at9:
                pairs += [(a,b) for a in pos_at12
                          for b in neg_at9]
            if neg_at12 and pos_at9:
                pairs += [(a,b) for a in neg_at12
                          for b in pos_at9]
            for an1, an2 in pairs:
                try:
                    d = r12[an1] - r9[an2]
                    if d <= SB_CUTOFF:
                        l12 = pdb_to_local(
                            rn12, ch_nsp12, struct_id)
                        l9  = pdb_to_local(
                            rn9,  ch_nsp9,  struct_id)
                        key = (l12, l9)
                        if key not in sb_seen:
                            sb_seen.add(key)
                            salt_bridges.append((
                                f"{rname12}{l12}",
                                f"{rname9}{l9}",
                                float(round(d,2))))
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+3
                            contact_score[rn9]  = \
                                contact_score.get(rn9, 0)+3
                except Exception:
                    pass

    # H-bonds + hydrophobic
    for rn12, r12 in res_12.items():
        rname12 = r12.resname.strip()
        for rn9, r9 in res_9.items():
            rname9 = r9.resname.strip()
            for a1 in r12.get_atoms():
                for a2 in r9.get_atoms():
                    try:
                        d = a1 - a2
                    except Exception:
                        continue
                    l12 = pdb_to_local(
                        rn12, ch_nsp12, struct_id)
                    l9  = pdb_to_local(
                        rn9,  ch_nsp9,  struct_id)
                    if (d <= 3.5 and
                        a1.name in HBOND_ATOMS and
                        a2.name in HBOND_ATOMS):
                        key = (l12, l9)
                        if key not in hb_seen:
                            hb_seen.add(key)
                            hbonds.append((
                                f"{rname12}{l12}",
                                f"{rname9}{l9}",
                                float(round(d,2))))
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+2
                            contact_score[rn9]  = \
                                contact_score.get(rn9, 0)+2
                    if (d <= 4.5 and
                        rname12 in HYDROPHOBIC and
                        rname9  in HYDROPHOBIC and
                        a1.element == "C" and
                        a2.element == "C"):
                        key = (l12, l9)
                        if key not in hy_seen:
                            hy_seen.add(key)
                            hydrophobic.append((
                                f"{rname12}{l12}",
                                f"{rname9}{l9}",
                                float(round(d,2))))
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+1
                            contact_score[rn9]  = \
                                contact_score.get(rn9, 0)+1

    nsp12_scores = {
        pdb_to_local(k, ch_nsp12, struct_id): v
        for k,v in contact_score.items()
        if k in res_12}
    nsp9_scores = {
        pdb_to_local(k, ch_nsp9, struct_id): v
        for k,v in contact_score.items()
        if k in res_9}

    top20_12 = sorted(nsp12_scores.items(),
                      key=lambda x:-x[1])[:20]
    top20_9  = sorted(nsp9_scores.items(),
                      key=lambda x:-x[1])[:20]
    total    = (len(salt_bridges)*3 +
                len(hbonds) +
                len(hydrophobic))

    print(f"  {struct_id}: SB={len(salt_bridges)} "
          f"HB={len(hbonds)} "
          f"HY={len(hydrophobic)} total={total}")
    for a,b,d in salt_bridges:
        print(f"    SB: {a} — {b} ({d} Å)")

    return {
        "struct":        struct_id,
        "salt_bridges":  salt_bridges,
        "n_hbonds":      len(hbonds),
        "n_hydrophobic": len(hydrophobic),
        "total_score":   total,
        "top20_NSP12":   top20_12,
        "top20_NSP9":    top20_9,
    }


def consensus_hotspots(results):
    from collections import Counter
    c12 = Counter()
    c9  = Counter()
    for r in results:
        c12.update([p for p,_ in r["top20_NSP12"]])
        c9.update( [p for p,_ in r["top20_NSP9"]])
    # With only 2 structures lower threshold to ≥1
    hot12 = sorted([p for p,cnt in c12.items()
                    if cnt >= 1],
                   key=lambda x:-c12[x])
    hot9  = sorted([p for p,cnt in c9.items()
                    if cnt >= 1],
                   key=lambda x:-c9[x])
    # Consensus ≥2 (present in both)
    con12 = sorted([p for p,cnt in c12.items()
                    if cnt >= 2],
                   key=lambda x:-c12[x])
    con9  = sorted([p for p,cnt in c9.items()
                    if cnt >= 2],
                   key=lambda x:-c9[x])
    return hot12, hot9, con12, con9


def main():
    print("\n" + "="*57)
    print("  Script 05_5: Interface Analysis — NSP9-NSP12")
    print("  Note: single PDB (8SQK) + AF3")
    print("="*57 + "\n")

    parser = PDB.PDBParser(QUIET=True)
    af3_s  = parser.get_structure(
        "AF3", AF3_DIR/"NSP9_NSP12_best_model.pdb")
    pdb_s  = parser.get_structure(
        "8SQK", PDB_DIR/"8SQK_NSP9-NSP12.pdb")

    print("  Building alignment maps...")
    build_align_map(pdb_s, af3_s)

    results = []
    for sid,(pdb_f,ch12,ch9) in STRUCTURES.items():
        results.append(
            analyze_interface(sid,pdb_f,ch12,ch9))

    hot12,hot9,con12,con9 = consensus_hotspots(results)
    print(f"\n  All hotspots (8SQK + AF3 union):")
    print(f"  NSP12 ({len(hot12)}): {hot12[:20]}")
    print(f"  NSP9  ({len(hot9)}):  {hot9[:20]}")
    print(f"\n  Consensus hotspots (both structures):")
    print(f"  NSP12 ({len(con12)}): {con12[:15]}")
    print(f"  NSP9  ({len(con9)}):  {con9[:15]}")

    print(f"\n  Salt bridge summary:")
    all_sbs = {}
    for r in results:
        for a,b,d in r["salt_bridges"]:
            key = f"{a}—{b}"
            all_sbs.setdefault(key,[]).append(
                (r["struct"],d))
    if all_sbs:
        for pair,occ in sorted(
                all_sbs.items(),
                key=lambda x:-len(x[1])):
            structs  = [s for s,_ in occ]
            best_d   = min(d for _,d in occ)
            print(f"    {pair:<35} "
                  f"[{', '.join(structs)}] "
                  f"min={best_d:.2f}Å")
    else:
        print("    No salt bridges detected")

    class NpEncoder(json.JSONEncoder):
        def default(self, o):
            if isinstance(o, (np.floating,
                               np.float32,
                               np.float64)):
                return float(o)
            if isinstance(o, np.integer):
                return int(o)
            return super().default(o)

    out = RES_DIR / "interface_analysis_5.json"
    with open(out,"w") as f:
        json.dump({
            "complex":         "NSP9-NSP12",
            "primary_pdb":     "8SQK",
            "note":            "single PDB + AF3 only",
            "hotspots_NSP12":  con12 if con12 else hot12,
            "hotspots_NSP9":   con9  if con9  else hot9,
            "salt_bridges_all":all_sbs,
            "per_structure": {
                r["struct"]: {
                    "salt_bridges":  r["salt_bridges"],
                    "n_hbonds":      r["n_hbonds"],
                    "n_hydrophobic": r["n_hydrophobic"],
                    "total_score":   r["total_score"],
                    "top20_NSP12":   r["top20_NSP12"],
                    "top20_NSP9":    r["top20_NSP9"],
                } for r in results},
        }, f, indent=2, cls=NpEncoder)

    print(f"\n  Saved: 02-validation/NSP9-NSP12/"
          f"interface_analysis_5.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
