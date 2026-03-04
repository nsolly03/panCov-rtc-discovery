"""
Script 05_4: Interface Analysis — NSP12-NSP8
============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB: 7BV2 (primary), 6NUR, 7C2K + AF3
Chain A = NSP12, Chain B = NSP8 (all structures)

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/05_interface_NSP12-NSP8_4.py
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
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP12-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP12-NSP8"
RES_DIR.mkdir(parents=True, exist_ok=True)

CHAIN_NSP12_PDB = "A"
CHAIN_NSP8_PDB  = "B"
CHAIN_NSP12_AF3 = "A"
CHAIN_NSP8_AF3  = "B"

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
    "7BV2": (PDB_DIR/"7BV2_NSP12-NSP8.pdb",
             CHAIN_NSP12_PDB, CHAIN_NSP8_PDB),
    "6NUR": (PDB_DIR/"6NUR_NSP12-NSP8.pdb",
             CHAIN_NSP12_PDB, CHAIN_NSP8_PDB),
    "7C2K": (PDB_DIR/"7C2K_NSP12-NSP8.pdb",
             CHAIN_NSP12_PDB, CHAIN_NSP8_PDB),
    "AF3":  (AF3_DIR/"NSP12_NSP8_best_model.pdb",
             CHAIN_NSP12_AF3, CHAIN_NSP8_AF3),
}

_ALIGN_MAPS = {}


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


def build_align_maps(struct_id, pdb_file,
                      ch12, ch8, af3_s):
    parser2  = PDB.PDBParser(QUIET=True)
    pdb_s    = parser2.get_structure(struct_id, pdb_file)

    af3_seq12, af3_rns12 = get_seq_resnums(
        af3_s, CHAIN_NSP12_AF3)
    af3_seq8,  af3_rns8  = get_seq_resnums(
        af3_s, CHAIN_NSP8_AF3)
    pdb_seq12, pdb_rns12 = get_seq_resnums(pdb_s, ch12)
    pdb_seq8,  pdb_rns8  = get_seq_resnums(pdb_s, ch8)

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

    _ALIGN_MAPS[struct_id] = {
        "nsp12": amap(pdb_seq12, pdb_rns12,
                      af3_seq12, af3_rns12),
        "nsp8":  amap(pdb_seq8,  pdb_rns8,
                      af3_seq8,  af3_rns8),
    }


def pdb_to_local(resnum, chain, struct_id):
    if struct_id == "AF3":
        return resnum
    maps = _ALIGN_MAPS.get(struct_id, {})
    if chain == CHAIN_NSP12_PDB:
        return maps.get("nsp12", {}).get(resnum, resnum)
    else:
        return maps.get("nsp8",  {}).get(resnum, resnum)


def analyze_interface(struct_id, pdb_file,
                       chain_nsp12, chain_nsp8):
    parser = PDB.PDBParser(QUIET=True)
    s      = parser.get_structure(struct_id, pdb_file)

    res_12 = {r.id[1]: r
              for m in s for c in m
              if c.id == chain_nsp12
              for r in c if r.id[0] == " "}
    res_8  = {r.id[1]: r
              for m in s for c in m
              if c.id == chain_nsp8
              for r in c if r.id[0] == " "}

    salt_bridges  = []
    hydrophobic   = []
    hbonds        = []
    contact_score = {}
    sb_seen       = set()

    # Pass 1: salt bridges
    for rn12, r12 in res_12.items():
        rname12 = r12.resname.strip()
        for rn8, r8 in res_8.items():
            rname8 = r8.resname.strip()
            pos_at12 = CHARGED_POS_ATOMS.get(rname12,[])
            neg_at12 = CHARGED_NEG_ATOMS.get(rname12,[])
            pos_at8  = CHARGED_POS_ATOMS.get(rname8, [])
            neg_at8  = CHARGED_NEG_ATOMS.get(rname8, [])
            pairs = []
            if pos_at12 and neg_at8:
                pairs += [(a,b) for a in pos_at12
                          for b in neg_at8]
            if neg_at12 and pos_at8:
                pairs += [(a,b) for a in neg_at12
                          for b in pos_at8]
            for an1, an2 in pairs:
                try:
                    d = r12[an1] - r8[an2]
                    if d <= SB_CUTOFF:
                        l12 = pdb_to_local(
                            rn12, chain_nsp12, struct_id)
                        l8  = pdb_to_local(
                            rn8,  chain_nsp8,  struct_id)
                        key = (l12, l8)
                        if key not in sb_seen:
                            sb_seen.add(key)
                            salt_bridges.append((
                                f"{rname12}{l12}",
                                f"{rname8}{l8}",
                                float(round(d,2))))
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+3
                            contact_score[rn8]  = \
                                contact_score.get(rn8, 0)+3
                except Exception:
                    pass

    # Pass 2: H-bonds + hydrophobic
    hb_seen = set()
    hy_seen = set()
    for rn12, r12 in res_12.items():
        rname12 = r12.resname.strip()
        for rn8, r8 in res_8.items():
            rname8 = r8.resname.strip()
            for a1 in r12.get_atoms():
                for a2 in r8.get_atoms():
                    try:
                        d = a1 - a2
                    except Exception:
                        continue
                    l12 = pdb_to_local(
                        rn12, chain_nsp12, struct_id)
                    l8  = pdb_to_local(
                        rn8,  chain_nsp8,  struct_id)
                    if (d <= 3.5 and
                        a1.name in HBOND_ATOMS and
                        a2.name in HBOND_ATOMS):
                        key = (l12, l8)
                        if key not in hb_seen:
                            hb_seen.add(key)
                            hbonds.append((
                                f"{rname12}{l12}",
                                f"{rname8}{l8}",
                                float(round(d,2))))
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+2
                            contact_score[rn8]  = \
                                contact_score.get(rn8, 0)+2
                    if (d <= 4.5 and
                        rname12 in HYDROPHOBIC and
                        rname8  in HYDROPHOBIC and
                        a1.element == "C" and
                        a2.element == "C"):
                        key = (l12, l8)
                        if key not in hy_seen:
                            hy_seen.add(key)
                            hydrophobic.append((
                                f"{rname12}{l12}",
                                f"{rname8}{l8}",
                                float(round(d,2))))
                            contact_score[rn12] = \
                                contact_score.get(rn12,0)+1
                            contact_score[rn8]  = \
                                contact_score.get(rn8, 0)+1

    nsp12_scores = {
        pdb_to_local(k, chain_nsp12, struct_id): v
        for k,v in contact_score.items()
        if k in res_12}
    nsp8_scores  = {
        pdb_to_local(k, chain_nsp8, struct_id): v
        for k,v in contact_score.items()
        if k in res_8}

    top20_12 = sorted(nsp12_scores.items(),
                      key=lambda x:-x[1])[:20]
    top20_8  = sorted(nsp8_scores.items(),
                      key=lambda x:-x[1])[:20]

    total = (len(salt_bridges)*3 +
             len(hbonds) +
             len(hydrophobic))

    print(f"  {struct_id}: SB={len(salt_bridges)} "
          f"HB={len(hbonds)} "
          f"HY={len(hydrophobic)} total={total}")
    for a,b,d in salt_bridges[:6]:
        print(f"    SB: {a} — {b} ({d} Å)")

    return {
        "struct":       struct_id,
        "salt_bridges": salt_bridges,
        "n_hbonds":     len(hbonds),
        "n_hydrophobic":len(hydrophobic),
        "total_score":  total,
        "top20_NSP12":  top20_12,
        "top20_NSP8":   top20_8,
    }


def consensus_hotspots(results):
    from collections import Counter
    nsp12_counts = Counter()
    nsp8_counts  = Counter()
    for r in results:
        nsp12_counts.update(
            [pos for pos,_ in r["top20_NSP12"]])
        nsp8_counts.update(
            [pos for pos,_ in r["top20_NSP8"]])
    hot12 = sorted(
        [p for p,c in nsp12_counts.items() if c >= 2],
        key=lambda x: -nsp12_counts[x])
    hot8  = sorted(
        [p for p,c in nsp8_counts.items()  if c >= 2],
        key=lambda x: -nsp8_counts[x])
    return hot12, hot8


def main():
    print("\n" + "="*57)
    print("  Script 05_4: Interface Analysis — NSP12-NSP8")
    print("="*57 + "\n")

    parser = PDB.PDBParser(QUIET=True)
    af3_s  = parser.get_structure(
        "AF3", AF3_DIR / "NSP12_NSP8_best_model.pdb")

    results = []
    for sid, (pdb_f, ch12, ch8) in STRUCTURES.items():
        if sid != "AF3":
            build_align_maps(sid, pdb_f,
                             ch12, ch8, af3_s)
        results.append(
            analyze_interface(sid, pdb_f, ch12, ch8))

    hot12, hot8 = consensus_hotspots(results)
    print(f"\n  Consensus hotspots (≥ 2 structures):")
    print(f"  NSP12 ({len(hot12)}): {hot12[:15]}")
    print(f"  NSP8  ({len(hot8)}):  {hot8[:15]}")

    print(f"\n  Salt bridge summary:")
    all_sbs = {}
    for r in results:
        for a,b,d in r["salt_bridges"]:
            key = f"{a}—{b}"
            all_sbs.setdefault(key,[]).append(
                (r["struct"],d))
    for pair, occ in sorted(all_sbs.items(),
                             key=lambda x:-len(x[1])):
        structs = [s for s,_ in occ]
        best_d  = min(d for _,d in occ)
        print(f"    {pair:<35} "
              f"[{', '.join(structs)}] "
              f"min={best_d:.2f}Å")

    out = RES_DIR / "interface_analysis_4.json"
    import numpy as np
    class NpEncoder(json.JSONEncoder):
        def default(self, o):
            if isinstance(o, (np.floating, np.float32,
                               np.float64)):
                return float(o)
            if isinstance(o, np.integer):
                return int(o)
            return super().default(o)

    with open(out,"w") as f:
        json.dump({
            "complex":          "NSP12-NSP8",
            "primary_pdb":      "7BV2",
            "hotspots_NSP12":   hot12,
            "hotspots_NSP8":    hot8,
            "salt_bridges_all": all_sbs,
            "per_structure": {
                r["struct"]: {
                    "salt_bridges":  r["salt_bridges"],
                    "n_hbonds":      r["n_hbonds"],
                    "n_hydrophobic": r["n_hydrophobic"],
                    "total_score":   r["total_score"],
                    "top20_NSP12":   r["top20_NSP12"],
                    "top20_NSP8":    r["top20_NSP8"],
                } for r in results},
        }, f, indent=2, cls=NpEncoder)

    print(f"\n  Saved: 02-validation/NSP12-NSP8/"
          f"interface_analysis_4.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
