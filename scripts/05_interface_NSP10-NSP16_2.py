"""
Script 05_2: Interface Analysis — NSP10-NSP16
=============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

PDB structures: 6W4H (primary, 1.80Å), 6WVN (1.95Å), 6WKQ (2.05Å)
AF3 model     : iptm=0.79, ranking=0.84
Chain layout  : PDB Chain A=NSP16, Chain B=NSP10
               AF3 Chain A=NSP10,  Chain B=NSP16

Contacts analyzed: H-bonds, salt bridges, hydrophobic, VdW
Output: hotspot residues + consensus across all 4 structures

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/05_interface_NSP10-NSP16_2.py
"""

import json
import numpy as np
from pathlib import Path
from Bio import PDB, Align
from Bio.PDB.Polypeptide import PPBuilder

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3" / "NSP10-NSP16"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
RES_DIR.mkdir(parents=True, exist_ok=True)

PDB_STRUCTURES = [
    {"file": PDB_DIR/"6W4H.pdb", "label":"6W4H", "res":"1.80Å",
     "chain_nsp10":"B", "chain_nsp16":"A"},
    {"file": PDB_DIR/"6WVN.pdb", "label":"6WVN", "res":"1.95Å",
     "chain_nsp10":"B", "chain_nsp16":"A"},
    {"file": PDB_DIR/"6WKQ.pdb", "label":"6WKQ", "res":"2.05Å",
     "chain_nsp10":"B", "chain_nsp16":"A"},
]
AF3_FILE        = AF3_DIR / "NSP10_NSP16_best_model.pdb"
CHAIN_NSP10_AF3 = "A"
CHAIN_NSP16_AF3 = "B"

HBOND_DIST   = 3.5
SALTBR_DIST  = 4.0
HYDRO_DIST   = 4.5
VDW_DIST     = 5.0
CHARGED_POS  = {"LYS","ARG","HIS"}
CHARGED_NEG  = {"ASP","GLU"}
HYDRO_AA     = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","TYR","PRO"}
HBOND_DONORS    = {"N","O","NZ","NH1","NH2","NE","OG","OG1",
                   "OE1","OE2","ND1","ND2","NE2"}
HBOND_ACCEPTORS = {"O","OD1","OD2","OE1","OE2","OG","OG1","OH","N","NZ"}


# ── Sequence utilities ─────────────────────────────────────

def get_sequence(structure, chain_id):
    resnums, seq = [], []
    for m in structure:
        for c in m:
            if c.id != chain_id:
                continue
            for r in sorted(c.get_residues(), key=lambda x: x.id[1]):
                if r.id[0] != " ":
                    continue
                try:
                    aa = PDB.Polypeptide.index_to_one(
                         PDB.Polypeptide.three_to_index(r.resname))
                    resnums.append(r.id[1])
                    seq.append(aa)
                except Exception:
                    pass
    return resnums, "".join(seq)


def build_resmap(pdb_resnums, pdb_seq, af3_seq):
    """Map PDB residue numbers → AF3 local positions via pairwise alignment."""
    aligner = Align.PairwiseAligner()
    aligner.mode             = "global"
    aligner.match_score      = 2
    aligner.mismatch_score   = -1
    aligner.open_gap_score   = -3
    aligner.extend_gap_score = -0.5
    aln     = next(iter(aligner.align(pdb_seq, af3_seq)))
    mapping = {}
    pi, ai  = 0, 0
    for col in zip(*aln):
        pa, aa = col
        if pa != "-" and aa != "-":
            if pi < len(pdb_resnums):
                mapping[pdb_resnums[pi]] = ai + 1
            pi += 1; ai += 1
        elif pa == "-":
            ai += 1
        else:
            pi += 1
    return mapping


# ── Contact analysis ───────────────────────────────────────

def analyze_contacts(structure,
                     chain_nsp10, chain_nsp16,
                     map_nsp10=None, map_nsp16=None,
                     label=""):
    """
    Analyze all contact types across the NSP10-NSP16 interface.
    Returns per-residue scores and salt bridge list.
    """
    scores_nsp10  = {}
    scores_nsp16  = {}
    salt_bridges  = []
    contact_counts = {"HBOND":0,"SALT_BRIDGE":0,
                      "HYDROPHOBIC":0,"VDW":0}

    # Collect residues
    res_nsp10 = {r.id[1]: r for m in structure
                 for c in m if c.id == chain_nsp10
                 for r in c if r.id[0] == " "}
    res_nsp16 = {r.id[1]: r for m in structure
                 for c in m if c.id == chain_nsp16
                 for r in c if r.id[0] == " "}

    atoms_nsp10 = [(rn, a) for rn, r in res_nsp10.items()
                   for a in r.get_atoms()]
    atoms_nsp16 = [(rn, a) for rn, r in res_nsp16.items()
                   for a in r.get_atoms()]

    for rn10, a10 in atoms_nsp10:
        for rn16, a16 in atoms_nsp16:
            try:
                dist = a10 - a16
            except Exception:
                continue
            if dist > VDW_DIST:
                continue

            res10_name = res_nsp10[rn10].resname.strip()
            res16_name = res_nsp16[rn16].resname.strip()

            # Map to AF3 local positions
            pos10 = (map_nsp10.get(rn10, rn10)
                     if map_nsp10 else rn10)
            pos16 = (map_nsp16.get(rn16, rn16)
                     if map_nsp16 else rn16)
            key10 = f"{res10_name}{pos10}"
            key16 = f"{res16_name}{pos16}"

            score = 0

            # Salt bridge
            if dist <= SALTBR_DIST:
                if ((res10_name in CHARGED_POS and
                     res16_name in CHARGED_NEG) or
                    (res10_name in CHARGED_NEG and
                     res16_name in CHARGED_POS)):
                    score += 3
                    contact_counts["SALT_BRIDGE"] += 1
                    pair = tuple(sorted([key10, key16]))
                    if not any(s["res_a"] == pair[0] and
                                s["res_b"] == pair[1]
                                for s in salt_bridges):
                        salt_bridges.append({
                            "res_a": pair[0],
                            "res_b": pair[1],
                            "type": "SALT_BRIDGE",
                            "distance": round(dist, 2)})

            # H-bond
            if (dist <= HBOND_DIST and
                a10.name in HBOND_DONORS and
                    a16.name in HBOND_ACCEPTORS):
                score += 2
                contact_counts["HBOND"] += 1

            # Hydrophobic
            if (dist <= HYDRO_DIST and
                res10_name in HYDRO_AA and
                    res16_name in HYDRO_AA and
                    a10.element == "C" and
                    a16.element == "C"):
                score += 1.5
                contact_counts["HYDROPHOBIC"] += 1

            # VdW
            if dist <= VDW_DIST:
                score += 0.5
                contact_counts["VDW"] += 1

            if score > 0:
                scores_nsp10[key10] = (scores_nsp10.get(key10, 0)
                                       + score)
                scores_nsp16[key16] = (scores_nsp16.get(key16, 0)
                                       + score)

    top_nsp10 = sorted(scores_nsp10.items(),
                       key=lambda x: x[1], reverse=True)[:20]
    top_nsp16 = sorted(scores_nsp16.items(),
                       key=lambda x: x[1], reverse=True)[:20]
    total     = sum(contact_counts.values())

    return {
        "source":         label,
        "total_contacts": total,
        "contact_counts": contact_counts,
        "top20_NSP10":    top_nsp10,
        "top20_NSP16":    top_nsp16,
        "salt_bridges":   salt_bridges,
    }


def consensus_hotspots(all_results, top_n=10):
    """Residues appearing in top20 of ≥2 structures."""
    from collections import defaultdict
    votes10 = defaultdict(int)
    votes16 = defaultdict(int)
    for r in all_results.values():
        for res, _ in r["top20_NSP10"]:
            votes10[res] += 1
        for res, _ in r["top20_NSP16"]:
            votes16[res] += 1
    cons10 = sorted([k for k,v in votes10.items() if v>=2],
                    key=lambda x: -votes10[x])[:top_n]
    cons16 = sorted([k for k,v in votes16.items() if v>=2],
                    key=lambda x: -votes16[x])[:top_n]
    return cons10, cons16


def main():
    print("\n" + "="*57)
    print("  Script 05_2: Interface Analysis — NSP10-NSP16")
    print("="*57 + "\n")

    parser_pdb = PDB.PDBParser(QUIET=True)
    all_results = {}

    # ── Load AF3 + get sequences ───────────────────────────
    af3_struct = parser_pdb.get_structure("af3", AF3_FILE)
    _, af3_seq_nsp10 = get_sequence(af3_struct, CHAIN_NSP10_AF3)
    _, af3_seq_nsp16 = get_sequence(af3_struct, CHAIN_NSP16_AF3)

    # ── PDB structures ─────────────────────────────────────
    for pdb_info in PDB_STRUCTURES:
        lbl = pdb_info["label"]
        if not pdb_info["file"].exists():
            print(f"  SKIP {lbl} — file not found")
            continue
        print(f"  Analyzing {lbl} ({pdb_info['res']})...")
        s = parser_pdb.get_structure(lbl, pdb_info["file"])
        rn10, seq10 = get_sequence(s, pdb_info["chain_nsp10"])
        rn16, seq16 = get_sequence(s, pdb_info["chain_nsp16"])
        map10 = build_resmap(rn10, seq10, af3_seq_nsp10)
        map16 = build_resmap(rn16, seq16, af3_seq_nsp16)
        res = analyze_contacts(s,
                               pdb_info["chain_nsp10"],
                               pdb_info["chain_nsp16"],
                               map_nsp10=map10,
                               map_nsp16=map16,
                               label=lbl)
        all_results[f"PDB_{lbl}"] = res
        print(f"    Total contacts : {res['total_contacts']}")
        print(f"    Salt bridges   : {len(res['salt_bridges'])}")
        for sb in res["salt_bridges"]:
            print(f"      {sb['res_a']} — {sb['res_b']} "
                  f"({sb['distance']} Å)")
        print(f"    Top 5 NSP10    : "
              f"{[x[0] for x in res['top20_NSP10'][:5]]}")
        print(f"    Top 5 NSP16    : "
              f"{[x[0] for x in res['top20_NSP16'][:5]]}")
        print()

    # ── AF3 ───────────────────────────────────────────────
    print("  Analyzing AF3 model...")
    res_af3 = analyze_contacts(af3_struct,
                               CHAIN_NSP10_AF3,
                               CHAIN_NSP16_AF3,
                               label="AF3_model")
    all_results["AF3_model"] = res_af3
    print(f"    Total contacts : {res_af3['total_contacts']}")
    print(f"    Salt bridges   : {len(res_af3['salt_bridges'])}")
    for sb in res_af3["salt_bridges"]:
        print(f"      {sb['res_a']} — {sb['res_b']} "
              f"({sb['distance']} Å)")
    print(f"    Top 5 NSP10    : "
          f"{[x[0] for x in res_af3['top20_NSP10'][:5]]}")
    print(f"    Top 5 NSP16    : "
          f"{[x[0] for x in res_af3['top20_NSP16'][:5]]}")
    print()

    # ── Consensus ─────────────────────────────────────────
    cons10, cons16 = consensus_hotspots(all_results)
    all_results["consensus"] = {
        "NSP10_hotspots": cons10,
        "NSP16_hotspots": cons16,
    }

    print("  ── Consensus hotspots (≥2 structures) ──")
    print(f"  NSP10: {cons10}")
    print(f"  NSP16: {cons16}")
    print()

    # ── All salt bridges ───────────────────────────────────
    all_sb = {}
    for k, r in all_results.items():
        if "salt_bridges" in r:
            for sb in r["salt_bridges"]:
                pair = f"{sb['res_a']}|{sb['res_b']}"
                all_sb.setdefault(pair, []).append(k)
    print("  ── Salt bridges across structures ──")
    for pair, sources in all_sb.items():
        print(f"  {pair:<35} present in: "
              f"{[s.replace('PDB_','') for s in sources]}")

    # ── Save ──────────────────────────────────────────────
    import numpy as np

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.floating, np.float32, np.float64)):
                return float(obj)
            if isinstance(obj, (np.integer,)):
                return int(obj)
            return super().default(obj)

    out = RES_DIR / "interface_analysis.json"
    with open(out, "w") as f:
        json.dump(all_results, f, indent=2, cls=NumpyEncoder)
    print(f"\n  Saved: 02-validation/NSP10-NSP16/interface_analysis.json")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
