"""
Script 05_6: Interface Analysis — NSP7-NSP8 (Dual Binding Mode)
===============================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

TWO NON-OVERLAPPING INTERFACES IDENTIFIED:
  Mode A (Crystal): NSP8 C-terminal (162-182) — NSP7 loop (24-27)
                    Observed in 7BV2, 6NUR crystal structures
  Mode B (AF3):     NSP8 N-terminal helix (80-121) — NSP7 body (2-79)
                    Predicted by AF3 (iptm=0.87, PAE=0.90 A)

Hypothesis: conformational flexibility allows two binding orientations
            Mode A = crystal packing artifact OR low-affinity state
            Mode B = physiological interface (larger, druggable)

PDB: 7BV2 (Chain B=NSP8, Chain C=NSP7) primary
     6NUR (Chain B=NSP8, Chain C=NSP7) secondary
AF3: Chain A=NSP8, Chain B=NSP7

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/05_interface_NSP7-NSP8_6.py
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
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP7-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP7-NSP8"
RES_DIR.mkdir(parents=True, exist_ok=True)

CHAIN_NSP8_PDB = "B"
CHAIN_NSP7_PDB = "C"
CHAIN_NSP8_AF3 = "A"
CHAIN_NSP7_AF3 = "B"

CHARGED_POS      = {"LYS","ARG","HIS"}
CHARGED_NEG      = {"ASP","GLU"}
HYDROPHOBIC      = {"ALA","VAL","ILE","LEU","MET",
                    "PHE","TRP","TYR","PRO"}
HBOND_ATOMS      = {"N","O","NZ","NH1","NH2","NE",
                    "OG","OG1","OE1","OE2","ND1",
                    "ND2","NE2","OH","OD1","OD2"}
CHARGED_POS_ATOMS = {
    "LYS":["NZ"],"ARG":["NH1","NH2","NE"],
    "HIS":["ND1","NE2"]}
CHARGED_NEG_ATOMS = {
    "ASP":["OD1","OD2"],"GLU":["OE1","OE2"]}

STRUCTURES = {
    "7BV2": (PDB_DIR/"7BV2_NSP7-NSP8.pdb",
              CHAIN_NSP8_PDB, CHAIN_NSP7_PDB,
              "crystal"),
    "6NUR": (PDB_DIR/"6NUR_NSP7-NSP8.pdb",
              CHAIN_NSP8_PDB, CHAIN_NSP7_PDB,
              "crystal"),
    "AF3":  (AF3_DIR/"NSP7_NSP8_best_model.pdb",
              CHAIN_NSP8_AF3, CHAIN_NSP7_AF3,
              "af3"),
}


def analyze_interface(struct_id, pdb_file,
                       ch_nsp8, ch_nsp7, mode):
    parser = PDB.PDBParser(QUIET=True)
    s      = parser.get_structure(struct_id, pdb_file)

    res_8 = {r.id[1]: r
             for m in s for c in m
             if c.id == ch_nsp8
             for r in c if r.id[0] == " "}
    res_7 = {r.id[1]: r
             for m in s for c in m
             if c.id == ch_nsp7
             for r in c if r.id[0] == " "}

    salt_bridges = []
    hbonds       = []
    hydrophobic  = []
    sb_seen = set(); hb_seen = set(); hy_seen = set()
    contact_8 = {}; contact_7 = {}

    for rn8, r8 in res_8.items():
        rname8 = r8.resname.strip()
        for rn7, r7 in res_7.items():
            rname7  = r7.resname.strip()
            # Salt bridges
            pos_at8 = CHARGED_POS_ATOMS.get(rname8,[])
            neg_at8 = CHARGED_NEG_ATOMS.get(rname8,[])
            pos_at7 = CHARGED_POS_ATOMS.get(rname7,[])
            neg_at7 = CHARGED_NEG_ATOMS.get(rname7,[])
            pairs = []
            if pos_at8 and neg_at7:
                pairs += [(a,b) for a in pos_at8
                          for b in neg_at7]
            if neg_at8 and pos_at7:
                pairs += [(a,b) for a in neg_at8
                          for b in pos_at7]
            for an1,an2 in pairs:
                try:
                    d = r8[an1] - r7[an2]
                    if d <= 5.0:
                        key = (rn8,rn7)
                        if key not in sb_seen:
                            sb_seen.add(key)
                            salt_bridges.append((
                                f"{rname8}{rn8}",
                                f"{rname7}{rn7}",
                                float(round(d,2))))
                            contact_8[rn8] = \
                                contact_8.get(rn8,0)+3
                            contact_7[rn7] = \
                                contact_7.get(rn7,0)+3
                except Exception:
                    pass
            # H-bonds + hydrophobic
            for a1 in r8.get_atoms():
                for a2 in r7.get_atoms():
                    try:
                        d = a1 - a2
                    except Exception:
                        continue
                    if (d <= 3.5 and
                        a1.name in HBOND_ATOMS and
                        a2.name in HBOND_ATOMS):
                        key = (rn8,rn7)
                        if key not in hb_seen:
                            hb_seen.add(key)
                            hbonds.append((
                                f"{rname8}{rn8}",
                                f"{rname7}{rn7}",
                                float(round(d,2))))
                            contact_8[rn8] = \
                                contact_8.get(rn8,0)+2
                            contact_7[rn7] = \
                                contact_7.get(rn7,0)+2
                    if (d <= 4.5 and
                        rname8 in HYDROPHOBIC and
                        rname7 in HYDROPHOBIC and
                        a1.element == "C" and
                        a2.element == "C"):
                        key = (rn8,rn7)
                        if key not in hy_seen:
                            hy_seen.add(key)
                            hydrophobic.append((
                                f"{rname8}{rn8}",
                                f"{rname7}{rn7}",
                                float(round(d,2))))
                            contact_8[rn8] = \
                                contact_8.get(rn8,0)+1
                            contact_7[rn7] = \
                                contact_7.get(rn7,0)+1

    total = (len(salt_bridges)*3 +
             len(hbonds) + len(hydrophobic))
    top_8 = sorted(contact_8.items(),
                   key=lambda x:-x[1])[:20]
    top_7 = sorted(contact_7.items(),
                   key=lambda x:-x[1])[:20]

    iface_8 = sorted(contact_8.keys())
    iface_7 = sorted(contact_7.keys())

    print(f"  {struct_id} [{mode}]: "
          f"SB={len(salt_bridges)} "
          f"HB={len(hbonds)} "
          f"HY={len(hydrophobic)} total={total}")
    print(f"    NSP8 interface: {iface_8}")
    print(f"    NSP7 interface: {iface_7}")
    for a,b,d in salt_bridges:
        print(f"    SB: {a} — {b} ({d} Å)")

    return {
        "struct":        struct_id,
        "mode":          mode,
        "salt_bridges":  salt_bridges,
        "n_hbonds":      len(hbonds),
        "n_hydrophobic": len(hydrophobic),
        "total_score":   total,
        "iface_NSP8":    iface_8,
        "iface_NSP7":    iface_7,
        "top20_NSP8":    top_8,
        "top20_NSP7":    top_7,
    }


def main():
    print("\n" + "="*60)
    print("  Script 05_6: Interface Analysis — NSP7-NSP8")
    print("  DUAL BINDING MODE HYPOTHESIS")
    print("="*60 + "\n")

    results      = []
    crystal_res  = []
    af3_res      = []

    for sid,(pf,c8,c7,mode) in STRUCTURES.items():
        r = analyze_interface(sid,pf,c8,c7,mode)
        results.append(r)
        if mode == "crystal":
            crystal_res.append(r)
        else:
            af3_res.append(r)

    # Mode A consensus (crystal)
    from collections import Counter
    c8_A = Counter()
    c7_A = Counter()
    for r in crystal_res:
        c8_A.update(r["iface_NSP8"])
        c7_A.update(r["iface_NSP7"])
    hot8_A = sorted([p for p,n in c8_A.items()
                     if n >= 1],
                    key=lambda x:-c8_A[x])
    hot7_A = sorted([p for p,n in c7_A.items()
                     if n >= 1],
                    key=lambda x:-c7_A[x])

    # Mode B (AF3 only)
    af3_r   = af3_res[0]
    hot8_B  = [p for p,_ in af3_r["top20_NSP8"]]
    hot7_B  = [p for p,_ in af3_r["top20_NSP7"]]

    print(f"\n  ── Mode A: Crystal interface ──")
    print(f"  NSP8 hotspots ({len(hot8_A)}): {hot8_A}")
    print(f"  NSP7 hotspots ({len(hot7_A)}): {hot7_A}")

    print(f"\n  ── Mode B: AF3 interface ──")
    print(f"  NSP8 hotspots ({len(hot8_B)}): "
          f"{hot8_B[:15]}")
    print(f"  NSP7 hotspots ({len(hot7_B)}): "
          f"{hot7_B[:15]}")

    print(f"\n  ── Salt bridge summary ──")
    all_sbs = {}
    for r in results:
        for a,b,d in r["salt_bridges"]:
            key = f"{a}—{b}"
            all_sbs.setdefault(key,[]).append(
                (r["struct"],r["mode"],d))
    if all_sbs:
        for pair,occ in sorted(
                all_sbs.items(),
                key=lambda x:-len(x[1])):
            structs = [s for s,_,_ in occ]
            best_d  = min(d for _,_,d in occ)
            modes   = list(set(m for _,m,_ in occ))
            print(f"    {pair:<35} "
                  f"[{', '.join(structs)}] "
                  f"mode={modes} "
                  f"min={best_d:.2f}Å")
    else:
        print("    No salt bridges detected")

    class NpEncoder(json.JSONEncoder):
        def default(self, o):
            if isinstance(o,(np.floating,
                              np.float32,
                              np.float64)):
                return float(o)
            if isinstance(o, np.integer):
                return int(o)
            return super().default(o)

    out = RES_DIR / "interface_analysis_6.json"
    with open(out,"w") as f:
        json.dump({
            "complex":       "NSP7-NSP8",
            "dual_mode":     True,
            "mode_A_crystal":{
                "description":  "C-terminal NSP8 — "
                                "NSP7 loop",
                "pdbs":         ["7BV2","6NUR"],
                "hotspots_NSP8":hot8_A,
                "hotspots_NSP7":hot7_A,
            },
            "mode_B_af3":{
                "description":  "N-terminal NSP8 helix"
                                " — NSP7 body",
                "source":       "AF3 (iptm=0.87)",
                "hotspots_NSP8":hot8_B,
                "hotspots_NSP7":hot7_B,
            },
            "salt_bridges_all": all_sbs,
            "per_structure": {
                r["struct"]:{
                    "mode":          r["mode"],
                    "salt_bridges":  r["salt_bridges"],
                    "n_hbonds":      r["n_hbonds"],
                    "n_hydrophobic": r["n_hydrophobic"],
                    "total_score":   r["total_score"],
                    "iface_NSP8":    r["iface_NSP8"],
                    "iface_NSP7":    r["iface_NSP7"],
                } for r in results},
        }, f, indent=2, cls=NpEncoder)

    print(f"\n  Saved: 02-validation/NSP7-NSP8/"
          f"interface_analysis_6.json")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
