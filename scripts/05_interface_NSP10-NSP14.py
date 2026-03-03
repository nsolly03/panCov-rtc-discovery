"""
Script 05: Interface Analysis — NSP10-NSP14
===========================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Analyzes the NSP10-NSP14 interface across ALL available structures:
  - 7DIY  : SARS-CoV-2 crystal structure (2.69 A) — primary
  - 5C8T  : SARS-CoV-1 crystal structure (3.20 A) — cross-validation
  - AF3   : AlphaFold3 predicted model            — working model

  Hotspots found in ALL THREE structures are the most reliable
  drug targets and will anchor the docking box in Script 07.

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/05_interface_NSP10-NSP14.py
"""

import json
from pathlib import Path
from Bio import PDB
from Bio.PDB import NeighborSearch

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT / "01-alphafold3"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
RES_DIR.mkdir(parents=True, exist_ok=True)

# ── All structures to analyze ──────────────────────────────
STRUCTURES = [
    {
        "label":    "PDB_7DIY",
        "file":     PDB_DIR / "7DIY.pdb",
        "chain_a":  "A",
        "chain_b":  "B",
        "label_a":  "NSP10",
        "label_b":  "NSP14",
        "source":   "SARS-CoV-2 crystal 2.69A",
    },
    {
        "label":    "PDB_5C8T",
        "file":     PDB_DIR / "5C8T.pdb",
        "chain_a":  "A",
        "chain_b":  "B",
        "label_a":  "NSP10",
        "label_b":  "NSP14",
        "source":   "SARS-CoV-1 crystal 3.20A",
    },
    {
        "label":    "AF3_model",
        "file":     AF3_DIR / "NSP10-NSP14" / "NSP10_NSP14_best_model.pdb",
        "chain_a":  "A",
        "chain_b":  "B",
        "label_a":  "NSP10",
        "label_b":  "NSP14",
        "source":   "AlphaFold3 iptm=0.89",
    },
]

# ── Residue classifications ────────────────────────────────
HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO","TYR"}
POSITIVE    = {"ARG","LYS","HIS"}
NEGATIVE    = {"ASP","GLU"}
HBOND_DONOR = {"ARG","LYS","SER","THR","TYR","TRP","ASN","GLN","HIS"}
HBOND_ACC   = {"ASP","GLU","SER","THR","TYR","ASN","GLN","HIS"}


def get_atoms(structure, chain_id):
    atoms = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                atoms.extend(list(chain.get_atoms()))
    return atoms


def get_residues(structure, chain_id):
    residues = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for res in chain:
                    if res.id[0] == " ":
                        residues.append(res)
    return residues


def analyze_interface(structure, chain_a, chain_b):
    """Analyze contacts — one entry per residue pair"""
    res_a  = get_residues(structure, chain_a)
    atoms_a = get_atoms(structure, chain_a)
    ns_b   = NeighborSearch(get_atoms(structure, chain_b))
    ns_a   = NeighborSearch(atoms_a)

    contacts   = []
    seen_pairs = set()
    hotspot_a  = {}
    hotspot_b  = {}

    for res in res_a:
        neighbors = ns_b.search(
            res.center_of_mass(), 8.0, "R")

        for nb_res in neighbors:
            if nb_res.id[0] != " ":
                continue
            rnum_a = res.id[1]
            rnum_b = nb_res.id[1]
            pair   = (rnum_a, rnum_b)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)

            rname_a = res.resname
            rname_b = nb_res.resname

            # Minimum atom-atom distance
            min_dist = 999.0
            for at_a in res.get_atoms():
                for at_b in nb_res.get_atoms():
                    d = float(at_a - at_b)
                    if d < min_dist:
                        min_dist = d

            if min_dist > 5.0:
                continue

            # Classify
            contact_type = "VDW"
            if (rname_a in POSITIVE and rname_b in NEGATIVE or
                rname_a in NEGATIVE and rname_b in POSITIVE):
                if min_dist < 4.0:
                    contact_type = "SALT_BRIDGE"
            elif (rname_a in HBOND_DONOR and rname_b in HBOND_ACC or
                  rname_a in HBOND_ACC   and rname_b in HBOND_DONOR):
                if min_dist < 3.5:
                    contact_type = "HBOND"
            elif rname_a in HYDROPHOBIC and rname_b in HYDROPHOBIC:
                contact_type = "HYDROPHOBIC"

            score = {"SALT_BRIDGE":3,
                     "HBOND":2,
                     "HYDROPHOBIC":1,
                     "VDW":0.5}[contact_type]

            key_a = f"{rname_a}{rnum_a}"
            key_b = f"{rname_b}{rnum_b}"
            hotspot_a[key_a] = hotspot_a.get(key_a, 0) + score
            hotspot_b[key_b] = hotspot_b.get(key_b, 0) + score

            contacts.append({
                "res_a":    key_a,
                "res_b":    key_b,
                "type":     contact_type,
                "distance": round(min_dist, 2),
            })

    top_a = sorted(hotspot_a.items(),
                   key=lambda x: x[1], reverse=True)[:20]
    top_b = sorted(hotspot_b.items(),
                   key=lambda x: x[1], reverse=True)[:20]

    return contacts, top_a, top_b


def print_summary(contacts, top_a, top_b,
                  label_a, label_b, source):
    type_counts = {}
    for c in contacts:
        type_counts[c["type"]] = type_counts.get(
            c["type"], 0) + 1

    print(f"\n  Source : {source}")
    print(f"  Total contacts : {len(contacts)}")
    for t, n in sorted(type_counts.items()):
        print(f"    {t:<15} : {n}")

    print(f"\n  Top 10 {label_a} hotspots:")
    print(f"  {'Residue':<12} {'Score':>8}")
    print(f"  {'-'*12} {'-'*8}")
    for res, score in top_a[:10]:
        print(f"  {res:<12} {score:>8.1f}")

    print(f"\n  Top 10 {label_b} hotspots:")
    print(f"  {'Residue':<12} {'Score':>8}")
    print(f"  {'-'*12} {'-'*8}")
    for res, score in top_b[:10]:
        print(f"  {res:<12} {score:>8.1f}")

    salt = [c for c in contacts
            if c["type"] == "SALT_BRIDGE"]
    if salt:
        print(f"\n  Salt bridges:")
        for sb in salt:
            print(f"  {label_a}:{sb['res_a']} -- "
                  f"{label_b}:{sb['res_b']} "
                  f"({sb['distance']} A)")


def main():
    print("\n" + "="*55)
    print("  Script 05: Interface Analysis — NSP10-NSP14")
    print("  Structures: 7DIY (CoV-2) + 5C8T (CoV-1) + AF3")
    print("="*55)

    parser  = PDB.PDBParser(QUIET=True)
    results = {}
    all_top_a = {}
    all_top_b = {}

    for s in STRUCTURES:
        print(f"\n{'='*55}")
        print(f"  Analyzing: {s['label']}  ({s['source']})")
        print(f"{'='*55}")

        if not s["file"].exists():
            print(f"  SKIP — file not found: {s['file'].name}")
            continue

        structure = parser.get_structure(
            s["label"], s["file"])
        contacts, top_a, top_b = analyze_interface(
            structure, s["chain_a"], s["chain_b"])

        print_summary(contacts, top_a, top_b,
                      s["label_a"], s["label_b"],
                      s["source"])

        results[s["label"]] = {
            "source":        s["source"],
            "total_contacts": len(contacts),
            "top20_NSP10":   top_a,
            "top20_NSP14":   top_b,
            "salt_bridges":  [c for c in contacts
                              if c["type"] == "SALT_BRIDGE"],
        }

        # Collect top residues per structure
        all_top_a[s["label"]] = {r[0] for r in top_a}
        all_top_b[s["label"]] = {r[0] for r in top_b}

    # ── Consensus across all structures ───────────────────
    print(f"\n{'='*55}")
    print("  CONSENSUS HOTSPOTS")
    print("  Residues in top 20 of ALL structures")
    print(f"{'='*55}")

    labels = list(all_top_a.keys())

    if len(labels) >= 2:
        consensus_a = set.intersection(*all_top_a.values())
        consensus_b = set.intersection(*all_top_b.values())
    else:
        consensus_a = set()
        consensus_b = set()

    print(f"\n  NSP10 consensus hotspots : {sorted(consensus_a)}")
    print(f"  NSP14 consensus hotspots : {sorted(consensus_b)}")
    print(f"\n  These residues are your highest confidence")
    print(f"  drug targets for the docking box.")

    results["consensus"] = {
        "NSP10": sorted(consensus_a),
        "NSP14": sorted(consensus_b),
    }

    # ── Save ──────────────────────────────────────────────
    out = RES_DIR / "interface_analysis.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n  Saved: 02-validation/NSP10-NSP14/"
          f"interface_analysis.json")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
