#!/usr/bin/env python3
"""
Script 07_pocket_NSP13-Helicase_7.py

fpocket pocket detection on NSP13 homodimer interface.
Runs on 7NIO (chains A+E), 6XEZ (chains E+F), and AF3 monomer.

Primary target: LYS414-ASP580/ASP583 dual salt bridge region (SARS-selective)
Conserved backbone anchors: GLY415, GLY478, THR552, ALA553

Docking box defined from hotspot Cα coordinates with 6.0 A padding.

Output: 02-validation/NSP13-Helicase/pocket_analysis_7.json
"""

import os, json, subprocess, shutil
import numpy as np
from Bio import PDB

# ── Paths ──────────────────────────────────────────────────────────────────
BASE    = os.path.expanduser("~/projects/rtc-pan-coronavirus")
PDB_DIR = f"{BASE}/00-reference/pdb_structures"
AF3_DIR = f"{BASE}/01-alphafold3/NSP13-Helicase"
VAL_DIR = f"{BASE}/02-validation/NSP13-Helicase"
TMP_DIR = f"{BASE}/tmp_fpocket_NSP13"
os.makedirs(VAL_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

# ── Structures to analyze ──────────────────────────────────────────────────
STRUCTURES = {
    "7NIO_dimer" : {
        "pdb"    : f"{PDB_DIR}/7NIO.pdb",
        "chains" : ["A", "E"],
        "label"  : "7NIO A+E (primary dimer, 2.80 A)"
    },
    "6XEZ_dimer" : {
        "pdb"    : f"{PDB_DIR}/6XEZ.pdb",
        "chains" : ["E", "F"],
        "label"  : "6XEZ E+F (RdRp complex, 3.50 A)"
    },
    "AF3_monomer" : {
        "pdb"    : f"{AF3_DIR}/NSP13_Helicase_best_model.pdb",
        "chains" : ["A"],
        "label"  : "AF3 monomer (ptm=0.910)"
    },
}

# ── Hotspot residues (7NIO chain A numbering) ──────────────────────────────
ALL_HOTSPOTS = [116,413,414,415,477,478,479,480,481,482,
                551,552,553,579,580,583,584]

# Primary pharmacophore residues (salt bridge anchors + top contacts)
PRIMARY_HOTSPOTS = {414, 477, 480, 482, 579, 580, 583, 584}

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

parser = PDB.PDBParser(QUIET=True)

# ── Helper: extract chains to temp PDB ────────────────────────────────────
def extract_chains_to_pdb(source_pdb, chains, out_path):
    struct = parser.get_structure("tmp", source_pdb)
    io     = PDB.PDBIO()

    class ChainSelector(PDB.Select):
        def accept_chain(self, chain):
            return chain.get_id() in chains
        def accept_residue(self, residue):
            return residue.get_id()[0] == " "  # ATOM only, no HETATM

    io.set_structure(struct)
    io.save(out_path, ChainSelector())

# ── Helper: run fpocket ────────────────────────────────────────────────────
def run_fpocket(pdb_path, label):
    out_dir = pdb_path.replace(".pdb", "_out")
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    cmd    = f"fpocket -f {pdb_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                            cwd=TMP_DIR)
    if result.returncode != 0:
        print(f"  fpocket error for {label}: {result.stderr[:300]}")
        return None
    if not os.path.exists(out_dir):
        # fpocket sometimes puts output in cwd
        basename = os.path.basename(pdb_path).replace(".pdb","_out")
        alt_dir  = os.path.join(TMP_DIR, basename)
        if os.path.exists(alt_dir):
            return alt_dir
        print(f"  WARNING: fpocket output dir not found for {label}")
        return None
    return out_dir

# ── Helper: parse fpocket pocket info file ────────────────────────────────
def parse_fpocket_info(info_path):
    pockets = []
    current = {}
    with open(info_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                pocket_num = int(line.split()[1].rstrip(":"))
                current = {"pocket_num": pocket_num}
            elif ":" in line:
                key, _, val = line.partition(":")
                key = key.strip().lower().replace(" ","_")
                try:
                    current[key] = float(val.strip().split()[0])
                except:
                    current[key] = val.strip()
    if current:
        pockets.append(current)
    return pockets

# ── Helper: find residues in pocket ───────────────────────────────────────
def get_pocket_residues(pocket_pdb_path):
    residues = set()
    try:
        with open(pocket_pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        res_num = int(line[22:26].strip())
                        residues.add(res_num)
                    except:
                        pass
    except:
        pass
    return residues

# ── Helper: get Cα coordinates for residues ───────────────────────────────
def get_ca_coords(struct, chain_ids, res_nums):
    coords = []
    for chain_id in chain_ids:
        try:
            chain = struct[0][chain_id]
        except:
            continue
        for res in chain:
            if not PDB.is_aa(res):
                continue
            if res.get_id()[1] in res_nums and "CA" in res:
                coords.append(res["CA"].get_vector().get_array())
    return coords

# ── Main analysis ──────────────────────────────────────────────────────────
print("=" * 60)
print("POCKET DETECTION — NSP13 HOMODIMER")
print("=" * 60)

results = {}

for struct_key, info in STRUCTURES.items():
    print(f"\n{'─'*60}")
    print(f"Structure: {info['label']}")
    print(f"{'─'*60}")

    # Extract relevant chains to temp PDB
    tmp_pdb = f"{TMP_DIR}/{struct_key}.pdb"
    extract_chains_to_pdb(info["pdb"], info["chains"], tmp_pdb)
    print(f"  Extracted chains {info['chains']} -> {tmp_pdb}")

    # Run fpocket
    out_dir = run_fpocket(tmp_pdb, struct_key)
    if out_dir is None:
        print(f"  SKIP: fpocket failed")
        continue

    # Find info file
    info_files = [f for f in os.listdir(out_dir) if f.endswith("_info.txt")]
    if not info_files:
        print(f"  WARNING: no info file in {out_dir}")
        continue
    info_file = os.path.join(out_dir, info_files[0])
    pockets   = parse_fpocket_info(info_file)
    print(f"  Total pockets found: {len(pockets)}")

    if not pockets:
        continue

    # Score each pocket by overlap with hotspot residues
    pockets_dir = os.path.join(out_dir, "pockets")
    scored = []

    for pocket in pockets:
        pnum     = pocket["pocket_num"]
        drug_key = [k for k in pocket if "druggab" in k]
        drug_score = float(pocket[drug_key[0]]) if drug_key else 0.0

        # Find pocket residues from pocket PDB file
        pocket_pdb = os.path.join(pockets_dir, f"pocket{pnum}_atm.pdb")
        pocket_res = get_pocket_residues(pocket_pdb) if os.path.exists(pocket_pdb) else set()

        # Count overlap with hotspots (weighted)
        overlap_all  = len(pocket_res & set(ALL_HOTSPOTS))
        overlap_prim = len(pocket_res & PRIMARY_HOTSPOTS)
        has_lys414   = 414 in pocket_res
        has_asp580   = 580 in pocket_res
        has_asp583   = 583 in pocket_res
        has_ile480   = 480 in pocket_res
        has_his482   = 482 in pocket_res

        # Composite overlap score (primary residues weighted x3)
        overlap_score = overlap_all + (overlap_prim * 2) + \
                        (3 if has_lys414 else 0) + \
                        (3 if (has_asp580 or has_asp583) else 0)

        # Volume key
        vol_key = [k for k in pocket if "volume" in k and "monte" not in k]
        volume  = float(pocket[vol_key[0]]) if vol_key else 0.0

        scored.append({
            "pocket_num"   : pnum,
            "drug_score"   : round(drug_score, 4),
            "volume"       : round(volume, 1),
            "overlap_all"  : overlap_all,
            "overlap_prim" : overlap_prim,
            "overlap_score": overlap_score,
            "has_LYS414"   : has_lys414,
            "has_ASP580"   : has_asp580,
            "has_ASP583"   : has_asp583,
            "has_ILE480"   : has_ile480,
            "has_HIS482"   : has_his482,
            "pocket_residues": sorted(pocket_res)
        })

    # Sort by overlap score first, then druggability
    scored.sort(key=lambda x: (-x["overlap_score"], -x["drug_score"]))

    best = scored[0] if scored else None

    print(f"\n  Top 5 pockets by interface overlap:")
    print(f"  {'#':>3}  {'Drug':>6}  {'Vol':>7}  {'Ovlp':>5}  "
          f"{'Prim':>5}  {'LYS414':>7}  {'ASP580':>7}  {'ASP583':>7}  {'ILE480':>7}  {'HIS482':>7}")
    for p in scored[:5]:
        print(f"  {p['pocket_num']:>3}  {p['drug_score']:>6.3f}  "
              f"{p['volume']:>7.1f}  {p['overlap_all']:>5}  "
              f"{p['overlap_prim']:>5}  "
              f"{'✅' if p['has_LYS414'] else '❌':>7}  "
              f"{'✅' if p['has_ASP580'] else '❌':>7}  "
              f"{'✅' if p['has_ASP583'] else '❌':>7}  "
              f"{'✅' if p['has_ILE480'] else '❌':>7}  "
              f"{'✅' if p['has_HIS482'] else '❌':>7}")

    if best:
        print(f"\n  Best pocket selected: Pocket {best['pocket_num']}")
        print(f"    Druggability  : {best['drug_score']}")
        print(f"    Volume        : {best['volume']} A3")
        print(f"    Hotspot overlap: {best['overlap_all']}/{len(ALL_HOTSPOTS)} total, "
              f"{best['overlap_prim']}/{len(PRIMARY_HOTSPOTS)} primary")
        print(f"    Contains LYS414: {'YES' if best['has_LYS414'] else 'NO'}")
        print(f"    Contains ASP580: {'YES' if best['has_ASP580'] else 'NO'}")
        print(f"    Contains ASP583: {'YES' if best['has_ASP583'] else 'NO'}")
        print(f"    Contains ILE480: {'YES' if best['has_ILE480'] else 'NO'}")
        print(f"    Contains HIS482: {'YES' if best['has_HIS482'] else 'NO'}")

    results[struct_key] = {
        "label"        : info["label"],
        "n_pockets"    : len(pockets),
        "best_pocket"  : best,
        "top5_pockets" : scored[:5]
    }

# ── Docking box from 7NIO hotspot Cα coords ───────────────────────────────
print(f"\n{'='*60}")
print("DOCKING BOX — from 7NIO primary hotspot Cα coordinates")
print(f"{'='*60}")

struct_7NIO = parser.get_structure("7NIO", f"{PDB_DIR}/7NIO.pdb")

# Use all 17 hotspot residues for box, anchored on primary
ca_all     = get_ca_coords(struct_7NIO, ["A"], ALL_HOTSPOTS)
ca_primary = get_ca_coords(struct_7NIO, ["A"], list(PRIMARY_HOTSPOTS))

if ca_all:
    arr  = np.array(ca_all)
    cmin = arr.min(axis=0)
    cmax = arr.max(axis=0)
    PAD  = 6.0
    box_min    = cmin - PAD
    box_max    = cmax + PAD
    box_center = (box_min + box_max) / 2
    box_size   = box_max - box_min
    box_volume = float(np.prod(box_size))

    print(f"\n  Hotspot Cα atoms used : {len(ca_all)} (all hotspots)")
    print(f"  Primary Cα atoms      : {len(ca_primary)}")
    print(f"  Padding               : {PAD} A")
    print(f"\n  Docking box:")
    print(f"    Center : ({box_center[0]:.3f}, {box_center[1]:.3f}, {box_center[2]:.3f})")
    print(f"    Size   : {box_size[0]:.3f} x {box_size[1]:.3f} x {box_size[2]:.3f} A")
    print(f"    Volume : {box_volume:,.0f} A3")

    # Verify key residues in box
    print(f"\n  Key residue Cα positions:")
    for chain in struct_7NIO[0]:
        if chain.get_id() != "A":
            continue
        for res in chain:
            if PDB.is_aa(res) and res.get_id()[1] in PRIMARY_HOTSPOTS:
                if "CA" in res:
                    ca  = res["CA"].get_vector().get_array()
                    ins = all(box_min[i] <= ca[i] <= box_max[i] for i in range(3))
                    print(f"    {res.get_resname()}{res.get_id()[1]:4d}  "
                          f"({ca[0]:.1f}, {ca[1]:.1f}, {ca[2]:.1f})  "
                          f"{'in box ✅' if ins else 'OUT OF BOX ❌'}")
else:
    print("  ERROR: no Cα coordinates found for hotspot residues")
    box_center = [0,0,0]
    box_size   = [20,20,20]
    box_volume = 0.0

# ── Save ───────────────────────────────────────────────────────────────────
docking_box = {
    "center_x" : round(float(box_center[0]), 3),
    "center_y" : round(float(box_center[1]), 3),
    "center_z" : round(float(box_center[2]), 3),
    "size_x"   : round(float(box_size[0]), 3),
    "size_y"   : round(float(box_size[1]), 3),
    "size_z"   : round(float(box_size[2]), 3),
    "volume_A3": round(box_volume, 0),
    "padding_A": PAD,
    "anchor_residues": sorted(PRIMARY_HOTSPOTS),
    "note": "LYS414-ASP580/ASP583 dual salt bridge region — SARS-CoV-1/2 selective"
}

output = {
    "complex"     : "NSP13-Helicase",
    "script"      : "07_pocket_NSP13-Helicase_7.py",
    "structures"  : results,
    "docking_box" : docking_box,
    "scientific_note": (
        "NSP13 homodimer interface is SARS-CoV-1/2 selective. "
        "ASP580 loses charge in MERS (D->A) and HCoV (D->T). "
        "Low fpocket druggability expected for PPI interface — "
        "docking box defined from hotspot coordinates regardless."
    )
}

out_path = f"{VAL_DIR}/pocket_analysis_7.json"
with open(out_path, "w") as f:
    json.dump(output, f, indent=2)
print(f"\n  Saved: {out_path}")

# Cleanup temp dir
shutil.rmtree(TMP_DIR, ignore_errors=True)
print("  Temp files cleaned up")
