#!/usr/bin/env python3
"""
Script 06_conservation_NSP12-NSP13_8.py

Conservation analysis of NSP12-NSP13 interface hotspots
across 5 coronaviruses.

NSP12 coordinates (polyprotein):
  SARS-CoV-2  P0DTD1  4393-5324  (932 aa)
  SARS-CoV-1  P0C6X7  4370-5301  (932 aa)
  MERS-CoV    K9N7C7  4378-5310  (933 aa)
  HCoV-229E   P0C6X1  4069-4995  (927 aa)
  HCoV-NL63   P0C6X5  4044-4970  (927 aa)

NSP13 coordinates (polyprotein):
  SARS-CoV-2  P0DTD1  5325-5925  (601 aa)
  SARS-CoV-1  P0C6X7  5302-5902  (601 aa)
  MERS-CoV    K9N7C7  5311-5908  (598 aa)
  HCoV-229E   P0C6X1  4996-5592  (597 aa)
  HCoV-NL63   P0C6X5  4971-5567  (597 aa)

Hotspot residues (6XEZ PDB numbering):
  NSP12: LEU900, ASP901, MET902, TYR903, SER904  (C-terminal tail)
  NSP13: PHE90, GLY91, LEU92, TYR93, LYS94, ASN95, THR96 (N-terminal)

Output: 02-validation/NSP12-NSP13/
  conservation_NSP12.csv
  conservation_NSP13.csv
  conservation_summary_8.json
"""

import os, json, math, subprocess, requests, time
import pandas as pd
from Bio import PDB, Align

# ── Paths ──────────────────────────────────────────────────────────────────
BASE     = os.path.expanduser("~/projects/rtc-pan-coronavirus")
SEQ_DIR  = f"{BASE}/00-reference/sequences/conservation"
PDB_DIR  = f"{BASE}/00-reference/pdb_structures"
VAL_DIR  = f"{BASE}/02-validation/NSP12-NSP13"
os.makedirs(SEQ_DIR, exist_ok=True)
os.makedirs(VAL_DIR, exist_ok=True)

# ── Coordinates ────────────────────────────────────────────────────────────
NSP12_COORDS = {
    "SARS-CoV-2" : ("P0DTD1", 4393, 5324),
    "SARS-CoV-1" : ("P0C6X7", 4370, 5301),
    "MERS-CoV"   : ("K9N7C7", 4378, 5310),
    "HCoV-229E"  : ("P0C6X1", 4069, 4995),
    "HCoV-NL63"  : ("P0C6X5", 4044, 4970),
}
NSP13_COORDS = {
    "SARS-CoV-2" : ("P0DTD1", 5325, 5925),
    "SARS-CoV-1" : ("P0C6X7", 5302, 5902),
    "MERS-CoV"   : ("K9N7C7", 5311, 5908),
    "HCoV-229E"  : ("P0C6X1", 4996, 5592),
    "HCoV-NL63"  : ("P0C6X5", 4971, 5567),
}

# ── Hotspot residues (6XEZ PDB numbering) ─────────────────────────────────
NSP12_HOTSPOTS = [
    (900, "LEU"), (901, "ASP"), (902, "MET"), (903, "TYR"), (904, "SER"),
]
NSP13_HOTSPOTS = [
    (90, "PHE"), (91, "GLY"), (92, "LEU"), (93, "TYR"),
    (94, "LYS"), (95, "ASN"), (96, "THR"),
]

# Primary pharmacophore residues
NSP12_PRIMARY = {901, 902}   # ASP901 (SB), MET902 (dominant hydrophobic)
NSP13_PRIMARY = {94}         # LYS94  (SB partner)

AA_MAP_3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}
AA_MAP_1TO3 = {v:k for k,v in AA_MAP_3TO1.items()}

# ── Helper: download and slice ─────────────────────────────────────────────
def fetch_nsp(acc, start, end, org, nsp_name):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            break
        except Exception as e:
            print(f"    Attempt {attempt+1} failed: {e}")
            time.sleep(5)
    else:
        return None
    lines = r.text.strip().split("\n")
    seq   = "".join(l for l in lines if not l.startswith(">"))
    sliced = seq[start-1:end]
    expected = end - start + 1
    print(f"    {org}: polyprotein={len(seq)} aa | "
          f"{nsp_name} [{start}-{end}] = {len(sliced)} aa (expected {expected})")
    if len(sliced) < expected * 0.9:
        print(f"    WARNING: sequence too short — check coordinates")
    time.sleep(1)
    return sliced

# ── Helper: Shannon conservation ──────────────────────────────────────────
def shannon_conservation(col_residues):
    residues = [aa for aa in col_residues if aa != "-"]
    if len(residues) < 2:
        return 1.0
    counts = {}
    for aa in residues:
        counts[aa] = counts.get(aa, 0) + 1
    total   = len(residues)
    entropy = -sum((c/total)*math.log2(c/total) for c in counts.values())
    max_ent = math.log2(min(total, 20))
    return round(1.0 - entropy/max_ent, 4) if max_ent > 0 else 1.0

# ── Helper: MUSCLE align + map positions ──────────────────────────────────
def run_conservation(nsp_name, coords_dict, hotspots, primary_set,
                     pdb_chain, pdb_file, sb_label=""):

    print(f"\n{'='*60}")
    print(f"{nsp_name} CONSERVATION")
    print(f"{'='*60}")

    # Step 1: Download sequences
    print("\nDownloading sequences ...")
    fastas = {}
    for org, (acc, start, end) in coords_dict.items():
        seq = fetch_nsp(acc, start, end, org, nsp_name)
        if seq:
            fastas[org] = seq

    if len(fastas) < 3:
        print("ERROR: not enough sequences")
        return []

    # Step 2: Write combined FASTA
    combined = f"{SEQ_DIR}/{nsp_name.replace('-','_')}_NSP12_NSP13_combined.fasta"
    with open(combined, "w") as f:
        for org, seq in fastas.items():
            acc = coords_dict[org][0]
            f.write(f">{org}|{acc}\n{seq}\n")

    # Step 3: MUSCLE alignment
    aligned_path = f"{SEQ_DIR}/{nsp_name.replace('-','_')}_NSP12_NSP13_aligned.fasta"
    cmd = f"muscle -align {combined} -output {aligned_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MUSCLE error: {result.stderr[:300]}")
        return []
    print(f"Alignment complete")

    # Step 4: Parse alignment
    aligned = {}
    current = None
    with open(aligned_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split("|")[0]
                aligned[current] = ""
            elif current:
                aligned[current] += line

    orgs    = list(aligned.keys())
    sars2   = aligned.get("SARS-CoV-2", "")
    if not sars2:
        print("ERROR: SARS-CoV-2 not in alignment")
        return []

    # Step 5: Map PDB local positions to alignment columns
    # Need to align PDB chain sequence against SARS-CoV-2 NSP sequence
    parser    = PDB.PDBParser(QUIET=True)
    struct    = parser.get_structure("pdb", pdb_file)
    pdb_reslist = [(r.get_id()[1], r.get_resname())
                   for r in struct[0][pdb_chain]
                   if PDB.is_aa(r) and r.get_resname() in AA_MAP_3TO1]
    pdb_seq   = "".join(AA_MAP_3TO1[r[1]] for r in pdb_reslist)
    sars2_seq = fastas.get("SARS-CoV-2", "")

    aligner = Align.PairwiseAligner()
    aligner.mode              = "global"
    aligner.match_score       =  2
    aligner.mismatch_score    = -1
    aligner.open_gap_score    = -5
    aligner.extend_gap_score  = -0.5

    aln = aligner.align(sars2_seq, pdb_seq)[0]
    sars2_to_pdb = {}
    for (rs, re), (qs, qe) in zip(aln.aligned[0], aln.aligned[1]):
        for rp, qp in zip(range(rs, re), range(qs, qe)):
            sars2_to_pdb[rp] = qp
    pdb_to_sars2 = {v: k for k, v in sars2_to_pdb.items()}

    # PDB residue number -> seq index in PDB seq
    pdb_num_to_seqidx = {pdb_num: i
                         for i, (pdb_num, _) in enumerate(pdb_reslist)}

    # SARS-CoV-2 local seq index -> alignment column
    sars2_pos_to_col = {}
    local_pos = 0
    for col, aa in enumerate(sars2):
        if aa != "-":
            sars2_pos_to_col[local_pos] = col
            local_pos += 1

    print(f"PDB sequence: {len(pdb_seq)} aa | "
          f"SARS-CoV-2 sequence: {len(sars2_seq)} aa | "
          f"Alignment columns: {len(sars2)}")

    # Step 6: Conservation at hotspot positions
    print(f"\n  {'PDB':>5}  {'Res':>4}  {'Score':>6}  {'Prim':>5}  "
          + "  ".join(f"{o[:8]:>8}" for o in orgs))
    print("  " + "-"*(30 + 10*len(orgs)))

    conservation_data = []
    for pdb_rn, resname_3 in hotspots:
        if pdb_rn not in pdb_num_to_seqidx:
            print(f"  {pdb_rn:>5}  {resname_3:>4}  NOT IN PDB")
            continue

        pdb_seqidx = pdb_num_to_seqidx[pdb_rn]
        if pdb_seqidx not in pdb_to_sars2:
            print(f"  {pdb_rn:>5}  {resname_3:>4}  NOT ALIGNED")
            continue

        sars2_seqidx = pdb_to_sars2[pdb_seqidx]
        if sars2_seqidx not in sars2_pos_to_col:
            print(f"  {pdb_rn:>5}  {resname_3:>4}  NO COLUMN")
            continue

        col      = sars2_pos_to_col[sars2_seqidx]
        col_res  = [aligned[org][col] for org in orgs]
        score    = shannon_conservation(col_res)
        is_prim  = "★" if pdb_rn in primary_set else ""

        print(f"  {pdb_rn:>5}  {resname_3:>4}  {score:>6.3f}  {is_prim:>5}  "
              + "  ".join(f"{aa:>8}" for aa in col_res))

        conservation_data.append({
            "pdb_res_num"  : pdb_rn,
            "res_name_3"   : resname_3,
            "res_name_1"   : AA_MAP_3TO1.get(resname_3, "?"),
            "conservation" : score,
            "is_primary"   : pdb_rn in primary_set,
            **{f"aa_{org.replace('-','_').replace('/','_')}": col_res[i]
               for i, org in enumerate(orgs)}
        })

    conserved = [d for d in conservation_data if d["conservation"] >= 0.8]
    variable  = [d for d in conservation_data if d["conservation"] <  0.8]

    print(f"\n  Conserved (>=0.8): {len(conserved)}/{len(conservation_data)}")
    for d in sorted(conserved, key=lambda x: -x["conservation"]):
        prim = " ★" if d["is_primary"] else ""
        print(f"    {d['res_name_3']}{d['pdb_res_num']:4d}  "
              f"score={d['conservation']:.3f}{prim}")

    if variable:
        print(f"  Variable (<0.8):")
        for d in sorted(variable, key=lambda x: x["conservation"]):
            prim = " ★" if d["is_primary"] else ""
            print(f"    {d['res_name_3']}{d['pdb_res_num']:4d}  "
                  f"score={d['conservation']:.3f}{prim}")

    print(f"\n  PRIMARY pharmacophore conservation:")
    for d in conservation_data:
        if d["is_primary"]:
            print(f"    {d['res_name_3']}{d['pdb_res_num']:4d}  "
                  f"score={d['conservation']:.3f}  {sb_label}")

    return conservation_data

# ── Run conservation for both NSPs ─────────────────────────────────────────
pdb_6XEZ = f"{PDB_DIR}/6XEZ.pdb"

nsp12_cons = run_conservation(
    "NSP12", NSP12_COORDS, NSP12_HOTSPOTS, NSP12_PRIMARY,
    pdb_chain="A", pdb_file=pdb_6XEZ,
    sb_label="[SB: ASP901-LYS94]"
)

nsp13_cons = run_conservation(
    "NSP13", NSP13_COORDS, NSP13_HOTSPOTS, NSP13_PRIMARY,
    pdb_chain="E", pdb_file=pdb_6XEZ,
    sb_label="[SB: LYS94-ASP901]"
)

# ── Save CSVs ──────────────────────────────────────────────────────────────
if nsp12_cons:
    df12 = pd.DataFrame(nsp12_cons)
    p12  = f"{VAL_DIR}/conservation_NSP12.csv"
    df12.to_csv(p12, index=False)
    print(f"\nSaved: {p12}")

if nsp13_cons:
    df13 = pd.DataFrame(nsp13_cons)
    p13  = f"{VAL_DIR}/conservation_NSP13.csv"
    df13.to_csv(p13, index=False)
    print(f"Saved: {p13}")

# ── Summary ────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("CONSERVATION SUMMARY — NSP12-NSP13")
print(f"{'='*60}")

for nsp_name, cons_data in [("NSP12", nsp12_cons), ("NSP13", nsp13_cons)]:
    if not cons_data:
        continue
    conserved = [d for d in cons_data if d["conservation"] >= 0.8]
    print(f"\n  {nsp_name}: {len(conserved)}/{len(cons_data)} conserved >=0.8")
    for d in cons_data:
        prim = " ★ PRIMARY" if d["is_primary"] else ""
        print(f"    {d['res_name_3']}{d['pdb_res_num']:4d}  "
              f"cons={d['conservation']:.3f}{prim}")

# ── Save JSON ──────────────────────────────────────────────────────────────
summary = {
    "complex"              : "NSP12-NSP13",
    "script"               : "06_conservation_NSP12-NSP13_8.py",
    "coronaviruses"        : list(NSP12_COORDS.keys()),
    "nsp12_conservation"   : nsp12_cons,
    "nsp13_conservation"   : nsp13_cons,
    "salt_bridge_anchors"  : {
        "ASP901_NSP12": next((d["conservation"] for d in nsp12_cons
                              if d["pdb_res_num"]==901), None),
        "LYS94_NSP13" : next((d["conservation"] for d in nsp13_cons
                              if d["pdb_res_num"]==94),  None),
        "MET902_NSP12": next((d["conservation"] for d in nsp12_cons
                              if d["pdb_res_num"]==902), None),
    }
}

json_path = f"{VAL_DIR}/conservation_summary_8.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)
print(f"\nSaved: {json_path}")
