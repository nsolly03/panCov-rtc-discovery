#!/usr/bin/env python3
"""
Script 06_conservation_NSP13-Helicase_7.py

Conservation analysis of NSP13 homodimer interface hotspots
across 5 coronaviruses using MUSCLE alignment + Shannon entropy.

UniProt coordinates (polyprotein positions):
  SARS-CoV-2  P0DTD1  5325-5925
  SARS-CoV-1  P0C6X7  5302-5902
  MERS-CoV    K9N7C7  5311-5908
  HCoV-229E   P0C6X1  4996-5592
  HCoV-NL63   P0C6X5  4971-5567

Hotspot residues (7NIO chain A local numbering — primary interface):
  ASN116, THR413, LYS414, GLY415, LYS477, GLY478, VAL479,
  ILE480, THR481, HIS482, GLU551, THR552, ALA553,
  ARG579, ASP580, ASP583, LYS584

Output: 02-validation/NSP13-Helicase/conservation_NSP13.csv
        02-validation/NSP13-Helicase/conservation_summary_7.json
"""

import os, json, math, subprocess, requests, time
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────
BASE     = os.path.expanduser("~/projects/rtc-pan-coronavirus")
SEQ_DIR  = f"{BASE}/00-reference/sequences/conservation"
VAL_DIR  = f"{BASE}/02-validation/NSP13-Helicase"
os.makedirs(SEQ_DIR, exist_ok=True)
os.makedirs(VAL_DIR, exist_ok=True)

# ── NSP13 UniProt coordinates ──────────────────────────────────────────────
CORONAVIRUSES = {
    "SARS-CoV-2" : ("P0DTD1", 5325, 5925),
    "SARS-CoV-1" : ("P0C6X7", 5302, 5902),
    "MERS-CoV"   : ("K9N7C7", 5311, 5908),
    "HCoV-229E"  : ("P0C6X1", 4996, 5592),
    "HCoV-NL63"  : ("P0C6X5", 4971, 5567),
}

# ── Hotspot residues (7NIO chain A local numbering) ───────────────────────
HOTSPOT_RESIDUES = [
    (116, "ASN"), (413, "THR"), (414, "LYS"), (415, "GLY"),
    (477, "LYS"), (478, "GLY"), (479, "VAL"), (480, "ILE"),
    (481, "THR"), (482, "HIS"), (551, "GLU"), (552, "THR"),
    (553, "ALA"), (579, "ARG"), (580, "ASP"), (583, "ASP"),
    (584, "LYS"),
]

# Salt bridge anchors — extra attention
PRIMARY_HOTSPOTS = {414, 477, 480, 482, 579, 580, 583, 584}

AA_MAP_1TO3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS",
    "Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE",
    "L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
    "S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"
}
AA_MAP_3TO1 = {v:k for k,v in AA_MAP_1TO3.items()}

# ── Step 1: Download and slice NSP13 sequences ─────────────────────────────
print("=" * 60)
print("STEP 1: Download NSP13 sequences from UniProt")
print("=" * 60)

nsp13_fastas = {}

for org, (acc, start, end) in CORONAVIRUSES.items():
    print(f"\n  Fetching {org} ({acc}) positions {start}-{end} ...")
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
        print(f"    ERROR: could not fetch {acc}")
        continue

    lines = r.text.strip().split("\n")
    seq   = "".join(l for l in lines if not l.startswith(">"))
    nsp13 = seq[start-1:end]  # 1-based slice

    expected_len = end - start + 1
    print(f"    Full polyprotein: {len(seq)} aa")
    print(f"    NSP13 slice [{start}-{end}]: {len(nsp13)} aa  (expected {expected_len})")

    if len(nsp13) < 400:
        print(f"    WARNING: NSP13 too short — check coordinates")

    nsp13_fastas[org] = nsp13

    fasta_path = f"{SEQ_DIR}/NSP13_{org.replace('-','_').replace('/','_')}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{org}|{acc}|NSP13|{start}-{end}\n{nsp13}\n")
    print(f"    Saved: {fasta_path}")
    time.sleep(1)

# ── Step 2: Write combined FASTA for MUSCLE ────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Write combined FASTA")
print("=" * 60)

combined_fasta = f"{SEQ_DIR}/NSP13_all_coronaviruses.fasta"
with open(combined_fasta, "w") as f:
    for org, seq in nsp13_fastas.items():
        acc = CORONAVIRUSES[org][0]
        f.write(f">{org}|{acc}|NSP13\n{seq}\n")
print(f"  Written: {combined_fasta}")
print(f"  Sequences: {len(nsp13_fastas)}")

# ── Step 3: MUSCLE alignment ───────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: MUSCLE multiple sequence alignment")
print("=" * 60)

aligned_fasta = f"{SEQ_DIR}/NSP13_aligned.fasta"
cmd = f"muscle -align {combined_fasta} -output {aligned_fasta}"
print(f"  Running: {cmd}")
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
if result.returncode != 0:
    print(f"  MUSCLE stderr: {result.stderr[:500]}")
else:
    print(f"  Alignment complete")

# ── Step 4: Parse alignment ────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Parse alignment")
print("=" * 60)

def parse_fasta(path):
    seqs = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split("|")[0]
                seqs[current] = ""
            elif current:
                seqs[current] += line
    return seqs

aligned = parse_fasta(aligned_fasta)
print(f"  Organisms aligned: {list(aligned.keys())}")
aln_lengths = set(len(s) for s in aligned.values())
print(f"  Alignment length: {aln_lengths}")

# ── Step 5: Map SARS-CoV-2 local positions to alignment columns ────────────
print("\n" + "=" * 60)
print("STEP 5: Map SARS-CoV-2 positions to alignment columns")
print("=" * 60)

sars2_aln = aligned.get("SARS-CoV-2", "")
if not sars2_aln:
    print("  ERROR: SARS-CoV-2 not found in alignment")
    exit(1)

# Build map: local_pos (1-based, gaps excluded) -> alignment_col (0-based)
sars2_pos_to_col = {}
local_pos = 0
for col, aa in enumerate(sars2_aln):
    if aa != "-":
        local_pos += 1
        sars2_pos_to_col[local_pos] = col

print(f"  SARS-CoV-2 aligned length (no gaps): {local_pos} aa")
print(f"  Total alignment columns: {len(sars2_aln)}")

# Verify hotspot positions are reachable
missing_hotspots = [pos for pos, _ in HOTSPOT_RESIDUES if pos not in sars2_pos_to_col]
if missing_hotspots:
    print(f"  WARNING: hotspot positions not found in alignment: {missing_hotspots}")
else:
    print(f"  All {len(HOTSPOT_RESIDUES)} hotspot positions mapped successfully")

# ── Step 6: Shannon entropy conservation ──────────────────────────────────
print("\n" + "=" * 60)
print("STEP 6: Shannon entropy per alignment column")
print("=" * 60)

orgs = list(aligned.keys())
n    = len(orgs)
aln_len = len(list(aligned.values())[0])

def shannon_conservation(col_residues):
    """1 - normalized Shannon entropy. 1.0 = identical, 0.0 = max diverse."""
    residues = [aa for aa in col_residues if aa != "-"]
    if len(residues) < 2:
        return 1.0
    counts = {}
    for aa in residues:
        counts[aa] = counts.get(aa, 0) + 1
    total = len(residues)
    entropy = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            entropy -= p * math.log2(p)
    max_entropy = math.log2(min(len(residues), 20))
    if max_entropy == 0:
        return 1.0
    return round(1.0 - entropy / max_entropy, 4)

def get_col_residues(col):
    return [aligned[org][col] for org in orgs]

# ── Step 7: Conservation at hotspot positions ──────────────────────────────
print("\n" + "=" * 60)
print("STEP 7: Conservation scores at hotspot residues")
print("=" * 60)

print(f"\n  {'Pos':>4}  {'Res':>3}  {'Score':>6}  {'Primary':>8}  "
      + "  ".join(f"{o[:8]:>8}" for o in orgs))
print("  " + "-" * (40 + 10 * len(orgs)))

conservation_data = []

for local_pos, resname_3 in HOTSPOT_RESIDUES:
    if local_pos not in sars2_pos_to_col:
        print(f"  {local_pos:>4}  {resname_3:>3}  NOT FOUND IN ALIGNMENT")
        continue

    col     = sars2_pos_to_col[local_pos]
    col_res = get_col_residues(col)
    score   = shannon_conservation(col_res)
    is_prim = "★" if local_pos in PRIMARY_HOTSPOTS else ""

    # AA at each organism
    aa_per_org = []
    for org in orgs:
        aa = aligned[org][col] if col < len(aligned[org]) else "-"
        aa_per_org.append(aa)

    print(f"  {local_pos:>4}  {resname_3:>3}  {score:>6.3f}  {is_prim:>8}  "
          + "  ".join(f"{aa:>8}" for aa in aa_per_org))

    conservation_data.append({
        "local_pos"    : local_pos,
        "res_name_3"   : resname_3,
        "res_name_1"   : AA_MAP_3TO1.get(resname_3, "?"),
        "conservation" : score,
        "is_primary"   : local_pos in PRIMARY_HOTSPOTS,
        **{f"aa_{org.replace('-','_').replace('/','_')}": aa_per_org[i]
           for i, org in enumerate(orgs)}
    })

# ── Step 8: Summary ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 8: Conservation summary")
print("=" * 60)

conserved     = [d for d in conservation_data if d["conservation"] >= 0.8]
not_conserved = [d for d in conservation_data if d["conservation"] <  0.8]

print(f"\n  Total hotspots analyzed : {len(conservation_data)}")
print(f"  Conserved (>= 0.8)      : {len(conserved)}")
print(f"  Variable  (<  0.8)      : {len(not_conserved)}")

print(f"\n  CONSERVED hotspots:")
for d in sorted(conserved, key=lambda x: -x["conservation"]):
    prim = " ★ PRIMARY" if d["is_primary"] else ""
    print(f"    {d['res_name_3']}{d['local_pos']:4d}  score={d['conservation']:.3f}{prim}")

if not_conserved:
    print(f"\n  VARIABLE hotspots (<0.8):")
    for d in sorted(not_conserved, key=lambda x: x["conservation"]):
        prim = " ★ PRIMARY" if d["is_primary"] else ""
        print(f"    {d['res_name_3']}{d['local_pos']:4d}  score={d['conservation']:.3f}{prim}")

# Primary hotspot summary
print(f"\n  PRIMARY hotspot conservation (salt bridge anchors + top contacts):")
for d in conservation_data:
    if d["is_primary"]:
        note = ""
        if d["local_pos"] in {414, 580, 583}:
            note = " [salt bridge]"
        elif d["local_pos"] in {480, 482}:
            note = " [top contact, 60 each]"
        elif d["local_pos"] in {477, 579, 584}:
            note = " [electrostatic anchor]"
        print(f"    {d['res_name_3']}{d['local_pos']:4d}  "
              f"score={d['conservation']:.3f}{note}")

# ── Save CSV ───────────────────────────────────────────────────────────────
df = pd.DataFrame(conservation_data)
csv_path = f"{VAL_DIR}/conservation_NSP13.csv"
df.to_csv(csv_path, index=False)
print(f"\n  Saved CSV: {csv_path}")

# ── Save JSON ──────────────────────────────────────────────────────────────
summary = {
    "complex"                : "NSP13-Helicase",
    "script"                 : "06_conservation_NSP13-Helicase_7.py",
    "coronaviruses"          : list(CORONAVIRUSES.keys()),
    "hotspots_total"         : len(conservation_data),
    "hotspots_conserved_08"  : len(conserved),
    "hotspots_variable"      : len(not_conserved),
    "primary_hotspots"       : sorted(PRIMARY_HOTSPOTS),
    "conservation_data"      : conservation_data,
    "salt_bridge_anchors"    : {
        "LYS414" : next((d["conservation"] for d in conservation_data
                         if d["local_pos"]==414), None),
        "ASP580" : next((d["conservation"] for d in conservation_data
                         if d["local_pos"]==580), None),
        "ASP583" : next((d["conservation"] for d in conservation_data
                         if d["local_pos"]==583), None),
    }
}

json_path = f"{VAL_DIR}/conservation_summary_7.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Saved JSON: {json_path}")
