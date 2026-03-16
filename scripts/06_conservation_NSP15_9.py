"""
Script 06_9 — Conservation Analysis NSP15
Complex  : NSP15-NSP15 (NendoU homodimer)
Suffix   : _9
Note     : Homodimer — conservation of NSP15 sequence itself across 5 CoV
           Key question: are dimer interface hotspots MORE conserved
           than the rest of the protein surface?
"""

import json
import os
import numpy as np
import requests
from Bio import pairwise2
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

OUT_DIR  = "02-validation/NSP15"
OUT_CSV  = f"{OUT_DIR}/conservation_NSP15.csv"
OUT_JSON = f"{OUT_DIR}/conservation_summary_9.json"
os.makedirs(OUT_DIR, exist_ok=True)

# ── UniProt accessions and NSP15 coordinates ──────────────────
ACCESSIONS = {
    "SARS-CoV-2" : ("P0DTD1", 6453, 6798),
    "SARS-CoV-1" : ("P0C6X7", 6430, 6775),
    "MERS-CoV"   : ("K9N7C7", 6433, 6775),
    "HCoV-229E"  : ("P0C6X1", 6111, 6458),
    "HCoV-NL63"  : ("P0C6X5", 6086, 6429),
}

print("=" * 60)
print("Script 06_9 — Conservation Analysis NSP15")
print("=" * 60)

# ── Download sequences ─────────────────────────────────────────
print("\nDownloading NSP15 sequences from UniProt...")
sequences = {}

for org, (acc, start, end) in ACCESSIONS.items():
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    r = requests.get(url, timeout=30).json()
    seq = r["sequence"]["value"]
    nsp15 = seq[start-1:end]
    sequences[org] = nsp15
    print(f"  {org:<16} {acc}  {start}-{end}  {len(nsp15)} aa")

ref_seq  = sequences["SARS-CoV-2"]
ref_len  = len(ref_seq)
print(f"\nReference (SARS-CoV-2): {ref_len} aa")

# ── Load interface hotspots from Script 05_9 output ───────────
with open(f"{OUT_DIR}/interface_analysis_9.json") as f:
    iface = json.load(f)

# Combine A and B — homodimer so positions are from same sequence
consensus_A = set(iface["consensus_A"])
consensus_B = set(iface["consensus_B"])
# Convert PDB positions to UniProt positions (offset = PDB - 1)
OFFSET = 1
all_hotspots_pdb = sorted(consensus_A | consensus_B)
all_hotspots = sorted([p - OFFSET for p in all_hotspots_pdb])
print(f"PDB hotspot positions converted: offset=-{OFFSET}")
print(f"Interface hotspots loaded: {len(all_hotspots)} positions")

# ── Align all sequences to SARS-CoV-2 reference ───────────────
print("\nAligning sequences to SARS-CoV-2 reference...")
alignments = {"SARS-CoV-2": {i+1: aa for i,aa in enumerate(ref_seq)}}

for org, seq in sequences.items():
    if org == "SARS-CoV-2":
        continue
    alns = pairwise2.align.globalms(
        ref_seq, seq, 2, -1, -5, -0.5, one_alignment_only=True)
    aln_ref, aln_qry = alns[0][0], alns[0][1]

    pos_map = {}
    i_ref = 0
    for r_aa, q_aa in zip(aln_ref, aln_qry):
        if r_aa != '-':
            i_ref += 1
            pos_map[i_ref] = q_aa if q_aa != '-' else '-'
    alignments[org] = pos_map
    identity = sum(1 for p,aa in pos_map.items()
                   if aa == ref_seq[p-1]) / ref_len * 100
    print(f"  {org:<16}: {identity:.1f}% identity to SARS-CoV-2")

# ── Compute conservation score per position ────────────────────
print("\nComputing conservation scores...")
conservation = {}

for pos in range(1, ref_len + 1):
    ref_aa = ref_seq[pos-1]
    matches = sum(1 for org, aln in alignments.items()
                  if aln.get(pos, '-') == ref_aa)
    score = matches / len(alignments)
    conservation[pos] = {
        "score"  : round(score, 3),
        "ref_aa" : ref_aa,
        "in_interface": pos in all_hotspots,
        "per_org": {org: aln.get(pos,'-')
                    for org, aln in alignments.items()},
    }

# ── Interface vs non-interface conservation ────────────────────
iface_scores    = [conservation[p]["score"] for p in all_hotspots
                   if p in conservation]
surface_scores  = [conservation[p]["score"] for p in range(1, ref_len+1)
                   if p not in all_hotspots]

print(f"\nConservation summary:")
print(f"  Interface hotspots  ({len(iface_scores):3d} pos): "
      f"mean={np.mean(iface_scores):.3f}  "
      f"min={np.min(iface_scores):.3f}  "
      f"max={np.max(iface_scores):.3f}")
print(f"  Non-interface       ({len(surface_scores):3d} pos): "
      f"mean={np.mean(surface_scores):.3f}  "
      f"min={np.min(surface_scores):.3f}  "
      f"max={np.max(surface_scores):.3f}")

# ── Hotspot conservation detail ────────────────────────────────
print("\nHotspot conservation (interface positions, sorted by score):")
print("{:<6} {:<6} {:<8} {}".format("Pos","AA","Score","Per organism"))
print("-" * 70)

hotspot_results = []
for pos in sorted(all_hotspots):
    if pos not in conservation:
        continue
    c = conservation[pos]
    per_org = "  ".join(
        f"{org[:4]}:{aa}"
        for org, aa in c["per_org"].items()
    )
    conserved = "✓" if c["score"] >= 0.8 else " "
    print(f"{pos:<6} {c['ref_aa']:<6} {c['score']:<8.3f} {per_org}  {conserved}")
    hotspot_results.append({
        "position"    : pos,
        "ref_aa"      : c["ref_aa"],
        "score"       : c["score"],
        "conserved"   : c["score"] >= 0.8,
        "per_org"     : c["per_org"],
    })

n_conserved = sum(1 for h in hotspot_results if h["conserved"])
print(f"\nConserved hotspots (score >= 0.8): {n_conserved}/{len(hotspot_results)}")

# Primary SB residues specifically
print("\nPrimary salt bridge residue conservation:")
for pos, name in [(39,"ASP40_PDB40"),(90,"ARG91_PDB91"),(61,"ARG62_PDB62"),(266,"GLU267_PDB267")]:
    if pos in conservation:
        c = conservation[pos]
        per = "  ".join(f"{o[:4]}:{a}" for o,a in c["per_org"].items())
        print(f"  {name:<8} score={c['score']:.3f}  {per}")

# ── Save CSV ───────────────────────────────────────────────────
with open(OUT_CSV, "w") as f:
    orgs = list(ACCESSIONS.keys())
    f.write("position,ref_aa,conservation_score,in_interface,"
            + ",".join(orgs) + "\n")
    for pos in range(1, ref_len+1):
        c = conservation[pos]
        row = [str(pos), c["ref_aa"], str(c["score"]),
               str(c["in_interface"])]
        row += [c["per_org"].get(org,"-") for org in orgs]
        f.write(",".join(row) + "\n")

# ── Save JSON summary ──────────────────────────────────────────
summary = {
    "complex"              : "NSP15-homodimer",
    "suffix"               : "_9",
    "n_hotspots"           : len(hotspot_results),
    "n_conserved_hotspots" : n_conserved,
    "mean_iface_conservation"   : round(float(np.mean(iface_scores)),3),
    "mean_surface_conservation" : round(float(np.mean(surface_scores)),3),
    "primary_sb_conservation"   : {
        "ASP40" : conservation.get(40,{}).get("score"),
        "ARG91" : conservation.get(91,{}).get("score"),
        "ARG62" : conservation.get(62,{}).get("score"),
        "GLU267": conservation.get(267,{}).get("score"),
    },
    "hotspot_details"      : hotspot_results,
    "note"                 : (
        "Homodimer: interface positions from both chain A and B "
        "map to same NSP15 sequence. Conservation compared between "
        "interface hotspots and rest of protein surface."
    )
}

with open(OUT_JSON, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nSaved: {OUT_CSV}")
print(f"Saved: {OUT_JSON}")
print("Next : Script 07_9 (pocket detection)")
