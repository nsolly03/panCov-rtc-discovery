"""
Script 06: Conservation Analysis — NSP10-NSP14
===============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Downloads pp1ab polyproteins from UniProt for 5 coronaviruses,
  slices NSP10 and NSP14 using verified coordinates,
  runs MUSCLE alignment, calculates conservation scores,
  and maps scores onto consensus hotspot residues from Script 05.

Coronaviruses: SARS-CoV-2, SARS-CoV-1, MERS-CoV, HCoV-229E, HCoV-NL63

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/06_conservation_NSP10-NSP14.py
"""

import json
import math
import subprocess
import time
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP14"
SEQ_DIR.mkdir(parents=True, exist_ok=True)
RES_DIR.mkdir(parents=True, exist_ok=True)

# ── Verified UniProt accessions and NSP coordinates ───────
CORONAVIRUSES = {
    "SARS-CoV-2": {
        "accession": "P0DTD1",
        "NSP10": (4254, 4392),
        "NSP14": (5926, 6452),
    },
    "SARS-CoV-1": {
        "accession": "P0C6X7",
        "NSP10": (4231, 4369),
        "NSP14": (5903, 6429),
    },
    "MERS-CoV": {
        "accession": "K9N7C7",
        "NSP10": (4238, 4377),
        "NSP14": (5909, 6432),
    },
    "HCoV-229E": {
        "accession": "P0C6X1",
        "NSP10": (3934, 4068),
        "NSP14": (5593, 6110),
    },
    "HCoV-NL63": {
        "accession": "P0C6X5",
        "NSP10": (3909, 4043),
        "NSP14": (5568, 6085),
    },
}

# ── Consensus hotspots from Script 05 ─────────────────────
HOTSPOTS = {
    "NSP10": [3,4,5,12,16,18,19,21,33,40,42,44,45,80,93,96],
    "NSP14": [4,5,7,8,9,10,20,21,24,25,26,27,41,69,
              126,127,128,131,201],
}


def fetch_polyprotein(cov, accession):
    """Download pp1ab from UniProt"""
    cached = SEQ_DIR / f"{cov.replace(' ','_')}_{accession}.fasta"
    if cached.exists():
        records = list(SeqIO.parse(cached, "fasta"))
        if records and len(records[0].seq) > 1000:
            print(f"  SKIP  {cov} (cached, "
                  f"{len(records[0].seq)} aa)")
            return str(records[0].seq)

    print(f"  Downloading {cov} ({accession})...")
    url = (f"https://rest.uniprot.org/uniprotkb/"
           f"{accession}.fasta")
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r:
            fasta = r.read().decode()
        cached.write_text(fasta)
        lines = fasta.strip().split("\n")
        seq   = "".join(lines[1:])
        print(f"  OK    {cov} — {len(seq)} aa")
        time.sleep(0.4)
        return seq
    except Exception as e:
        print(f"  FAIL  {cov}: {e}")
        return None


def run_muscle(fasta_in, fasta_out):
    """Run MUSCLE 5 alignment"""
    cmd = ["muscle",
           "-align",  str(fasta_in),
           "-output", str(fasta_out)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  MUSCLE stderr: {result.stderr[:300]}")
    return result.returncode == 0


def shannon_entropy(column):
    counts = {}
    for aa in column:
        if aa == "-":
            continue
        counts[aa] = counts.get(aa, 0) + 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def conservation_score(entropy):
    return round(1.0 - (entropy / math.log2(20)), 4)


def analyze_conservation(aln_file, hotspots):
    """Compute per-position conservation scores"""
    records = list(SeqIO.parse(aln_file, "fasta"))
    if not records:
        return None

    # Use SARS-CoV-2 as reference
    ref = next(
        (r for r in records if "SARS-CoV-2" in r.id),
        records[0])

    aln_len = len(records[0].seq)
    matrix  = [[str(r.seq[i]) for r in records]
               for i in range(aln_len)]

    ref_pos = 0
    rows    = []
    for i, col in enumerate(matrix):
        ref_aa = str(ref.seq[i])
        if ref_aa != "-":
            ref_pos += 1
            entropy = shannon_entropy(col)
            score   = conservation_score(entropy)
            aas     = [aa for aa in col if aa != "-"]
            rows.append({
                "position":     ref_pos,
                "aa_SARS2":     ref_aa,
                "conservation": score,
                "entropy":      entropy,
                "is_hotspot":   ref_pos in hotspots,
                "all_aas":      "".join(aas),
                "n_sequences":  len(aas),
            })

    return pd.DataFrame(rows)


def print_hotspot_table(df, nsp_name):
    hdf = df[df["is_hotspot"]].sort_values(
        "conservation", ascending=False)
    print(f"\n  {nsp_name} hotspot conservation "
          f"({len(hdf)} residues, "
          f"{df['n_sequences'].max()} coronaviruses):")
    print(f"  {'Pos':<5} {'AA':<4} {'Score':>8} "
          f"{'Entropy':>8}  All AAs")
    print(f"  {'-'*5} {'-'*4} {'-'*8} {'-'*8}  {'-'*15}")
    for _, r in hdf.iterrows():
        flag = " ⭐" if r["conservation"] >= 0.8 else ""
        print(f"  {int(r['position']):<5} "
              f"{r['aa_SARS2']:<4} "
              f"{r['conservation']:>8.3f} "
              f"{r['entropy']:>8.4f}  "
              f"{r['all_aas']}{flag}")
    conserved = hdf[
        hdf["conservation"] >= 0.8]["position"].tolist()
    print(f"\n  ⭐ Conserved hotspots (score ≥0.8): "
          f"{[int(p) for p in conserved]}")
    return [int(p) for p in conserved]


def main():
    print("\n" + "="*55)
    print("  Script 06: Conservation — NSP10-NSP14")
    print("  5 coronaviruses | Verified UniProt coordinates")
    print("="*55)

    # ── Step 1: Download polyproteins ─────────────────────
    print("\n  Step 1: Downloading polyproteins from UniProt\n")
    polyproteins = {}
    for cov, info in CORONAVIRUSES.items():
        seq = fetch_polyprotein(cov, info["accession"])
        if seq:
            polyproteins[cov] = seq

    print(f"\n  Downloaded: {len(polyproteins)}/5")

    # ── Step 2: Slice NSPs ────────────────────────────────
    print("\n  Step 2: Slicing NSP sequences\n")
    for nsp in ["NSP10", "NSP14"]:
        records = []
        for cov, poly_seq in polyproteins.items():
            start, end = CORONAVIRUSES[cov][nsp]
            seq    = poly_seq[start-1:end]
            exp_len = end - start + 1
            status = "✅" if abs(len(seq)-exp_len) < 5 else "⚠️"
            print(f"  {status} {cov:<20} {nsp}: "
                  f"{len(seq)} aa "
                  f"(expected ~{exp_len})")
            safe_id = cov.replace(" ","_")
            records.append(SeqRecord(
                Seq(seq), id=safe_id,
                description=f"{nsp}|{cov}"))

        out = SEQ_DIR / f"{nsp}_all_coronaviruses.fasta"
        SeqIO.write(records, out, "fasta")
        print(f"  Saved: {out.name}\n")

    # ── Step 3: MUSCLE alignment ──────────────────────────
    print("\n  Step 3: Running MUSCLE alignments\n")
    alignments = {}
    for nsp in ["NSP10", "NSP14"]:
        fasta_in  = SEQ_DIR / f"{nsp}_all_coronaviruses.fasta"
        fasta_out = SEQ_DIR / f"{nsp}_aligned.fasta"
        if fasta_out.exists():
            fasta_out.unlink()
        print(f"  Aligning {nsp} ({len(polyproteins)} sequences)...")
        ok = run_muscle(fasta_in, fasta_out)
        if ok:
            print(f"  OK    {fasta_out.name}")
            alignments[nsp] = fasta_out
        else:
            print(f"  FAIL  {nsp} alignment")

    # ── Step 4: Conservation scores ───────────────────────
    print("\n  Step 4: Calculating conservation scores")
    summary = {}

    for nsp, aln_file in alignments.items():
        print(f"\n{'='*55}")
        print(f"  {nsp} conservation results")
        print(f"{'='*55}")

        df = analyze_conservation(aln_file, HOTSPOTS[nsp])
        if df is None:
            continue

        csv_out = RES_DIR / f"conservation_{nsp}.csv"
        df.to_csv(csv_out, index=False)

        conserved = print_hotspot_table(df, nsp)
        summary[nsp] = {
            "n_coronaviruses":    len(polyproteins),
            "coronaviruses":      list(polyproteins.keys()),
            "conserved_hotspots": conserved,
            "hotspot_scores":     df[df["is_hotspot"]][
                ["position","aa_SARS2","conservation","all_aas"]
            ].to_dict("records"),
        }

    # ── Save summary ──────────────────────────────────────
    out = RES_DIR / "conservation_summary.json"
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*55}")
    print(f"  Summary saved: conservation_summary.json")
    print(f"  CSVs saved:    conservation_NSP10.csv")
    print(f"                 conservation_NSP14.csv")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()
