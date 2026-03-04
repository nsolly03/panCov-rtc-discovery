"""
Script 06_3: Conservation Analysis — NSP12-NSP7
================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Sequences: NSP12_aligned.fasta, NSP7_aligned.fasta
Method: Shannon entropy on MUSCLE-aligned FASTA (5 coronaviruses)
Hotspots from Script 05_3 interface analysis

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/06_conservation_NSP12-NSP7_3.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from Bio import AlignIO

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "02-validation" / "NSP12-NSP7"
RES_DIR.mkdir(parents=True, exist_ok=True)

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

# Consensus hotspots from Script 05_3
HOTSPOTS_NSP12 = [440,412,442,443,420,843,409,
                   40,33,41,37,413,415,14,23]
HOTSPOTS_NSP7  = [40,14,33,41,37,11,23,5,15,29,12,4,1]

# Primary salt bridge residues
PRIMARY_NSP12 = {431}   # GLU431 — LYS2(NSP7)
PRIMARY_NSP7  = {2}     # LYS2   — GLU431(NSP12)


def shannon_entropy(column):
    """Shannon entropy of an alignment column."""
    col   = [aa for aa in column if aa != "-"]
    if not col:
        return 0.0
    freqs = {}
    for aa in col:
        freqs[aa] = freqs.get(aa, 0) + 1
    n     = len(col)
    H     = -sum((c/n)*np.log2(c/n)
                 for c in freqs.values())
    max_H = np.log2(min(n, 20))
    return 1.0 - H/max_H if max_H > 0 else 1.0


def analyze_conservation(aln_file, hotspots,
                           primary_set, nsp_label):
    aln     = AlignIO.read(aln_file, "fasta")
    ref_seq = next(
        r for r in aln
        if "SARS-CoV-2" in r.id.replace("_"," "))

    # Build position map: local 1-based -> alignment col
    pos_map  = {}
    ref_pos  = 0
    for col_i, aa in enumerate(str(ref_seq.seq)):
        if aa != "-":
            ref_pos += 1
            pos_map[ref_pos] = col_i

    rows = []
    for pos in sorted(set(hotspots) | primary_set):
        if pos not in pos_map:
            continue
        col_i = pos_map[pos]
        col   = str(ref_seq.seq)[col_i]
        cons  = shannon_entropy(
            [str(r.seq)[col_i] for r in aln])

        # Per-species AA
        aa_species = {}
        for r in aln:
            rid = r.id.replace("_"," ")
            for cov in CORONAVIRUSES:
                if cov in rid:
                    if cov not in aa_species:
                        aa_species[cov] = \
                            str(r.seq)[col_i]

        sars2 = aa_species.get("SARS-CoV-2","?")
        sars1 = aa_species.get("SARS-CoV-1","?")
        mers  = aa_species.get("MERS-CoV",  "?")
        hcov1 = aa_species.get("HCoV-229E", "?")
        hcov2 = aa_species.get("HCoV-NL63", "?")

        sb_flag = "★" if pos in primary_set else " "
        gate    = "✅" if cons >= 0.8 else "  "

        rows.append({
            "position":     pos,
            "aa_SARS2":     sars2,
            "conservation": round(cons, 3),
            "is_hotspot":   pos in hotspots,
            "primary_sb":   pos in primary_set,
            "SARS-CoV-2":   sars2,
            "SARS-CoV-1":   sars1,
            "MERS-CoV":     mers,
            "HCoV-229E":    hcov1,
            "HCoV-NL63":    hcov2,
        })

        print(f"  {pos:<6} {sars2}  {cons:.3f}  "
              f"      {sb_flag}  "
              f"SARS:{sars2} SARS:{sars1} "
              f"MERS:{mers} HCoV:{hcov1} "
              f"HCoV:{hcov2} {gate}")

    df      = pd.DataFrame(rows)
    n_cons  = len(df[df["conservation"] >= 0.8])
    n_total = len(df[df["is_hotspot"]])
    print(f"\n  {nsp_label}: {n_cons}/{n_total} "
          f"hotspots conserved ≥ 0.8")
    return df


def main():
    print("\n" + "="*57)
    print("  Script 06_3: Conservation — NSP12-NSP7")
    print("="*57 + "\n")

    print("  NSP12 hotspot conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6}  "
          f"     {'SB'}  Species breakdown")
    print(f"  {'-'*55}")
    df12 = analyze_conservation(
        SEQ_DIR / "NSP12_aligned.fasta",
        HOTSPOTS_NSP12, PRIMARY_NSP12, "NSP12")

    print(f"\n  NSP7 hotspot conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6}  "
          f"     {'SB'}  Species breakdown")
    print(f"  {'-'*55}")
    df7  = analyze_conservation(
        SEQ_DIR / "NSP7_aligned.fasta",
        HOTSPOTS_NSP7, PRIMARY_NSP7, "NSP7")

    # Primary salt bridge residues summary
    print(f"\n  ── Primary salt bridge residues ──")
    for pos in sorted(PRIMARY_NSP12):
        row = df12[df12["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            print(f"  {'GLU'+str(pos):<10} (NSP12) "
                  f"cons={c:.3f}  [salt bridge anchor]")
    for pos in sorted(PRIMARY_NSP7):
        row = df7[df7["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            print(f"  {'LYS'+str(pos):<10} (NSP7)  "
                  f"cons={c:.3f}  [salt bridge anchor]")

    # Save
    df12.to_csv(RES_DIR / "conservation_NSP12.csv",
                index=False)
    df7.to_csv( RES_DIR / "conservation_NSP7.csv",
                index=False)
    print(f"\n  Saved: conservation_NSP12.csv")
    print(f"  Saved: conservation_NSP7.csv")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
