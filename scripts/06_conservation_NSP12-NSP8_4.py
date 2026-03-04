"""
Script 06_4: Conservation Analysis — NSP12-NSP8
================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/06_conservation_NSP12-NSP8_4.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from Bio import AlignIO

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "02-validation" / "NSP12-NSP8"
RES_DIR.mkdir(parents=True, exist_ok=True)

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

HOTSPOTS_NSP12 = [387,129,389,271,330,131,380,523,
                   91,87,332,95,117,517,99]
HOTSPOTS_NSP8  = [117,129,80,115,131,112,91,87,
                   116,95,98,83,128,90,121]

PRIMARY_NSP12  = {523, 332, 517}
PRIMARY_NSP8   = {80, 99, 79}


def shannon_entropy(column):
    col   = [aa for aa in column if aa != "-"]
    if not col:
        return 0.0
    freqs = {}
    for aa in col:
        freqs[aa] = freqs.get(aa,0) + 1
    n     = len(col)
    H     = -sum((c/n)*np.log2(c/n)
                 for c in freqs.values())
    max_H = np.log2(min(n,20))
    return 1.0 - H/max_H if max_H > 0 else 1.0


def analyze_conservation(aln_file, hotspots,
                          primary_set, nsp_label):
    aln     = AlignIO.read(aln_file, "fasta")
    ref_seq = next(
        r for r in aln
        if "SARS-CoV-2" in r.id.replace("_"," "))

    pos_map, ref_pos = {}, 0
    for col_i, aa in enumerate(str(ref_seq.seq)):
        if aa != "-":
            ref_pos += 1
            pos_map[ref_pos] = col_i

    rows = []
    for pos in sorted(set(hotspots) | primary_set):
        if pos not in pos_map:
            continue
        col_i = pos_map[pos]
        cons  = shannon_entropy(
            [str(r.seq)[col_i] for r in aln])

        aa_species = {}
        for r in aln:
            rid = r.id.replace("_"," ")
            for cov in CORONAVIRUSES:
                if cov in rid and cov not in aa_species:
                    aa_species[cov] = str(r.seq)[col_i]

        sars2 = aa_species.get("SARS-CoV-2","?")
        sars1 = aa_species.get("SARS-CoV-1","?")
        mers  = aa_species.get("MERS-CoV",  "?")
        hcov1 = aa_species.get("HCoV-229E", "?")
        hcov2 = aa_species.get("HCoV-NL63", "?")
        sb_flag = "★" if pos in primary_set else " "
        gate    = "✅" if cons >= 0.8 else "  "

        print(f"  {pos:<6} {sars2}  {cons:.3f}  "
              f"{sb_flag}  "
              f"SARS2:{sars2} SARS1:{sars1} "
              f"MERS:{mers} 229E:{hcov1} "
              f"NL63:{hcov2} {gate}")

        rows.append({
            "position":     pos,
            "aa_SARS2":     sars2,
            "conservation": round(cons,3),
            "is_hotspot":   pos in hotspots,
            "primary_sb":   pos in primary_set,
            "SARS-CoV-2":   sars2,
            "SARS-CoV-1":   sars1,
            "MERS-CoV":     mers,
            "HCoV-229E":    hcov1,
            "HCoV-NL63":    hcov2,
        })

    df      = pd.DataFrame(rows)
    n_cons  = len(df[df["conservation"] >= 0.8])
    n_total = len(df[df["is_hotspot"]])
    print(f"\n  {nsp_label}: {n_cons}/{n_total} "
          f"hotspots conserved ≥ 0.8")
    return df


def main():
    print("\n" + "="*57)
    print("  Script 06_4: Conservation — NSP12-NSP8")
    print("="*57 + "\n")

    print("  NSP12 hotspot conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6}  "
          f"{'SB'}  Species breakdown")
    print(f"  {'-'*55}")
    df12 = analyze_conservation(
        SEQ_DIR / "NSP12_aligned.fasta",
        HOTSPOTS_NSP12, PRIMARY_NSP12, "NSP12")

    print(f"\n  NSP8 hotspot conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6}  "
          f"{'SB'}  Species breakdown")
    print(f"  {'-'*55}")
    df8  = analyze_conservation(
        SEQ_DIR / "NSP8_aligned.fasta",
        HOTSPOTS_NSP8, PRIMARY_NSP8, "NSP8")

    print(f"\n  ── Primary salt bridge residues ──")
    for pos in sorted(PRIMARY_NSP12):
        row = df12[df12["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            a = row.iloc[0]["aa_SARS2"]
            print(f"  {a}{pos:<8} (NSP12) "
                  f"cons={c:.3f}  [salt bridge anchor]")
    for pos in sorted(PRIMARY_NSP8):
        row = df8[df8["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            a = row.iloc[0]["aa_SARS2"]
            print(f"  {a}{pos:<8} (NSP8)  "
                  f"cons={c:.3f}  [salt bridge anchor]")

    df12.to_csv(RES_DIR / "conservation_NSP12.csv",
                index=False)
    df8.to_csv( RES_DIR / "conservation_NSP8.csv",
                index=False)
    print(f"\n  Saved: conservation_NSP12.csv")
    print(f"  Saved: conservation_NSP8.csv")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
