"""
Script 06_5: Conservation Analysis — NSP9-NSP12
================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/06_conservation_NSP9-NSP12_5.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from Bio import AlignIO

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "02-validation" / "NSP9-NSP12"
RES_DIR.mkdir(parents=True, exist_ok=True)

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

# Consensus hotspots from Script 05_5
HOTSPOTS_NSP12 = [38,1,3,4,96,733,202,103,
                   221,233,291,2,223]
HOTSPOTS_NSP9  = [38,1,3,4,96,103,2,97]

# AF3-predicted SB anchors
PRIMARY_NSP12  = {740, 744}
PRIMARY_NSP9   = {36}


def shannon_entropy(column):
    col  = [aa for aa in column if aa != "-"]
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
    all_pos = sorted(set(hotspots) | primary_set)
    for pos in all_pos:
        if pos not in pos_map:
            print(f"  ⚠️  {nsp_label} pos {pos} "
                  f"not in alignment — skipping")
            continue
        col_i = pos_map[pos]
        cons  = shannon_entropy(
            [str(r.seq)[col_i] for r in aln])

        aa_species = {}
        for r in aln:
            rid = r.id.replace("_"," ")
            for cov in CORONAVIRUSES:
                if cov in rid and \
                        cov not in aa_species:
                    aa_species[cov] = \
                        str(r.seq)[col_i]

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
    hs_df   = df[df["is_hotspot"]]
    n_cons  = len(hs_df[hs_df["conservation"] >= 0.8])
    n_total = len(hs_df)
    print(f"\n  {nsp_label}: {n_cons}/{n_total} "
          f"hotspots conserved ≥ 0.8")
    return df


def main():
    print("\n" + "="*57)
    print("  Script 06_5: Conservation — NSP9-NSP12")
    print("="*57 + "\n")

    print("  NSP12 hotspot conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6}  "
          f"{'SB'}  Species breakdown")
    print(f"  {'-'*56}")
    df12 = analyze_conservation(
        SEQ_DIR / "NSP12_aligned.fasta",
        HOTSPOTS_NSP12, PRIMARY_NSP12, "NSP12")

    print(f"\n  NSP9 hotspot conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6}  "
          f"{'SB'}  Species breakdown")
    print(f"  {'-'*56}")
    df9  = analyze_conservation(
        SEQ_DIR / "NSP9_aligned.fasta",
        HOTSPOTS_NSP9, PRIMARY_NSP9, "NSP9")

    print(f"\n  ── AF3-predicted salt bridge residues ──")
    for pos in sorted(PRIMARY_NSP12):
        row = df12[df12["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            a = row.iloc[0]["aa_SARS2"]
            print(f"  {a}{pos:<8} (NSP12) "
                  f"cons={c:.3f}  "
                  f"[AF3-predicted SB anchor]")
    for pos in sorted(PRIMARY_NSP9):
        row = df9[df9["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            a = row.iloc[0]["aa_SARS2"]
            print(f"  {a}{pos:<8} (NSP9)  "
                  f"cons={c:.3f}  "
                  f"[AF3-predicted SB anchor]")

    df12.to_csv(RES_DIR / "conservation_NSP12.csv",
                index=False)
    df9.to_csv( RES_DIR / "conservation_NSP9.csv",
                index=False)
    print(f"\n  Saved: conservation_NSP12.csv")
    print(f"  Saved: conservation_NSP9.csv")
    print("="*57 + "\n")


if __name__ == "__main__":
    main()
