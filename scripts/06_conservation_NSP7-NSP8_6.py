"""
Script 06_6: Conservation Analysis — NSP7-NSP8
===============================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Analyzes conservation for BOTH binding modes:
  Mode A (crystal): NSP8 res 163,178-180 / NSP7 res 24,26,27
  Mode B (AF3):     NSP8 res 84-120,150,190 / NSP7 res 2-76

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/06_conservation_NSP7-NSP8_6.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from Bio import AlignIO

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "02-validation" / "NSP7-NSP8"
RES_DIR.mkdir(parents=True, exist_ok=True)

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

# Mode A — crystal hotspots
HOTSPOTS_NSP8_A = [163, 178, 179, 180]
HOTSPOTS_NSP7_A = [24, 26, 27]

# Mode B — AF3 hotspots (top 15 each)
HOTSPOTS_NSP8_B = [84,87,89,90,91,92,94,95,96,
                    98,102,103,106,107,110,111,
                    116,119,120,150,190]
HOTSPOTS_NSP7_B = [2,6,12,13,16,19,28,35,49,50,
                    52,53,54,56,57,58,59,60,66,
                    68,69,71,74,75,76]

# Primary SB — AF3 Mode B only
PRIMARY_NSP8    = {190}   # ARG190
PRIMARY_NSP7    = {50}    # GLU50


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


def analyze(aln_file, hotspots_A, hotspots_B,
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

    all_pos = sorted(set(hotspots_A) |
                     set(hotspots_B) | primary_set)
    rows    = []
    for pos in all_pos:
        if pos not in pos_map:
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
        in_A  = pos in hotspots_A
        in_B  = pos in hotspots_B
        sb    = pos in primary_set
        mode  = ("AB" if in_A and in_B
                 else "A" if in_A else "B")
        gate  = "✅" if cons >= 0.8 else "  "
        sbf   = "★" if sb else " "

        rows.append({
            "position":     pos,
            "aa_SARS2":     sars2,
            "conservation": round(cons,3),
            "mode":         mode,
            "is_hotspot":   True,
            "primary_sb":   sb,
            "SARS-CoV-2":   sars2,
            "SARS-CoV-1":   aa_species.get(
                             "SARS-CoV-1","?"),
            "MERS-CoV":     aa_species.get(
                             "MERS-CoV","?"),
            "HCoV-229E":    aa_species.get(
                             "HCoV-229E","?"),
            "HCoV-NL63":    aa_species.get(
                             "HCoV-NL63","?"),
        })
        print(f"  {pos:<6} {sars2}  "
              f"{cons:.3f} {sbf} [{mode}]  "
              f"SARS2:{sars2} SARS1:"
              f"{aa_species.get('SARS-CoV-1','?')} "
              f"MERS:"
              f"{aa_species.get('MERS-CoV','?')} "
              f"229E:"
              f"{aa_species.get('HCoV-229E','?')} "
              f"NL63:"
              f"{aa_species.get('HCoV-NL63','?')} "
              f"{gate}")

    df     = pd.DataFrame(rows)
    # Stats per mode
    for mode_lbl, hs in [("A", hotspots_A),
                          ("B", hotspots_B)]:
        sub = df[df["position"].isin(hs)]
        n_c = len(sub[sub["conservation"] >= 0.8])
        print(f"\n  {nsp_label} Mode {mode_lbl}: "
              f"{n_c}/{len(sub)} hotspots ≥ 0.8")
    return df


def main():
    print("\n" + "="*60)
    print("  Script 06_6: Conservation — NSP7-NSP8")
    print("  Dual mode: A=crystal, B=AF3")
    print("="*60 + "\n")

    print("  NSP8 conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6} "
          f"{'SB'} [Mode]  Species")
    print(f"  {'-'*58}")
    df8 = analyze(
        SEQ_DIR/"NSP8_aligned.fasta",
        HOTSPOTS_NSP8_A, HOTSPOTS_NSP8_B,
        PRIMARY_NSP8, "NSP8")

    print(f"\n  NSP7 conservation:")
    print(f"  {'Pos':<6} {'AA'} {'Cons':>6} "
          f"{'SB'} [Mode]  Species")
    print(f"  {'-'*58}")
    df7 = analyze(
        SEQ_DIR/"NSP7_aligned.fasta",
        HOTSPOTS_NSP7_A, HOTSPOTS_NSP7_B,
        PRIMARY_NSP7, "NSP7")

    print(f"\n  ── AF3 SB residues ──")
    for pos in sorted(PRIMARY_NSP8):
        row = df8[df8["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            a = row.iloc[0]["aa_SARS2"]
            print(f"  {a}{pos} (NSP8) cons={c:.3f}")
    for pos in sorted(PRIMARY_NSP7):
        row = df7[df7["position"]==pos]
        if not row.empty:
            c = row.iloc[0]["conservation"]
            a = row.iloc[0]["aa_SARS2"]
            print(f"  {a}{pos} (NSP7) cons={c:.3f}")

    df8.to_csv(RES_DIR/"conservation_NSP8.csv",
               index=False)
    df7.to_csv(RES_DIR/"conservation_NSP7.csv",
               index=False)
    print(f"\n  Saved: conservation_NSP8.csv")
    print(f"  Saved: conservation_NSP7.csv")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
