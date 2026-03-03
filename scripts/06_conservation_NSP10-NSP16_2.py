"""
Script 06_2: Conservation Analysis — NSP10-NSP16
=================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

Computes per-residue conservation across 5 coronaviruses
for all hotspot residues identified in Script 05_2.
Uses Shannon entropy on MUSCLE-aligned FASTA files.

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/06_conservation_NSP10-NSP16_2.py
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import AlignIO

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
RES_DIR = PROJECT / "02-validation" / "NSP10-NSP16"
RES_DIR.mkdir(parents=True, exist_ok=True)

CORONAVIRUSES = ["SARS-CoV-2","SARS-CoV-1","MERS-CoV",
                 "HCoV-229E","HCoV-NL63"]

# Hotspots from Script 05_2 consensus
HOTSPOTS_NSP10 = [5,40,42,43,44,45,71,76,78,80,93,94,95,96]
HOTSPOTS_NSP16 = [40,41,44,76,83,87,102,104,106,244,247]

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def shannon_conservation(column):
    """
    Conservation = 1 - normalized Shannon entropy.
    1.0 = identical in all sequences.
    0.0 = maximum diversity.
    """
    aas = [aa for aa in column if aa not in ("-","X","?")]
    if not aas:
        return 0.0
    counts = {}
    for aa in aas:
        counts[aa] = counts.get(aa, 0) + 1
    n = len(aas)
    probs = [c/n for c in counts.values()]
    entropy = -sum(p * np.log2(p) for p in probs if p > 0)
    max_entropy = np.log2(min(n, 20))
    if max_entropy == 0:
        return 1.0
    return round(1.0 - entropy / max_entropy, 4)


def analyze_conservation(nsp, aln_file, hotspot_positions):
    """
    Compute conservation for all positions, flag hotspots.
    Returns DataFrame with columns:
      position, aa_SARS2, conservation, is_hotspot, aa_per_species
    """
    aln = AlignIO.read(aln_file, "fasta")

    # Find SARS-CoV-2 reference
    ref = next((r for r in aln
                if "SARS-CoV-2" in r.id or
                "SARS_CoV_2" in r.id), aln[0])

    # Map alignment column → sequence position
    ref_pos, pos_map = 0, {}
    for col_i, aa in enumerate(str(ref.seq)):
        if aa != "-":
            ref_pos += 1
            pos_map[ref_pos] = col_i

    rows = []
    for seq_pos in range(1, ref_pos + 1):
        if seq_pos not in pos_map:
            continue
        col_i   = pos_map[seq_pos]
        column  = [str(r.seq)[col_i] for r in aln]
        ref_aa  = str(ref.seq)[col_i]
        cons    = shannon_conservation(column)

        # AA per species — exact match to avoid SARS-CoV-2/1 collision
        aa_species = {}
        for r in aln:
            rid = r.id.replace("_"," ")
            for cov in CORONAVIRUSES:
                if cov in rid or cov.replace("-","_") in r.id:
                    # Don't overwrite if already assigned
                    if cov not in aa_species:
                        aa_species[cov] = str(r.seq)[col_i]

        rows.append({
            "position":    seq_pos,
            "aa_SARS2":    ref_aa,
            "conservation": cons,
            "is_hotspot":  seq_pos in hotspot_positions,
            **{f"aa_{c.replace('-','_').replace(' ','_')}":
               aa_species.get(c, "?")
               for c in CORONAVIRUSES},
        })

    return pd.DataFrame(rows)


def main():
    print("\n" + "="*57)
    print("  Script 06_2: Conservation Analysis — NSP10-NSP16")
    print("="*57 + "\n")

    results = {}

    for nsp, aln_file, hotspots in [
        ("NSP10",
         SEQ_DIR / "NSP10_aligned.fasta",
         HOTSPOTS_NSP10),
        ("NSP16",
         SEQ_DIR / "NSP16_aligned.fasta",
         HOTSPOTS_NSP16),
    ]:
        print(f"  Analyzing {nsp}...")
        df = analyze_conservation(nsp, aln_file, hotspots)

        hotspot_df = df[df["is_hotspot"]].copy()
        conserved  = hotspot_df[
            hotspot_df["conservation"] >= 0.8]

        print(f"  {'Pos':<6} {'AA':<4} {'Cons':>6}  "
              f"{'Hotspot':>8}  Species AAs")
        print(f"  {'-'*6} {'-'*4} {'-'*6}  "
              f"{'-'*8}  {'-'*30}")
        for _, row in hotspot_df.sort_values(
                "conservation", ascending=False).iterrows():
            flag  = " ✅" if row["conservation"] >= 0.8 else ""
            aas   = " ".join(
                f"{c[:4]}:{row['aa_'+c.replace('-','_').replace(' ','_')]}"
                for c in CORONAVIRUSES)
            print(f"  {int(row['position']):<6} "
                  f"{row['aa_SARS2']:<4} "
                  f"{row['conservation']:>6.3f}  "
                  f"{'★' if row['position'] in {93,95,106,102} else ' ':>8}  "
                  f"{aas}{flag}")

        print(f"\n  {nsp}: {len(conserved)}/{len(hotspot_df)} "
              f"hotspots conserved ≥ 0.8\n")

        # Save CSV
        out_csv = RES_DIR / f"conservation_{nsp}.csv"
        df.to_csv(out_csv, index=False)
        print(f"  Saved: conservation_{nsp}.csv")
        results[nsp] = df

    # Summary
    print("\n  ── Primary drug target residues ──")
    primaries = {
        "NSP10": [(93,"LYS93","salt bridge anchor"),
                  (95,"LYS95","salt bridge anchor"),
                  (80,"HIS80","salt bridge with ASP102")],
        "NSP16": [(106,"ASP106","salt bridge anchor x2"),
                  (102,"ASP102","salt bridge with HIS80")],
    }
    for nsp, plist in primaries.items():
        df = results[nsp]
        for pos, name, role in plist:
            row = df[df["position"]==pos]
            if not row.empty:
                cons = row.iloc[0]["conservation"]
                flag = "✅" if cons >= 0.8 else "⚠️"
                print(f"  {flag} {name:<10} ({nsp}) "
                      f"cons={cons:.3f}  [{role}]")

    print("\n" + "="*57 + "\n")


if __name__ == "__main__":
    main()
