"""
Script 03: Download SARS-CoV-2 NSP sequences from UniProt
==========================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Downloads the full SARS-CoV-2 polyprotein (P0DTD1) from UniProt
  then slices it into individual NSP sequences using known
  amino acid positions. This gives complete sequences for AF3.

  Output: one FASTA file per NSP in 00-reference/sequences/
  These are ready to paste directly into AlphaFold3 server.

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/03_download_sequences.py
"""

import requests
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences"
SEQ_DIR.mkdir(parents=True, exist_ok=True)

# ── NSP positions in P0DTD1 polyprotein (UniProt) ─────────
NSP_INFO = {
    "NSP7":  {"start": 3860, "end": 3942,
               "function": "Cofactor for RdRp"},
    "NSP8":  {"start": 3943, "end": 4140,
               "function": "Primase cofactor"},
    "NSP9":  {"start": 4141, "end": 4253,
               "function": "RNA-binding protein"},
    "NSP10": {"start": 4254, "end": 4392,
               "function": "Exonuclease cofactor"},
    "NSP12": {"start": 4393, "end": 5324,
               "function": "RNA-dependent RNA polymerase"},
    "NSP13": {"start": 5325, "end": 5925,
               "function": "Helicase"},
    "NSP14": {"start": 5926, "end": 6452,
               "function": "Exonuclease / NMTase"},
    "NSP16": {"start": 6799, "end": 7096,
               "function": "2-O-Methyltransferase"},
}


def fetch_polyprotein():
    """Download full P0DTD1 polyprotein from UniProt"""
    cached = SEQ_DIR / "P0DTD1_polyprotein.fasta"
    if cached.exists():
        print("  SKIP  P0DTD1 polyprotein (already downloaded)")
        lines = cached.read_text().strip().split("\n")
        return "".join(lines[1:])

    print("  Downloading P0DTD1 polyprotein from UniProt...")
    url = "https://rest.uniprot.org/uniprotkb/P0DTD1.fasta"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        cached.write_text(r.text)
        lines = r.text.strip().split("\n")
        seq = "".join(lines[1:])
        print(f"  OK    P0DTD1 downloaded ({len(seq)} aa total)")
        return seq
    except Exception as e:
        print(f"  FAIL  {e}")
        return None


def main():
    print("\n" + "="*55)
    print("  Script 03: Downloading NSP sequences from UniProt")
    print("  Saving to: 00-reference/sequences/")
    print("="*55 + "\n")

    # Step 1 — get full polyprotein
    full_seq = fetch_polyprotein()
    if not full_seq:
        print("Cannot continue without polyprotein sequence.")
        return

    print(f"\n  Slicing into individual NSP sequences...\n")

    # Step 2 — slice each NSP
    records = []
    for name, info in NSP_INFO.items():
        start = info["start"] - 1   # convert to 0-indexed
        end   = info["end"]
        seq   = full_seq[start:end]
        length = len(seq)

        record = SeqRecord(
            Seq(seq),
            id=name,
            description=f"SARS-CoV-2 | {info['function']} | "
                        f"pos {info['start']}-{info['end']}"
        )
        records.append(record)

        # Save individual FASTA
        out = SEQ_DIR / f"{name}.fasta"
        SeqIO.write(record, out, "fasta")
        print(f"  OK    {name}.fasta  ({length} aa | "
              f"pos {info['start']}-{info['end']})")

    # Step 3 — save combined FASTA
    combined = SEQ_DIR / "all_nsps_combined.fasta"
    SeqIO.write(records, combined, "fasta")
    print(f"\n  Combined file: all_nsps_combined.fasta")

    # Step 4 — print AF3 submission guide
    print("\n" + "="*55)
    print("  AF3 SUBMISSION GUIDE")
    print("  Go to: https://alphafoldserver.com")
    print("  Submit 8 separate jobs — one per binary pair")
    print("="*55)
    jobs = [
        ("1", "NSP10-NSP16",    "NSP10.fasta + NSP16.fasta"),
        ("2", "NSP12-NSP7",     "NSP12.fasta + NSP7.fasta"),
        ("3", "NSP12-NSP8",     "NSP12.fasta + NSP8.fasta"),
        ("4", "NSP7-NSP8",      "NSP7.fasta  + NSP8.fasta"),
        ("5", "NSP9-NSP12",     "NSP9.fasta  + NSP12.fasta"),
        ("6", "NSP10-NSP14",    "NSP10.fasta + NSP14.fasta"),
        ("7", "NSP13-Helicase", "NSP13.fasta (monomer)"),
        ("8", "NSP12-NSP13",    "NSP12.fasta + NSP13.fasta"),
    ]
    print(f"\n  {'#':<4} {'Complex':<20} {'Sequences'}")
    print(f"  {'-'*4} {'-'*20} {'-'*30}")
    for num, job, seqs in jobs:
        print(f"  {num:<4} {job:<20} {seqs}")
    print(f"\n  After download save each AF3 zip to:")
    print(f"  01-alphafold3/<complex>/")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
