"""
Fetch NSP12 and NSP7 sequences for 5 coronaviruses
for conservation analysis.
"""
import time
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = "olivier.nsekuye@uliege.be"

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"
SEQ_DIR = PROJECT / "00-reference" / "sequences" / "conservation"
SEQ_DIR.mkdir(parents=True, exist_ok=True)

# UniProt IDs for NSP12 (RdRp) across 5 coronaviruses
# NSP12 = pp1ab residues ~4393-5324 (SARS-CoV-2)
NSP12_IDS = {
    "SARS-CoV-2": ("P0DTD1", 4393, 5324),
    "SARS-CoV-1": ("P0C6X7", 4374, 5305),
    "MERS-CoV":   ("K9N7C7", 4386, 5317),
    "HCoV-229E":  ("P0C6X1", 4071, 4989),
    "HCoV-NL63":  ("P0C6X5", 4093, 5014),
}

# NSP7 = pp1ab (SARS-CoV-2 residues ~3859-3942)
NSP7_IDS = {
    "SARS-CoV-2": ("P0DTD1", 3859, 3942),
    "SARS-CoV-1": ("P0C6X7", 3840, 3923),
    "MERS-CoV":   ("K9N7C7", 3852, 3935),
    "HCoV-229E":  ("P0C6X1", 3574, 3651),
    "HCoV-NL63":  ("P0C6X5", 3596, 3673),
}


def fetch_uniprot_segment(uniprot_id, start, end, label):
    """Fetch segment from UniProt via Entrez."""
    try:
        handle = Entrez.efetch(
            db="protein", id=uniprot_id,
            rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        seq = str(record.seq)[start-1:end]
        print(f"  {label}: {len(seq)} aa "
              f"(UniProt {uniprot_id} [{start}-{end}])")
        return seq
    except Exception as e:
        print(f"  {label}: ERROR — {e}")
        return None


def write_fasta(sequences, out_file):
    with open(out_file, "w") as f:
        for label, seq in sequences.items():
            if seq:
                f.write(f">{label}\n{seq}\n")
    print(f"  Written: {out_file.name}")


print("Fetching NSP12 sequences...")
nsp12_seqs = {}
for cov, (uid, s, e) in NSP12_IDS.items():
    nsp12_seqs[cov] = fetch_uniprot_segment(uid, s, e, cov)
    time.sleep(0.5)

write_fasta(nsp12_seqs,
            SEQ_DIR / "NSP12_unaligned.fasta")

print("\nFetching NSP7 sequences...")
nsp7_seqs = {}
for cov, (uid, s, e) in NSP7_IDS.items():
    nsp7_seqs[cov] = fetch_uniprot_segment(uid, s, e, cov)
    time.sleep(0.5)

write_fasta(nsp7_seqs,
            SEQ_DIR / "NSP7_unaligned.fasta")

print("\nDone.")
