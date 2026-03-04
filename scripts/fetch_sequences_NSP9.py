"""Fetch NSP9 sequences for 5 coronaviruses."""
import time
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = "olivier.nsekuye@uliege.be"

SEQ_DIR = (Path.home() / "projects" / "rtc-pan-coronavirus"
           / "00-reference" / "sequences" / "conservation")
SEQ_DIR.mkdir(parents=True, exist_ok=True)

# NSP9 = pp1ab segments
NSP9_IDS = {
    "SARS-CoV-2": ("P0DTD1", 4141, 4253),
    "SARS-CoV-1": ("P0C6X7", 4122, 4234),
    "MERS-CoV":   ("K9N7C7", 4134, 4246),
    "HCoV-229E":  ("P0C6X1", 3845, 3957),
    "HCoV-NL63":  ("P0C6X5", 3867, 3979),
}

def fetch_segment(uid, start, end, label):
    try:
        handle = Entrez.efetch(
            db="protein", id=uid,
            rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        seq = str(record.seq)[start-1:end]
        print(f"  {label}: {len(seq)} aa "
              f"(UniProt {uid} [{start}-{end}])")
        return seq
    except Exception as e:
        print(f"  {label}: ERROR — {e}")
        return None

print("Fetching NSP9 sequences...")
seqs = {}
for cov, (uid, s, e) in NSP9_IDS.items():
    seqs[cov] = fetch_segment(uid, s, e, cov)
    time.sleep(0.5)

out = SEQ_DIR / "NSP9_unaligned.fasta"
with open(out, "w") as f:
    for label, seq in seqs.items():
        if seq:
            f.write(f">{label}\n{seq}\n")
print(f"Written: {out.name}")
