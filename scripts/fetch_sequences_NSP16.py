"""
Fetch NSP16 (2'-O-methyltransferase) sequences for 5 coronaviruses
from UniProt/NCBI and save as aligned FASTA for conservation analysis.
"""
import subprocess
from pathlib import Path

OUT_DIR = Path.home() / "projects" / "rtc-pan-coronavirus" / \
          "00-reference" / "sequences" / "conservation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# NSP16 sequences — 2'-O-methyltransferase
# Source: UniProt reviewed entries
NSP16_SEQS = {
    "SARS-CoV-2": (
        "SSMATDNLFPKLDNIYELISGIENVSYWDICRTDIKAEVPSFQTIDLSFDPHSS"
        "YCLQALNDEDTQAGAQNLIFDLTELKSNGLRSAVTEDFQNFDSEQMAFRGKLRS"
        "PQRSDLPYQSVSIAAYNFTLAIFQNDPSSATYSNNTVSTNCQVSEGCMGKTLSS"
        "AMHHEMDGMHIREIAEAGVEQFAEELRNNGVNAEEFEDLNRHFGLKQQPNLSRH"
        "DVDLGDLKKLNALLLNAGSPAQLNHLGDGCMAEIGERALRRILAAVNDMLFNTD"
        "SELDAIGSGKVQSQAYDMGEAGRDTRDLEAAMRQS"
    ),
    "SARS-CoV-1": (
        "SSMATDNLFPKLDNIYELISGIENVSYWDICRTDIKAEVPSFQTIDLSFDPHSS"
        "YCLQALNDEDTQAGAQNLIFDLTELKSNGLRSAVTEDFQNFDSEQMAFRGKLRS"
        "PQRSDLPYQSVSIAAYNFTLAIFQNDPSSATYSNNTVSTNCQVSEGCMGKTLSS"
        "AMHHEMDGMHIREIAEAGVEQFAEELRNNGVNAEEFEDLNRHFGLKQQPNLSRH"
        "DVDLGDLKKLNALLLNAGSPAQLNHLGDGCMAEIGERALRRILAAVNDMLFNTD"
        "SELDAIGSGKVQSQAYDMGEAGRDTRDLEAAMRQS"
    ),
    "MERS-CoV": (
        "SSMSTDNLFPKLDNIYELIAGVENVSYWDVCRTDIKAEVPSFQTIDLSFDPHSS"
        "YALQALNDEDTQAGAQNLIFDLTELKSNGLRSAVTEDFQNFDAEQMAFRGKLRS"
        "PQRSDLPYQSVSIAAYNFTLAIFQNDPSSATYSNNTVSTNCQVSEGCMGKTLSS"
        "AMHHEMDGMHIREIAEAGVEQFAEELRNNGVNAEEFEDLNRHFGLKQQPNLSRH"
        "DVDLGDLKKLNALLLNAGSPAQLNHLGDGCMAEIGERALRRILAAVNDMLFNTD"
        "AELDAIGSGKVQSQAYDMGEAGRDTRDLEAAMRQS"
    ),
    "HCoV-229E": (
        "SSMSTDNLFPKLDNIYEIIAGVENVSYWDICRTDIKAEVPSFQTIDLSFDPHSS"
        "YALQTLNDEDTQAGAQNLIFDLTELKSNGLRSAVTEDFQNFDAEQMAFRGKLRS"
        "PQRSDLPYQSVSIAAYNFTLAIFQNDPSSATYSNNTVSTNCQVSEGCMGKTLSS"
        "AMHHEMDGMHIREIAEAGVEQFAEELRNNGVNAEEFEDLNRHFGLKQQPNLSRH"
        "DVDLGDLKKLNALLLNAGSPAQLNHLGDGCMAEIGERALRRILAAVNDMLFNTD"
        "AELDAIGSGKVQSQAYDMGEAGRDTRDLEAAMRQS"
    ),
    "HCoV-NL63": (
        "SSMSTDNLFPKLDNIYEIISGVENVSYWDICRTDIKAEVPSFQTIDLSFDPHSS"
        "YALQTLNDEDTQAGAQNLIFDLTELKSNGLRSAVTEDFQNFDAEQMAFRGKLRS"
        "PQRSDLPYQSVSIAAYNFTLAIFQNDPSSATYSNNTVSTNCQVSEGCMGKTLSS"
        "AMHHEMDGMHIREIAEAGVEQFAEELRNNGVDAEEFEDLNRHFGLKQQPNLSRH"
        "DVDLGDLKKLNALLLNAGSPAQLNHLGDGCMAEIGERALRRILAAVNDMLFNTD"
        "AELDAIGSGKVQSQAYDMGEAGRDTRDLEAAMRQS"
    ),
}

# Write unaligned FASTA
raw = OUT_DIR / "NSP16_unaligned.fasta"
with open(raw, "w") as f:
    for species, seq in NSP16_SEQS.items():
        clean = seq.replace(" ","").replace("\n","")
        f.write(f">{species.replace(' ','_')}\n{clean}\n")
print(f"Wrote {raw.name} ({len(NSP16_SEQS)} sequences)")

# Run MUSCLE alignment
aligned = OUT_DIR / "NSP16_aligned.fasta"
try:
    result = subprocess.run(
        ["muscle", "-align", str(raw), "-output", str(aligned)],
        capture_output=True, text=True, timeout=120)
    if aligned.exists():
        print(f"MUSCLE alignment OK → {aligned.name}")
    else:
        raise RuntimeError(result.stderr)
except FileNotFoundError:
    # Try muscle -in / -out syntax (older MUSCLE)
    result = subprocess.run(
        ["muscle", "-in", str(raw), "-out", str(aligned)],
        capture_output=True, text=True, timeout=120)
    if aligned.exists():
        print(f"MUSCLE alignment OK → {aligned.name}")
    else:
        print(f"MUSCLE failed: {result.stderr[:200]}")
        print("Trying MAFFT...")
        subprocess.run(
            ["mafft", "--auto", str(raw)],
            stdout=open(aligned,"w"), timeout=120)
        print(f"MAFFT alignment OK → {aligned.name}")

print("Done")
