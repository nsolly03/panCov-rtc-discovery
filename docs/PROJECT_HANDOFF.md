# Project Handoff Document
## RTC Pan-Coronavirus Inhibitor Discovery
### NSP10-NSP14 Interface — Full Pipeline

**Candidate:** Olivier Nsekuye  
**Institution:** University of Liège — GIGA-VIN Lab  
**Supervisor:** Prof. Jean-Claude Twizere  
**Funding:** FRIA-B1 Fellowship  
**GitHub:** https://github.com/nsolly03/panCov-rtc-discovery  
**Last updated:** 2026-03-03  

---

## How to use this document
This document captures everything done in this project so far.
If you are continuing in a new chat or a new collaborator is
joining, start here. Follow the steps in order. Every command
that was run is recorded. Every problem encountered is documented
with its fix.

---

## Environment setup

**Operating system:** Ubuntu 24.04 (WSL2 on Windows)  
**Conda environment:** rtc-discovery (Python 3.10)  

**Activate environment:**
```
conda activate rtc-discovery
cd ~/projects/rtc-pan-coronavirus
```

**Key packages installed:**
```
biopython    — PDB parsing, sequence handling
numpy        — numerical operations
pandas       — data tables
matplotlib   — plots
seaborn      — heatmaps
requests     — HTTP downloads
tqdm         — progress bars
jupyter      — notebooks
muscle 5.3   — multiple sequence alignment (conda bioconda)
py3Dmol      — 3D visualization in Jupyter (pip)
nglview      — 3D viewer (pip)
fpocket 4.0  — pocket detection (conda bioconda)
```

**Recreate environment from scratch:**
```
conda env create -f environment.yml
conda activate rtc-discovery
conda install -c bioconda muscle fpocket -y
pip install py3Dmol nglview
```

---

## Project structure
```
rtc-pan-coronavirus/
├── 00-reference/
│   ├── pdb_structures/      — 13 PDB files + 14 extracted binaries
│   ├── sequences/           — NSP FASTA files from UniProt
│   │   └── conservation/    — Per-coronavirus sequences + alignments
│   └── known_interfaces/
│       └── chain_map.tsv    — Verified chain assignments for all PDBs
├── 01-alphafold3/
│   └── NSP10-NSP14/         — AF3 prediction files
│       ├── NSP10_NSP14_best_model.cif
│       ├── NSP10_NSP14_best_model.pdb  (converted)
│       └── NSP10_NSP14_confidence.json
├── 02-validation/
│   └── NSP10-NSP14/
│       ├── validation_result.json
│       ├── interface_analysis.json
│       ├── conservation_NSP10.csv
│       ├── conservation_NSP14.csv
│       └── conservation_summary.json
├── 03-virtual-screening/    — (pending)
├── 04-hits/                 — (pending)
├── 05-experimental/         — (pending)
├── scripts/                 — all Python scripts
├── notebooks/               — Jupyter notebooks (pending)
├── results/                 — figures (pending)
├── docs/                    — this document + literature
├── WORKLOG.md               — detailed step-by-step log
├── README.md                — project overview
└── environment.yml          — conda environment
```

---

## 8 Binary interface targets

| # | Complex | Priority | PDB files | Status |
|---|---------|----------|-----------|--------|
| 1 | NSP10-NSP16 | Critical | 6W4H, 6WVN, 6WKQ | ⏳ pending |
| 2 | NSP12-NSP7 | Critical | 7BV2, 6NUR, 7C2K | ⏳ pending |
| 3 | NSP12-NSP8 | Critical | 7BV2, 6NUR, 7C2K | ⏳ pending |
| 4 | NSP7-NSP8 | Low-Med | 7BV2, 6NUR | ⏳ pending |
| 5 | NSP9-NSP12 | Critical | 8SQK | ⏳ pending |
| 6 | NSP10-NSP14 | High | 7DIY, 5C8T | ✅ active |
| 7 | NSP13-Helicase | High | 6XEZ, 7NIO | ⏳ pending |
| 8 | NSP12-NSP13 | Medium | 6XEZ, 7CXM, 7RDY | ⏳ pending |

---

## Scripts — what each one does

### Script 01 — download_structures.py
Downloads all 13 PDB files from RCSB.
```
python scripts/01_download_structures.py
```
Output: 00-reference/pdb_structures/*.pdb

### Script 02 — extract_chains.py
Extracts binary interface chains from multi-chain PDB files.
```
python scripts/02_extract_chains.py
```
Output: 14 binary PDB files e.g. 7BV2_NSP12-NSP7.pdb

### Script 03 — download_sequences.py
Downloads all NSP sequences from UniProt P0DTD1 polyprotein.
```
python scripts/03_download_sequences.py
```
Output: 00-reference/sequences/NSP7.fasta ... NSP16.fasta

### Script 04 — validate_NSP10-NSP14.py
Validates AF3 predicted interface against PDB crystal structure.
Computes Precision, Recall, F1. Gate: F1 >= 0.70 to proceed.
```
python scripts/04_validate_NSP10-NSP14.py
```
Output: 02-validation/NSP10-NSP14/validation_result.json
Result: F1=0.952 ✅ PASS

### Script 05 — interface_NSP10-NSP14.py
Analyzes interface contacts across 3 structures:
7DIY (SARS-CoV-2), 5C8T (SARS-CoV-1), AF3 model.
Identifies hydrogen bonds, salt bridges, hydrophobic contacts.
Ranks hotspot residues. Finds consensus across all structures.
```
python scripts/05_interface_NSP10-NSP14.py
```
Output: 02-validation/NSP10-NSP14/interface_analysis.json

**Key finding:**
HIS80(NSP10) — ASP126(NSP14) salt bridge at 3.65A
Present in SARS-CoV-2, SARS-CoV-1 AND AF3 model.
Primary drug target for docking box.

**Consensus hotspots NSP10 (17 residues):**
ALA4, ALA18, ASN3, ASN40, GLU6, HIS80, LEU45, LYS93,
MET44, PHE16, PHE19, SER33, THR5, THR12, TYR96, VAL21, VAL42

**Consensus hotspots NSP14 (17 residues):**
ASP10, ASP126, ASP41, HIS26, ILE201, LEU7, LEU27, LYS9,
PHE8, PRO20, PRO24, THR5, THR21, THR25, TYR69, VAL4, VAL66

### Script 06 — conservation_NSP10-NSP14.py
Downloads NSP10 and NSP14 from 5 coronaviruses via UniProt.
Runs MUSCLE alignment. Calculates Shannon entropy conservation.
Maps scores onto consensus hotspot residues.
```
python scripts/06_conservation_NSP10-NSP14.py
```
Output:
  02-validation/NSP10-NSP14/conservation_NSP10.csv
  02-validation/NSP10-NSP14/conservation_NSP14.csv
  02-validation/NSP10-NSP14/conservation_summary.json
Status: ⏳ running

### Script 07 — pocket_NSP10-NSP14.py (PENDING)
Runs fpocket on 7DIY and AF3 model.
Ranks pockets by druggability score.
Defines docking box centered on conserved hotspots.
```
python scripts/07_pocket_NSP10-NSP14.py
```

---

## AlphaFold3 predictions

**How to get AF3 predictions:**
1. Go to https://alphafoldserver.com
2. Sign in with Google
3. New prediction → add protein sequences
4. Submit one job per binary pair
5. Download zip → extract into 01-alphafold3/<complex>/
6. Convert .cif to .pdb:
```
python3 -c "
from Bio import PDB
parser = PDB.MMCIFParser(QUIET=True)
structure = parser.get_structure('name', 'model.cif')
io = PDB.PDBIO()
io.set_structure(structure)
io.save('model.pdb')
"
```

**Minimum quality thresholds:**
- iptm > 0.6 (interface confidence)
- ptm  > 0.5 (overall fold)
- has_clash = 0.0 (no clashes)

**NSP10-NSP14 AF3 scores:**
- iptm = 0.89 ✅ excellent
- ptm  = 0.88 ✅ excellent
- ranking_score = 0.91 ✅ excellent

---

## Verified PDB chain assignments

| PDB | Chain A | Chain B | Chain C | Chain E/F | Chain G |
|-----|---------|---------|---------|-----------|---------|
| 6W4H | NSP16 | NSP10 | — | — | — |
| 6WVN | NSP16 | NSP10 | — | — | — |
| 6WKQ | NSP16 | NSP10 | — | — | — |
| 7BV2 | NSP12 | NSP8 | NSP7 | — | — |
| 6NUR | NSP12 | NSP8 | NSP7 | — | — |
| 7C2K | NSP12 | NSP8 | NSP7 | — | — |
| 8SQK | NSP12 | NSP8 | NSP7 | — | NSP9 |
| 7DIY | NSP10 | NSP14 | — | — | — |
| 5C8T | NSP10 | NSP14 | — | — | — |
| 6XEZ | NSP12 | NSP8 | NSP7 | NSP13 | — |
| 7CXM | NSP12 | NSP8 | NSP7 | NSP13 | — |
| 7RDY | NSP12 | NSP8 | NSP7 | NSP13 | — |
| 7NIO | NSP13 | — | — | — | — |

---

## Known issues and fixes

### Issue 1 — NCBI Entrez returns wrong proteins
**Symptom:** Sequences of 82 aa, 306 aa, 1356 aa for NSP10
**Cause:** Generic search query returned wrong organisms
**Fix:** Use verified UniProt accessions with slice coordinates

### Issue 2 — NSP14 annotation varies by coronavirus
**Symptom:** No NSP14 found for HCoV-229E, HCoV-NL63
**Cause:** These annotate it as "Exoribonuclease" not "NSP14"
**Fix:** Use coordinate-based slicing from feature table

### Issue 3 — BatCoV-RaTG13 unavailable
**Symptom:** HTTP 400 for QHR63298/QHR63299
**Cause:** No reviewed UniProt entry available
**Fix:** Excluded from conservation analysis (5 coronaviruses used)

### Issue 4 — HCoV-OC43 and HKU1 pp1ab missing
**Symptom:** Only pp1a (4383 aa) available — NSP14 in pp1b region
**Fix:** Excluded from conservation analysis

### Issue 5 — JSON serialization fails with BioPython float32
**Symptom:** TypeError: Object of type float32 is not JSON serializable
**Fix:** Use float() to cast before json.dump()

### Issue 6 — Large AF3 MSA files exceed GitHub 50MB limit
**Symptom:** GitHub warning on git push
**Fix:** Added to .gitignore:
  01-alphafold3/**/msas/
  01-alphafold3/**/*full_data*.json
  01-alphafold3/**/*.a3m

---

## Next steps (in order)

1. ⏳ Finish Script 06 — conservation analysis
2. ⏳ Script 07 — fpocket pocket detection on 7DIY + AF3
3. ⏳ Script 08 — define docking box from conserved hotspots
4. ⏳ Notebook 01 — visualization (contact map, conservation
                     heatmap, 3D hotspot viewer with py3Dmol)
5. ⏳ Repeat Scripts 04-08 for all 7 remaining complexes
6. ⏳ Virtual screening setup with VirtualFlow

---

## Pipeline per complex (repeat for each of 8 targets)
```
Script 04 — AF3 validation vs PDB     Gate: F1 >= 0.70
Script 05 — Interface analysis        3 structures: PDB1+PDB2+AF3
Script 06 — Conservation analysis     5 coronaviruses
Script 07 — Pocket detection          fpocket on PDB + AF3
Script 08 — Docking box definition    Conserved hotspots only
Notebook  — Visualization             Contact map + heatmap + 3D
```

---

## GitHub repository
https://github.com/nsolly03/panCov-rtc-discovery

**Clone:**
```
git clone https://github.com/nsolly03/panCov-rtc-discovery.git
cd panCov-rtc-discovery
conda env create -f environment.yml
conda activate rtc-discovery
```

