# PAN-CORONAVIRUS RTC INHIBITOR DISCOVERY
## Comprehensive Project Handoff & Continuation Guide

**Candidate:** Olivier Nsekuye  
**Institution:** University of Liège — GIGA-VIN Lab  
**Supervisor:** Prof. Jean-Claude Twizere  
**Funding:** FRIA-B1 Fellowship (2025–2029)  
**GitHub:** https://github.com/nsolly03/panCov-rtc-discovery  
**Last updated:** 2026-03-03  

---

## How to use this document

This document is the single source of truth for continuing this project. It is written for both human collaborators and AI assistants. Every step that was run is recorded in order. Every problem encountered is documented with its solution. If starting a new AI chat session, paste this document at the beginning of the conversation along with the current `WORKLOG.md` from the repository.

---

## 1. Quick Start — Resume the Project

```bash
# Activate environment
conda activate rtc-discovery
cd ~/projects/rtc-pan-coronavirus

# Check where we are
cat WORKLOG.md | tail -60

# Next pending script
python scripts/07_pocket_NSP10-NSP14.py   # does not exist yet — write it
```

---

## 2. Environment Setup

**System:** Ubuntu 24.04 (WSL2 on Windows 11) | Conda 25.11.1 | Git 2.43.0

**Conda environment:** `rtc-discovery` (Python 3.10)

| Package | Version | How installed | Purpose |
|---------|---------|---------------|---------|
| biopython | latest | conda | PDB parsing, sequence handling |
| numpy | latest | conda | Numerical operations |
| pandas | latest | conda | Data tables and CSV output |
| matplotlib | latest | conda | Plots and figures |
| seaborn | latest | conda | Heatmaps and statistical plots |
| requests | latest | conda | HTTP downloads |
| tqdm | latest | conda | Progress bars |
| jupyter | latest | conda | Interactive notebooks |
| ipykernel | latest | conda | Jupyter kernel |
| muscle | 5.3 | conda bioconda | Multiple sequence alignment |
| fpocket | 4.0 | conda bioconda | Pocket detection and druggability |
| py3Dmol | latest | pip | 3D structure visualization in Jupyter |
| nglview | latest | pip | Alternative 3D viewer in Jupyter |

**Recreate environment from scratch:**
```bash
conda env create -f environment.yml
conda activate rtc-discovery
conda install -c bioconda muscle fpocket -y
pip install py3Dmol nglview
python scripts/01_download_structures.py   # re-download PDB files
```

---

## 3. Project Structure

```
rtc-pan-coronavirus/
├── 00-reference/
│   ├── pdb_structures/          # 13 original PDB files + 14 extracted binary PDBs
│   ├── sequences/               # 8 NSP FASTA files from UniProt P0DTD1
│   │   └── conservation/        # Per-coronavirus FASTAs + MUSCLE alignments
│   └── known_interfaces/
│       └── chain_map.tsv        # Verified chain assignments for all 13 PDB files
├── 01-alphafold3/
│   └── NSP10-NSP14/             # AF3 files: best_model.cif/.pdb, confidence.json
│   └── <other complexes>/       # AF3 files present, not yet processed
├── 02-validation/
│   └── NSP10-NSP14/             # validation_result.json, interface_analysis.json,
│                                #   conservation_NSP10/NSP14.csv, conservation_summary.json
│   └── <other complexes>/       # Empty — pending
├── 03-virtual-screening/        # Empty — pending
├── 04-hits/                     # Empty — pending
├── 05-experimental/             # Empty — pending
├── scripts/                     # All Python scripts (01-06 complete)
├── notebooks/                   # Empty — visualization notebooks pending
├── results/                     # Empty — figures pending
├── docs/                        # This document + literature
├── WORKLOG.md                   # Detailed step-by-step lab notebook
├── README.md                    # Project overview
├── .gitignore                   # Excludes PDB files, AF3 MSA folders, large JSONs
└── environment.yml              # Conda environment definition
```

---

## 4. The 8 Binary Interface Targets

| # | Complex | Priority | PDB files | Chain A | Chain B | Pipeline status |
|---|---------|----------|-----------|---------|---------|-----------------|
| 1 | NSP10-NSP16 | Critical | 6W4H, 6WVN, 6WKQ | NSP10 (B) | NSP16 (A) | ⏳ Pending |
| 2 | NSP12-NSP7 | Critical | 7BV2, 6NUR, 7C2K | NSP12 (A) | NSP7 (C) | ⏳ Pending |
| 3 | NSP12-NSP8 | Critical | 7BV2, 6NUR, 7C2K | NSP12 (A) | NSP8 (B) | ⏳ Pending |
| 4 | NSP7-NSP8 | Low-Med | 7BV2, 6NUR | NSP7 (C) | NSP8 (B) | ⏳ Pending + fpocket gate |
| 5 | NSP9-NSP12 | Critical | 8SQK | NSP12 (A) | NSP9 (G) | ⏳ Pending |
| 6 | NSP10-NSP14 | High | 7DIY, 5C8T | NSP10 (A) | NSP14 (B) | ✅ Scripts 01-06 complete |
| 7 | NSP13-Helicase | High | 6XEZ, 7NIO | NSP13 (E) | — | ⏳ Pending |
| 8 | NSP12-NSP13 | Medium | 6XEZ, 7CXM, 7RDY | NSP12 (A) | NSP13 (E) | ⏳ Pending |

> ⚠️ **NSP7-NSP8** requires an additional fpocket druggability gate (score > 0.5) before proceeding to docking, due to its small interface area (~800 Å²).

---

## 5. Scripts — Complete Record

### RULE: Update WORKLOG.md after every script run

```bash
cat >> WORKLOG.md << 'EOF'

## Entry XXX — <what you did>
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/XX_scriptname.py
**Result:** <paste key output values>
**Status:** ✅ Done

---
EOF

git add .
git commit -m "<ComplexName>: <what was done>"
git push
```

---

### Script 01 — `download_structures.py` ✅

Downloads all 13 PDB files from RCSB PDB.

```bash
python scripts/01_download_structures.py
```

**Output:** `00-reference/pdb_structures/*.pdb` (13 files)  
**WORKLOG entry:** 005

---

### Script 02 — `extract_chains.py` ✅

Extracts binary interface chains from multi-chain PDB files using BioPython PDBIO and a ChainSelector class.

```bash
python scripts/02_extract_chains.py
```

**Output:** 14 binary PDB files e.g. `7BV2_NSP12-NSP7.pdb`  
**WORKLOG entry:** 006

---

### Script 03 — `download_sequences.py` ✅

Downloads all NSP sequences by slicing the SARS-CoV-2 pp1ab polyprotein (P0DTD1) from UniProt using verified amino acid position coordinates.

```bash
python scripts/03_download_sequences.py
```

**Output:** `00-reference/sequences/NSP7.fasta` through `NSP16.fasta`  
**Output:** `00-reference/sequences/all_nsps_combined.fasta`  
**WORKLOG entry:** 007

> ⚠️ Do NOT use NCBI Entrez to search for NSP sequences — it returns wrong proteins. Always use UniProt P0DTD1 with slice coordinates. See Known Issues section.

---

### Script 04 — `validate_NSP10-NSP14.py` ✅

Validates AF3 predicted interface against PDB crystal structure. Computes Precision, Recall, F1 at 5.0 Å cutoff.  
**Gate: F1 >= 0.70 required to proceed.**

```bash
python scripts/04_validate_NSP10-NSP14.py
```

**Input:** `7DIY.pdb` + `NSP10_NSP14_best_model.pdb`  
**Output:** `02-validation/NSP10-NSP14/validation_result.json`  
**Result:** NSP10 F1=0.952, NSP14 F1=0.953, Overall F1=0.952 ✅ PASS  
**WORKLOG entry:** 009

**AF3 confidence scores for NSP10-NSP14:**
- iptm = 0.89 ✅ (threshold > 0.6)
- ptm = 0.88 ✅ (threshold > 0.5)
- ranking_score = 0.91 ✅
- has_clash = 0.0 ✅
- fraction_disordered = 0.04 ✅

---

### Script 05 — `interface_NSP10-NSP14.py` ✅

Analyzes interface contacts across ALL available structures: 7DIY (SARS-CoV-2, 2.69Å), 5C8T (SARS-CoV-1, 3.20Å), and AF3 model. Identifies H-bonds, salt bridges, hydrophobic contacts. Ranks hotspot residues. Finds consensus across all 3 structures.

```bash
python scripts/05_interface_NSP10-NSP14.py
```

**Input:** `7DIY.pdb`, `5C8T.pdb`, `NSP10_NSP14_best_model.pdb`  
**Output:** `02-validation/NSP10-NSP14/interface_analysis.json`  
**WORKLOG entry:** 010

**Key findings:**

| Finding | Detail |
|---------|--------|
| PRIMARY DRUG TARGET | HIS80(NSP10) — ASP126(NSP14) salt bridge: 3.65Å (SARS-CoV-2), 2.59Å (SARS-CoV-1), 2.91Å (AF3) — present in ALL 3 structures |
| Secondary salt bridge | LYS93(NSP10) — GLU128(NSP14) at 3.11Å (SARS-CoV-1 only) |
| Total contacts (7DIY) | 174: 17 H-bonds, 41 hydrophobic, 1 salt bridge, 115 VDW |

**NSP10 consensus hotspots (17 residues, in all 3 structures):**  
`ALA4, ALA18, ASN3, ASN40, GLU6, HIS80, LEU45, LYS93, MET44, PHE16, PHE19, SER33, THR5, THR12, TYR96, VAL21, VAL42`

**NSP14 consensus hotspots (17 residues, in all 3 structures):**  
`ASP10, ASP126, ASP41, HIS26, ILE201, LEU7, LEU27, LYS9, PHE8, PRO20, PRO24, THR5, THR21, THR25, TYR69, VAL4, VAL66`

---

### Script 06 — `conservation_NSP10-NSP14.py` ✅

Downloads NSP10 and NSP14 polyprotein sequences from UniProt for 5 coronaviruses. Slices NSPs using verified coordinates. Runs MUSCLE alignment. Calculates Shannon entropy conservation score per residue. Maps onto consensus hotspots.

```bash
python scripts/06_conservation_NSP10-NSP14.py
```

**Input:** UniProt accessions P0DTD1, P0C6X7, K9N7C7, P0C6X1, P0C6X5  
**Output:** `conservation/NSP10_aligned.fasta`, `conservation/NSP14_aligned.fasta`  
**Output:** `02-validation/NSP10-NSP14/conservation_NSP10.csv`  
**Output:** `02-validation/NSP10-NSP14/conservation_NSP14.csv`  
**Output:** `02-validation/NSP10-NSP14/conservation_summary.json`  
**WORKLOG entry:** 013

**NSP10 conserved hotspots (score ≥ 0.8, across 5 coronaviruses):**

| Residue | Score | Notes |
|---------|-------|-------|
| THR5 | 1.000 | Identical in all 5 |
| PHE19 | 1.000 | Identical in all 5 |
| VAL21 | 1.000 | Identical in all 5 |
| ASN40 | 1.000 | Identical in all 5 |
| MET44 | 1.000 | Identical in all 5 |
| LEU45 | 1.000 | Identical in all 5 |
| **HIS80** | **1.000** | **PRIMARY TARGET — salt bridge anchor** |
| LYS93 | 1.000 | Identical in all 5 |
| VAL42 | 0.833 | Nearly conserved |

**NSP14 conserved hotspots (score ≥ 0.8):**

| Residue | Score | Notes |
|---------|-------|-------|
| LEU7 | 1.000 | Identical in all 5 |
| PHE8 | 1.000 | Identical in all 5 |
| LYS9 | 1.000 | Identical in all 5 |
| PRO20 | 1.000 | Identical in all 5 |
| THR25 | 1.000 | Identical in all 5 |
| THR127 | 1.000 | Identical in all 5 |
| ILE201 | 1.000 | Identical in all 5 |
| ASP10 | 0.833 | Nearly conserved |
| VAL4 | 0.833 | Nearly conserved |
| LEU27 | 0.833 | Nearly conserved |

> ✅ **Critical insight:** HIS80(NSP10) is fully conserved (1.000) across all 5 coronaviruses AND forms the primary salt bridge with ASP126(NSP14). This is the core scientific argument for pan-coronavirus activity in the FRIA proposal.

---

### Script 07_2 — `pocket_NSP10-NSP14_2.py` ✅

Runs fpocket on 7DIY, 5C8T, and AF3 model. Parses all pocket scores. Identifies which pocket overlaps with conserved hotspot residues. Scores overlap using primary residues (HIS80, LYS93, ASP126, THR127) weighted ×3. Computes docking box center and dimensions from hotspot Cα coordinates with 6.0 Å padding.

```bash
python scripts/07_pocket_NSP10-NSP14_2.py
```

**Input:** `7DIY.pdb`, `5C8T.pdb`, `NSP10_NSP14_best_model.pdb`  
**Output:** `02-validation/NSP10-NSP14/pocket_analysis_2.json`  
**WORKLOG entry:** 015

**Results:**

| Structure | Best pocket | Druggability | Spans interface | Contains HIS80 |
|-----------|-------------|-------------|-----------------|----------------|
| 7DIY | Pocket 3 | 0.000 | ✅ | ✅ |
| 5C8T | Pocket 90 | 0.001 | ✅ | ✅ |
| AF3 | Pocket 9 | 0.000 | ✅ | ❌ (slightly offset) |

> ⚠️ Low druggability scores (0.000) are **expected and normal** for PPI interfaces — fpocket is optimised for enzyme active sites. The docking box is defined from hotspot coordinates, not from druggability score.

> ⚠️ 5C8T returned 99 pockets because it contains 2 copies of the complex (chains A/B + C/D). This inflates the pocket count but does not affect the analysis.

**Selected docking box (7DIY, Pocket 3 + conserved hotspot Cα coords):**
- Center: (-4.776, 7.298, -25.886)
- Size: 31.685 × 34.288 × 48.982 Å
- Volume: 53,215 Å³
- Anchor residue: HIS80(NSP10)

---

### Script 08_2 — `docking_prep_NSP10-NSP14_2.py` ✅

Prepares the receptor file for virtual screening. Strips waters and heteroatoms from 7DIY. Keeps chains A+B only. Writes AutoDock Vina config and VirtualFlow JSON config. Verifies all hotspot residues are present in the receptor.

```bash
python scripts/08_docking_prep_NSP10-NSP14_2.py
```

**Input:** `7DIY.pdb`  
**Output:** `03-virtual-screening/NSP10-NSP14_2/receptor_NSP10-NSP14_2.pdb`  
**Output:** `03-virtual-screening/NSP10-NSP14_2/vina_config_NSP10-NSP14_2.txt`  
**Output:** `03-virtual-screening/NSP10-NSP14_2/virtualflow_config_NSP10-NSP14_2.json`  
**WORKLOG entry:** 016

**Receptor:** 417 residues, chains A+B, clean  
**Hotspot verification:** 9/9 NSP10 + 10/10 NSP14 conserved residues present ✅

---

## 6. Full Pipeline Per Complex

Repeat these steps for each of the 8 binary complexes. NSP10-NSP14 is the completed reference.

| Step | Script naming pattern | What it does | Gate / Key output |
|------|-----------------------|--------------|-------------------|
| 1 | `04_validate_<COMPLEX>.py` | AF3 vs PDB interface F1 at 5.0 Å cutoff | F1 ≥ 0.70 required |
| 2 | `05_interface_<COMPLEX>.py` | Contacts across ALL PDB structures + AF3. Consensus hotspots. | Ranked hotspot list + salt bridges |
| 3 | `06_conservation_<COMPLEX>.py` | 5 coronavirus sequences, MUSCLE, Shannon entropy per residue | Conserved hotspot residues |
| 4 | `07_pocket_<COMPLEX>.py` | fpocket on PDB + AF3. Druggability score. | Docking box coords. NSP7-NSP8: score > 0.5 required |
| 5 | `08_docking_prep_<COMPLEX>.py` | Prepare receptor, define docking box | VirtualFlow config file |
| 6 | `notebooks/<COMPLEX>.ipynb` | Contact map, conservation heatmap, 3D py3Dmol viewer | Publication figures |

**Naming convention:**
- Scripts: `04_validate_NSP10-NSP16.py`, `05_interface_NSP10-NSP16.py` etc.
- Validation: `02-validation/NSP10-NSP16/`
- AF3: `01-alphafold3/NSP10-NSP16/`

---

## 7. AlphaFold3 Guide

### Submit jobs
1. Go to https://alphafoldserver.com
2. Sign in with Google
3. New prediction → type job name e.g. `NSP10-NSP16`
4. Add Protein 1: paste sequence from `00-reference/sequences/NSP10.fasta`
5. Add Protein 2: paste sequence from `00-reference/sequences/NSP16.fasta`
6. Submit — takes 30–60 minutes per job
7. Download zip → extract into `01-alphafold3/<ComplexName>/`

> ⚠️ NSP13-Helicase is a monomer — submit NSP13.fasta only (single chain).

### Convert .cif to .pdb after download
```bash
python3 -c "
from Bio import PDB
parser = PDB.MMCIFParser(QUIET=True)
structure = parser.get_structure('NAME', '01-alphafold3/NSP10-NSP16/NSP10_NSP16_best_model.cif')
io = PDB.PDBIO()
io.set_structure(structure)
io.save('01-alphafold3/NSP10-NSP16/NSP10_NSP16_best_model.pdb')
print('Converted OK')
"
```

### Minimum quality thresholds

| Score | Minimum | Meaning |
|-------|---------|---------|
| iptm | > 0.6 | Interface prediction confidence |
| ptm | > 0.5 | Overall fold confidence |
| has_clash | = 0.0 | No structural clashes |
| fraction_disordered | < 0.2 | Less than 20% disordered |

---

## 8. Verified PDB Chain Assignments

> ⚠️ Always verify chain assignments by running `grep "^COMPND" structure.pdb` — never trust literature chain labels.

| PDB | Resolution | Chain A | Chain B | Chain C | Chain E | Chain G | Used for |
|-----|-----------|---------|---------|---------|---------|---------|----------|
| 6W4H | 1.80Å | NSP16 | NSP10 | — | — | — | NSP10-NSP16 primary |
| 6WVN | 1.95Å | NSP16 | NSP10 | — | — | — | NSP10-NSP16 + SAM cofactor |
| 6WKQ | 2.05Å | NSP16 | NSP10 | — | — | — | NSP10-NSP16 third structure |
| 7BV2 | 2.50Å | NSP12 | NSP8 | NSP7 | — | — | NSP12-NSP7/8, NSP7-NSP8 |
| 6NUR | 3.10Å | NSP12 | NSP8 | NSP7 | — | — | NSP12-NSP7/8, NSP7-NSP8 |
| 7C2K | 2.90Å | NSP12 | NSP8 | NSP7 | — | — | NSP12-NSP7/8 |
| 8SQK | 3.01Å | NSP12 | NSP8 | NSP7 | — | NSP9 | NSP9-NSP12 |
| 7DIY | 2.69Å | NSP10 | NSP14 | — | — | — | NSP10-NSP14 primary (SARS-CoV-2) |
| 5C8T | 3.20Å | NSP10 | NSP14 | — | — | — | NSP10-NSP14 (SARS-CoV-1) |
| 6XEZ | 3.50Å | NSP12 | NSP8 | NSP7 | NSP13 | — | NSP13-Helicase, NSP12-NSP13 |
| 7CXM | 3.20Å | NSP12 | NSP8 | NSP7 | NSP13 | — | NSP12-NSP13 |
| 7RDY | 3.10Å | NSP12 | NSP8 | NSP7 | NSP13 | — | NSP12-NSP13 |
| 7NIO | 2.80Å | NSP13 | — | — | — | — | NSP13-Helicase (dimer) |

**PDB files removed (wrong structures):**
- `7DFG` — contained RdRp complex, not NSP10-NSP16
- `9FW2` — contained NSP11+NSP14, not NSP9-NSP12

---

## 9. UniProt Accessions for Conservation Analysis

> ⚠️ Use UniProt REST API only. Do NOT use NCBI Entrez — it returns wrong proteins.

| Coronavirus | Accession | Total aa | NSP10 positions | NSP14 positions |
|-------------|-----------|----------|-----------------|-----------------|
| SARS-CoV-2 | P0DTD1 | 7096 | 4254–4392 | 5926–6452 |
| SARS-CoV-1 | P0C6X7 | 7073 | 4231–4369 | 5903–6429 |
| MERS-CoV | K9N7C7 | 7078 | 4238–4377 | 5909–6432 |
| HCoV-229E | P0C6X1 | 6758 | 3934–4068 | 5593–6110 (annotated as "Exoribonuclease") |
| HCoV-NL63 | P0C6X5 | 6729 | 3909–4043 | 5568–6085 (annotated as "Exoribonuclease") |

**Excluded coronaviruses:**
- HCoV-OC43, HCoV-HKU1: Only pp1a available — NSP14 is in the pp1b region, unavailable
- BatCoV-RaTG13: No reviewed UniProt entry available

> For other complexes (NSP10-NSP16, NSP12-NSP7 etc.) fetch NSP coordinates fresh from UniProt feature annotations using `type == 'Chain'` and filtering by description.

---

## 10. Visualization Plan

All visualizations go in `notebooks/<ComplexName>.ipynb`. Use py3Dmol for 3D, matplotlib/seaborn for 2D.

| Visualization | Tool | Status |
|---------------|------|--------|
| Interface contact map (residue vs residue) | seaborn heatmap | ⏳ Pending |
| Hotspot bar chart (score per residue) | matplotlib barh | ⏳ Pending |
| Conservation heatmap (residue × coronavirus) | seaborn | ⏳ Pending |
| 3D interface viewer — hotspots colored by conservation | py3Dmol | ⏳ Pending |
| 3D salt bridge viewer — HIS80-ASP126 highlighted | py3Dmol | ⏳ Pending |
| Pocket druggability plot | matplotlib | ⏳ Pending (in Jupyter notebook) |
| Summary figure across all 8 complexes | matplotlib | ⏳ Pending (final) |

**Launch Jupyter:**
```bash
conda activate rtc-discovery
cd ~/projects/rtc-pan-coronavirus
jupyter lab
```

---

## 11. Known Issues and Fixes

### Issue 1 — NCBI Entrez returns wrong proteins
- **Symptom:** Sequences of 82 aa, 306 aa, or 1356 aa instead of expected NSP lengths
- **Cause:** Generic Entrez search queries match wrong organisms (e.g. Porcine reproductive syndrome virus)
- **Fix:** Use UniProt REST API with verified accession numbers and slice coordinates. Never use Entrez esearch for NSP sequences

### Issue 2 — NSP14 annotation varies by coronavirus
- **Symptom:** No NSP14 found for HCoV-229E and HCoV-NL63 in UniProt feature search
- **Cause:** Annotated as "Exoribonuclease" not "Non-structural protein 14"
- **Fix:** Use coordinate-based slicing. HCoV-229E: 5593–6110. HCoV-NL63: 5568–6085

### Issue 3 — Large AF3 MSA files exceed GitHub 50MB limit
- **Symptom:** GitHub warning about files > 50MB on git push
- **Fix:** Already in `.gitignore`: `01-alphafold3/**/msas/`, `01-alphafold3/**/*full_data*.json`, `01-alphafold3/**/*.a3m`

### Issue 4 — BioPython float32 not JSON serializable
- **Symptom:** `TypeError: Object of type float32 is not JSON serializable`
- **Fix:** Cast with `float()` before `json.dump()`, or use a custom `NpEncoder(json.JSONEncoder)` class

### Issue 5 — Salt bridge duplicates in interface analysis
- **Symptom:** Same salt bridge printed multiple times
- **Fix:** Use `seen_pairs = set()` keyed by `(residue_number_A, residue_number_B)` to deduplicate before appending contacts

### Issue 6 — sed command corrupts script files
- **Symptom:** `ValueError: I/O operation on closed file` after using `sed -i` to patch a script
- **Fix:** Always rewrite entire script with `cat > script.py << 'EOF'` instead of patching with sed

### Issue 7 — PDB files not in GitHub
- **Note:** PDB files are in `.gitignore` by design — repository too large otherwise
- **Fix:** Run `python scripts/01_download_structures.py` after cloning to get all PDB files

### Issue 8 — Wrong structures downloaded initially
- **7DFG:** Contained RdRp complex (NSP12/NSP7/NSP8), not NSP10-NSP16 — removed
- **9FW2:** Contained NSP11+NSP14, not NSP9-NSP12 — removed
- **7LYJ:** RNA only, no proteins — removed
- **Fix:** 6WKQ added as correct third NSP10-NSP16 structure

---

## 12. Git Workflow

```bash
# Standard commit after every step
git add .
git commit -m "<ComplexName>: <what was done>"
git push

# Clone fresh
git clone https://github.com/nsolly03/panCov-rtc-discovery.git
cd panCov-rtc-discovery
```

**Credentials:** Stored via `git config --global credential.helper store` — no password needed.

**Commit message conventions:**
```
NSP10-NSP14: AF3 validation complete, F1=0.952
NSP10-NSP14: interface analysis, HIS80-ASP126 primary hotspot
NSP10-NSP14: conservation complete, 9 NSP10 residues conserved
Add script 07: fpocket pocket detection NSP10-NSP14
```

---

## 13. NSP10-NSP14 Pipeline — Completed Milestones

| Step | Script | Result | WORKLOG |
|------|--------|--------|---------|
| AF3 validation | 04_validate_NSP10-NSP14.py | F1=0.952 ✅ | 009 |
| Interface analysis | 05_interface_NSP10-NSP14.py | HIS80-ASP126 salt bridge identified ✅ | 010 |
| Conservation | 06_conservation_NSP10-NSP14.py | 9 NSP10 + 10 NSP14 residues conserved ✅ | 013 |
| Pocket detection | 07_pocket_NSP10-NSP14_2.py | 7DIY Pocket 3 selected, box defined ✅ | 015 |
| Docking prep | 08_docking_prep_NSP10-NSP14_2.py | Receptor 417 aa, all hotspots verified ✅ | 016 |
| Visualization | notebooks/NSP10-NSP14_2.ipynb | ⏳ Next step |  |

---

## 14. Immediate Next Steps (in order)

| Priority | Action | Notes |
|----------|--------|-------|
| 1 — NOW | Build `notebooks/NSP10-NSP14_2.ipynb` | Contact map, conservation heatmap, py3Dmol 3D with HIS80-ASP126, pocket visualization |
| 2 | Repeat scripts 04_2–08_2 for NSP10-NSP16 | 6W4H primary, 6WVN + 6WKQ secondary. All new scripts get _2 suffix |
| 3 | Repeat for NSP12-NSP7 | 7BV2 primary PDB |
| 4 | Repeat for NSP12-NSP8 | 7BV2 primary PDB |
| 5 | Repeat for NSP9-NSP12 | 8SQK chains A+G |
| 6 | Repeat for NSP13-Helicase | 6XEZ chain E |
| 7 | Repeat for NSP12-NSP13 | 6XEZ primary |
| 8 | Repeat for NSP7-NSP8 | 7BV2 — additional fpocket gate score > 0.5 required |
| 9 | Virtual screening setup | VirtualFlow after all 8 complexes have docking boxes |

> ⚠️ **Naming convention from this point:** All new scripts, output files, and notebooks use `_2` suffix to avoid overwriting any existing files.

---

*Document last updated: March 3, 2026 | Olivier Nsekuye | University of Liège GIGA-VIN Lab*
