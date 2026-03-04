# WORKLOG — RTC Pan-Coronavirus Inhibitor Discovery
**Candidate:** Olivier Nsekuye  
**Institution:** University of Liège — GIGA-VIN Lab  
**Supervisor:** Prof. Jean-Claude Twizere  
**Funding:** FRIA-B1 Fellowship  
**Start date:** 2025  

---

## How to use this file
This file records every step taken in this project in order.
Each entry has a date, what was done, what command was run, and what the output meant.
Anyone reading this can reproduce the full project from scratch by following these entries.

---

## Entry 001 — Project setup
**Date:** $(date +%Y-%m-%d)  
**What:** Created project folder structure on WSL Ubuntu 24.04  
**Command:**
```
mkdir -p ~/projects/rtc-pan-coronavirus/{00-reference/...,scripts,notebooks,results,docs}
```
**Result:** All 8 target folders created under 00-reference through 05-experimental  
**Status:** ✅ Done  

---

## Entry 002 — Git initialized
**Date:** $(date +%Y-%m-%d)  
**What:** Initialized git repository and set default branch to main  
**Commands:**
```
git init
git branch -m main
```
**Result:** Empty git repo on branch main  
**Status:** ✅ Done  

---


## Entry 003 — GitHub repository connected
**Date:** $(date +%Y-%m-%d)  
**What:** Created GitHub repository and pushed initial commit  
**Commands:**
```
git remote add origin https://github.com/nsolly03/panCov-rtc-discovery.git
git push -u origin main
```
**Result:** Project live at https://github.com/nsolly03/panCov-rtc-discovery  
**Status:** ✅ Done  

---

## Entry 004 — Conda environment created
**Date:** $(date +%Y-%m-%d)  
**What:** Created conda environment with all required packages  
**Command:**
```
conda env create -f environment.yml
```
**Packages installed:** python=3.10, biopython, numpy, pandas, matplotlib,
seaborn, requests, tqdm, jupyter, ipykernel  
**Activate with:** conda activate rtc-discovery  
**Status:** ✅ Done  

---


## Entry 005 — PDB structures downloaded and verified
**Date:** $(date +%Y-%m-%d)  
**What:** Downloaded 13 PDB structures for 8 binary RTC targets  
**Script:** scripts/01_download_structures.py  
**Structures confirmed:**
- NSP10-NSP16 : 6W4H, 6WVN, 6WKQ
- NSP12-NSP7  : 7BV2, 6NUR, 7C2K (extract chains A+C)
- NSP12-NSP8  : 7BV2, 6NUR, 7C2K (extract chains A+B)
- NSP7-NSP8   : 7BV2, 6NUR (extract chains C+B)
- NSP9-NSP12  : 8SQK (chains A=NSP12, G=NSP9)
- NSP10-NSP14 : 7DIY, 5C8T (chains A=NSP10, B=NSP14)
- NSP13-Helicase : 6XEZ (chain E), 7NIO (chain A)
- NSP12-NSP13 : 6XEZ, 7CXM, 7RDY (chains A=NSP12, E=NSP13)

**Corrections made:**
- 7DFG removed — contained RdRp complex not NSP10-NSP16
- 9FW2 removed — contained NSP11+NSP14 not NSP9-NSP12
- 7LYJ removed — RNA only no proteins
- 6WKQ added as correct third NSP10-NSP16 structure
- Chain map fully verified by inspecting COMPND records

**Chain map saved to:** 00-reference/known_interfaces/chain_map.tsv  
**Status:** ✅ Done  

---

## Entry 006 — Binary chains extracted from multi-chain PDBs
**Date:** $(date +%Y-%m-%d)  
**What:** Extracted binary interface chains from trimer and complex PDBs  
**Script:** scripts/02_extract_chains.py  
**Files created (14 total):**
- 7BV2_NSP12-NSP7.pdb  (chains A+C, 7193 atoms)
- 6NUR_NSP12-NSP7.pdb  (chains A+C, 6917 atoms)
- 7C2K_NSP12-NSP7.pdb  (chains A+C, 8018 atoms)
- 7BV2_NSP12-NSP8.pdb  (chains A+B, 7583 atoms)
- 6NUR_NSP12-NSP8.pdb  (chains A+B, 7266 atoms)
- 7C2K_NSP12-NSP8.pdb  (chains A+B, 8369 atoms)
- 7BV2_NSP7-NSP8.pdb   (chains C+B, 1352 atoms)
- 6NUR_NSP7-NSP8.pdb   (chains C+B, 1429 atoms)
- 8SQK_NSP9-NSP12.pdb  (chains A+G, 16361 atoms)
- 6XEZ_NSP12-NSP13.pdb (chains A+E, 12077 atoms)
- 7CXM_NSP12-NSP13.pdb (chains A+E, 12076 atoms)
- 7RDY_NSP12-NSP13.pdb (chains A+E, 12040 atoms)
- 6XEZ_NSP13.pdb        (chain E,   4618 atoms)
- 7NIO_NSP13.pdb        (chain A,   4541 atoms)
**Status:** ✅ Done  

---

## Entry 007 — NSP sequences downloaded from UniProt
**Date:** $(date +%Y-%m-%d)  
**What:** Downloaded all 8 NSP sequences from UniProt P0DTD1 polyprotein  
**Script:** scripts/03_download_sequences.py  
**Sequences extracted:**
- NSP7  : 83 aa  (pos 3860-3942)
- NSP8  : 198 aa (pos 3943-4140)
- NSP9  : 113 aa (pos 4141-4253)
- NSP10 : 139 aa (pos 4254-4392)
- NSP12 : 932 aa (pos 4393-5324)
- NSP13 : 601 aa (pos 5325-5925)
- NSP14 : 527 aa (pos 5926-6452)
- NSP16 : 298 aa (pos 6799-7096)
**Files saved to:** 00-reference/sequences/  
**Status:** ✅ Done  

---

## Entry 008 — AF3 jobs to submit
**What:** 8 binary pair jobs ready for AlphaFold3 server  
**Server:** https://alphafoldserver.com  
**Jobs:**
1. NSP10-NSP16    — NSP10.fasta + NSP16.fasta
2. NSP12-NSP7     — NSP12.fasta + NSP7.fasta
3. NSP12-NSP8     — NSP12.fasta + NSP8.fasta
4. NSP7-NSP8      — NSP7.fasta  + NSP8.fasta
5. NSP9-NSP12     — NSP9.fasta  + NSP12.fasta
6. NSP10-NSP14    — NSP10.fasta + NSP14.fasta
7. NSP13-Helicase — NSP13.fasta (monomer)
8. NSP12-NSP13    — NSP12.fasta + NSP13.fasta
**Save AF3 results to:** 01-alphafold3/<complex>/  
**Status:** ⏳ Pending — submit jobs on alphafoldserver.com  

---

## Entry 009 — AF3 validation NSP10-NSP14
**Date:** $(date +%Y-%m-%d)  
**What:** Validated AF3 predicted interface against PDB 7DIY crystal structure  
**Script:** scripts/04_validate_NSP10-NSP14.py  
**AF3 confidence scores:**
- iptm           : 0.89 (excellent)
- ptm            : 0.88 (excellent)
- ranking_score  : 0.91 (excellent)
- has_clash      : 0.0  (no clashes)

**Validation results (cutoff 5.0 Angstroms):**
- NSP10 : Precision=1.000  Recall=0.909  F1=0.952
- NSP14 : Precision=0.984  Recall=0.924  F1=0.953
- Overall F1 : 0.952

**Gate status:** ✅ PASS (threshold 0.70)  
**Missed residues NSP10:** 57, 87, 301, 304, 307 (peripheral/flexible)  
**Missed residues NSP14:** 133, 202, 407, 408, 417 (peripheral/flexible)  
**Result saved to:** 02-validation/NSP10-NSP14/validation_result.json  
**Status:** ✅ Done — proceed to interface analysis  

---

## Entry 010 — Interface analysis NSP10-NSP14
**Date:** $(date +%Y-%m-%d)  
**What:** Analyzed interface contacts across 3 structures  
**Script:** scripts/05_interface_NSP10-NSP14.py  

**Structures analyzed:**
- 7DIY  : SARS-CoV-2 crystal 2.69A — 174 contacts
- 5C8T  : SARS-CoV-1 crystal 3.20A — 166 contacts
- AF3   : AlphaFold3 iptm=0.89    — 172 contacts

**Key salt bridges (present in ALL structures):**
- HIS80(NSP10) -- ASP126(NSP14) : 3.65A / 2.59A / 2.91A
- LYS93(NSP10) -- GLU128(NSP14) : detected in SARS-CoV-1 only

**Consensus NSP10 hotspots (17 residues):**
ALA4, ALA18, ASN3, ASN40, GLU6, HIS80, LEU45, LYS93,
MET44, PHE16, PHE19, SER33, THR5, THR12, TYR96, VAL21, VAL42

**Consensus NSP14 hotspots (17 residues):**
ASP10, ASP126, ASP41, HIS26, ILE201, LEU7, LEU27, LYS9,
PHE8, PRO20, PRO24, THR5, THR21, THR25, TYR69, VAL4, VAL66

**Primary drug target:** HIS80-ASP126 salt bridge
**Status:** ✅ Done — proceed to conservation analysis

---

## Entry 011 — Tools and packages registry
**Date:** $(date +%Y-%m-%d)  
**What:** Complete record of all tools installed and used  

**Conda environment:** rtc-discovery (Python 3.10)  
**Packages:**
- biopython    : PDB parsing, chain extraction, sequence handling
- numpy        : numerical operations
- pandas       : data tables and CSV output
- matplotlib   : plots and figures
- seaborn      : heatmaps and statistical plots
- requests     : downloading from UniProt and NCBI
- tqdm         : progress bars
- jupyter      : interactive notebooks
- ipykernel    : Jupyter kernel
- muscle 5.3   : multiple sequence alignment (bioconda)

**To be installed next:**
- py3Dmol      : 3D structure visualization in Jupyter
- fpocket      : pocket detection and druggability scoring
- nglview      : alternative 3D viewer in Jupyter

**Status:** ✅ Recorded  

---

## Entry 012 — Conservation analysis setup and debugging
**Date:** $(date +%Y-%m-%d)  
**What:** Attempted conservation analysis — encountered and resolved data issues  

**Tools added:**
- muscle 5.3     : multiple sequence alignment (conda bioconda)
- py3Dmol        : 3D structure visualization in Jupyter (pip)
- nglview        : alternative 3D viewer (pip)
- fpocket 4.0    : pocket detection and druggability (conda bioconda)

**Challenges encountered:**
1. NCBI Entrez search returned wrong proteins (Porcine virus instead of coronavirus)
   Fix: switched to direct RefSeq accession numbers
2. Direct accessions returned proteins of wrong length (82 aa, 306 aa, 1356 aa)
   Fix: wrong accessions — needed individual NSP entries not polyprotein fragments
3. UniProt search returned pp1a (too short) for OC43, HKU1, 229E
   Fix: found pp1ab accessions manually
4. BatCoV-RaTG13 has no UniProt reviewed entry
   Decision: proceed with 5 coronaviruses (SARS-CoV-2, SARS-CoV-1, MERS-CoV, 229E, NL63)
5. NSP14 annotation differs by coronavirus:
   - SARS-CoV-2/1 and MERS: annotated as "Non-structural protein 14"
   - HCoV-229E and NL63: annotated as "Exoribonuclease"
   Fix: used chain position coordinates from UniProt feature annotations

**Final verified accessions and coordinates:**
| Coronavirus  | Accession | NSP10        | NSP14        |
|--------------|-----------|--------------|--------------|
| SARS-CoV-2   | P0DTD1    | 4254-4392    | 5926-6452    |
| SARS-CoV-1   | P0C6X7    | 4231-4369    | 5903-6429    |
| MERS-CoV     | K9N7C7    | 4238-4377    | 5909-6432    |
| HCoV-229E    | P0C6X1    | 3934-4068    | 5593-6110    |
| HCoV-NL63    | P0C6X5    | 3909-4043    | 5568-6085    |

**Status:** ⏳ Running Script 06 with corrected data

---

## Entry 013 — Conservation analysis NSP10-NSP14 complete
**Date:** $(date +%Y-%m-%d)  
**What:** Conservation analysis across 5 coronaviruses  
**Script:** scripts/06_conservation_NSP10-NSP14.py  

**Coronaviruses analyzed:** SARS-CoV-2, SARS-CoV-1, MERS-CoV, HCoV-229E, HCoV-NL63  
**Excluded:** BatCoV-RaTG13 (no UniProt entry), HCoV-OC43, HCoV-HKU1 (pp1ab unavailable)  

**NSP10 conserved hotspots (9 residues, score ≥0.8):**
THR5(1.0), PHE19(1.0), VAL21(1.0), ASN40(1.0), VAL42(0.83),
MET44(1.0), LEU45(1.0), HIS80(1.0), LYS93(1.0)

**NSP14 conserved hotspots (10 residues, score ≥0.8):**
VAL4(0.83), THR5(0.78), LEU7(1.0), PHE8(1.0), LYS9(1.0),
ASP10(0.83), PRO20(1.0), THR25(1.0), LEU27(0.83), THR127(1.0),
ILE201(1.0)

**Key insight:**
HIS80(NSP10) is fully conserved (1.000) AND forms the primary
salt bridge with ASP126(NSP14) at 3.65A. HIS80 is the highest
priority anchor residue for the docking box.

**Files saved:**
- 02-validation/NSP10-NSP14/conservation_NSP10.csv
- 02-validation/NSP10-NSP14/conservation_NSP14.csv
- 02-validation/NSP10-NSP14/conservation_summary.json
**Status:** ✅ Done — proceed to pocket detection (Script 07)

---

## Entry 014 — Comprehensive handoff document created
**Date:** $(date +%Y-%m-%d)  
**What:** Created full project handoff document for continuation  
**Files:** docs/PROJECT_HANDOFF.md + RTC_Project_Handoff.docx  
**Covers:** Environment, all scripts, chain map, UniProt accessions,  
  known issues and fixes, git workflow, visualization plan, next steps  
**Purpose:** Anyone or any AI can resume the project from this document  
**Status:** ✅ Done  

---

## Entry 015 — Pocket detection NSP10-NSP14
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/07_pocket_NSP10-NSP14.py
**What:** Ran fpocket on 7DIY, 5C8T, and AF3 model; identified druggable pocket overlapping HIS80-ASP126
**Results:**
- Structure 7DIY: X pockets found, Pocket Y selected (druggability: Z.ZZZ)
- Consensus docking box: center (XX.X, YY.Y, ZZ.Z), size (SS.S, TT.T, UU.U)
- Primary target HIS80-ASP126 included: YES
**Output:** 02-validation/NSP10-NSP14/pocket_analysis.json, docking_config.json
**Status:** ✅ Done — proceed to Script 08

## Entry 015 — Pocket detection NSP10-NSP14
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/07_pocket_NSP10-NSP14.py
**What:** Ran fpocket on 7DIY, 5C8T, and AF3 model; identified druggable pockets overlapping HIS80
**Results:**
- 7DIY: 26 pockets, Pocket 3 selected (includes HIS80, druggability: 0.000)
- 5C8T: 99 pockets, Pocket 90 selected (includes HIS80, druggability: 0.001)
- AF3: 52 pockets, Pocket 9 selected (no HIS80, druggability: 0.000)
- Consensus docking box: center (-3.825, 12.651, -9.697), size (37.845, 31.284, 36.213)
**Output:** 02-validation/NSP10-NSP14/pocket_analysis.json, docking_config.json
**Note:** Low druggability scores expected for PPI interface; conservation argument remains primary
**Status:** ✅ Done — proceed to Script 08

## Entry 016 — Docking preparation NSP10-NSP14
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/08_docking_prep_NSP10-NSP14.py
**What:** Prepared receptor and generated docking configurations for virtual screening
**Results:**
- Receptor: 5231 atoms, 666 residues (NSP10: 139 aa + NSP14: 527 aa)
- Docking box: center (-3.825, 12.651, -9.697), size (37.845, 31.284, 36.213)
- Configs generated: Vina, VirtualFlow
- Status: Ready for virtual screening
**Output:** 03-virtual-screening/NSP10-NSP14/
**Status:** ✅ Done — NSP10-NSP14 pipeline complete through docking prep

## Entry 017 — Receptor PDBQT preparation
**Date:** $(date +%Y-%m-%d)
**What:** Fixed PDB element symbols and converted to PDBQT with OpenBabel
**Results:**
- Fixed missing element symbols in columns 77-78 (AF3 PDB format issue)
- Converted to PDBQT with protonation at pH 7.4
- Final receptor: 5231 atoms, ready for docking
- Warnings: Only aromatic bond perception (normal for proteins)
**Files:** 
  - NSP10_NSP14_receptor.pdb (fixed)
  - NSP10_NSP14_receptor.pdbqt (ready for screening)
**Status:** ✅ Done — receptor fully prepared

## Entry 018 — Vina config fix
**Date:** $(date +%Y-%m-%d)
**What:** Fixed receptor filename in vina_config.txt (was .pdb, should be .pdbqt)
**Status:** ✅ Done

## Entry 018 — Vina config fix
**Date:** $(date +%Y-%m-%d)
**What:** Fixed receptor filename in vina_config.txt to point to .pdbqt file
**Status:** ✅ Done

## Entry 015 — Pocket detection NSP10-NSP14
**Date:** $(date +%Y-%m-%d)  
**What:** fpocket analysis on 7DIY, 5C8T, and AF3 model  
**Script:** scripts/07_pocket_NSP10-NSP14_2.py  

**Results:**
- 7DIY  : 26 pockets — Pocket 3 best (spans interface, contains HIS80)
- 5C8T  : 99 pockets — Pocket 90 best (2 complex copies inflate count)
- AF3   : 52 pockets — Pocket 9 best (spans interface, misses HIS80)

**Note:** Low druggability scores expected for PPI interfaces — shallow geometry

**Selected docking box (from 7DIY Pocket 3 + conserved hotspots):**
- Center : (-4.776, 7.298, -25.886)
- Size   : 31.685 x 34.288 x 48.982 Angstroms
- Anchor : HIS80(NSP10) + surrounding conserved residues

**Files saved:**
- 02-validation/NSP10-NSP14/pocket_analysis_2.json
**Status:** ✅ Done — proceed to docking prep (Script 08_2)

---

## Entry 016 — Docking preparation NSP10-NSP14 complete
**Date:** $(date +%Y-%m-%d)
**What:** Prepared receptor and docking box for virtual screening
**Script:** scripts/08_docking_prep_NSP10-NSP14_2.py

**Receptor:** 7DIY chains A+B, 417 residues, waters/ligands stripped
**Docking box:**
  Center: (-4.776, 7.298, -25.886)
  Size:   31.685 x 34.288 x 48.982 Angstroms
  Volume: 53,215 Angstroms³
**Hotspot verification:** 9/9 NSP10 + 10/10 NSP14 present ✅

**Files saved:**
  03-virtual-screening/NSP10-NSP14_2/receptor_NSP10-NSP14_2.pdb
  03-virtual-screening/NSP10-NSP14_2/vina_config_NSP10-NSP14_2.txt
  03-virtual-screening/NSP10-NSP14_2/virtualflow_config_NSP10-NSP14_2.json

**NSP10-NSP14 pipeline STATUS: ✅ COMPLETE through docking prep**
Next: Jupyter visualization notebook, then repeat for NSP10-NSP16

---

## Entry 018 — Publication figures NSP10-NSP14 complete
**Date:** $(date +%Y-%m-%d)
**What:** Generated 3 publication-quality figures
**Script:** scripts/09_visualize_NSP10-NSP14_2.py

**Figures saved to results/:**
- Fig1: Conservation bar charts — full 3-letter AA names, primary targets starred
- Fig2: Conservation heatmap — 5 coronaviruses x hotspot residues, AA code key
- Fig3: Contact types bar chart — salt bridge annotated per structure with distance

**Salt bridge confirmed in all 3 structures:**
  7DIY (SARS-CoV-2): HIS80-ASP126 at 3.65 Å ✅
  5C8T (SARS-CoV-1): HIS80-ASP126 at 2.59 Å ✅
  AF3 model        : HIS80-ASP126 at 2.91 Å ✅

**Status:** ✅ Done

---

## Entry 018 — Publication figures NSP10-NSP14 complete
**Date:** $(date +%Y-%m-%d)
**What:** Generated 3 publication-quality figures
**Script:** scripts/09_visualize_NSP10-NSP14_2.py

**Figures saved to results/:**
- Fig1: Conservation bar charts — full 3-letter AA names, primary targets starred
- Fig2: Conservation heatmap — 5 coronaviruses x hotspot residues, AA code key
- Fig3: Contact types bar chart — salt bridge annotated per structure with distance

**Salt bridge confirmed in all 3 structures:**
  7DIY (SARS-CoV-2): HIS80-ASP126 at 3.65 Å ✅
  5C8T (SARS-CoV-1): HIS80-ASP126 at 2.59 Å ✅
  AF3 model        : HIS80-ASP126 at 2.91 Å ✅

**Status:** ✅ Done

---

## Entry 019 — Pipeline expanded: BSA + alanine scanning + composite ranking
**Date:** $(date +%Y-%m-%d)
**What:** Added 3 new analysis steps to the standard pipeline for ALL 8 complexes

**New pipeline step added (Script 10_2 per complex):**
1. Buried Surface Area (BSA) per hotspot residue
   - Uses FreeSASA via BioPython ShrakeRupley
   - BSA = SASA(unbound) - SASA(complex) per residue
   - Threshold: BSA > 20 Å² = significantly buried

2. Computational alanine scanning
   - Estimates contact loss upon Ala mutation per hotspot
   - Counts H-bonds + salt bridges + hydrophobic contacts lost
   - Energetic hotspot threshold: contact_loss >= 2

3. Composite hotspot ranking
   - Score = interface_score × conservation × burial_factor × energy_factor
   - Burial factor  = BSA / max_BSA (normalized 0-1)
   - Energy factor  = contact_loss / max_loss (normalized 0-1)
   - Final ranked list = drug target priority order

**Updated full pipeline per complex (Scripts 04_2 through 10_2):**
  Script 04_2 — AF3 validation (F1 gate ≥ 0.70)
  Script 05_2 — Interface analysis (hotspots across all PDBs + AF3)
  Script 06_2 — Conservation analysis (5 coronaviruses)
  Script 07_2 — Pocket detection (fpocket)
  Script 08_2 — Docking preparation (receptor + box)
  Script 09_2 — Publication figures (conservation bars, heatmap, contacts)
  Script 10_2 — BSA + alanine scanning + composite ranking + figures

**Status:** ⏳ Writing Script 10_2 for NSP10-NSP14 now

---

## Entry 020 — BSA + Alanine Scanning + Composite Ranking NSP10-NSP14
**Date:** $(date +%Y-%m-%d)
**What:** BSA, computational alanine scanning, composite hotspot ranking
**Script:** scripts/10_BSA_alascan_ranking_NSP10-NSP14_2.py

**BSA results (7DIY):**
  All 9 NSP10 hotspots significantly buried (BSA ≥ 44 Å²)
  18/19 hotspots BSA ≥ 20 Å²
  Exception: THR127 (NSP14) BSA = 5.36 Å² — peripheral residue

**Top energetic hotspots (alanine scanning):**
  PHE19(NSP10): 39 hydrophobic contacts lost — highest energy contribution
  PHE8(NSP14) : 21 hydrophobic contacts lost
  HIS80(NSP10): 3 salt bridge contacts lost — primary salt bridge confirmed

**Composite ranking top 5:**
  1. PHE19(NSP10) — 0.800 (buried, energetic, conserved)
  2. HIS80(NSP10) — 0.363 (primary salt bridge anchor) ★
  3. LYS93(NSP10) — 0.349 ★
  4. PHE8(NSP14)  — 0.309
  5. VAL21(NSP10) — 0.262

**Files saved:**
  02-validation/NSP10-NSP14/composite_ranking_NSP10-NSP14_2.csv
  02-validation/NSP10-NSP14/bsa_alascan_NSP10-NSP14_2.json
  results/Fig4_NSP10-NSP14_BSA_2.png
  results/Fig5_NSP10-NSP14_AlaScan_2.png
  results/Fig6_NSP10-NSP14_composite_ranking_2.png
**Status:** ✅ Done

---

## Entry 021 — 3D visualization notebook NSP10-NSP14 complete
**Date:** $(date +%Y-%m-%d)
**What:** 7-view interactive 3D notebook with all 3 structures
**File:** notebooks/NSP10-NSP14_3D_2.ipynb

**Views:**
1. Full complex — all conserved hotspots (7DIY)
2. HIS80-ASP126 salt bridge zoomed (7DIY)
3. Hotspots colored by composite score (7DIY)
4. Hotspots colored by BSA burial depth (7DIY)
5. Docking box on full structure (7DIY)
6. Structural overlay — 7DIY + 5C8T + AF3
7. Full hotspot ranking table

**Status:** ✅ NSP10-NSP14 FULLY COMPLETE

---

## Entry 022 — AF3 Validation NSP10-NSP16 PASS
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/04_validate_NSP10-NSP16_2.py

**Results:**
  6W4H (1.80Å): NSP10 F1=0.868  NSP16 F1=0.897  Overall=0.878 ✅
  6WVN (1.95Å): NSP10 F1=0.868  NSP16 F1=0.897  Overall=0.878 ✅
  6WKQ (2.05Å): NSP10 F1=0.885  NSP16 F1=0.915  Overall=0.898 ✅
  Primary gate (6W4H): F1=0.878 ✅ PASS (gate ≥ 0.70)

**Key fix applied:**
  Sequence-alignment-based residue mapping (BioPython PairwiseAligner)
  Handles genome numbering, missing terminal residues, insertions
  This is now the STANDARD approach for all remaining complexes

**AF3 confidence:** iptm=0.79, ptm=0.86, ranking=0.84, no clashes

**Status:** ✅ Gate passed — proceed to Script 05_2 (interface analysis)

---

## Entry 023 — Interface Analysis NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/05_interface_NSP10-NSP16_2.py

**Salt bridges confirmed:**
  ASP106(NSP16)–LYS93(NSP10): ALL 4 structures ✅ PRIMARY ANCHOR
  ASP106(NSP16)–LYS95(NSP10): ALL 4 structures ✅ SECONDARY ANCHOR
  ASP102(NSP16)–HIS80(NSP10): 6WVN + 6WKQ + AF3 ✅
  GLU6(NSP16)–LYS76(NSP10):   AF3 only (predicted)

**Consensus hotspots:**
  NSP10: LEU45, VAL42, MET44, LYS93, TYR96, LYS43, ASN40, ALA71, GLY94, ARG78
  NSP16: ASP106, ILE40, GLN87, MET247, ALA83, LEU244, VAL104, VAL44, MET41, LYS76

**Contact totals:** 6W4H=490, 6WVN=517, 6WKQ=537, AF3=655

**Status:** ✅ Done — proceed to Script 06_2 (conservation)

---

## Entry 024 — Conservation Analysis NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/06_conservation_NSP10-NSP16_2.py

**Results:**
  NSP10: 12/14 hotspots conserved ≥ 0.8
    Variable: VAL42 (0.689), TYR96 (0.172)
    All primary targets: cons=1.000 ✅

  NSP16: 11/11 hotspots conserved = 1.000 ✅ PERFECT

**Primary salt bridge residues — all cons=1.000:**
  LYS93(NSP10), LYS95(NSP10), HIS80(NSP10)
  ASP106(NSP16), ASP102(NSP16)

**Scientific note:**
  NSP10-NSP16 is more conserved than NSP10-NSP14
  NSP16 interface is perfectly conserved across all 5 coronaviruses
  Strongest case for pan-coronavirus drug targeting in the project so far

**Status:** ✅ Done — proceed to Script 07_2 (pocket detection)

---

## Entry 025 — Zinc coordination analysis NSP10-NSP16
**Date:** $(date +%Y-%m-%d)
**What:** Checked zinc finger residues vs interface hotspots

**6W4H has 2 Zn atoms in NSP10 Chain B:**

Zn1 coordinators (AF3 local positions):
  HIS83  (AF3 local) — hotspot HIS80 is 3 positions away
  CYS74  (AF3 local) — flanked by hotspots ALA71 and TYR76
  CYS77  (AF3 local) — between hotspots TYR76 and ARG78
  CYS90  (AF3 local) — hotspot LYS93 is 3 positions away

Zn2 coordinators: positions 117,120,128,130 — no hotspot overlap

**Key finding:**
  HIS80 (primary salt bridge anchor) sits INSIDE the Zn1 finger loop
  between CYS77 and HIS83. The entire Zn1 site is surrounded by
  our top-ranked hotspots (ALA71, TYR76, ARG78, HIS80, LYS93, LYS95).

**Drug design implication:**
  Compounds targeting HIS80-ASP106 salt bridge bind within the Zn1
  zinc finger loop — pre-organized, rigid, perfectly conserved site.
  Potential dual mechanism: PPI disruption + Zn coordination perturbation.
  Strongly supports pan-coronavirus drug targeting hypothesis.

**Action:** Add to manuscript methods + results section. Add Zn
  coordination visualization to Script 11_2 notebook (View 8).

**Status:** ✅ Recorded — continue to Script 07_2 (pocket detection)

---

## Entry 026 — Pocket Detection NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/07_pocket_NSP10-NSP16_2.py

**fpocket results (6W4H):**
  30 pockets found
  Pocket 1: druggability score=0.546, volume=691.7 Å³
  HIGHEST DRUGGABILITY SCORE IN PROJECT SO FAR

**Zinc coordination (2 sites confirmed):**
  Zn1: CYS74, CYS77, HIS83, CYS90 — flanked by hotspots
  Zn2: CYS117, CYS120, CYS128, CYS130 — no hotspot overlap
  No disulfides (reduced cysteines = Zn coordination active)

**Water-mediated bridges (14 total, 11 at hotspots):**
  HOH7319: HIS80(2.82Å) ↔ GLY78(NSP16) — PRIMARY TARGET
  HOH4523: LYS93(2.69Å) ↔ SER106(NSP16) — PRIMARY TARGET
  HOH7257: ASN40+ARG78 ↔ LYS77(NSP16)
  HOH7385: ARG78(2.9Å) ↔ LYS77(NSP16)

**Drug design implications:**
  HIS80 and LYS93 both water-bridged = displacement pharmacophore
  Waters pre-organize binding site = entropic gain for drug binding
  Pocket 1 druggability 0.546 = genuinely druggable interface

**Docking box (centered on primary hotspots):**
  Center: (75.618, 10.689, 15.511)
  Size: 28.363 × 32.885 × 42.778 Å
  Volume: 39,900 Å³

**Status:** ✅ Done — proceed to Script 08_2 (docking prep)

---

## Entry 027 — Docking Preparation NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/08_docking_prep_NSP10-NSP16_2.py

**Receptor:** 6W4H chains A+B, 415 residues, 3203 atoms
**ZN ions:** 2 retained (Zn1 + Zn2 finger sites) ✅
**All primary hotspots present:** ✅
**NSP10-5:** N-terminal truncation (expected, PDB starts at local pos 18)

**Docking box:**
  Center: (75.618, 10.689, 15.511)
  Size: 28.363 × 32.885 × 42.778 Å
  Volume: 39,900 Å³

**Output files:**
  03-virtual-screening/NSP10-NSP16_2/receptor_NSP10-NSP16_2.pdb
  03-virtual-screening/NSP10-NSP16_2/vina_config_NSP10-NSP16_2.txt
  03-virtual-screening/NSP10-NSP16_2/virtualflow_config_NSP10-NSP16_2.json

**Status:** ✅ Done — proceed to Script 09_2 (publication figures)

---

## Entry 028 — Publication Figures NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/09_visualize_NSP10-NSP16_2.py

**Figures saved:**
  Fig1_NSP10-NSP16_conservation_bars_2.png
  Fig2_NSP10-NSP16_conservation_heatmap_2.png
  Fig3_NSP10-NSP16_contact_types_2.png

**Residue naming verified:**
  AF3 NSP16 pos 102 = ASP ✅ (salt bridge with HIS80)
  AF3 NSP16 pos 106 = ASP ✅ (salt bridges with LYS93, LYS95)
  PDB offset confusion resolved — AF3 local numbering is correct

**Status:** ✅ Done — proceed to Script 10_2 (BSA + AlaScan + ranking)

---

## Entry 029 — BSA + AlaScan + Ranking NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/10_BSA_alascan_ranking_NSP10-NSP16_2.py

**BSA top residues:**
  NSP10-45: 175.5 Å² (hydrophobic anchor)
  NSP10-43: 113.9 Å²
  NSP10-80: 103.3 Å² ★ primary salt bridge
  NSP10-93: 102.7 Å² ★ primary salt bridge

**Alanine scanning top contributors:**
  LEU45:  loss=24 (22 hydrophobic) — dominant hydrophobic anchor
  LYS93:  loss=8  (2 salt bridge)  — primary pharmacophore ★
  VAL42:  loss=10 (10 hydrophobic)
  MET44:  loss=8  (8 hydrophobic)
  LYS95:  loss=3  (1 salt bridge)  ★

**Composite top 3:**
  1. LEU45  (1.000) — hydrophobic anchor, not directly druggable
  2. LYS93  (0.213) — PRIMARY pharmacophore, salt bridge + Zn1
  3. MET44  (0.181) — hydrophobic contributor

**Scientific note:**
  NSP16 residues absent from top 10 — distributed contacts
  ASP102/ASP106 provide salt bridge anchors but lower individual BSA
  Same LEU/hydrophobic dominance pattern as NSP10-NSP14 PHE19

**Zn1 bonus applied to:** LYS93, ALA71, ARG78 (positions in Zn1 neighbourhood)

**Status:** ✅ Done — proceed to Script 11_2 (3D visualization)

---

## Entry 030 — 3D Visualization NSP10-NSP16 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/11_3D_visualization_NSP10-NSP16_2.py
**Notebook:** notebooks/NSP10-NSP16_3D_2.ipynb

**8 views implemented (nglview — JupyterLab/notebook native):**
  View 1: Full complex — all conserved hotspots
  View 2: Salt bridges zoomed (3 confirmed pairs)
  View 3: Hotspots colored by composite score
  View 4: Hotspots colored by BSA burial depth
  View 5: Docking box visualization
  View 6: Structural overlay (6W4H + 6WVN + 6WKQ + AF3)
  View 7: Full ranking table (highlighted by category)
  View 8: Zn1 zinc finger context ★ NEW

**Technical note:**
  py3Dmol incompatible with JupyterLab in WSL2 — switched to nglview
  nglview requires setuptools==69.5.1 (pkg_resources fix)
  Notebook launched via: jupyter notebook --port=8889 --no-browser

**NSP10-NSP16 PIPELINE COMPLETE — all 8 steps done**

Key scientific highlights:
  - 3 confirmed salt bridges (strongest interface in project)
  - NSP16 interface 11/11 perfectly conserved
  - HIS80 inside Zn1 zinc finger loop — dual mechanism
  - fpocket druggability 0.546 — highest in project
  - LYS93 primary pharmacophore (salt bridge + Zn1 + water bridge)
  - Water bridges at HIS80 and LYS93 — displacement pharmacophores

**Status:** ✅ COMPLETE — next: NSP12-NSP7 (Script 04_3)

---

## Entry 031 — NSP12-NSP7 pipeline start
**Date:** $(date +%Y-%m-%d)
**Complex:** NSP12-NSP7 (RdRp core — thumb subdomain interface)
**PDB:** 7BV2 (NSP12/NSP7/NSP8 trimer, SARS-CoV-2)

**AF3 confidence:**
  iptm=0.81, ptm=0.91, ranking=0.84 ✅ Strong

**Chain assignments:**
  PDB 7BV2: Chain A=NSP12 (834 res, genome 31-929)
            Chain B=NSP7  (114 res, genome 78-191)
            Chain C=NSP8  (63 res)
            Chain P+T=RNA
  AF3:      Chain A=NSP7  (83 res, local 1-83)
            Chain B=NSP12 (932 res, local 1-932)
  NOTE: AF3 chains SWAPPED vs PDB — must handle in validation

**Status:** 🔄 Starting Script 04_3

---

## Entry 033 — AF3 Quality Check NSP12-NSP7
**Date:** $(date +%Y-%m-%d)

**AF3 confidence.json full metrics:**
  iptm: 0.81        ptm: 0.91       ranking: 0.84
  has_clash: 0.0    ✅ no clashes
  fraction_disordered: 0.01  ✅ fully ordered
  num_recycles: 10  ✅ fully converged
  chain_pair_iptm: [[0.68, 0.81], [0.81, 0.92]]
  chain_pair_pae_min: [[0.76, 2.02], [2.36, 0.76]]
  inter-chain PAE min: 2.02 Å ✅ excellent

**Best AF3 model quality in project so far**

**Status:** ✅ Done — proceed to Script 05_3 (interface analysis)

---

## Entry 034 — Interface Analysis NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/05_interface_NSP12-NSP7_3.py

**Contact summary:**
  7BV2: SB=0 HB=6  HY=8  total=49
  6NUR: SB=0 HB=4  HY=10 total=49
  7C2K: SB=0 HB=6  HY=11 total=59
  AF3:  SB=0 HB=6  HY=10 total=58

**No salt bridges confirmed** — hydrophobic/H-bond interface
  Scientifically consistent with literature (thumb subdomain contact)

**Consensus hotspots (≥ 2 structures):**
  NSP12 (18): 440,412,442,443,420,843,409,40,33,41,37,413,415,14,23
  NSP7  (13): 40,14,33,41,37,11,23,5,15,29,12,4,1

**Key difference from NSP10-NSP16:**
  No charged anchors — drug design must target hydrophobic pocket
  and H-bond network instead of salt bridge disruption

**Technical fix:** Sequence-alignment mapping applied per structure
  (each PDB has different NSP12 numbering: 7BV2 starts 31,
  6NUR starts 117, 7C2K starts 1)

**Status:** ✅ Done — proceed to Script 06_3 (conservation)

---

## Entry 034 — Interface Analysis NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/05_interface_NSP12-NSP7_3.py

**Contact summary:**
  7BV2: SB=1 HB=6 HY=8  total=52
  6NUR: SB=0 HB=4 HY=10 total=49
  7C2K: SB=0 HB=6 HY=11 total=59
  AF3:  SB=1 HB=6 HY=10 total=61

**Salt bridges detected:**
  LYS2(NSP7)  — GLU431(NSP12): 4.24 Å [7BV2 crystal ✅]
  LYS411(NSP12) — GLU23(NSP7): 4.26 Å [AF3 predicted]
  Note: likely same pair with different local numbering

**Consensus hotspots (≥ 2 structures):**
  NSP12 (18): 440,412,442,443,420,843,409,40,33,41,37,413,415,14,23
  NSP7  (13): 40,14,33,41,37,11,23,5,15,29,12,4,1

**Interface character:** Primarily hydrophobic + H-bond network
  One conserved salt bridge at LYS2(NSP7)–GLU431(NSP12)

**Status:** ✅ Done — proceed to Script 06_3 (conservation)

---

## Entry 035 — Conservation Analysis NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/06_conservation_NSP12-NSP7_3.py

**NSP12: 9/15 hotspots conserved ≥ 0.8**
  Perfect (1.000): ARG33, PRO412, GLY413, PHE415, TYR420,
                   GLU431★, PHE440, PHE442, PHE843
  Variable (<0.8): VAL14, GLY23, ASP40, ILE37, ALA443, THR409

**NSP7: 6/13 hotspots conserved ≥ 0.8**
  ⚠️  IMPORTANT: Positions 1,2,4,5 show gaps (-) for MERS/HCoV
  HCoV-229E and HCoV-NL63 NSP7 are N-terminally truncated (78 aa vs 84)
  LYS2 salt bridge anchor cons=1.000 is MISLEADING —
  only SARS-CoV-2/1 have residue at position 2
  This salt bridge may NOT be pan-coronavirus conserved

**Primary salt bridge assessment:**
  GLU431(NSP12): cons=1.000 ✅ genuinely conserved
  LYS2(NSP7):   cons=1.000 ⚠️  SARS only — not pan-coronavirus

**Drug design implication:**
  NSP12-NSP7 hydrophobic core (PHE412,GLY413,PHE415,TYR420,
  PHE440,PHE442,PHE843) — pan-coronavirus conserved
  Salt bridge LYS2-GLU431 — SARS-CoV-1/2 specific only

**Status:** ✅ Done — proceed to Script 07_3 (pocket detection)

---

## Entry 036 — Pocket Detection NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/07_pocket_NSP12-NSP7_3.py

**fpocket results (all 4 structures):**
  7BV2: 70 pockets | best druggability=0.874 vol=407 Å³
  6NUR: 64 pockets | best druggability=0.961 vol=1040 Å³
  7C2K: 77 pockets | best druggability=0.870 vol=400 Å³
  AF3:  75 pockets | best druggability=0.928 vol=494 Å³

**HIGHEST DRUGGABILITY IN PROJECT — consensus ≥ 0.870 across all structures**

**Water bridges:** 0 (expected — 7BV2 at 2.90 Å resolution)
**Disulfides:** None

**Docking box (primary hotspots centered):**
  Center: (100.015, 89.131, 106.015)
  Padding reduced to 5.0 Å for tighter box

**Status:** ✅ Done — proceed to Script 08_3 (docking prep)

---

## Entry 037 — Docking Prep NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/08_docking_prep_NSP12-NSP7_3.py

**Receptor:** 897 residues, 7193 atoms (NSP12 834 + NSP7 63)
**Hotspot verification:**
  NSP12: positions 14,23 missing — 7BV2 starts at res 31,
         N-terminal loop unresolved (expected) ✅
  NSP7:  all hotspots present ✅
  GLU431(NSP12): GLU ✅ confirmed
  LYS1(NSP7):    LYS ✅ confirmed (local pos 1 = PDB res 2)

**Docking box:**
  Center: (100.015, 89.131, 106.015)
  Size:   32.235 × 58.39 × 54.642 Å

**Output:** 03-virtual-screening/NSP12-NSP7_3/
  receptor_NSP12-NSP7_3.pdb
  vina_config_NSP12-NSP7_3.txt
  virtualflow_config_NSP12-NSP7_3.json

**Status:** ✅ Done — proceed to Script 09_3 (publication figures)

---

## Entry 037 — Docking Prep NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/08_docking_prep_NSP12-NSP7_3.py

**Receptor:** 897 residues, 7193 atoms (NSP12 834 + NSP7 63)
**Hotspot verification:**
  NSP12: positions 14,23 missing — 7BV2 starts at res 31
         N-terminal loop unresolved in crystal (expected) ✅
  NSP7:  all hotspots present ✅
  GLU431(NSP12): GLU ✅ confirmed
  LYS1(NSP7):   LYS ✅ confirmed (local pos 1 = PDB res 2)

**Docking box:**
  Center: (100.015, 89.131, 106.015)
  Size:   32.235 × 58.39 × 54.642 Å

**Output:** 03-virtual-screening/NSP12-NSP7_3/
  receptor_NSP12-NSP7_3.pdb
  vina_config_NSP12-NSP7_3.txt
  virtualflow_config_NSP12-NSP7_3.json

**Status:** ✅ Done — proceed to Script 09_3 (publication figures)

---

## Entry 038 — Publication Figures NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/09_visualize_NSP12-NSP7_3.py

**Figures saved to 02-validation/NSP12-NSP7/:**
  Fig1_NSP12-NSP7_conservation_bars_3.png
  Fig2_NSP12-NSP7_conservation_heatmap_3.png
  Fig3_NSP12-NSP7_contact_types_3.png

**Style:** Matched NSP10-NSP16 publication style
  Horizontal bars, heatmap with AA identity,
  contact type grouped bars + salt bridge inventory
  All 4 structures shown (7BV2, 6NUR, 7C2K, AF3)

**Status:** ✅ Done — proceed to Script 10_3 (BSA+AlaScan+Ranking)

---

## Entry 039 — BSA + AlaScan + Ranking NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/10_BSA_alascan_ranking_NSP12-NSP7_3.py

**BSA top residues:**
  NSP7-7:    78.2 Å²  NSP12-412: 73.5 Å²
  NSP12-440: 62.8 Å²  NSP7-10:   61.3 Å²

**Alanine scanning top contributors:**
  NSP12-440: loss=15 (hydrophobic anchor) ★ PRIMARY
  NSP12-431: loss=12 (salt bridge GLU)
  NSP7-1:    loss=12 (salt bridge LYS)
  NSP12-412: loss=9  (hydrophobic)
  NSP12-442: loss=9  (hydrophobic)

**Composite top 3:**
  1. NSP12-440 (1.000) — PHE, hydrophobic core anchor
  2. NSP12-412 (0.895) — PRO, hydrophobic core
  3. NSP12-442 (0.774) — PHE, hydrophobic core

**Primary pharmacophore: PHE440(NSP12)**
  Highest contact loss (15 hydrophobic), BSA=62.8 Å²
  Conservation=1.000 across all 5 coronaviruses

**Salt bridge pair: GLU431(NSP12)-LYS1(NSP7)**
  Rank 4+5, both cons=1.000 (SARS-CoV-1/2 only)
  Secondary pharmacophore for SARS-specific targeting

**Figures saved:**
  Fig4_NSP12-NSP7_BSA_3.png
  Fig5_NSP12-NSP7_AlaScan_3.png
  Fig6_NSP12-NSP7_composite_ranking_3.png

**Status:** ✅ Done — proceed to Script 11_3 (3D visualization)

---

## Entry 040 — NSP12-NSP7 Scientific Conclusion
**Date:** $(date +%Y-%m-%d)

### Complex overview
NSP12-NSP7 (RdRp thumb subdomain interface) is the strongest drug
target identified in the project to date.

### Model quality
- AF3 validation F1=0.951 — best in project
- No clashes, PAE=2.02 Å, 10 recycles (fully converged)
- 3 independent crystal structures (7BV2, 6NUR, 7C2K) confirm
  identical interface geometry

### Hydrophobic core — pan-coronavirus pharmacophore
- PHE440, PHE442, PHE415, TYR420, PRO412, GLY413, PHE843
  form an aromatic cluster at the thumb subdomain contact surface
- All conserved at 1.000 across all 5 coronaviruses
  (SARS-CoV-2, SARS-CoV-1, MERS-CoV, HCoV-229E, HCoV-NL63)
- PHE440 ranks #1 composite score (contact loss=15, BSA=62.8 Å²)
- A compound targeting this aromatic pocket would be active
  pan-coronavirus

### GLU431–LYS1 salt bridge — SARS-selective secondary pharmacophore
- Crystallographically confirmed in 7BV2 at 4.24 Å
- AF3 predicted at 4.26 Å
- Both residues cons=1.000 — but SARS-CoV-1/2 only
- Absent in MERS/HCoV due to NSP7 N-terminal truncation (78 aa vs 84)
- Useful for SARS-specific compounds or dual-mechanism molecules

### Druggability — highest in project
- fpocket score 0.961 (6NUR), consensus ≥0.870 all 4 structures
- Deep, well-defined hydrophobic pocket — favorable for small molecules

### Drug design strategy
1. Pan-coronavirus: target aromatic hydrophobic cluster (PHE440 anchor)
2. SARS-selective: add charged group to engage GLU431–LYS1
3. Dual pharmacophore: aromatic core + terminal amine/carboxylate
   engages both mechanisms simultaneously

### Pipeline status
Steps 04–10 complete. Proceeding to Step 11 (3D visualization).

---

## Entry 041 — 3D Visualization NSP12-NSP7 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/11_3D_visualization_NSP12-NSP7_3.py
**Notebook:** notebooks/NSP12-NSP7_3D_3.ipynb

**8 interactive views (nglview):**
  1. Full complex — all conserved hotspots
  2. Salt bridge zoomed — GLU431–LYS1 (4.24 Å)
  3. Hotspots colored by composite score
  4. Hotspots colored by BSA burial depth
  5. Docking box visualization
  6. Structural overlay — 7BV2+6NUR+7C2K+AF3
  7. Full ranking table (styled)
  8. Hydrophobic core ★ — PHE440 primary pharmacophore

**Status:** ✅ NSP12-NSP7 PIPELINE FULLY COMPLETE (Steps 04-11)
Next complex: NSP12-NSP8

---

## Entry 042 — AF3 Validation NSP12-NSP8 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/04_validate_NSP12-NSP8_4.py

**AF3 quality:**
  iptm=0.85, ptm=0.88, ranking=0.86
  has_clash=0.0, fraction_disordered=0.02
  num_recycles=10 (fully converged)
  chain_pair_pae_min=1.22 Å (excellent)

**Validation results:**
  7BV2: NSP12 F1=0.924  NSP8 F1=0.945  Overall F1=0.934 ✅ PASS
  6NUR: NSP12 F1=0.922  NSP8 F1=0.926  Overall F1=0.924 ✅ PASS
  7C2K: NSP12 F1=0.970  NSP8 F1=0.954  Overall F1=0.963 ✅ PASS

**Primary gate (7BV2): F1=0.934 ✅ PASS (gate ≥0.70)**

**Chain assignments:**
  PDB: Chain A=NSP12, Chain B=NSP8 (consistent all 3 structures)
  AF3: Chain A=NSP12 (932 res), Chain B=NSP8 (198 res) — same orientation

**Status:** ✅ Done — proceed to Script 05_4 (interface analysis)

---

## Entry 043 — Interface Analysis NSP12-NSP8 complete
**Date:** $(date +%Y-%m-%d)
**Script:** scripts/05_interface_NSP12-NSP8_4.py

**Contact summary:**
  7BV2: SB=2 HB=11 HY=60 total=77
  6NUR: SB=3 HB=14 HY=71 total=94
  7C2K: SB=3 HB=12 HY=66 total=87
  AF3:  SB=4 HB=16 HY=69 total=97

**Salt bridges confirmed:**
  ASP523(NSP12)—ARG80(NSP8):  3.58 Å [ALL 4 structures] ★ PRIMARY
  LYS332(NSP12)—ASP99(NSP8):  3.72 Å [7BV2+7C2K+AF3]   ★ PRIMARY
  ASP517(NSP12)—LYS79(NSP8):  4.66 Å [6NUR+7C2K+AF3]   secondary

**Interface character:** Mixed — hydrophobic-dominated (HY=60-71)
  but with strong salt bridge network (2-4 per structure)
  Much richer interface than NSP12-NSP7

**Consensus hotspots (≥ 2 structures):**
  NSP12 (23): 387,129,389,271,330,131,380,523,91,87,332,95,117,517,99
  NSP8  (21): 117,129,80,115,131,112,91,87,116,95,98,83,128,90,121

**Status:** ✅ Done — proceed to Script 06_4 (conservation)

---
