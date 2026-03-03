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
