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
