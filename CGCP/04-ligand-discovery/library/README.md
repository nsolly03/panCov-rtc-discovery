# Screening Library

This directory contains the compound library for Phase 4 virtual screening.
Large files are excluded from Git via .gitignore.

## Enamine REAL HAC 22-23 (1.1B compounds)

**File:** 2025.02_Enamine_REAL_HAC_22_23_1.1B_CXSMILES.cxsmiles.bz2  
**Size:** 24 GB compressed  
**Format:** Tab-separated CXSMILES with pre-computed properties  
**Source:** https://enamine.net/compound-collections/real-compounds  

### Key columns
- smiles, id, MW, HAC, sLogP, HBA, HBD, RotBonds, QED, TPSA
- PPI_modulators — Enamine pre-flagged PPI-suitable compounds
- InChIKey — for deduplication

### To obtain this file
Download from Enamine website or request from project PI.
Place in: CGCP/04-ligand-discovery/library/enamine-real/
