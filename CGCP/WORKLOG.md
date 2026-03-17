# WORKLOG — CGCP: Conservation-Guided Cluster Pharmacophore

**Project:** Systematic Discovery of Pan-Coronavirus RTC Inhibitors
**Candidate:** Olivier Nsekuye
**Institution:** University of Liège — GIGA-VIN Lab
**Supervisor:** Prof. Jean-Claude Twizere
**Funding:** FRIA-B1 Fellowship
**Start Date:** 2026-03-17
**Working Directory:** ~/projects/rtc-pan-coronavirus/CGCP/

---

## Project Background

### The Problem

Coronaviruses continue to pose pandemic threats. SARS-CoV-2 variants evade existing vaccines and therapeutics. Broad-spectrum antivirals targeting conserved viral machinery are urgently needed.

### Why Protein-Protein Interfaces?

| Advantage | Rationale |
|-----------|-----------|
| Conserved across species | Interface residues under functional constraint |
| Lower resistance risk | Multiple contact points required for function |
| Novel mechanism | Distinct from existing polymerase inhibitors |
| Pan-coronavirus potential | Target conserved interfaces, not variable active sites |

### Nine Target Interfaces

| # | Complex | Function | Status |
|---|---------|----------|--------|
| 1 | NSP10-NSP14 | Exoribonuclease regulation | Pending |
| 2 | NSP10-NSP16 | Methyltransferase activation | Pending |
| 3 | NSP12-NSP7 | RdRp cofactor binding | Pending |
| 4 | NSP12-NSP8 | RdRp processivity | Pending |
| 5 | NSP9-NSP12 | NiRAN domain interaction | Pending |
| 6 | NSP7-NSP8 | Cofactor assembly | Pending |
| 7 | NSP13-Helicase | Helicase dimerization | Pending |
| 8 | NSP12-NSP13 | Polymerase-helicase junction | Pending |
| 9 | NSP15 | Endoribonuclease dimerization | Pending |

---

## Entry 001 — CGCP Project Initialization

**Date:** 2026-03-17
**Phase:** 0 — Ground Truth Calibration
**Step:** Project setup

**What:** Created complete CGCP directory structure

**Status:** Complete

**Next:** Phase 0.1 — BCL-2/BAX positive control setup

## Entry 003 — Phase 0.1 Complete: BCL-2/BAX Contact Mapping Passed

**Date:** 2026-03-17
**Phase:** 0 — Ground Truth Calibration
**Step:** 0.1 — BCL-2/BAX contact mapping validation

**What:** Ran contact mapping pipeline on BCL-XL/BIM complex. All known hotspots detected.

**Issue encountered and resolved:**
- 2YXJ was incorrectly a BCL-XL homodimer — no BIM peptide present
- Replaced with 3FDL (BCL-XL chain A, 158 res + BIM chain B, 26 res)
- Diagnostic worked correctly: pipeline detected symmetric contacts as warning

**Results:**
- Total contacts: 65
- Unique chain A residues at interface: 27
- Unique chain B residues at interface: 20
- PHE105: 2 contacts (LEU94, ILE97) ✓
- VAL126: 4 contacts (ALA91, GLU87, LEU94) ✓
- PHE146: 2 contacts (LEU94, ILE90) ✓

**Validation:** PASSED — pipeline correctly identifies known hotspots at 5.0A cutoff

**Status:** ✅ Complete
**Next:** Phase 0.2 — MDM2/p53 positive control
