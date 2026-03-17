# Phase 2 Step 1 — Structural Audit: NSP12-NSP7

**Date:** 2026-03-17  
**PDB file:** `receptor_NSP12-NSP7_3.pdb`  
**Chain A:** NSP12 (834 residues)  
**Chain C:** NSP7  (63 residues)  
**Overall result:** PASS  

## Checklist

| Check | Status | Detail |
|-------|--------|--------|
| C1: Both chains present | ✅ PASS | Chains found: ['A', 'C'] |
| C2: Chain sizes correct | ✅ PASS | Chain A: 834 res (expected ~834)  Chain C: 63 res (expected ~63) |
| C3: Chains adjacent (<10 Å) | ✅ PASS | Min inter-chain distance (sample): 3.00 Å |
| C4: Backbone complete | ✅ PASS | Chain A missing: 0 residues  Chain C missing: 0 residues |
| C5: PHE440 present and correct | ✅ PASS | Present: True  Identity: PHE (expected PHE) |
| C6 [+ctrl]: PHE440 near NSP7 (<15.0 Å) | ✅ PASS | PHE440 to NSP7 min dist: 3.69 Å |
| C7 [-ctrl]: GLY200 far from NSP7 (>25.0 Å) | ✅ PASS | GLY200 to NSP7 min dist: 66.97 Å |

## Controls

**Positive control:** PHE440 (primary anchor) to NSP7 = 3.690000057220459 Å  
Expected: < 15.0 Å (at interface)  

**Negative control:** RES200 (surface) to NSP7 = 66.97000122070312 Å  
Expected: > 25.0 Å (away from interface)  

## Decision

Structure **accepted** for Phase 2 Step 2 (contact mapping).
