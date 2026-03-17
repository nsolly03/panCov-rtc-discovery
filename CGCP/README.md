# CGCP: Conservation-Guided Cluster Pharmacophore

**Location:** ~/projects/rtc-pan-coronavirus/CGCP/
**Started:** 2026-03-17
**Status:** Phase 0 — Ground Truth Calibration

## Quick Start

Phase 0: Calibrate on known systems
    cd 00-ground-truth/positive-controls/bcl2-bax
    # Download structures
    curl -O https://files.rcsb.org/download/2YXJ.pdb

Phase 1: Evaluate 9 interfaces
    cd 01-interface-selection

Phase 2: Deep dive
    cd 02-deep-dive/[SELECTED-INTERFACE]
