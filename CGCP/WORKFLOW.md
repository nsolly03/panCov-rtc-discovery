# CGCP End-to-End Workflow
**Project:** Conservation-Guided Cluster Pharmacophore for Pan-Coronavirus RTC Inhibitor Discovery  
**Version:** 1.0  
**Date:** 2026-03-17  
**Working Directory:** `~/projects/rtc-pan-coronavirus/CGCP/`

---

## 1. Project Rationale

### 1.1 The Problem

Coronaviruses continue to pose pandemic threats. SARS-CoV-2 variants evade existing vaccines and therapeutics. Broad-spectrum antivirals targeting conserved viral machinery are urgently needed.

The **Replication-Transcription Complex (RTC)** is the viral machinery responsible for genome replication and mRNA synthesis. It comprises multiple non-structural proteins (NSPs) that assemble through protein-protein interactions (PPIs).

### 1.2 Why Target Protein-Protein Interfaces?

| Advantage | Rationale |
|-----------|-----------|
| Conserved across species | Interface residues under functional constraint |
| Lower resistance risk | Multiple contact points required for function |
| Novel mechanism | Distinct from existing polymerase inhibitors |
| Pan-coronavirus potential | Target conserved interfaces, not variable active sites |

### 1.3 The Challenge

PPI interfaces are traditionally "undruggable":
- Flat, extended surfaces (not deep pockets)
- Large interaction areas (1500-3000 Å²)
- Dynamic, conformationally variable
- Low fpocket druggability scores (often 0.000)

### 1.4 Our Solution: CGCP

**Conservation-Guided Cluster Pharmacophore** combines:
1. **Targeting residue clusters, not single hotspots** — Small molecules bridge 3-6 complementary features
2. **Using conservation as primary filter** — Pan-coronavirus activity built into design
3. **Defining abstract pharmacophores** — Not chasing pre-formed pockets
4. **Constraint-based screening** — Dock to clusters, not entire interfaces

---

## 2. Core Concepts

### 2.1 Definitions

| Term | Definition | Example |
|------|-----------|---------|
| **Anchor** | Single conserved, buried residue with high energetic contribution. It serves as the geometric and energetic center for cluster design. | A PHE or LYS residue conserved across all 5 coronaviruses, with high BSA |
| **Cluster** | Spatial grouping of 3-6 residues within 8-12 Å of each other, containing complementary chemical features that a small molecule can simultaneously engage | Aromatic anchor + charged residue + hydrophobic patch within 10 Å |
| **Pharmacophore** | The abstract geometric arrangement of chemical features (aromatic ring, positive charge, H-bond donor, etc.) required for biological activity. Not the residues themselves, but their spatial relationships | "Aromatic ring 5 Å from H-bond donor, 7 Å from hydrophobic patch" |
| **Conservation-guided** | Using evolutionary conservation across coronaviruses to prioritize features that will yield broad-spectrum activity | Selecting clusters where ≥3 residues have cons ≥0.8 across SARS-CoV-2/1, MERS, 229E, NL63 |
| **Induced fit** | Conformational change in the protein upon ligand binding that creates a binding site not present in the apo structure | Flat interface becomes druggable pocket when fragment binds |
| **Cryptic pocket** | A pocket invisible in apo structures but revealed by molecular dynamics or ligand binding | NSP7-NSP8 Mode B: N-terminal helix creates pocket not seen in crystal structures |

### 2.2 Why Not Other Approaches?

| Approach | Why Rejected | CGCP Alternative |
|----------|-------------|------------------|
| Single hotspot targeting | Small molecule binds one residue — insufficient affinity/specificity | Cluster with multiple anchors |
| Pure pocket detection (fpocket) | PPI interfaces often score 0.000 — no pocket to find | Conservation prioritizes clusters regardless of pocket score |
| Allosteric hunting | Requires unknown sites, months of MD, not pan-coronavirus | Interface is known, measurable, targetable |
| Traditional docking to whole interface | Too large, non-specific binding, no pharmacophore constraint | Cluster defines precise search space |
| Fragment-based without guidance | Random exploration, no pan-coronavirus guarantee | Conservation ensures broad-spectrum potential from design stage |

### 2.3 The Nine Target Interfaces

| # | Complex | Function | Druggability | Status |
|---|---------|----------|-------------|--------|
| 1 | NSP10-NSP14 | Exoribonuclease regulation | 0.000 | Pending |
| 2 | NSP10-NSP16 | Methyltransferase activation | 0.546 | Pending |
| 3 | NSP12-NSP7 | RdRp cofactor binding | **0.961** | Pending |
| 4 | NSP12-NSP8 | RdRp processivity | **0.874** | Pending |
| 5 | NSP9-NSP12 | NiRAN domain interaction | **0.895** | Pending |
| 6 | NSP7-NSP8 | Cofactor assembly | 0.531 | Pending |
| 7 | NSP13-Helicase | Helicase dimerization | 0.001 | Pending |
| 8 | NSP12-NSP13 | Polymerase-helicase junction | 0.000 | Pending |
| 9 | NSP15 | Endoribonuclease dimerization | 0.900 | Pending |

**Note:** Druggability scores above are preliminary and must be independently validated per interface in Phase 2. They inform but do not determine tier classification.

**Tier 1 (>0.80):** NSP12-NSP7, NSP9-NSP12, NSP12-NSP8, NSP15  
**Tier 2 (0.40-0.80):** NSP10-NSP16, NSP7-NSP8  
**Tier 3 (<0.10):** NSP10-NSP14, NSP13-Helicase, NSP12-NSP13

---

## 3. Workflow Overview

```
PHASE 0: Ground Truth Calibration (2-3 days)
│   Validate methods on known, published PPI inhibitor systems
│   Establish pass/fail thresholds before touching your 9 interfaces
│
│   0.1 Positive Control: BCL-2/BAX (Venetoclax — FDA-approved)
│   0.2 Positive Control: MDM2/p53 (Nutlins — classic hotspot)
│   0.3 Negative Control: Undruggable PPI (confirm rejection)
│
PHASE 1: Interface Selection (1 day)
│   Blind, objective evaluation of all 9 interfaces
│   No reference to any previous pipeline results
│
│   1.1 Structural quality audit (all 9)
│   1.2 Tier classification S/T/D
│   1.3 Selection of first interface to analyze
│
PHASE 2: Deep Dive — Single Interface (2-3 weeks per interface)
│   Complete all 8 steps before moving to next interface
│
│   Step 1: Structural Verification
│   Step 2: Residue-Level Contact Mapping (blind, raw)
│   Step 3: Feature Classification (no scoring yet)
│   Step 4: Spatial Clustering (geometry only)
│   Step 5: Conservation Overlay (independent data layer)
│   Step 6: Integrated Assessment (your manual decision)
│   Step 7: Pharmacophore Hypothesis (abstract, no molecule)
│   Step 8: Control Experiments (validate hypothesis)
│
PHASE 3: Validation and Expansion (1 week)
│   Cross-interface analysis after ≥2 interfaces complete
│
PHASE 4: Ligand Discovery (2-3 weeks per interface)
│   Step 9:  Library preparation and filtering
│   Step 10: Constraint docking setup
│   Step 11: Local test screen (100-500 compounds)
│   Step 12: HPC production screen (NIC5)
│
PHASE 5: Hit Analysis (1-2 weeks per interface)
│   Step 13: Hit classification
│   Step 14: Multi-target analysis
│   Step 15: Visual verification and pose inspection
│
PHASE 6: Documentation (ongoing)
    Daily WORKLOG entries
    Weekly summaries
    Manuscript drafting
    Reproducibility bundle
```

**Core philosophy:** One interface at a time. Complete all steps before moving on.

---

## 4. Phase 0: Ground Truth Calibration

### Purpose

Before running our pipeline on unknown interfaces, we validate it on systems where the answers are already known. If we cannot recover published results on BCL-2 and MDM2, we cannot trust results on new interfaces.

### 4.1 Positive Control 1: BCL-2/BAX — Venetoclax System

**Why this system:**
- Venetoclax is FDA-approved — ground truth exists
- PHE105 is a well-characterized anchor residue
- Hydrophobic groove is clearly documented
- Multiple high-resolution structures available

**Structures to use:**

| PDB | System | Resolution | Purpose |
|-----|--------|------------|---------|
| 2YXJ | BCL-XL / BIM peptide | 2.2 Å | Primary validation |
| 4LVT | BCL-2 / Venetoclax analog | 2.05 Å | Inhibitor-bound |
| 6O0F | MCL-1 / S63845 | 2.12 Å | Related system |

**Pipeline success criteria (must ALL pass):**

| Criterion | Expected | If Failed |
|-----------|----------|-----------|
| PHE105 in top 5 interface residues | YES | Check chain assignment |
| Hydrophobic cluster identified | YES | Adjust distance cutoff |
| Druggability score >0.5 | YES | Re-run fpocket |
| Known inhibitor matches cluster | YES | Revise cluster definition |

**Scripts:**
- `scripts/phase-0/setup_bcl2.sh` — Download and validate structures
- `scripts/phase-0/audit_bcl2.pml` — PyMOL visual inspection
- `scripts/phase-0/map_contacts_bcl2.py` — Contact mapping with hotspot validation

### 4.2 Positive Control 2: MDM2/p53 — Nutlin System

**Why this system:**
- Nutlins are classic "hotspot" PPI inhibitors
- Three p53 residues (TRP23, PHE19, LEU26) are the known anchors
- Confirms pipeline works on a different interface type

**Structures to use:**

| PDB | System | Resolution | Purpose |
|-----|--------|------------|---------|
| 1YCR | MDM2 / p53 peptide | 2.6 Å | Primary validation |
| 4HG7 | MDM2 / Nutlin-3a | 2.2 Å | Inhibitor-bound |

**Success criteria:**
- TRP23, PHE19, LEU26 identified as top anchors
- MDM2 hydrophobic cleft detected as cluster
- Pharmacophore matches three-point nutlin binding model

**Scripts:**
- `scripts/phase-0/setup_mdm2.sh`
- `scripts/phase-0/map_contacts_mdm2.py`

### 4.3 Negative Control: Undruggable PPI

**Purpose:** Confirm that our pipeline correctly identifies when a PPI is NOT a good drug target.

**Success criteria:**
- Druggability score <0.1
- No compact clusters with diverse features
- Contact map shows extended, featureless interface

---

## 5. Phase 1: Interface Selection

### 5.1 Evaluation Criteria (Applied Blindly)

Evaluate all 9 interfaces against objective criteria **before** any detailed analysis and **without reference to prior pipeline results**.

| Criterion | How Measured | Threshold | Points |
|-----------|-------------|-----------|--------|
| Structural resolution | PDB REMARK 2 header | ≤3.0 Å | 1 |
| Number of PDB structures available | PDB search | ≥2 | 1 |
| Interface size | FreeSASA or NACCESS | BSA ≥500 Å² | 1 |
| Sequence coverage | Resolved residues / full length | ≥80% | 1 |
| Species for MSA | PDB taxonomy / UniProt | ≥3 CoV species | 1 |
| Preliminary druggability | fpocket score | ≥0.3 | 1 |

### 5.2 Tier Classification

| Tier | Score | Decision |
|------|-------|----------|
| **S** | ≥4/6 | Full Phase 2 deep dive |
| **T** | 2-3/6 | Abbreviated analysis |
| **D** | 0-1/6 | Deprioritize |

### 5.3 Interface Selection Rule

- Start with randomly selected Tier S interface
- Complete Phase 2 in full before selecting next
- If no Tier S available, escalate Tier T

---

## 6. Phase 2: Deep Dive — Single Interface

All 8 steps apply identically to every interface. No skipping.

### Step 1: Structural Verification

**Goal:** Confirm the structure is real, complete, and correctly assigned before any analysis.

**Manual audit (PyMOL):**
- [ ] Both chains visible and adjacent (< 10 Å across interface)
- [ ] No crystal packing contacts masquerading as interface
- [ ] Missing residues < 20% of interface region
- [ ] B-factors in interface region are reasonable (not red)
- [ ] Resolution ≤ 3.0 Å

**Positive control check:** Run same audit on BCL-2 (2YXJ) — should pass all.  
**Negative control check:** Apply to a non-adjacent chain pair — should fail adjacency check.

**Output files:**
```
02-deep-dive/[INTERFACE]/step-01-structure/
├── audit-report.md
├── [INTERFACE]_audit.pse       (PyMOL session)
└── visuals/
    └── step01_overview.png     (300 DPI, ray-traced)
```

### Step 2: Residue-Level Contact Mapping

**Goal:** Generate the complete, unfiltered list of residues in contact across the interface.

**Method:**
- Distance cutoff: 5.0 Å between any heavy atoms
- Include both chains
- Exclude: HETATM, water molecules
- NO filtering by conservation, BSA, or any score yet

**Positive control:** Compare output with PDBsum or PISA for same structure.  
**Negative control:** Apply cutoff to non-adjacent chains — should return zero contacts.

**Output files:**
```
02-deep-dive/[INTERFACE]/step-02-contacts/
├── raw-contact-map.tsv         (all contacts, no filtering)
├── contact-map-summary.md
└── visuals/
    ├── step02_2D_contact_matrix.png
    └── step02_3D_contacts.png
```

### Step 3: Feature Classification

**Goal:** Assign chemical features to every interface residue. No scoring.

| Feature Type | Residues | Chemical Rationale |
|-------------|----------|--------------------|
| Aromatic | PHE, TYR, TRP | π-π stacking, CH-π interactions |
| Positive charge | LYS, ARG | Salt bridges, H-bonds |
| Negative charge | ASP, GLU | Salt bridges, H-bonds |
| H-bond donor | SER, THR, TYR, ASN, GLN, LYS, ARG, HIS, TRP | Directional interactions |
| H-bond acceptor | SER, THR, TYR, ASN, GLN, ASP, GLU, HIS | Directional interactions |
| Hydrophobic | ALA, VAL, LEU, ILE, PRO, MET | Desolvation, vdW |
| Sulfur | CYS, MET | Metal coordination, H-bonds |
| Special | HIS (dual), CYS (metal) | Context-dependent |

**Output files:**
```
02-deep-dive/[INTERFACE]/step-03-features/
├── feature-classification.tsv
└── visuals/
    └── step03_feature_map.png    (residues as spheres, colored by type)
```

### Step 4: Spatial Clustering

**Goal:** Group interface residues into geometrically compact clusters using coordinates alone — no chemical or conservation weighting at this stage.

**Method:**
- Input: Cα coordinates of all interface residues
- Algorithm: DBSCAN (ε=12 Å, min_samples=3)
- Cluster quality = number of members × feature diversity

**Positive control:** BCL-2 hydrophobic groove should appear as a single compact cluster.  
**Negative control:** Shuffled coordinates → no clusters (random placement).

**Output files:**
```
02-deep-dive/[INTERFACE]/step-04-clusters/
├── cluster-assignments.tsv
├── cluster-summary.md
└── visuals/
    ├── step04_clusters_3D.png
    └── step04_cluster_radii.png
```

### Step 5: Conservation Overlay

**Goal:** Bring conservation data in as an independent layer — do not modify cluster boundaries based on conservation at this stage.

**Coronaviruses analyzed:**
- SARS-CoV-2 (reference)
- SARS-CoV-1
- MERS-CoV
- HCoV-229E
- HCoV-NL63

**Conservation score:** 0.0 (fully variable) → 1.0 (identical in all 5 species)

**Conservation classification:**

| Score | Class | Drug Design Implication |
|-------|-------|------------------------|
| 1.000 | Identical | Pan-coronavirus target |
| 0.800-0.999 | Conserved | Broad-spectrum potential |
| 0.500-0.799 | Semi-conserved | Selective activity possible |
| <0.500 | Variable | SARS-CoV-2 specific only |

**Positive control:** Known pan-coronavirus residues must score ≥0.8.  
**Negative control:** Known variable positions (e.g., receptor binding domain) must score <0.5.

**Output files:**
```
02-deep-dive/[INTERFACE]/step-05-conservation/
├── msa/
│   ├── [PROTEIN]_alignment.fasta
│   └── conservation-scores.tsv
├── conservation-overlay.tsv       (cluster members + conservation)
└── visuals/
    ├── step05_conservation_heatmap.png
    └── step05_cluster_conservation.png
```

### Step 6: Integrated Assessment

**Goal:** YOU make a documented decision about which clusters are worth pursuing and why.

**Decision table (complete for every cluster):**

| Cluster | Members | Feature Diversity | Mean Conservation | BSA Contribution | Assessment | Rationale |
|---------|---------|-------------------|-------------------|-----------------|------------|-----------|
| 1 | 5 res | Aromatic+charged+HBD | 0.92 | High | **Pursue** | 3 conserved anchors, compact |
| 2 | 4 res | Hydrophobic only | 0.61 | Low | Skip | Insufficient diversity |
| 3 | 3 res | Charged only | 0.45 | Medium | Defer | Variable, single feature type |

**This is a manual step.** Write your rationale in complete sentences, not just numbers.

**Output files:**
```
02-deep-dive/[INTERFACE]/step-06-assessment/
├── assessment-report.md        (your written rationale)
└── summary-visuals/
    └── step06_integrated_view.png
```

### Step 7: Pharmacophore Hypothesis

**Goal:** Define what a ligand must look like in abstract terms — without picking a specific molecule.

**Pharmacophore definition format:**

```
Cluster [X] Pharmacophore:
  ANCHOR:   [Residue] — [interaction type]
  REQUIRED: [Feature 1] — [distance from anchor ± tolerance]
  REQUIRED: [Feature 2] — [distance from anchor ± tolerance]
  DESIRED:  [Feature 3] — [distance from anchor ± tolerance]
  TOLERATED: [Feature 4] — [region description]
  FORBIDDEN: [Charge type] — [reason: e.g. clashes with LYS20]

  Estimated ligand span: X–Y Å
  Estimated MW: 300–500 Da
  Expected rotatable bonds: ≤8
```

**Output files:**
```
02-deep-dive/[INTERFACE]/step-07-pharmacophore/
├── pharmacophore-definition.md
├── 2d-sketches/
│   └── pharmacophore_sketch.png   (hand-drawn or ChemDraw)
└── 3d-definitions/
    └── pharmacophore_3D.pse       (PyMOL spheres with tolerances)
```

### Step 8: Control Experiments

**Goal:** Validate the pharmacophore hypothesis using controls before any real screening.

| Control | Method | Expected Result | If Failed |
|---------|--------|-----------------|-----------|
| Published inhibitor | Dock known PPI inhibitor to cluster | Matches pharmacophore features | Re-examine cluster definition |
| Benzene (inert) | Dock benzene to cluster | Poor score, no feature match | Check scoring function |
| Shifted pharmacophore | Move center 20 Å away | No binding, poor scores | Validate search space setup |
| Active site residues (if distinct) | Add catalytic triad to analysis | Must NOT appear in cluster | Confirms clean interface/active site separation |

**Output files:**
```
02-deep-dive/[INTERFACE]/step-08-controls/
├── control-results.md
├── positive-control-results/
└── negative-control-results/
```

---

## 7. Phase 3: Validation and Expansion

After completing Phase 2 for at least 2 interfaces:

- **Cross-interface conservation:** Residues shared across interfaces (e.g., NSP12 appears in 4 complexes) — verify consistent classification
- **Cluster overlap analysis:** Do any clusters from different interfaces converge on the same pharmacophore type?
- **Druggability reassessment:** Does cluster-based druggability correlate with fpocket score?

---

## 8. Phase 4: Ligand Discovery

### Step 9: Library Preparation

**Sources:**
- ZINC20 (primary)
- Enamine REAL (if HPC download succeeds)

**Filters applied (in order):**

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| Molecular weight | 200-500 Da | Drug-like range |
| LogP | -1 to 5 | Membrane permeability |
| H-bond donors | ≤5 | Lipinski |
| H-bond acceptors | ≤10 | Lipinski |
| Rotatable bonds | ≤8 | Conformational entropy |
| PAINS filter | Remove flagged | Eliminate frequent hitters |
| Pharmacophore pre-filter | Must match anchor feature | Reduce search space |

### Step 10: Constraint Docking Setup

**Tool:** AutoDock Vina

**Docking box:** 15 Å cube centered on cluster anchor (not whole interface)

**Pharmacophore constraints:**
- Molecule must interact with anchor residue
- Must engage ≥2 additional cluster features
- Scored by cluster complementarity, not just docking score

### Step 11: Local Test Screen

**Purpose:** Validate setup before committing to HPC.

- 100-500 compounds from filtered library
- Verify top hits engage expected cluster features
- Check for artifacts (too-high scores, clashing poses)

### Step 12: HPC Production Screen (NIC5)

```bash
# Transfer to NIC5
scp -r 04-ligand-discovery/[INTERFACE]/ nic5:~/cgcp-screening/

# Submit SLURM array
sbatch scripts/phase-4/screen_array.sh

# Retrieve results
scp -r nic5:~/cgcp-screening/results/ 04-ligand-discovery/[INTERFACE]/results/
```

---

## 9. Phase 5: Hit Analysis

### Step 13: Hit Classification

| Class | Criteria | Action |
|-------|----------|--------|
| **A** | Score ≤-7.0 + pose engages anchor + ≥2 cluster features | Priority for Step 14 |
| **B** | Score ≤-6.0 + good pose | Secondary analysis |
| **C** | Borderline score or single-feature engagement | Document, deprioritize |
| **Rejected** | Poor pose, clashing, PAINS, artifacts | Remove |

### Step 14: Multi-Target Analysis

For each Class A hit:
- Re-dock to all Tier 1 interfaces
- Compare poses across interfaces
- Flag hits active in ≥2 interfaces as **multi-target candidates**

### Step 15: Visual Verification and Pose Inspection

Every Class A hit must be manually inspected in PyMOL:
- Verify pose is chemically reasonable
- Confirm no steric clashes
- Confirm anchor interaction is genuine
- Generate publication-quality figure for each confirmed hit

---

## 10. Phase 6: Documentation

### Daily WORKLOG Entry Format

```
## Entry XXX — [Brief Description]

**Date:** YYYY-MM-DD
**Phase:** [0-6]
**Step:** [Step name]
**Interface:** [If applicable]

**What:** [What you did in 1-2 sentences]

**Commands:**
[paste key commands]

**Observations:**
- [Observation 1]
- [Observation 2]

**Decisions:**
- [Decision with explicit rationale]

**Controls:**
- Positive: [Result]
- Negative: [Result]

**Status:** [✅ Complete / ⏳ In Progress / ❌ Failed]
**Next:** [Next action]
```

---

## 11. Visual Standards

All figures must meet these standards before inclusion in any report or manuscript.

### Technical Requirements

| Element | Standard | Command |
|---------|----------|---------|
| Resolution | 300 DPI minimum | `ray 2400, 1800` |
| Output size | 2400 × 1800 pixels | `png fig.png, dpi=300` |
| Background | White, transparent option | `bg_color white` |
| Camera | Orthographic (no perspective distortion) | `set orthoscopic, 1` |
| Anti-aliasing | On | `set antialias, 2` |
| Shadows | Off | `set ray_shadows, 0` |
| Fog | Off | `set ray_trace_fog, 0` |

### PyMOL Render Template

```python
# Standard render settings — apply before every ray trace
bg_color white
set ray_opaque_background, off
set orthoscopic, 1
set antialias, 2
set ray_shadows, 0
set ray_trace_fog, 0
set hash_max, 300

# Render
ray 2400, 1800
png output_figure.png, dpi=300
```

### Color Conventions (Colorblind-Safe)

| Element | Color | Hex |
|---------|-------|-----|
| Protein chain A | Steel blue | #4682B4 |
| Protein chain B | Forest green | #228B22 |
| Anchor residue | Crimson | #DC143C |
| Cluster members | Orange | #FF8C00 |
| Conserved (1.0) | Dark green | #006400 |
| Variable (<0.5) | Dark red | #8B0000 |
| H-bond | Yellow dashed | — |
| Hydrophobic | Gray | #808080 |

### Figure Panels per Interface

| Figure | Content |
|--------|---------|
| Panel A | Full complex overview, both chains |
| Panel B | Interface zoom, 20 Å sphere |
| Panel C | Cluster view, colored by feature type |
| Panel D | Conservation overlay (green→red gradient) |
| Panel E | Pharmacophore hypothesis (3D spheres + tolerances) |
| Panel F | Top docking hit pose (after Phase 4) |

---

## 12. Key Decisions and Rationale

This section records the major methodological decisions made during project design.

### Why CGCP over allosteric site hunting

Allosteric sites require: (1) unknown target identification, (2) long MD simulations (≥100 ns × 3 replicates), (3) extensive biophysical validation, (4) no guarantee of pan-coronavirus conservation. CGCP targets known, measurable interfaces with established conservation data.

### Why clusters over single hotspots

A single hotspot gives insufficient affinity and specificity for a drug-like molecule. A cluster of 3-6 complementary residues within 10 Å provides multiple simultaneous contacts — mimicking how successful PPI inhibitors like Venetoclax actually bind.

### Why interface-by-interface

Prevents cherry-picking. Ensures every interface gets equal rigor. Builds methodological expertise progressively. Allows cross-interface comparison only after unbiased individual analysis.

### Why ground truth first

If BCL-2/MDM2 controls fail, no downstream result is trustworthy. Phase 0 is not optional.

---

## 13. Repository Structure

```
~/projects/rtc-pan-coronavirus/CGCP/
├── WORKFLOW.md                        ← This document
├── WORKLOG.md                         ← Daily execution log
├── README.md                          ← Quick start guide
│
├── 00-ground-truth/
│   ├── positive-controls/
│   │   ├── bcl2-bax/
│   │   └── mdm2-p53/
│   └── negative-controls/
│
├── 01-interface-selection/
│   ├── raw-evaluation-data/
│   ├── tier-classification/
│   └── selected-interface/
│
├── 02-deep-dive/
│   ├── TEMPLATE/                      ← Copy for each interface
│   │   ├── step-01-structure/
│   │   ├── step-02-contacts/
│   │   ├── step-03-features/
│   │   ├── step-04-clusters/
│   │   ├── step-05-conservation/
│   │   ├── step-06-assessment/
│   │   ├── step-07-pharmacophore/
│   │   └── step-08-controls/
│   └── [INTERFACE-NAME]/              ← One folder per interface
│
├── 03-validation/
│   └── cross-interface-analysis/
│
├── 04-ligand-discovery/
│   └── TEMPLATE/
│       ├── step-09-library/
│       ├── step-10-constraint-docking/
│       ├── step-11-local-test/
│       └── step-12-hpc-screen/
│
├── 05-hit-analysis/
│   └── TEMPLATE/
│       ├── step-13-classification/
│       ├── step-14-multi-target/
│       └── step-15-visual-verification/
│
├── 06-documentation/
│   ├── manuscript/
│   ├── supplementary/
│   └── reproducibility-bundle/
│
├── scripts/
│   ├── phase-0/
│   ├── phase-1/
│   ├── phase-2/
│   ├── phase-3/
│   ├── phase-4/
│   └── phase-5/
│
├── notebooks/
│   ├── manual-inspection/
│   └── visualization-templates/
│
└── logs/
    ├── 2026-03/
    ├── 2026-04/
    └── 2026-05/
```

---

*This workflow is a living document. Update version number and date when making structural changes.*  
*Last modified: 2026-03-17 | Version: 1.0*
