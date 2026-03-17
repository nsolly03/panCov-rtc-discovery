# PyMOL Audit Script for BCL-2/BAX (2YXJ)
# Run: pymol 2YXJ.pdb, then @audit_bcl2.pml

# Setup
bg_color white
set ray_opaque_background, off
set orthoscopic, 1
set antialias, 2

# Load structure (assumes already loaded)
# load 2YXJ.pdb

# Basic representation
show cartoon
color gray, all

# Define chains
select bclxl, chain A
select bim, chain B

color marine, bclxl
color forest, bim

# Known hotspot: PHE105 (BCL-XL) - critical for binding
select phe105, resi 105 and chain A
show sticks, phe105
color red, phe105
label phe105 and name CA, "PHE105"

# Additional known hotspots
select val126, resi 126 and chain A
select phe146, resi 146 and chain A
show sticks, val126 phe146
color orange, val126 phe146

# Measure interface distance
distance interface_5A, bclxl, bim, 5.0

# Zoom to interface
zoom phe105, 15

# Save session
save bcl2_audit_session.pse

# Render for documentation
ray 2400, 1800
png step01_bcl2_overview.png, dpi=300, width=2400, height=1800

# Print audit notes
print("========================================")
print("BCL-2/BAX Audit Complete")
print("========================================")
print("")
print("Manual verification checklist:")
print("  [ ] Chains A (BCL-XL) and B (BIM) are adjacent")
print("  [ ] PHE105 is exposed and not buried")
print("  [ ] Hydrophobic groove visible (LEU/VAL/PHE cluster)")
print("  [ ] No crystal packing contacts at interface")
print("  [ ] B-factors reasonable (no red in interface)")
print("")
print("Session saved: bcl2_audit_session.pse")
print("Image saved: step01_bcl2_overview.png")
