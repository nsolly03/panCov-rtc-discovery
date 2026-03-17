#!/bin/bash
# Phase 0.1: BCL-2/BAX Positive Control Setup
# Usage: bash setup_bcl2.sh

set -e  # Exit on error

echo "=========================================="
echo "Phase 0.1: BCL-2/BAX Positive Control"
echo "=========================================="

# Define working directory
WORK_DIR="$HOME/projects/rtc-pan-coronavirus/CGCP/00-ground-truth/positive-controls/bcl2-bax"
cd "$WORK_DIR"

echo ""
echo "Working directory: $WORK_DIR"
echo ""

# Define structures to download
declare -A STRUCTURES=(
    ["2YXJ"]="BCL-XL/BIM peptide (2.0A) - Primary validation"
    ["4LVT"]="BCL-2/Venetoclax analog (2.2A) - Inhibitor-bound"
    ["6O0F"]="MCL-1/S63845 (2.1A) - Related system"
)

echo "Step 1: Downloading PDB structures..."
echo ""

for PDB in "${!STRUCTURES[@]}"; do
    DESC="${STRUCTURES[$PDB]}"
    echo "  Downloading $PDB: $DESC"
    
    if curl -s -f -o "${PDB}.pdb" "https://files.rcsb.org/download/${PDB}.pdb"; then
        echo "    ✓ Success"
    else
        echo "    ✗ Failed to download $PDB"
        exit 1
    fi
done

echo ""
echo "Step 2: Validating downloads..."
echo ""

for PDB in "${!STRUCTURES[@]}"; do
    if [ -f "${PDB}.pdb" ]; then
        SIZE=$(du -h "${PDB}.pdb" | cut -f1)
        
        # Extract resolution from PDB header
        RESOLUTION=$(grep "REMARK   2 RESOLUTION" "${PDB}.pdb" | awk '{print $4}' || echo "unknown")
        
        # Count atoms
        ATOMS=$(grep -c "^ATOM" "${PDB}.pdb" || echo "0")
        
        echo "  $PDB:"
        echo "    Size: $SIZE"
        echo "    Resolution: ${RESOLUTION}A"
        echo "    Atoms: $ATOMS"
        
        # Validation checks
        if [ "$RESOLUTION" != "unknown" ] && (( $(echo "$RESOLUTION <= 3.0" | bc -l) )); then
            echo "    ✓ Resolution acceptable"
        else
            echo "    ⚠ Resolution may be too high"
        fi
        
        if [ "$ATOMS" -gt 1000 ]; then
            echo "    ✓ Structure has sufficient atoms"
        else
            echo "    ⚠ Structure seems incomplete"
        fi
    else
        echo "  $PDB: ✗ File not found"
        exit 1
    fi
    echo ""
done

echo "Step 3: Creating validation checklist..."
echo ""

cat > validation_checklist.md << 'CHECKLIST'
# BCL-2/BAX Validation Checklist

## Structural Quality
- [ ] 2YXJ: Resolution ≤ 3.0A
- [ ] 4LVT: Resolution ≤ 3.0A  
- [ ] 6O0F: Resolution ≤ 3.0A
- [ ] All structures: No missing backbone atoms in interface region

## Interface Verification
- [ ] 2YXJ: Chains A and B physically adjacent (<5A)
- [ ] 2YXJ: PHE105 (BCL-XL) visible and solvent-exposed
- [ ] 2YXJ: Hydrophobic groove present at interface

## Control Validation
- [ ] PHE105 detected as hotspot by contact mapping
- [ ] Known inhibitor binding site identified
- [ ] Druggability score >0.5 for binding groove

## Manual Inspection
- [ ] Visual inspection in PyMOL completed
- [ ] No crystal packing artifacts observed
- [ ] B-factors reasonable in interface region

Status: ⏳ In Progress
Completed: 
CHECKLIST

echo "  Created: validation_checklist.md"

echo ""
echo "=========================================="
echo "Setup Complete"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Inspect structures: pymol 2YXJ.pdb"
echo "  2. Run audit: @path/to/audit_bcl2.pml"
echo "  3. Map contacts: python3 map_contacts_bcl2.py"
echo "  4. Update checklist: validation_checklist.md"
echo ""
echo "Files in $WORK_DIR:"
ls -lh *.pdb *.md 2>/dev/null || ls -lh
