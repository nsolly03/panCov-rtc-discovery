#!/bin/bash
# Setup BCL-2/BAX positive control
# run: bash setup_bcl2.sh

cd ~/projects/rtc-pan-coronavirus/CGCP/00-ground-truth/positive-controls/bcl2-bax

echo "=== Phase 0.1: BCL-2/BAX Positive Control Setup ==="

# Download structures
echo "Downloading PDB structures..."
curl -s -O https://files.rcsb.org/download/2YXJ.pdb
curl -s -O https://files.rcsb.org/download/4LVT.pdb
curl -s -O https://files.rcsb.org/download/6O0F.pdb

# Verify downloads
echo ""
echo "Downloaded files:"
ls -lh *.pdb

# Basic validation
echo ""
echo "Validating structures..."
for pdb in 2YXJ 4LVT 6O0F; do
    if grep -q "ATOM" ${pdb}.pdb; then
        resolution=$(grep "REMARK   2 RESOLUTION" ${pdb}.pdb | awk '{print $4}')
        echo "  ${pdb}: OK (Resolution: ${resolution}A)"
    else
        echo "  ${pdb}: FAILED"
    fi
done

echo ""
echo "Next steps:"
echo "  1. Inspect in PyMOL: pymol 2YXJ.pdb"
echo "  2. Run audit: @../../../scripts/phase-0/audit_bcl2.pml"
echo "  3. Map contacts: python3 ../../../scripts/phase-0/map_contacts_bcl2.py"
