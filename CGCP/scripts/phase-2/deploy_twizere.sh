#!/bin/bash
# ================================================================
# One-shot deploy: sets up Prof. Twizere's CGCP VFUparr from
# your laptop. Run AFTER he has added your SSH key.
# Usage: bash deploy_twizere.sh
# ================================================================

set -e

TWIZERE_HOST="nic5-twizere"
SETUP_SCRIPT="$(dirname $0)/setup_cgcp_vfuparr_twizere.sh"

echo "============================================================"
echo "  CGCP Deploy → Prof. Twizere NIC5 account"
echo "============================================================"

# Step 1: Test SSH connection
echo "[1/3] Testing SSH connection to twizere account..."
ssh -o ConnectTimeout=10 ${TWIZERE_HOST} "echo '  ✅ SSH OK — logged in as: \$(whoami)'"

# Step 2: Transfer setup script
echo "[2/3] Transferring setup script..."
scp ${SETUP_SCRIPT} ${TWIZERE_HOST}:/scratch/ulg/gigambd/twizere/
echo "  ✅ Script transferred"

# Step 3: Execute setup script remotely
echo "[3/3] Running setup on twizere account..."
ssh ${TWIZERE_HOST} "bash /scratch/ulg/gigambd/twizere/setup_cgcp_vfuparr_twizere.sh"

echo ""
echo "============================================================"
echo "  Deploy complete — Twizere's 4 interfaces are ready"
echo "  Launch parallel monitoring:"
echo "  bash $(dirname $0)/cgcp_tmux_parallel.sh"
echo "============================================================"
