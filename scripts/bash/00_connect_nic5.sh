#!/bin/bash
# Script 00 — Start SSH agent and connect to NIC5
# Usage: source scripts/bash/00_connect_nic5.sh
#   (use 'source' not 'bash' so the agent persists in current shell)

eval $(ssh-agent)
ssh-add ~/.ssh/id_rsa.ceci
echo "Agent ready. Connect with: ssh nic5"
echo "Transfer with: scp localfile nic5:/scratch/ulg/gigambd/onsekuye/"
