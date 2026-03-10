#!/bin/bash
# Script 17a — Download docked poses for triple-target hits from NIC5
# Pan-coronavirus RTC Inhibitor Discovery Pipeline
# Olivier Nsekuye | GIGA-VIN Lab, University of Liège | 2026-03-10
#
# Usage: bash scripts/bash/17a_download_poses_nic5.sh
# Requires: ssh agent running with NIC5 key loaded
#   eval $(ssh-agent) && ssh-add ~/.ssh/id_rsa.ceci

BASE=~/projects/rtc-pan-coronavirus
SCRATCH=/scratch/ulg/gigambd/onsekuye/rtc-screening

HITS="351017 13633807 13633805 351016 5024943 351018 3169307 5024944 351019 13662104 5024945"

mkdir -p $BASE/04-hits/poses/{NSP12-NSP7,NSP9-NSP12,NSP12-NSP8}

for TARGET in NSP12-NSP7 NSP9-NSP12 NSP12-NSP8; do
    echo "Downloading poses for $TARGET..."
    for ZID in $HITS; do
        scp nic5:$SCRATCH/results/${TARGET}/${ZID}_out.pdbqt \
            $BASE/04-hits/poses/${TARGET}/ 2>/dev/null
    done
    COUNT=$(ls $BASE/04-hits/poses/${TARGET}/ | wc -l)
    echo "  $COUNT files downloaded for $TARGET"
done

echo "Done."
