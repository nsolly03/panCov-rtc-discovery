#!/bin/bash
# ================================================================
# CGCP VFUparr Setup — Prof. Twizere NIC5 account
# Interfaces: NSP10-NSP14, NSP13-Helicase, NSP12-NSP13, NSP15
# Run from: Prof. Twizere's NIC5 account
# ================================================================

set -e

GLOBALSCRATCH="/scratch/ulg/gigambd/twizere"
VFUPARR_ROOT="${GLOBALSCRATCH}/cgcp-vfuparr"
OLIVIER_SCRATCH="/scratch/ulg/gigambd/onsekuye"
RECEPTOR_DIR="${OLIVIER_SCRATCH}/rtc-screening/receptors"
SMILES_FILE="${OLIVIER_SCRATCH}/cgcp-vfuparr/shared/smiles_tier2_vfuparr.smi"
VINA="/home/ulg/gigambd/onsekuye/.conda/envs/rtc-screening/bin/vina"
DATASET_CALC="${OLIVIER_SCRATCH}/cgcp-vfuparr/NSP12-NSP7/dataset_calc.py"

echo "============================================================"
echo "  CGCP VFUparr Setup — twizere account"
echo "============================================================"

# ── Step 1: Verify access to Olivier's files ─────────────────────
echo "[0/4] Verifying access to shared files..."
for f in "$SMILES_FILE" "$RECEPTOR_DIR" "$DATASET_CALC" "$VINA"; do
  if [ -r "$f" ]; then
    echo "  ✅ Accessible: $f"
  else
    echo "  ❌ NOT accessible: $f — contact onsekuye to fix permissions"
    exit 1
  fi
done

# ── Step 2: Create shared symlink directory ───────────────────────
echo "[1/4] Setting up shared directory..."
mkdir -p ${VFUPARR_ROOT}/shared
# Symlink to directory (NOT the file — learned from onsekuye's bug)
if [ ! -L "${VFUPARR_ROOT}/shared/smiles_tier2_vfuparr.smi" ]; then
  ln -s ${SMILES_FILE} ${VFUPARR_ROOT}/shared/smiles_tier2_vfuparr.smi
fi
echo "  ✅ SMILES linked"

# ── Step 3: Interface definitions ────────────────────────────────
declare -A CX CY CZ SX SY SZ RECEPTOR

CX["NSP10-NSP14"]="-7.499";   CY["NSP10-NSP14"]="4.093";    CZ["NSP10-NSP14"]="-19.631"
SX["NSP10-NSP14"]="58.0";     SY["NSP10-NSP14"]="56.0";     SZ["NSP10-NSP14"]="72.0"
RECEPTOR["NSP10-NSP14"]="receptor_NSP10-NSP14.pdbqt"

CX["NSP13-Helicase"]="-32.205"; CY["NSP13-Helicase"]="12.718";  CZ["NSP13-Helicase"]="-6.774"
SX["NSP13-Helicase"]="42.0";    SY["NSP13-Helicase"]="42.0";    SZ["NSP13-Helicase"]="42.0"
RECEPTOR["NSP13-Helicase"]="receptor_NSP13-Helicase.pdbqt"

CX["NSP12-NSP13"]="200.254";  CY["NSP12-NSP13"]="194.151";  CZ["NSP12-NSP13"]="156.838"
SX["NSP12-NSP13"]="30.0";     SY["NSP12-NSP13"]="28.0";     SZ["NSP12-NSP13"]="26.0"
RECEPTOR["NSP12-NSP13"]="receptor_NSP12-NSP13.pdbqt"

CX["NSP15"]="-59.490";  CY["NSP15"]="39.441";   CZ["NSP15"]="-21.127"
SX["NSP15"]="34.0";     SY["NSP15"]="52.0";     SZ["NSP15"]="52.0"
RECEPTOR["NSP15"]="receptor_NSP15.pdbqt"

# ── Step 4: Build interface directories ──────────────────────────
echo "[2/4] Building interface directories..."

for IFACE in NSP10-NSP14 NSP13-Helicase NSP12-NSP13 NSP15; do
  IFACE_DIR="${VFUPARR_ROOT}/${IFACE}"
  mkdir -p ${IFACE_DIR}/{DATA,OUTPUTS,logs}

  # Copy receptor
  cp ${RECEPTOR_DIR}/${RECEPTOR[$IFACE]} ${IFACE_DIR}/DATA/docking_receptor.pdbqt

  # Symlink vina
  ln -sf ${VINA} ${IFACE_DIR}/DATA/vina

  # Symlink shared directory (points to DIRECTORY not file)
  rm -f ${IFACE_DIR}/shared
  ln -s ${VFUPARR_ROOT}/shared ${IFACE_DIR}/shared

  # Copy dataset_calc.py from Olivier's validated setup
  cp ${DATASET_CALC} ${IFACE_DIR}/dataset_calc.py

  # Write all.ctrl
  cat > ${IFACE_DIR}/all.ctrl << CTRL
RECEPTOR_LOCATION=./DATA/docking_receptor.pdbqt
SMILES_FILES=./shared/smiles_tier2_vfuparr.smi
VINA_LOCATION=./DATA/vina
CENTER_X=${CX[$IFACE]}
CENTER_Y=${CY[$IFACE]}
CENTER_Z=${CZ[$IFACE]}
SIZE_X=${SX[$IFACE]}
SIZE_Y=${SY[$IFACE]}
SIZE_Z=${SZ[$IFACE]}
NUM_MODES=1
EXHAUSTIVENESS=8
SCORE_THRESHOLD=-7.5
NUM_CPUS=4
TOTAL_COMPOUNDS=225014734
NUM_TASKS=1000
CTRL

  # Write submit.sh
  cat > ${IFACE_DIR}/submit.sh << SLURM
#!/bin/bash
#SBATCH --job-name=CGCP_${IFACE}
#SBATCH --account=ceci
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4000M
#SBATCH --time=2-00:00:00
#SBATCH --array=1-1000
#SBATCH --output=logs/job_%A_task_%a.out
#SBATCH --error=logs/job_%A_task_%a.err
#SBATCH --open-mode=append

module --force purge
source /home/ulg/gigambd/onsekuye/.conda/envs/rtc-screening/../../etc/profile.d/conda.sh
conda activate /home/ulg/gigambd/onsekuye/.conda/envs/rtc-screening

echo "Vina: \$(which vina)"
echo "Python: \$(which python3)"
echo "Task: \$SLURM_ARRAY_TASK_ID"

cd \$SLURM_SUBMIT_DIR
python3 dataset_calc.py all.ctrl \$SLURM_ARRAY_TASK_ID
SLURM

  chmod +x ${IFACE_DIR}/submit.sh
  echo "  ✅ ${IFACE} — Box:(${CX[$IFACE]},${CY[$IFACE]},${CZ[$IFACE]}) size=(${SX[$IFACE]}×${SY[$IFACE]}×${SZ[$IFACE]})"
done

echo ""
echo "[3/4] Verifying setup..."
for IFACE in NSP10-NSP14 NSP13-Helicase NSP12-NSP13 NSP15; do
  IFACE_DIR="${VFUPARR_ROOT}/${IFACE}"
  echo "  ${IFACE}:"
  echo "    Receptor: $(ls -lh ${IFACE_DIR}/DATA/docking_receptor.pdbqt | awk '{print $5}')"
  echo "    SMILES:   $(ls -la ${IFACE_DIR}/shared/smiles_tier2_vfuparr.smi)"
done

echo ""
echo "[4/4] Summary"
echo "============================================================"
echo "  VFUparr root: ${VFUPARR_ROOT}"
echo "  Interfaces ready: NSP10-NSP14, NSP13-Helicase, NSP12-NSP13, NSP15"
echo ""
echo "  NEXT STEPS:"
echo "  1. Pilot test:"
echo "     cd ${VFUPARR_ROOT}/NSP10-NSP14"
echo "     sbatch --array=1-10 submit.sh"
echo "  2. Check pilot logs after ~30 min"
echo "  3. If OK: submit all 4 interfaces"
echo "============================================================"
