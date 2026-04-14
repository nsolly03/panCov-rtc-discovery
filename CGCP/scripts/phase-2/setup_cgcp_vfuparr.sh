#!/bin/bash
# =============================================================================
# CGCP VFUparr Setup Script — All 9 Interfaces
# =============================================================================
# Creates 9 independent VFUparr screening directories on NIC5, one per
# RTC interface, using v3 docking box coordinates.
#
# Prerequisites (already confirmed):
#   - VFUparr template cloned at $GLOBALSCRATCH/cgcp-vfuparr-template/
#   - Vina at ~/.conda/envs/rtc-screening/bin/vina
#   - Receptor PDBQTs at $GLOBALSCRATCH/rtc-screening/receptors/
#   - Filtered SMILES at $GLOBALSCRATCH/rtc-screening/enamine/filtered/
#
# Usage:
#   bash setup_cgcp_vfuparr.sh
#
# Author : CGCP Pipeline | Olivier Nsekuye | ULiège GIGA-VIN
# Date   : 2026-04-12
# =============================================================================

set -e

# ─────────────────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────────────────
SCRATCH="/scratch/ulg/gigambd/onsekuye"
TEMPLATE="$SCRATCH/cgcp-vfuparr-template"
VFUPARR_ROOT="$SCRATCH/cgcp-vfuparr"
RECEPTORS="$SCRATCH/rtc-screening/receptors"
SMILES_RAW="$SCRATCH/rtc-screening/enamine/filtered/filtered_tier2_ALL.smi"
SMILES_CONV="$SCRATCH/cgcp-vfuparr/shared/smiles_tier2_vfuparr.smi"
VINA_BIN="/home/ulg/gigambd/onsekuye/.conda/envs/rtc-screening/bin/vina"
CONDA_ENV="rtc-screening"

# SLURM settings for NIC5
SLURM_ACCOUNT="ulg"          # NIC5 account — adjust if needed
SLURM_PARTITION="batch"
SLURM_WALLTIME="2-00:00:00"
SLURM_MEM="4000M"
SLURM_NTASKS="1"
SLURM_CPUS="4"               # 4 CPUs per task for vina exhaustiveness

# Screening settings
EXHAUSTIVENESS=8
NUM_MOLS=225014734           # full tier2 filtered library
MAX_NUM_JOBS=1000            # 1000 array jobs = ~225k compounds per job
SCORE_THRESHOLD="-7.5"       # save poses only for score < -7.5 kcal/mol

# ─────────────────────────────────────────────────────────────────────────────
# V3 DOCKING BOXES (from verify_docking_boxes.py)
# Format: "CENTER_X CENTER_Y CENTER_Z SIZE_X SIZE_Y SIZE_Z"
# ─────────────────────────────────────────────────────────────────────────────
declare -A BOXES
BOXES["NSP12-NSP7"]="94.401 73.209 121.340 32.0 46.0 30.0"
BOXES["NSP12-NSP8"]="95.325 117.201 110.709 66.0 50.0 46.0"
BOXES["NSP9-NSP12"]="134.504 168.078 175.149 46.0 40.0 48.0"
BOXES["NSP10-NSP16"]="74.163 22.930 17.098 32.0 44.0 38.0"
BOXES["NSP7-NSP8"]="-9.252 -10.325 -13.930 34.0 40.0 32.0"
BOXES["NSP10-NSP14"]="-7.499 4.093 -19.631 58.0 56.0 72.0"
BOXES["NSP13-Helicase"]="-32.205 12.718 -6.774 42.0 42.0 42.0"
BOXES["NSP12-NSP13"]="200.254 194.151 156.838 30.0 28.0 26.0"
BOXES["NSP15"]="-59.490 39.441 -21.127 34.0 52.0 52.0"

# ─────────────────────────────────────────────────────────────────────────────
# RECEPTOR PDBQT FILES (in $SCRATCH/rtc-screening/receptors/)
# ─────────────────────────────────────────────────────────────────────────────
declare -A RECEPTOR_FILES
RECEPTOR_FILES["NSP12-NSP7"]="receptor_NSP12-NSP7.pdbqt"
RECEPTOR_FILES["NSP12-NSP8"]="receptor_NSP12-NSP8.pdbqt"
RECEPTOR_FILES["NSP9-NSP12"]="receptor_NSP9-NSP12.pdbqt"
RECEPTOR_FILES["NSP10-NSP16"]="receptor_NSP10-NSP16.pdbqt"
RECEPTOR_FILES["NSP7-NSP8"]="receptor_NSP7-NSP8.pdbqt"
RECEPTOR_FILES["NSP10-NSP14"]="receptor_NSP10-NSP14.pdbqt"
RECEPTOR_FILES["NSP13-Helicase"]="receptor_NSP13-Helicase.pdbqt"
RECEPTOR_FILES["NSP12-NSP13"]="receptor_NSP12-NSP13.pdbqt"
RECEPTOR_FILES["NSP15"]="receptor_NSP15.pdbqt"

INTERFACES=(
    "NSP12-NSP7"
    "NSP12-NSP8"
    "NSP9-NSP12"
    "NSP10-NSP16"
    "NSP7-NSP8"
    "NSP10-NSP14"
    "NSP13-Helicase"
    "NSP12-NSP13"
    "NSP15"
)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Convert SMILES format
# Raw:       SMILES ID E1=1 E2=1 E3=1  (space-separated)
# VFUparr:  SMILES,ID  (comma-separated, first two columns only)
# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "  CGCP VFUparr Setup"
echo "============================================================"

mkdir -p "$VFUPARR_ROOT/shared"

if [ -f "$SMILES_CONV" ]; then
    echo "[SKIP] SMILES already converted: $SMILES_CONV"
    echo "       Lines: $(wc -l < $SMILES_CONV)"
else
    echo "[1/4] Converting SMILES format..."
    echo "      Input:  $SMILES_RAW"
    echo "      Output: $SMILES_CONV"
    awk '{print $1","$2}' "$SMILES_RAW" > "$SMILES_CONV"
    N=$(wc -l < "$SMILES_CONV")
    echo "      Done: $N compounds"
    # Verify format
    echo "      Sample:"
    head -3 "$SMILES_CONV" | sed 's/^/        /'
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Check receptors
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "[2/4] Checking receptor PDBQT files..."
MISSING_RECEPTORS=0
for iface in "${INTERFACES[@]}"; do
    rfile="$RECEPTORS/${RECEPTOR_FILES[$iface]}"
    if [ -f "$rfile" ]; then
        size=$(du -sh "$rfile" | cut -f1)
        echo "      ✅ $iface — $size"
    else
        echo "      ❌ MISSING: $rfile"
        MISSING_RECEPTORS=$((MISSING_RECEPTORS + 1))
    fi
done

if [ $MISSING_RECEPTORS -gt 0 ]; then
    echo ""
    echo "  ⚠  $MISSING_RECEPTORS receptor(s) missing."
    echo "  These need to be prepared with obabel before continuing."
    echo "  See prepare_receptors.sh for conversion commands."
    echo "  Continuing setup for interfaces with available receptors..."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — Build per-interface directories
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "[3/4] Building interface directories..."

for iface in "${INTERFACES[@]}"; do
    IFACE_DIR="$VFUPARR_ROOT/$iface"
    rfile="$RECEPTORS/${RECEPTOR_FILES[$iface]}"

    echo ""
    echo "  ── $iface ──────────────────────────────────"

    # Skip if receptor missing
    if [ ! -f "$rfile" ]; then
        echo "  ⚠  Skipping — receptor not found"
        continue
    fi

    # Create directory structure
    mkdir -p "$IFACE_DIR/DATA"
    mkdir -p "$IFACE_DIR/OUTPUTS"

    # Parse box coordinates
    read CX CY CZ SX SY SZ <<< "${BOXES[$iface]}"

    # ── Copy receptor ──────────────────────────────────────────────────
    cp "$rfile" "$IFACE_DIR/DATA/docking_receptor.pdbqt"

    # ── Symlink vina executable ────────────────────────────────────────
    ln -sf "$VINA_BIN" "$IFACE_DIR/DATA/vina"
    chmod +x "$IFACE_DIR/DATA/vina" 2>/dev/null || true

    # ── Copy dataset_calc.py from template ────────────────────────────
    cp "$TEMPLATE/dataset_calc.py" "$IFACE_DIR/dataset_calc.py"

    # Patch dataset_calc.py to use vina instead of qvina
    sed -i 's|./DATA/qvina|./DATA/vina|g' "$IFACE_DIR/dataset_calc.py"
    sed -i 's|./DATA/smina|./DATA/vina|g' "$IFACE_DIR/dataset_calc.py"

    # ── Write all.ctrl ─────────────────────────────────────────────────
    cat > "$IFACE_DIR/all.ctrl" << CTRL
# CGCP VFUparr Configuration
# Interface: ${iface}
# Generated: $(date +%Y-%m-%d)
# Docking boxes: v3 (verify_docking_boxes.py, 2026-04-12)

# SMILES library (tier2 filtered, 225M compounds)
SMILES_FILES=./shared/smiles_tier2_vfuparr.smi
NUM_MOLS=${NUM_MOLS}

# Receptor
RECEPTOR_LOCATION=./DATA/docking_receptor.pdbqt

# Docking parameters
EXHAUSTIVENESS=${EXHAUSTIVENESS}

# Docking box — v3 coordinates
CENTER_X=${CX}
CENTER_Y=${CY}
CENTER_Z=${CZ}
SIZE_X=${SX}
SIZE_Y=${SY}
SIZE_Z=${SZ}

# Job array size
MAX_NUM_JOBS=${MAX_NUM_JOBS}

# Pose saving threshold (only save poses with score < this value)
DOCKING_SCORE_THRESHOLD=${SCORE_THRESHOLD}
CTRL

    # Symlink shared SMILES (avoid copying 432M file 9 times)
    ln -sf "$SMILES_CONV" "$IFACE_DIR/shared"

    # ── Write submit.sh (NIC5-adapted) ────────────────────────────────
    cat > "$IFACE_DIR/submit.sh" << SLURM
#!/bin/bash
#SBATCH --job-name=CGCP_${iface}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --ntasks=${SLURM_NTASKS}
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_WALLTIME}
#SBATCH --array=1-${MAX_NUM_JOBS}
#SBATCH --output=logs/job_%A_task_%a.out
#SBATCH --error=logs/job_%A_task_%a.err
#SBATCH --open-mode=append

# ── Environment ──────────────────────────────────────────────────────
module --force purge
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Confirm vina available
echo "Vina: \$(which vina)"
echo "Python: \$(which python3)"
echo "Task: \$SLURM_ARRAY_TASK_ID"

# ── Run docking ──────────────────────────────────────────────────────
cd \$SLURM_SUBMIT_DIR
python3 dataset_calc.py \$SLURM_ARRAY_TASK_ID
SLURM

    # Create logs directory
    mkdir -p "$IFACE_DIR/logs"

    echo "  ✅ Directory ready: $IFACE_DIR"
    echo "     Box: center=($CX,$CY,$CZ) size=($SX×$SY×$SZ)"
done

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Summary and next steps
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[4/4] Summary"
echo "============================================================"
echo ""
echo "  VFUparr root: $VFUPARR_ROOT"
echo "  Interfaces set up:"
for iface in "${INTERFACES[@]}"; do
    IFACE_DIR="$VFUPARR_ROOT/$iface"
    if [ -d "$IFACE_DIR" ]; then
        echo "    ✅ $iface"
    else
        echo "    ❌ $iface (skipped)"
    fi
done

echo ""
echo "  Chunk file check:"
beegfs-ctl --getquota --uid $(id -u) 2>/dev/null | tail -2

echo ""
echo "============================================================"
echo "  NEXT STEPS"
echo "============================================================"
echo ""
echo "  1. Submit a PILOT test (1 interface, 10 tasks) first:"
echo "     cd $VFUPARR_ROOT/NSP12-NSP7"
echo "     sbatch --array=1-10 submit.sh"
echo ""
echo "  2. Check pilot output after ~30 min:"
echo "     cat logs/job_*_task_*.out | tail -20"
echo "     ls OUTPUTS/ | head -20"
echo ""
echo "  3. If pilot OK, submit all 9 interfaces:"
echo "     for iface in NSP12-NSP7 NSP12-NSP8 NSP9-NSP12 NSP10-NSP16 \\"
echo "                  NSP7-NSP8 NSP10-NSP14 NSP13-Helicase NSP12-NSP13 NSP15; do"
echo "       cd $VFUPARR_ROOT/\$iface && sbatch submit.sh && cd -"
echo "     done"
echo ""
echo "  NOTE: Each interface screens 225M compounds across 1000 array jobs"
echo "        (~225k compounds per job). Estimated runtime: 24-48h per interface."
echo ""
