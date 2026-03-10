#!/bin/bash
#SBATCH --job-name=rtc_screen
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-97
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --time=04:00:00
#SBATCH --partition=batch

# ── Environment ──────────────────────────────────────────────────────────
module load Miniconda3/23.9.0-0
source activate screening

SCRATCH=/scratch/ulg/gigambd/onsekuye/rtc-screening
PDBQT_DIR=$SCRATCH/pdbqt
REC_DIR=$SCRATCH/receptors

# ── Target (set via --export=TARGET=NSP12-NSP7 at submission) ────────────
TARGET=${TARGET:-NSP12-NSP7}
REC=$REC_DIR/receptor_${TARGET}.pdbqt
CFG=$REC_DIR/config_${TARGET}.txt
OUT_DIR=$SCRATCH/results/${TARGET}
mkdir -p $OUT_DIR

# ── Ligand chunk for this array task ─────────────────────────────────────
LIGANDS=($(ls $PDBQT_DIR/*.pdbqt))
TOTAL=${#LIGANDS[@]}
CHUNK=100
START=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
END=$(( START + CHUNK ))
if [ $END -gt $TOTAL ]; then END=$TOTAL; fi

echo "Task $SLURM_ARRAY_TASK_ID: ligands $START-$END of $TOTAL | Target: $TARGET"
echo "Started: $(date)"

# ── Parse box from config ─────────────────────────────────────────────────
CX=$(grep center_x $CFG | awk -F= '{print $2}' | tr -d ' ')
CY=$(grep center_y $CFG | awk -F= '{print $2}' | tr -d ' ')
CZ=$(grep center_z $CFG | awk -F= '{print $2}' | tr -d ' ')
SX=$(grep size_x   $CFG | awk -F= '{print $2}' | tr -d ' ')
SY=$(grep size_y   $CFG | awk -F= '{print $2}' | tr -d ' ')
SZ=$(grep size_z   $CFG | awk -F= '{print $2}' | tr -d ' ')

echo "Box: center=($CX,$CY,$CZ) size=($SX,$SY,$SZ)"

# ── Docking loop ─────────────────────────────────────────────────────────
SCORES_FILE=$OUT_DIR/scores_task${SLURM_ARRAY_TASK_ID}.tsv
echo -e "zinc_id\tscore" > $SCORES_FILE

python3 << PYEOF
import os, sys
from vina import Vina
from pathlib import Path

rec     = "$REC"
out_dir = "$OUT_DIR"
cx, cy, cz = float("$CX"), float("$CY"), float("$CZ")
sx, sy, sz = float("$SX"), float("$SY"), float("$SZ")
scores_file = "$SCORES_FILE"

ligands_all = sorted(Path("$PDBQT_DIR").glob("*.pdbqt"))
start, end  = $START, $END
ligands     = ligands_all[start:end]

print(f"Docking {len(ligands)} ligands against $TARGET")

v = Vina(sf_name='vina', cpu=$SLURM_CPUS_PER_TASK, verbosity=0)
v.set_receptor(rec)
v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])

with open(scores_file, 'a') as sf:
    for i, lig in enumerate(ligands):
        zinc_id = lig.stem
        out_path = Path(out_dir) / f"{zinc_id}_out.pdbqt"
        if out_path.exists():
            continue
        try:
            v.set_ligand_from_file(str(lig))
            v.dock(exhaustiveness=16, n_poses=5)
            v.write_poses(str(out_path), n_poses=5, overwrite=True)
            energy = v.energies(n_poses=1)[0][0]
            sf.write(f"{zinc_id}\t{energy:.3f}\n")
            sf.flush()
        except Exception as e:
            sf.write(f"{zinc_id}\tERROR\n")
        if (i+1) % 10 == 0:
            print(f"  [{i+1}/{len(ligands)}] done")

print("Task complete.")
PYEOF

echo "Finished: $(date)"
