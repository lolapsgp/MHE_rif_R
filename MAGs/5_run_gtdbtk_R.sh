#!/bin/bash

#SBATCH --job-name=gtdbtk_R
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --error=%x_%j.err

source /home/lginerp/.bashrc
micromamba activate gtdbtk24


# === PATHS ===
BASE_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly"
MAGS_DIR_R="${BASE_DIR}/Segatella_copri/scripts/test_out/R/bins"
OUTDIR="${BASE_DIR}/Segatella_copri/results/gtdbtk/R"


# Temp folder
export TMPDIR="/fast/AG_Forslund/Lola/tmp/gtdbtk_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

# GTDB-Tk database
export GTDBTK_DATA_PATH="/fast/AG_Forslund/shared/GTDB-Tk/release220"

mkdir -p "$OUTDIR"

echo "Running GTDB-Tk on MAGs in:"
echo "$MAGS_DIR_R"

gtdbtk classify_wf \
    --genome_dir "$MAGS_DIR_R" \
    --extension fa.gz \
    --out_dir "$OUTDIR" \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --skip_ani_screen

echo "GTDB-Tk finished"

