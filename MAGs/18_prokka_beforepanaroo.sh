#!/bin/bash
#SBATCH --job-name=prokka_copri
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-7
#SBATCH --error=%x_%j.err

source /home/lginerp/.bashrc
micromamba activate prokka

# =====================
# PATHS
# =====================
INPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/drep_output/dereplicated_genomes"
PROKKA_OUT="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/prokka_dreplicated_genomes"

mkdir -p "$PROKKA_OUT" logs

# =====================
# GET MAG LIST (array-safe)
# =====================
MAG=$(ls "$INPUT_DIR"/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$MAG" ]; then
    echo "No MAG found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BASENAME=$(basename "$MAG" .fa)
OUT="${PROKKA_OUT}/${BASENAME}"

# Skip if already done
if [ -d "$OUT" ]; then
    echo "Skipping $BASENAME (already annotated)"
    exit 0
fi

echo "Annotating $BASENAME"

prokka \
    --outdir "$OUT" \
    --prefix "$BASENAME" \
    --locustag "$BASENAME" \
    --genus Prevotella \
    --usegenus \
    --metagenome \
    --kingdom Bacteria \
    --cpus 8 \
    "$MAG"

echo "Finished $BASENAME"
