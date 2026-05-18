#!/bin/bash
#SBATCH --job-name=eggnog_copri_NR
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=82G
#SBATCH --time=24:00:00
#SBATCH --array=1-7
#SBATCH --error=%x_%j.err
#SBATCH --error=%x_%j.log

# =====================
# Load environment
# =====================
source /home/lginerp/.bashrc
micromamba activate eggnog   # <-- eggNOG-mapper env

# =====================
# Paths
# =====================
PROKKA_OUT="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/prokka_dreplicated_genomes"
EGGNOG_OUT="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/eggnog_mapper_derep_gen"
EGGNOG_DB="/fast/AG_Forslund/shared/eggnog-mapper"   # eggNOG database dir
TMP_DIR="/fast/AG_Forslund/Lola/tmp"

mkdir -p "$EGGNOG_OUT"

# =====================
# Get MAG list safely
# =====================
MAG_DIR=$(find "$PROKKA_OUT" -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
BASENAME=$(basename "$MAG_DIR")
FAA="$MAG_DIR/$BASENAME.faa"
OUT_PREFIX="$BASENAME"

if [ ! -f "$FAA" ]; then
    echo "Protein FASTA not found for $BASENAME"
    exit 1
fi

if [ -f "${OUT_PREFIX}.annotations" ]; then
    echo "Skipping $BASENAME (already annotated)"
    exit 0
fi

echo "Running eggNOG-mapper for $BASENAME"

# =====================
# Run eggNOG-mapper
# =====================
emapper.py \
  -i "$FAA" \
  -o "$OUT_PREFIX" \
  --output_dir "$EGGNOG_OUT" \
  --itype proteins \
  -m diamond \
  --data_dir "$EGGNOG_DB" \
  --temp_dir "$TMP_DIR" \
  --override

echo "Finished $BASENAME"
