#!/bin/bash
#SBATCH --job-name=checkm2_pcopri
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=checkm2_%j.out
#SBATCH --error=checkm2_%j.err

# --------------------
# Load environment
# --------------------
source /home/lginerp/.bashrc
micromamba activate checkm2_env

# --------------------
# Define paths
# --------------------
BASE_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri"

R_DIR="${BASE_DIR}/R"
NR_DIR="${BASE_DIR}/NR"
ALL_MAGS="${BASE_DIR}/all_mags"

OUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/checkm2_out_all"

# --------------------
# Prepare input MAGs
# --------------------
mkdir -p "$ALL_MAGS"
mkdir -p "$OUT_DIR"

# Clean old symlinks
rm -f "$ALL_MAGS"/*.fa

# Symlink all MAGs (R + NR)
ln -s "$R_DIR"/*.fa "$ALL_MAGS"/
ln -s "$NR_DIR"/*.fa "$ALL_MAGS"/

echo "Number of MAGs:"
ls "$ALL_MAGS"/*.fa | wc -l

# --------------------
# Run CheckM2
# --------------------
checkm2 predict \
  --input "$ALL_MAGS" \
  --output-directory "$OUT_DIR" \
  --threads ${SLURM_CPUS_PER_TASK} \
  -x fa \
  --force


echo "CheckM2 completed successfully"
