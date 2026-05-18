#!/bin/bash
#SBATCH --job-name=phylophlan3
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

source /home/lginerp/.bashrc
micromamba activate phylophlan   # assuming phylophlan is here

# =====================
# PATHS
# =====================
INPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly"
OUTPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/phylophlan3"
DATABASE_DIR="/fast/AG_Forslund/Lola/Sofware/phylophlan"

CONFIG_FOLDER="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/phylophlan3"
CONFIG_FILE="supermatrix_aa.cfg"

mkdir -p "$OUTPUT_DIR"

# =====================
# RUN PHYLOPHLAN 3
# =====================

phylophlan \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  -d phylophlan \
  --databases_folder "$DATABASE_DIR" \
  --configs_folder "$CONFIG_FOLDER" \
  -f "$CONFIG_FILE" \
  --diversity low \
  --accurate \
  --nproc 8 \
  --genome_extension fa \
  --min_num_entries 2 \
  --verbose

