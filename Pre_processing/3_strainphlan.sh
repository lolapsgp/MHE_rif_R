#!/bin/bash

#SBATCH --job-name=strainphlan
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

# Load environment
source /home/lginerp/.bashrc
micromamba activate metaphlan_py310

# =========================
# Paths
# =========================
BASE_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/strainphlan"

CONS_MARKERS_DIR="${BASE_DIR}/consensus_markers"
CLADE_MARKERS="${BASE_DIR}/segatella_markers/t__SGB1626.fna"
REF_GENOME="/fast/AG_Forslund/Lola/Sofware/Prevotella_genome/Prevotella_copri_A_pan-genome.fna"
OUTPUT_DIR="${BASE_DIR}/strainphlan_out"

mkdir -p "$OUTPUT_DIR"

echo "Running StrainPhlAn for Segatella copri (SGB1626)"

# =========================
# Run StrainPhlAn
# =========================

strainphlan \
  -s "${CONS_MARKERS_DIR}"/*.json.bz2 \
  -m "$CLADE_MARKERS" \
  -r "$REF_GENOME" \
  -c t__SGB1626 \
  -o "$OUTPUT_DIR" \
  -n 6 \
  --phylophlan_mode accurate \
  --mutation_rates \
  --sample_with_n_markers 5 \
  --marker_in_n_samples_perc 20 \
  --sample_with_n_markers_perc 5 \
  --database /fast/AG_Forslund/shared/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl 

echo "Done! Results are in $OUTPUT_DIR"
