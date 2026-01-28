#!/bin/bash

#SBATCH --job-name=strainphlan_markers
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6
#SBATCH --array=1-35
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# Load environment
source /home/lginerp/.bashrc
micromamba activate metaphlan_py310

# =========================
# Paths
# =========================
BASE_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/strainphlan"
CONS_MARKERS_DIR="${BASE_DIR}/consensus_markers"
mkdir -p "$CONS_MARKERS_DIR"

# Sample list
SAMPLE_LIST="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/sample_list.txt"
SAMPLES=($(cat "$SAMPLE_LIST"))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}  # SLURM arrays are 1-indexed

echo "Processing sample: $SAMPLE"

# Path to input BAM / SAM file
SAMPLE_BAM="$BASE_DIR/bams/${SAMPLE}.sam.bz2"

# =========================
# Run sample2markers for this sample
# =========================
sample2markers.py \
    -i "$SAMPLE_BAM" \
    -o "$CONS_MARKERS_DIR" \
    -n 6 \
    --database /fast/AG_Forslund/shared/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl


echo "Consensus markers generated for $SAMPLE in $CONS_MARKERS_DIR"


