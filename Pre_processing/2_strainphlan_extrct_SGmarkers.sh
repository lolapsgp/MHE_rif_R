#!/bin/bash

#SBATCH --job-name=strainphlan_markers_SC
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# Load environment
source /home/lginerp/.bashrc
micromamba activate metaphlan_py310

# =========================
# Paths
# =========================
BASE_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/strainphlan"
SC_MARKERS_DIR="${BASE_DIR}/segatella_markers"
mkdir -p "$SC_MARKERS_DIR"


echo "Processing SC markers"


# =========================
# Run extract_markers for SC
# =========================

extract_markers.py \
    -c t__SGB1626 \
    -o "$SC_MARKERS_DIR" \
    --database /fast/AG_Forslund/shared/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl


echo "Segatella markers generated"


