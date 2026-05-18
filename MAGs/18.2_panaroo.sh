#!/bin/bash
#SBATCH --job-name=panaroo_MAGs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

# Load environment
source /home/lginerp/.bashrc
micromamba activate panaroo_env

# Directories
GFF_DIR=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/panaroo/gff_files
OUTPUT_DIR=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/panaroo/panaroo_output
mkdir -p /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/panaroo/panaroo_output
mkdir -p /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/panaroo/gff_files

# Run Panaroo
panaroo -i ${GFF_DIR}/*.gff \
        -o ${OUTPUT_DIR} \
        -t 16 \
        --clean-mode strict