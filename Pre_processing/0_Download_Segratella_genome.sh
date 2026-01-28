#!/bin/bash

#SBATCH --job-name=download_prevotella
#SBATCH --chdir=/fast/AG_Forslund/Lola/Sofware/Prevotella_genome
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH -n 1

# Target directory
TARGET_DIR="/fast/AG_Forslund/Lola/Sofware/Prevotella_genome"
mkdir -p "$TARGET_DIR"

# Pangenome download URL
URL="https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000003697/downloads/pan-genome.fna"

# Output file name
OUTPUT_FILE="${TARGET_DIR}/Prevotella_copri_A_pan-genome.fna"

cd "$TARGET_DIR"

echo "Downloading Prevotella copri_A pangenome..."
wget -c "$URL" -O "$OUTPUT_FILE"

if [ $? -eq 0 ]; then
    echo "Download complete: $OUTPUT_FILE"
else
    echo "Download failed. Check your network or URL."
fi

