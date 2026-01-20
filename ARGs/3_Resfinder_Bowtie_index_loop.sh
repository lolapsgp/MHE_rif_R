#!/bin/bash

#SBATCH --job-name=map_ARGs_to_contigs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024
#SBATCH --mem=10G
#SBATCH --time=00:25:00
#SBATCH -n 4
#SBATCH --array=1-19

source /home/lginerp/.bashrc
micromamba activate bowtie2

# Directories
FASTA_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Resfinder_results"
INDEX_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/bowtie_indexes_per_sample_resfinder"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_of_sample_ok.txt)

# Create output dir if it doesn't exist
mkdir -p "$INDEX_DIR"

# Loop over all FASTA files
for fasta in "$FASTA_DIR"/${SAMPLE}/ResFinder_Resistance_gene_seq_filtered.fsa; do
    # Extract base name (e.g., PC239_ARG_index)

    echo "Building Bowtie2 index for $SAMPLE..."
    
    # Run bowtie2-build
    bowtie2-build "$fasta" "$INDEX_DIR/$SAMPLE"
done

echo "All Bowtie2 indices built."
