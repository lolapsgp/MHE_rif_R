#!/bin/bash

#SBATCH --job-name=semibin_NR
#SBATCH --partition=gpu
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --error=%x_%j.err
#SBATCH --mem=40G
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=15
#SBATCH --ntasks=1

source /home/lginerp/.bashrc
micromamba activate semibin

GROUP_CONTIGS="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/NonR/concatenated.fa"
BAM_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/bams/NR"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/NonR"
THREADS="$SLURM_CPUS_PER_TASK"

mkdir -p "$OUTDIR"

# Collect ALL BAM files
BAM_FILES=$(find "$BAM_DIR" -maxdepth 1 -type f -name "*.sorted.bam" | sort)

if [ -z "$BAM_FILES" ]; then
    echo "No BAM files found in $BAM_DIR"
    exit 1
fi

echo "Running SemiBin2 multi_easy_bin"
echo "Using BAM files:"
echo "$BAM_FILES"

SemiBin2 multi_easy_bin \
    --input-fasta "$GROUP_CONTIGS" \
    --input-bam $BAM_FILES \
    --tmpdir /fast/AG_Forslund/Lola/tmp/ \
    --threads "$THREADS" \
    --output "$OUTDIR"

echo "Finished SemiBin2 multi_easy_bin"
echo "Bins likely in $OUTDIR/bins/"
