#!/bin/bash

#SBATCH --job-name=semibin_NR
#SBATCH --partition=gpu
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --error=%x_%A_%a.err
#SBATCH --mem=40G
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=15
#SBATCH --array=1
#SBATCH --ntasks-per-node=1

source /home/lginerp/.bashrc

micromamba activate semibin

GROUP_CONTIGS="${GROUP_CONTIGS:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/NonR/concatenated.fa}"
BAM_DIR="${BAM_DIR:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/bams}"
OUTDIR="${OUTDIR:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/NonR}"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "$OUTDIR"

# Read sample names from the Responders.txt file
SAMPLES=($(cat /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/NonR.txt))

# Get the sample name corresponding to the SLURM array task ID
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}"  # SLURM array starts at 1, arrays in bash start at 0

# Find BAM files corresponding to the current sample (space-separated)
BAM_FILES=$(find "$BAM_DIR" -maxdepth 1 -type f -name "${SAMPLE}.sorted.bam")

if [ -z "$BAM_FILES" ]; then
    echo "No BAM files found for sample $SAMPLE, skipping."
    exit 1
fi

echo "Running SemiBin for sample $SAMPLE with BAMs: $BAM_FILES"
SemiBin2 multi_easy_bin \
    --input-fasta "$GROUP_CONTIGS" \
    --input-bam $BAM_FILES \
    --threads 15\
    --output "$OUTDIR/$SAMPLE"

echo "Finished semibin for sample $SAMPLE, bins likely in $OUTDIR/$SAMPLE/bins/"

