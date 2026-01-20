#!/bin/bash

#SBATCH --job-name=map_ARGs_to_contigs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024
#SBATCH --mem=10G
#SBATCH --time=00:25:00
#SBATCH -n 4
#SBATCH --array=1-19

# Activate environment
source /home/lginerp/.bashrc
micromamba activate bowtie2

# Threads and directories
THREADS=4
IDX_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/bowtie_indexes_per_sample"
INDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/bowtie_out"
mkdir -p "$OUTDIR"

# Get the sample from the list
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_of_sample_ok.txt)

# Define input and index
R1="$INDIR/${SAMPLE}filtered.pair.1.fq.gz"
R2="$INDIR/${SAMPLE}filtered.pair.2.fq.gz"
INDEX="$IDX_DIR/${SAMPLE}"

# Check if index exists
if [[ ! -f "${INDEX}.1.bt2" ]]; then
    echo "[$(date)] ERROR: Index for $SAMPLE not found: ${INDEX}.1.bt2"
    exit 1
fi

echo "[$(date)] Mapping $SAMPLE to its own index..."

# Run bowtie2
bowtie2 -x "$INDEX" -1 "$R1" -2 "$R2" \
    --very-sensitive-local --no-unal -p "$THREADS" \
    -S "$OUTDIR/${SAMPLE}.sam"

# Convert and sort
samtools view -bS "$OUTDIR/${SAMPLE}.sam" | samtools sort -o "$OUTDIR/${SAMPLE}.sorted.bam"
samtools index "$OUTDIR/${SAMPLE}.sorted.bam"

# Generate idxstats
if [[ -s "$OUTDIR/${SAMPLE}.sorted.bam" && -s "$OUTDIR/${SAMPLE}.sorted.bam.bai" ]]; then
    TMPIDX="$OUTDIR/${SAMPLE}.idxstats.tmp"
    samtools idxstats "$OUTDIR/${SAMPLE}.sorted.bam" > "$TMPIDX"
    mv "$TMPIDX" "$OUTDIR/${SAMPLE}.idxstats.txt"
else
    echo "[$(date)] ERROR: BAM or index missing for $SAMPLE"
    exit 1
fi

# Cleanup
rm "$OUTDIR/${SAMPLE}.sam"

echo "[$(date)] Done with $SAMPLE"


