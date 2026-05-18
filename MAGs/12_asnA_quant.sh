#!/bin/bash

#SBATCH --job-name=asnA
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-20      # Total number of Illumina samples
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# Load environment
source /home/lginerp/.bashrc
micromamba activate bowtie2

# ------------------ USER INPUT ------------------
GENE_FA="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/asnA.fa"
ILLUMINA_SAMPLES_FILE="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/sample_list.txt"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/asnA"

# Shared Bowtie2 index (built once)
BT2_INDEX_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/bowtie2_index"
BT2_INDEX_PREFIX="${BT2_INDEX_DIR}/asnA"

mkdir -p $OUTDIR

# ------------------ Get this job's sample ------------------
# Only take lines with 3 columns (paired-end Illumina)
SAMPLE_LINE=$(awk 'NF==3' $ILLUMINA_SAMPLES_FILE | sed -n "${SLURM_ARRAY_TASK_ID}p")

SAMPLE_NAME=$(echo $SAMPLE_LINE | awk '{print $1}')
READ1=$(echo $SAMPLE_LINE | awk '{print $2}')
READ2=$(echo $SAMPLE_LINE | awk '{print $3}')

# Sample-specific output folder
SAMPLE_OUT="${OUTDIR}/${SAMPLE_NAME}"
mkdir -p $SAMPLE_OUT

# ------------------ Bowtie2 Mapping ------------------

# Build index (gene-specific)
mkdir -p "$BT2_INDEX_DIR"

# Build index only once (first job that reaches this)
if [ ! -f "${BT2_INDEX_PREFIX}.1.bt2" ]; then
    echo "Building Bowtie2 index for asnA..."
    bowtie2-build "$GENE_FA" "$BT2_INDEX_PREFIX"
else
    echo "Bowtie2 index already exists. Skipping build."
fi

echo "[1/3] Mapping Illumina reads for sample $SAMPLE_NAME ..."
bowtie2 -x "$BT2_INDEX_PREFIX" \
        -1 "$READ1" -2 "$READ2" \
        -S "${SAMPLE_OUT}/${SAMPLE_NAME}_mapped.sam" \
        --very-sensitive -p "$SLURM_CPUS_PER_TASK"


# Convert SAM → sorted BAM
echo "[2/3] Converting SAM to sorted BAM ..."
samtools view -bS ${SAMPLE_OUT}/${SAMPLE_NAME}_mapped.sam | samtools sort -o ${SAMPLE_OUT}/${SAMPLE_NAME}_mapped.sorted.bam

# Count mapped reads
echo "[3/3] Counting mapped reads ..."
MAPPED_READS=$(samtools view -c -F 4 ${SAMPLE_OUT}/${SAMPLE_NAME}_mapped.sorted.bam)

# Gene length in bp
GENE_LENGTH=$(grep -v ">" $GENE_FA | tr -d '\n' | wc -c)

# Total reads (Illumina paired)
TOTAL_READS=$(( $(zcat $READ1 | wc -l)/4 + $(zcat $READ2 | wc -l)/4 ))

# Calculate FPKM
FPKM=$(awk -v mapped=$MAPPED_READS -v gen_length=$GENE_LENGTH -v total=$TOTAL_READS \
        'BEGIN{print (mapped*1e9)/(gen_length*total)}')

# Save results
echo -e "Sample\tMapped_Reads\tGene_Length\tTotal_Reads\tFPKM" > ${SAMPLE_OUT}/${SAMPLE_NAME}_asnA_abundance.txt
echo -e "${SAMPLE_NAME}\t${MAPPED_READS}\t${GENE_LENGTH}\t${TOTAL_READS}\t${FPKM}" >> ${SAMPLE_OUT}/${SAMPLE_NAME}_asnA_abundance.txt

echo "Done for sample $SAMPLE_NAME."
