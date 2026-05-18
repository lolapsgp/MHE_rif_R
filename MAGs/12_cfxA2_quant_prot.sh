#!/bin/bash

#SBATCH --job-name=cfxA2_diamond
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-20
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# ------------------ Environment ------------------
source /home/lginerp/.bashrc
micromamba activate eggnog   

# ------------------ USER INPUT ------------------
ASN_PROT_FA="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/cfxA2.faa"
ILLUMINA_SAMPLES_FILE="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/sample_list_all.txt"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/cfxA2_diamond"

# DIAMOND database (built once)
DIAMOND_DB_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/diamond_db"
DIAMOND_DB="${DIAMOND_DB_DIR}/cfxA2"

mkdir -p "$OUTDIR" "$DIAMOND_DB_DIR"


# ------------------ Build DIAMOND DB (once) ------------------
if [ ! -f "${DIAMOND_DB}.dmnd" ]; then
    echo "Building DIAMOND database for cfxA2..."
    diamond makedb --in "$ASN_PROT_FA" -d "$DIAMOND_DB"
else
    echo "DIAMOND database already exists. Skipping."
fi

# ------------------ Get FASTQ pair for this array task ------------------

READ1=$(ls /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs/*pair.1.fq.gz | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$READ1" ]; then
    echo "No READ1 found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

READ2=${READ1/pair.1.fq.gz/pair.2.fq.gz}

SAMPLE_NAME=$(basename "$READ1" | sed 's/.pair.1.fq.gz//')

SAMPLE_OUT="${OUTDIR}/${SAMPLE_NAME}"
mkdir -p "$SAMPLE_OUT"

# ------------------ DIAMOND blastx ------------------
echo "[1/2] Running DIAMOND blastx for sample $SAMPLE_NAME ..."

diamond blastx \
    -d "$DIAMOND_DB" \
    -q "$READ1" \
    -o "${SAMPLE_OUT}/${SAMPLE_NAME}_cfxA2_hits.tsv" \
    --sensitive \
    --evalue 1e-5 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outfmt 6

# ------------------ Quantification ------------------
echo "[2/2] Quantifying cfxA2 hits ..."

# Number of reads with at least one hit
HIT_READS=$(cut -f1 "${SAMPLE_OUT}/${SAMPLE_NAME}_cfxA2_hits.tsv" | sort | uniq | wc -l)

# Total reads in R1
TOTAL_READS=$(zcat "$READ1" | wc -l)
TOTAL_READS=$((TOTAL_READS / 4))

# Gene length in amino acids
GENE_LEN_AA=$(grep -v ">" "$ASN_PROT_FA" | tr -d '\n' | wc -c)
GENE_LEN_AA=$((GENE_LEN_AA))

# RPK-like normalization
RPK=$(awk -v hits="$HIT_READS" -v len="$GENE_LEN_AA" \
          'BEGIN{print hits / (len/1000)}')

# Save results
OUTFILE="${SAMPLE_OUT}/${SAMPLE_NAME}_cfxA2_abundance.txt"

echo -e "Sample\tHit_Reads\tGene_Length_aa\tTotal_Reads\tRPK" >| "$OUTFILE"
echo -e "${SAMPLE_NAME}\t${HIT_READS}\t${GENE_LEN_AA}\t${TOTAL_READS}\t${RPK}" >> "$OUTFILE"

echo "Done for sample $SAMPLE_NAME."
