#!/bin/bash

#SBATCH --job-name=rpoB_blast_MAGs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-9
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# ------------------ Environment ------------------
source /home/lginerp/.bashrc
micromamba activate eggnog

# ------------------ USER INPUT ------------------
rpoB_PROT_FA="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/rpoB.clean.faa"
MAG_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/checkm2_out_all/protein_files"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/rpoB_blast_all_MAGs"

mkdir -p "$OUTDIR"

# ------------------ Get MAG ------------------
MAG=$(ls "$MAG_DIR"/*.faa | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -z "$MAG" ] && exit 1

MAG_NAME=$(basename "$MAG" .faa)
MAG_OUT="${OUTDIR}/${MAG_NAME}"
mkdir -p "$MAG_OUT"

echo "Processing $MAG_NAME"

# ------------------ Build BLAST DB (try regardless) ------------------
DB="${MAG_OUT}/${MAG_NAME}_blastdb"

echo "Building BLAST DB..."
makeblastdb -in "$MAG" -dbtype prot -out "$DB" 2> "${MAG_OUT}/makeblastdb.log"

# ------------------ Run BLASTP ------------------
OUTPUT_ALIGN="${MAG_OUT}/${MAG_NAME}_rpoB_alignment.txt"
OUTPUT_TAB="${MAG_OUT}/${MAG_NAME}_rpoB_hits.tsv"

echo "Running BLASTP..."

blastp \
    -query "$rpoB_PROT_FA" \
    -db "$DB" \
    -out "$OUTPUT_ALIGN" \
    -evalue 1e-3 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -max_target_seqs 10 \
    -outfmt 0 \
    2> "${MAG_OUT}/blastp_align.log"

blastp \
    -query "$rpoB_PROT_FA" \
    -db "$DB" \
    -out "$OUTPUT_TAB" \
    -evalue 1e-3 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -max_target_seqs 10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen \
                 qstart qend sstart send evalue bitscore \
                 qseq sseq" \
    2> "${MAG_OUT}/blastp_tab.log"

# ------------------ Presence (purely result-based) ------------------
PRESENCE_FILE="${MAG_OUT}/${MAG_NAME}_rpoB_presence.txt"

if [ -s "$OUTPUT_TAB" ]; then
    echo "$MAG_NAME: rpoB detected (BLAST)" > "$PRESENCE_FILE"
else
    echo "$MAG_NAME: rpoB NOT detected (BLAST)" > "$PRESENCE_FILE"
fi

echo "Done $MAG_NAME"