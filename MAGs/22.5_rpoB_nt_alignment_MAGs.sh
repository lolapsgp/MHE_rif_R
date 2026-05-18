#!/bin/bash
#SBATCH --job-name=rpoB_blast_align
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=rpoB_blast_%A_%a.out
#SBATCH --error=rpoB_blast_%A_%a.err

# ------------------ Environment ------------------
source /home/lginerp/.bashrc
micromamba activate eggnog

# ------------------ Paths ------------------
RPOB_SEQ_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/rpoB_tblastn"
COMBINED_FA="${RPOB_SEQ_DIR}/all_MAGs_rpoB.fasta"
OUTDIR="${RPOB_SEQ_DIR}/blast_alignments"
mkdir -p "$OUTDIR"

# ------------------ Combine all rpoB sequences ------------------
echo "Combining all rpoB sequences..."
cat ${RPOB_SEQ_DIR}/*/*_rpoB_nuc.fasta > "$COMBINED_FA"

# Make sure headers are unique
awk '/^>/{print ">"++i;next}{print}' "$COMBINED_FA" > "${COMBINED_FA%.fasta}_uniq.fasta"
COMBINED_FA="${COMBINED_FA%.fasta}_uniq.fasta"

# ------------------ Build BLAST DB ------------------
echo "Building BLAST database..."
makeblastdb -in "$COMBINED_FA" -dbtype nucl -out "${COMBINED_FA%.fasta}_db"

# ------------------ Run pairwise BLAST for each sequence ------------------
echo "Running pairwise BLAST..."
while read -r LINE; do
    # Skip FASTA header lines
    if [[ "$LINE" =~ ^">" ]]; then
        SEQ_ID=$(echo "$LINE" | cut -c2-)
        QUERY_FA="${OUTDIR}/${SEQ_ID}_query.fasta"
        echo "$LINE" > "$QUERY_FA"
        read -r SEQ
        echo "$SEQ" >> "$QUERY_FA"

        # Output alignment
        OUTPUT="${OUTDIR}/${SEQ_ID}_vs_all.tsv"
        blastn -query "$QUERY_FA" \
               -db "${COMBINED_FA%.fasta}_db" \
               -out "$OUTPUT" \
               -outfmt 0 \
               -max_target_seqs 10 \
               -num_threads $SLURM_CPUS_PER_TASK
    fi
done < "$COMBINED_FA"

echo "All pairwise BLAST alignments done! Outputs in $OUTDIR"

echo "All pairwise BLAST alignments done! Outputs in $OUTDIR"