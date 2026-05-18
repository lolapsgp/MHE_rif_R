#!/bin/bash
#SBATCH --job-name=tblastn_rpoB
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=1-9
#SBATCH --output=tblastn_rpoB_%A_%a.out
#SBATCH --error=tblastn_rpoB_%A_%a.err

# ------------------ Environment ------------------
source /home/lginerp/.bashrc
micromamba activate eggnog

# ------------------ Paths ------------------
RPOB_PROT="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/rpoB.clean.faa"
MAG_FA_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/rpoB_tblastn"
mkdir -p "$OUTDIR"

# ------------------ MAG list and BLAST hit IDs ------------------
# For completeness, store sseqid but tblastn works directly with protein
MAG_LIST=("PC304_3_SemiBin_21" "RS_03_d01_ST_SemiBin_1" "RS_03_d90_ST_SemiBin_1" \
          "RS_04_d01_ST_SemiBin_0" "RS_04_d90_ST_SemiBin_0" "RS_06_d01_ST_SemiBin_1" \
          "RS_25_d01_ST_SemiBin_0" "RS_25_d90_ST_SemiBin_0" "RS_37_d90_ST_SemiBin_0")

# SLURM_ARRAY_TASK_ID picks one MAG
MAG="${MAG_LIST[$SLURM_ARRAY_TASK_ID-1]}"
MAG_FA="${MAG_FA_DIR}/${MAG}.fa"
MAG_OUT="${OUTDIR}/${MAG}"
mkdir -p "$MAG_OUT"

# Skip if no MAG fasta
if [ ! -f "$MAG_FA" ]; then
    echo "$MAG: FASTA file not found, skipping."
    exit 1
fi

echo "Processing $MAG ..."

# ------------------ Make BLAST nucleotide DB ------------------
DB="${MAG_OUT}/${MAG}_db"
makeblastdb -in "$MAG_FA" -dbtype nucl -out "$DB" 2> "${MAG_OUT}/makeblastdb.log"

# ------------------ Run tblastn ------------------
TBLASTN_OUT="${MAG_OUT}/${MAG}_rpoB_tblastn.tsv"
tblastn -query "$RPOB_PROT" \
        -db "$DB" \
        -evalue 1e-3 \
        -num_threads $SLURM_CPUS_PER_TASK \
        -max_target_seqs 1 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "$TBLASTN_OUT" \
        2> "${MAG_OUT}/tblastn.log"

# ------------------ Extract nucleotide sequence ------------------
if [ -s "$TBLASTN_OUT" ]; then
    # Read first line (best hit)
    hit_line=$(head -n1 "$TBLASTN_OUT")
    contig=$(echo "$hit_line" | cut -f2)      # sseqid
    start=$(echo "$hit_line" | cut -f9)
    end=$(echo "$hit_line" | cut -f10)

    # Make sure start < end
    if [ "$start" -gt "$end" ]; then
        tmp=$start
        start=$end
        end=$tmp
        strand="-"
    else
        strand="+"
    fi

    OUT_FASTA="${MAG_OUT}/${MAG}_rpoB_nuc.fasta"

    # Extract sequence with samtools faidx (make sure faidx exists)
    samtools faidx "$MAG_FA" "${contig}:${start}-${end}" > "$OUT_FASTA"

    echo "$MAG: rpoB nucleotide sequence extracted -> $OUT_FASTA (strand $strand)"
else
    echo "$MAG: No tblastn hit found!"
fi

echo "Done $MAG"