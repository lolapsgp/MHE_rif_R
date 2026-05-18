#!/bin/bash
#SBATCH --job-name=rpoB_variant
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-35
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

source /home/lginerp/.bashrc
micromamba activate bowtie2

SEQ_A_FA="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/rpoB_A.fa"
SEQ_C_FA="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/rpoB.fa"
SAMPLE_FILE="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/sample_list.filtered.txt"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/rpoB"
mkdir -p "$OUTDIR"

SEQ_A=$(grep -v ">" "$SEQ_A_FA" | tr -d '\n')
SEQ_C=$(grep -v ">" "$SEQ_C_FA" | tr -d '\n')

SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_FILE")
IFS=',' read -r SAMPLE_NAME READ1 READ2 <<< "$SAMPLE_LINE"

if [ -z "$SAMPLE_NAME" ] || [ -z "$READ1" ]; then
    echo "ERROR parsing sample line: $SAMPLE_LINE"
    exit 1
fi

SAMPLE_OUT="${OUTDIR}/${SAMPLE_NAME}"
mkdir -p "$SAMPLE_OUT"

count_full_seq_gz() {
    local seq="$1"
    local file="$2"
    zcat "$file" | sed -n '2~4p' | grep -c "$seq"
}

COUNT_A=$(count_full_seq_gz "$SEQ_A" "$READ1")
COUNT_C=$(count_full_seq_gz "$SEQ_C" "$READ1")

if [ -n "$READ2" ] && [ -f "$READ2" ]; then
    COUNT_A=$((COUNT_A + $(count_full_seq_gz "$SEQ_A" "$READ2")))
    COUNT_C=$((COUNT_C + $(count_full_seq_gz "$SEQ_C" "$READ2")))
fi

TOTAL=$((COUNT_A + COUNT_C))
if [ "$TOTAL" -gt 0 ]; then
    FREQ_C=$(awk -v alt=$COUNT_C -v total=$TOTAL 'BEGIN{print alt/total}')
else
    FREQ_C=0
fi

OUTFILE="${SAMPLE_OUT}/${SAMPLE_NAME}_full_sequence_counts.txt"
echo -e "Sample\tSeq_A_Count\tSeq_C_Count\tTotal\tFreq_C" > "$OUTFILE"
echo -e "${SAMPLE_NAME}\t${COUNT_A}\t${COUNT_C}\t${TOTAL}\t${FREQ_C}" >> "$OUTFILE"

echo "Done for $SAMPLE_NAME: Seq_A=${COUNT_A}, Seq_C=${COUNT_C}, Freq_C=${FREQ_C}"