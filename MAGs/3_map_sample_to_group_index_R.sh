#!/bin/bash

#SBATCH --job-name=map2contigs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-12  # adjust to number of samples

source /home/lginerp/.bashrc


# --- Activate environment ---
micromamba activate semibin

# --- Variables (passed via --export or set defaults) ---
GROUP_CONTIGS="${GROUP_CONTIGS:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/Responders/concatenated.fa}"
NAMES_FILE="${NAMES_FILE:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Responders.txt}"
READS_SPAIN="${READS_SPAIN:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs}"
READS_UK="${READS_UK:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Patel_2021_09_10}"
OUT_BAM_DIR="${OUT_BAM_DIR:-/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/bams}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"


sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NAMES_FILE")
if [[ -z "$sample" ]]; then
  echo "No sample found for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

mkdir -p "$OUT_BAM_DIR"

index_base="${GROUP_CONTIGS%.*}"
index_lock="${index_base}.index.lock"
if [[ ! -e "${index_base}.1.bt2" ]]; then
  # Only one process should build the index
  if ( set -o noclobber; echo "$$" > "$index_lock" ) 2> /dev/null; then
    trap 'rm -f "$index_lock"; exit $?' INT TERM EXIT
    echo "Building Bowtie2 index for $GROUP_CONTIGS"
    bowtie2-build --threads "$THREADS" "$GROUP_CONTIGS" "$index_base"
    rm -f "$index_lock"
    trap - INT TERM EXIT
  else
    # Wait for index to be built by another process
    while [ -e "$index_lock" ]; do
      sleep 10
    done
  fi
fi

# --- Determine sequencing type ---
if [[ -f "${READS_SPAIN}/${sample}filtered.pair.1.fq.gz" && -f "${READS_SPAIN}/${sample}filtered.pair.2.fq.gz" ]]; then
    echo "$sample detected as paired-end (Illumina)"
    R1="${READS_SPAIN}/${sample}filtered.pair.1.fq.gz"
    R2="${READS_SPAIN}/${sample}filtered.pair.2.fq.gz"
    bowtie2 -x "$index_base" -1 "$R1" -2 "$R2" --threads "$THREADS" \
      | samtools view -@ "$THREADS" -bS -F 4 -q 1 - \
      | samtools sort -@ "$THREADS" -o "${OUT_BAM_DIR}/${sample}.sorted.bam"
elif [[ -f "${READS_UK}/${sample}filtered.fq.gz" ]]; then
    echo "$sample detected as single-end (Ion Torrent)"
    RSE="${READS_UK}/${sample}filtered.fq.gz"
    bowtie2 -x "$index_base" -U "$RSE" --threads "$THREADS" \
      | samtools view -@ "$THREADS" -bS -F 4 -q 1 - \
      | samtools sort -@ "$THREADS" -o "${OUT_BAM_DIR}/${sample}.sorted.bam"
else
    echo "ERROR: reads not found for $sample in either Spain or UK folders" >&2
    exit 1
fi

samtools index -@ "$THREADS" "${OUT_BAM_DIR}/${sample}.sorted.bam"
echo "Finished mapping $sample → ${OUT_BAM_DIR}/${sample}.sorted.bam"
