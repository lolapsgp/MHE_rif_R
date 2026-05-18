#!/bin/bash

#SBATCH --job-name=asnA_diamond_MAGs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-20      # Adjust to number of MAGs
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# ------------------ Environment ------------------
source /home/lginerp/.bashrc
micromamba activate eggnog   # Ensure diamond is available

# ------------------ USER INPUT ------------------
ASN_PROT_FA="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/asnA.clean.faa"
MAG_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/checkm2_out_all/protein_files"
OUTDIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/asnA_in_MAGs"

mkdir -p "$OUTDIR"

# ------------------ Check ASN_PROT_FA ------------------
if [ ! -s "$ASN_PROT_FA" ]; then
    echo "ERROR: ASN_PROT_FA file is missing or empty"
    exit 1
fi

if ! head -n1 "$ASN_PROT_FA" | grep -q "^>"; then
    echo "ERROR: ASN_PROT_FA file does not start with >"
    exit 1
fi

# ------------------ Get this job's MAG ------------------
MAG=$(ls "$MAG_DIR"/*.faa | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$MAG" ]; then
    echo "No MAG found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

MAG_NAME=$(basename "$MAG" .faa)
MAG_OUT="${OUTDIR}/${MAG_NAME}"
mkdir -p "$MAG_OUT"

echo "Processing MAG: $MAG_NAME"

# ------------------ Skip empty or malformed MAGs ------------------
if [ ! -s "$MAG" ]; then
    echo "WARNING: MAG file $MAG is empty. Skipping."
    echo "$MAG_NAME: SKIPPED (empty MAG)" > "${MAG_OUT}/${MAG_NAME}_asnA_presence.txt"
    exit 0
fi

if ! head -n1 "$MAG" | grep -q "^>"; then
    echo "WARNING: MAG file $MAG is malformed. Skipping."
    echo "$MAG_NAME: SKIPPED (malformed MAG)" > "${MAG_OUT}/${MAG_NAME}_asnA_presence.txt"
    exit 0
fi

# ------------------ Build DIAMOND DB (MAG proteins) ------------------
DIAMOND_DB="${MAG_OUT}/${MAG_NAME}_db"
if [ ! -f "${DIAMOND_DB}.dmnd" ]; then
    echo "Building DIAMOND database for $MAG_NAME..."
    diamond makedb --in "$MAG" -d "$DIAMOND_DB"
else
    echo "DIAMOND database already exists. Skipping."
fi

# ------------------ DIAMOND blastp (asnA vs MAG) ------------------
OUTPUT="${MAG_OUT}/${MAG_NAME}_asnA_hits.tsv"
echo "Searching for asnA in $MAG_NAME..."
diamond blastp \
    -q "$ASN_PROT_FA" \
    -d "$DIAMOND_DB" \
    -o "$OUTPUT" \
    --evalue 1e-5 \
    --sensitive \
    --threads $SLURM_CPUS_PER_TASK \
    --outfmt 6

# ------------------ Report presence/absence ------------------
PRESENCE_FILE="${MAG_OUT}/${MAG_NAME}_asnA_presence.txt"
if [ -s "$OUTPUT" ]; then
    echo "$MAG_NAME: asnA detected" > "$PRESENCE_FILE"
else
    echo "$MAG_NAME: asnA NOT detected" > "$PRESENCE_FILE"
fi

echo "Done for MAG $MAG_NAME"
