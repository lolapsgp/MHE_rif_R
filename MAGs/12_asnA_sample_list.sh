 #!/bin/bash

# ---------------- USER INPUT ----------------
ILLUMINA_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs"
IONTORRENT_DIR="/fast/AG_Forslund/Lola/Secuences_bibliography/Patel_2021_09_10"

OUTPUT_FILE="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/sample_list.txt"
# --------------------------------------------

echo "# SampleName ForwardRead ReverseRead" > $OUTPUT_FILE

# ---------------- Illumina paired-end ----------------
# Look for *.pair.1.fq.gz and match corresponding *.pair.2.fq.gz
for R1 in ${ILLUMINA_DIR}/*.pair.1.fq.gz; do
    # Remove ".pair.1.fq.gz" to get sample prefix
    SAMPLE=$(basename $R1 | sed 's/\.pair\.1\.fq\.gz//')
    R2="${ILLUMINA_DIR}/${SAMPLE}.pair.2.fq.gz"
    
    if [ -f "$R2" ]; then
        echo -e "${SAMPLE}\t${R1}\t${R2}" >> $OUTPUT_FILE
    else
        echo "Warning: R2 not found for $SAMPLE, skipping..."
    fi
done

# ---------------- IonTorrent single-end ----------------
# Make sure the glob expands to actual files
shopt -s nullglob
for SE in ${IONTORRENT_DIR}/*STfiltered.fq.gz; do
    SAMPLE=$(basename $SE | sed 's/.fq.gz//')
    echo -e "${SAMPLE}\t${SE}" >> $OUTPUT_FILE
done
shopt -u nullglob

echo "Sample list created at $OUTPUT_FILE"
