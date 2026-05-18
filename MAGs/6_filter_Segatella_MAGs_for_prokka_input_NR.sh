#!/bin/bash

BASE_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly"
GTDB_SUMMARY="${BASE_DIR}/Segatella_copri/results/gtdbtk/NR/gtdbtk.bac120.summary.tsv"
BINS_DIR="${BASE_DIR}/Segatella_copri/scripts/test_out/NR/bins"
OUT_DIR="${BASE_DIR}/Segatella_copri/MAGs/Prevotella_copri/NR"

mkdir -p "$OUT_DIR"

awk -F '\t' '
NR == 1 { next }                               # skip header
$2 ~ /s__Prevotella copri/ { print $1 }       # correct species name
' "$GTDB_SUMMARY" | while read -r MAG; do

    BIN_FILE="${BINS_DIR}/${MAG}.fa.gz"

    if [ -f "$BIN_FILE" ]; then
        cp "$BIN_FILE" "$OUT_DIR/"
    else
        echo "WARNING: $BIN_FILE not found"
    fi

done

echo "Done. Copied Prevotella copri MAGs to:"
echo "$OUT_DIR"