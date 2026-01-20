#!/bin/bash
#Filter contigs smaller than 1000bp

# Loop over all relevant contig files
for file in Assembly/assembled_contigs/*.contigs.fasta
do
    # Extract the basename, e.g., PC239filtered from PC239filtered.contigs.fasta
    basename=$(basename "$file" .contigs.fasta)
    # Construct output file path
    out="Assembly/assembled_filtered_contigs/${basename}.fcontigs.fasta"
    # Run seqkit to filter contigs ≥1000 bp
    seqkit seq -m 1000 "$file" > "$out"
done