#!/bin/bash
# Set input and output directories
INPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs"
OUTPUT_DIR="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/rgi_outputs"

# Load micromamba RGI environment
source ~/.bashrc
micromamba activate rgi  # replace with your actual environment name if different

cd /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024
rgi load --card_json /fast/AG_Forslund/shared/references/CARD/card.json --local


# Loop through all filtered contigs files
for contigs_file in "$INPUT_DIR"/*filtered.fcontigs.fasta; do
    # Extract the base sample name (e.g., PC239 from PC239filtered.fcontigs.fasta)
    filename=$(basename "$contigs_file")
    sample_name="${filename%filtered.fcontigs.fasta}"

    echo "Running RGI for sample: $sample_name"

    rgi main \
      --input_sequence "$contigs_file" \
      --output_file "$OUTPUT_DIR/${sample_name}rgi" \
      --input_type contig \
      --local \
      --clean
done

echo "RGI processing complete!"