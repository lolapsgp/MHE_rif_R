#!/bin/bash

source /home/lginerp/.bashrc

# Usage:
# ./prepare_group_contigs.sh /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/Responders /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Responders.txt /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/responders_contigs.fasta responders
# ./prepare_group_contigs.sh /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/NonR /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/NonR.txt /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/NonR_contigs.fasta NonR
contig_dir="$1"
names_file="$2"
out_fasta="$3"
group_name="$4"   # string used to prefix contig names, e.g. responders

mkdir -p "$(dirname "$out_fasta")"

> "$out_fasta"
while read -r sample; do
  sample_fasta="${contig_dir}/${sample}.fcontigs.fasta"
  if [[ ! -f "$sample_fasta" ]]; then
    echo "WARNING: $sample_fasta not found, skipping" >&2
    continue
  fi
  # prefix contig headers to retain sample origin:
  awk -v pref="${sample}:${group_name}_${sample}_" '/^>/{print ">"pref substr($0,2); next} {print}' "$sample_fasta" >> "$out_fasta"
done < "$names_file"

echo "Wrote combined contigs to $out_fasta"
