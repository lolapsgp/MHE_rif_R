#!/bin/bash

#SBATCH --job-name=tax_metaphlan
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024

#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH -n 6
#SBATCH --array=1-19


source /home/lginerp/.bashrc
conda activate metaphlan

#-------------------------------------------------------------------------------------------------------------------------------
#Run for multiple samples
#loop over all input files (make sure you have deleted previously generated Bowtie 2 output files to prevent errors):
#example wiki loop
#for i in SRS*.fq.gz; do metaphlan $i --input_type fasta --nproc 4 > ${i%.fq.gz}_profile.txt; done

SAMPLE_LIST="file_names.txt"
sample=$(sed "$SLURM_ARRAY_TASK_ID"'q;d' "$SAMPLE_LIST")

output_dir="./metaphlan_output/${sample}"
output_file="${output_dir}/${sample}_profile.txt"

    # Check if the output file already exists
    if [[ -f "${output_file}" ]]; then
        echo "Output file for ${sample} already exists. Skipping..."
        exit
    fi

    mkdir -p "${output_dir}"

    input1="./outputs/${sample}.pair.1.fq.gz"
    input2="./outputs/${sample}.pair.2.fq.gz"

    metaphlan "${input1}","${input2}" --input_type fastq --nproc 4 --index mpa_vJun23_CHOCOPhlAnSGB_202307 --force --bowtie2db /fast/AG_Forslund/shared/references/metaphlan_db_guix --bowtie2out "${output_dir}/${sample}.bowtie2.bz2" --samout "${output_dir}/${sample}.sam.bz2" --read_min_len 70 -o "${output_file}"

#Merge the output files (one per sample) in one table
merge_metaphlan_tables.py ./metaphlan_output/*/*_profile.txt > merged_abundance_table.txt
