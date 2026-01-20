#!/bin/bash

source /home/lginerp/.bashrc
conda activate biobakery

#$ -l os=any
#$ -N tax_metaphlan
#$ -cwd

#$ -l m_mem_free=6G 
#$ -l h_rt=24:00:00 
#$ -t 1-19 
#$ -tc 19 
#$ -pe smp 5

#-------------------------------------------------------------------------------------------------------------------------------
#Run for multiple samples
#loop over all input files (make sure you have deleted previously generated Bowtie 2 output files to prevent errors):
#example wiki loop 
#for i in SRS*.fq.gz; do metaphlan $i --input_type fasta --nproc 4 > ${i%.fq.gz}_profile.txt; done

#Loop for sequences with paired-end reads
for i in *.pair.1.fq.gz; do
    sample=${i%.pair.1.fq.gz}
    mkdir -p ./metaphlan_output/${sample}
    metaphlan ./outputs/${sample}_1.fq.gz,./outputs/${sample}_2.fq.gz \ 
        --nproc 5 \ 
        --input_type fastq \ 
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 \ 
        --force \ 
        --bowtie2db /fast/AG_Forslund/shared/references/metaphlan_db_guix \ 
        --bowtie2out ./metaphlan_output/${sample}/${sample}.bowtie2.bz2 \ 
        --samout ./metaphlan_output/${sample}/${sample}.sam.bz2 \ 
        --read_min_len 70 \ 
        -o ./metaphlan_output/${sample}/${sample}_profile.txt
done

#Merge the output files (one per sample) in one table
merge_metaphlan_tables.py ./metaphlan_output/*/*_profile.txt > merged_abundance_table.txt

