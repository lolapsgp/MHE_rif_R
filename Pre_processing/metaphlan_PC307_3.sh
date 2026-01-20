#Run for one sample

INPUT_DIR='/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs/'
metaphlan ${INPUT_DIR}PC307_3filtered.pair.1.fq.gz,${INPUT_DIR}PC307_3filtered.pair.2.fq.gz \ 
--nproc 5 --input_type fastq \ 
--index mpa_vJun23_CHOCOPhlAnSGB_202307 --force \ 
--bowtie2db /fast/AG_Forslund/shared/references/metaphlan_db_guix \ 
--bowtie2out ./bowtie/PC307_3_profiled_metaphlan.bowtie2.bz2 \ 
--samout ./metaplhan_out/PC307_3.sam.bz2 \ 
--read_min_len 70 -o ./metaplhan_out/PC307_3_profiled_metaphlan.tsv

