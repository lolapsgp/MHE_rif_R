#!/bin/bash 

#$ -l os=any
#$ -N ngless_fun_no_writing_bam_2
#$ -cwd

#$ -l m_mem_free=60G 
#$ -l h_rt=96:00:00 
#$ -t 1-19 
#$ -tc 19
#$ -pe smp 5

source /home/lginerp/.bashrc
conda activate ngless

ngless --create-report --trace --temporary-directory='/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/tmp' func_anno_4.ngl
