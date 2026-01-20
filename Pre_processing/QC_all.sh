#!/bin/bash 

source /home/lginerp/.bashrc
conda activate ngless


#$ -l os=centos7
#$ -N ngless_INCLIVA24
#$ -cwd

#$ -l m_mem_free=20G 
#$ -l h_rt=36:00:00 
#$ -t 1-19 
#$ -tc 19
#$ -pe smp 5

ngless --create-report --trace --temporary-directory='/scratch/AG_Forslund/Lola/Secuencias_INCLIVA_2024/tmp' QC_all.ngl
