#!/bin/bash 

source /home/lginerp/.bashrc
conda activate ngless


#$ -l os=centos7
#$ -N ngless_tax_INCLIVA24
#$ -cwd

#$ -l m_mem_free=60G 
#$ -l h_rt=72:00:00 
#$ -t 1-19 
#$ -tc 19
#$ -pe smp 5

ngless --create-report --trace --config-file '/fast/AG_Forslund/Lola/Sofware/ngless.conf' tax_anno.ngl
