#!/bin/bash

#SBATCH --job-name=metahit_assembly
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024

#SBATCH --mem=20G
#SBATCH --time=36:00:00
#SBATCH -n 6
#SBATCH --array=1-19


source /home/lginerp/.bashrc
conda activate ngless

#-------------------------------------------------------------------------------------------------------------------------------
ngless --create-report --trace --temporary-directory='/scratch/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/tmp' MEGAhit_assembly.ngl
