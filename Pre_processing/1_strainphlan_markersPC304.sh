#!/bin/bash

#SBATCH --job-name=strainphlan_markers
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

# Load environment
source /home/lginerp/.bashrc
micromamba activate metaphlan_py310

sample2markers.py     -i /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/metaphlan_output/PC304_3/PC304_3.sam.bz2     -o "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/strainphlan"     -n 8     --database /fast/AG_Forslund/shared/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl



