#!/bin/bash
#SBATCH --job-name=drep_MAGs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

source /home/lginerp/.bashrc
micromamba activate drep_env

# Run dRep
dRep dereplicate /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/drep_output \
  -g /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly/*.fa \
  -p 8 \
  -comp 10 \
  -con 5 \
  --S_algorithm fastANI \
  --multiround_primary_clustering \
  -sa 0.99 \
  --genomeInfo /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/checkm2_out_all/genomeInfo.csv
  
echo "dRep finished"