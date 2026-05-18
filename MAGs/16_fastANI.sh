#!/bin/bash
#SBATCH --job-name=fastANI_MAGs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

source /home/lginerp/.bashrc
micromamba activate fastani_env

# Run FastANI
fastANI \
--ql genome_list.txt \
--rl genome_list.txt \
-o ANI_all_vs_all.txt \
-t 4 \
--minFraction 0.05
echo "FastANI finished"