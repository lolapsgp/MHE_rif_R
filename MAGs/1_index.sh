#!/bin/bash

#SBATCH --job-name=map2contigs
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-1  # adjust to number of samples

source /home/lginerp/.bashrc


# --- Activate environment ---
micromamba activate semibin

bowtie2-build --threads 8   /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/NonR_contigs.fasta   /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/NonR_contigs
bowtie2-build --threads 8   /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/responders_contigs.fasta   /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/contigs_groups/responders_contigs
