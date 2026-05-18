#!/bin/bash
#SBATCH --job-name=download_sequences
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

source /home/lginerp/.bashrc

# Create destination directory 
mkdir -p /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Genomes_strains
cd /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Genomes_strains

echo "Starting download of Segatella copri genomes..."

# WGS assemblies (.fsa_nt.gz)
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/AC/BX/ACBX02/ACBX02.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/JA/HR/GH/JAHRGH01/JAHRGH01.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/VU/NF/VUNF01/VUNF01.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/JA/HR/GE/JAHRGE01/JAHRGE01.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/JA/HR/GB/JAHRGB01/JAHRGB01.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/JA/HR/GJ/JAHRGJ01/JAHRGJ01.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/JA/HR/GF/JAHRGF01/JAHRGF01.1.fsa_nt.gz

# Complete chromosomes (direct FASTA)
wget -O CP152351.1.fasta "https://www.ncbi.nlm.nih.gov/nuccore/CP152351.1?report=fasta"
wget -O CP137558.1.fasta "https://www.ncbi.nlm.nih.gov/nuccore/CP137558.1?report=fasta"

#MAGs bibliography
#GenBank: CAMKRF000000000.1
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/CA/MK/RF/CAMKRF01/CAMKRF01.1.fsa_nt.gz
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/CC/XN/CD/CCXNCD01/CCXNCD01.1.fsa_nt.gz

echo "Download finished!"