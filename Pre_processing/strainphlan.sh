#!/bin/bash

#SBATCH --job-name=tax_metaphlan
#SBATCH --chdir=/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024

#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH -n 6
#SBATCH --array=1-19


source /home/lginerp/.bashrc
conda activate metaphlan

#-------------------------------------------------------------------------------------------------------------------------------

#StrainPhlAn
#The first step is to run MetaPhlAn 4 to obtain the sam output files. The sam files contain the alignment information from mapping 
#the reads of each sample against the MetaPhlAn 4 marker database. For that, run the script at stranPhlAn wiki.
#Done in metaphlan.sh script

#Step 2: Run sample_to_markers
#Providing the sam output files, from the prior step, to the sample_to_markers script will generate a marker file for each sample.
# The marker files contain the consensus of unique marker genes for each species found in the sample, which will be used 
#for SNP profiling.


SAMPLE_LIST="file_names.txt"
sample=$(sed "$SLURM_ARRAY_TASK_ID"'q;d' "$SAMPLE_LIST")

input_dir="./metaphlan_output/${sample}"
output_dir1="./strainphlan_outputs/${sample}"
output_dir2="./strainphlan_outputs/"

mkdir -p ./strainphlan_outputs
mkdir -p ./strainphlan_outputs/consensus_markers
sample2markers.py -i "${input_dir}/${sample}.bowtie2.bz2" -o strainphlan_outputs/consensus_markers --nproc 4


#Step 3: Generate trees from alignments
output_file="${output_dir1}/${sample}_strain.txt"

    # Check if the output file already exists
    if [[ -f "${output_file}" ]]; then
        echo "Output file for ${sample} already exists. Skipping..."
        exit
    fi
mkdir -p "${output_dir2}/clade_markers"
mkdir -p "${output_dir2}/output"
extract_markers.py -c t__SGB4933 -o "${output_dir2}/clade_markers"
#Wehre are the reference genomes??
strainphlan -s "${output_dir2}/consensus_markers/*.json.bz2" -m "${output_dir2}/clade_markers/t__SGB4933.fna" -r reference_genomes/*.fna.bz2 -o "${output_file}" -c t__SGB4933 --phylophlan_mode fast --nproc 4
