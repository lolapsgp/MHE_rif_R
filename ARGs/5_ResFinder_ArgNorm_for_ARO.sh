#!/bin/bash


# Activate environment
source /home/lginerp/.bashrc
# Activate environment
micromamba activate argnorm

cd /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Resfinder_results

for d in *; do
  if [ -d "$d" ] && [ -f "$d/ResFinder_results_tab.txt" ]; then
    argnorm resfinder -i "$d/ResFinder_results_tab.txt" -o "$d/ResFinder_results_tab_argNorm.tsv"
  fi
done
