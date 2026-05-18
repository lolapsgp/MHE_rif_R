#!/bin/bash

cd /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/gtdbtk/R/identify/

#Get the list of Segatella copri bins
ls /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/R/*.fa \
  | xargs -n 1 basename \
  | sed 's/\.fa$//' \
  > segatella_bins_R.txt

ls /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/NR/*.fa \
  | xargs -n 1 basename \
  | sed 's/\.fa$//' \
  > segatella_bins_NR.txt

#Keep only markers from those bins
awk -F'\t' '
FNR==NR { a[$1]; next }          # read MAG names
FNR==1 { print; next }           # print header of GTDB file
($1 in a)
' \
segatella_bins_R.txt \
/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/gtdbtk/R/identify/gtdbtk.bac120.markers_summary.tsv \
> gtdbtk.bac120.markers_summary.segatella_R.tsv


awk -F'\t' '
FNR==NR { a[$1]; next }          # read MAG names
FNR==1 { print; next }           # print header of GTDB file
($1 in a)
' \
segatella_bins_NR.txt \
/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/gtdbtk/NR/identify/gtdbtk.bac120.markers_summary.tsv \
> gtdbtk.bac120.markers_summary.segatella_NR.tsv

#extract unique genes
cut -f1,6 gtdbtk.bac120.markers_summary.segatella_R.tsv \
| tail -n +2 \
> segatella_R.complete_markers.tsv

cut -f1,6 gtdbtk.bac120.markers_summary.segatella_NR.tsv \
| tail -n +2 \
> segatella_NR.complete_markers.tsv




