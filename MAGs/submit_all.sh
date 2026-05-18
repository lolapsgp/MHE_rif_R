#!/usr/bin/env bash
set -euo pipefail

# Spain
SPAIN_COUNT=$(wc -l < /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/names.txt)
sbatch --array=1-${SPAIN_COUNT} /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/scripts/run_gtdbtk.sh Spain

# UK
UK_COUNT=$(wc -l < /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/names_UK.txt)
sbatch --array=1-${UK_COUNT} /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/scripts/run_gtdbtk.sh UK
