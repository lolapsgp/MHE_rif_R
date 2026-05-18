library(omixerRpm)
#Import NGLess outputs
eggnog_all <- read_tsv(
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_abundance_normalized_strains.tsv",
  col_types = cols(.default = "c")
)

library(dplyr)

# Optional: convert numeric columns from character to numeric
eggnog_all <- eggnog_all %>%
  mutate(
    gene_count = as.numeric(gene_count),
    total_genes = as.numeric(total_genes),
    rel_abundance = as.numeric(rel_abundance)
  )

# Check structure
str(eggnog_all)
head(eggnog_all)

library(dplyr)
library(tidyr)

# Assume your merged table is called eggnog_all

library(tidyverse)

# Pivot to wide format: rows = KEGG_ko, columns = MAG
ko_matrix <- eggnog_all %>%
  select(MAG, KEGG_ko, rel_abundance) %>%   # use rel_abundance for omixerRPM
  pivot_wider(
    names_from = MAG,
    values_from = rel_abundance,
    values_fill = 0   # fill missing combinations with 0
  ) %>%
  column_to_rownames(var = "KEGG_ko")      # KEGG KO as rownames

# Check result
head(ko_matrix)

write.table(ko_matrix, file="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KEGG_input_GMM_strains.tsv", quote=FALSE, sep='\t', col.names = NA)

# Load the mapping database (v.1.09 created by me adding modules from 104 to 109)
listDB()
db <- loadDB("GMMs.v1.09")

# Run the module mapping on the loaded table.
mods <- rpm("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KEGG_input_GMM_strains.tsv", minimum.coverage=0.3, annotation = 1, module.db = db)

# get the name of the first predicted module
getNames(db, mods@annotation[1,]) 

GMMs<-data.frame(mods@abundance) 
names<-as.data.frame(mods@annotation)
rownames(GMMs)<-names$Module
long_names<-as.data.frame(mods@db@module.names)
GMMs$names = long_names$V2[match(rownames(GMMs),rownames(long_names))]
colnames(GMMs) <- gsub("filtered$", "", colnames(GMMs))

#Write the file to a csv to save it. 
write.csv(GMMs, file = "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/GMMs_strains.csv")



library(pheatmap)

# Keep 'names' column separately
gmm_names <- GMMs$names

# Convert all other columns to numeric
gmm_mat <- GMMs[, -which(names(GMMs) == "names")]  # remove 'names' column
gmm_mat <- apply(gmm_mat, 2, as.numeric)           # ensure numeric

# Assign row names
rownames(gmm_mat) <- gmm_names

# Define groups
annotation_col <- data.frame(
  Group = ifelse(colnames(gmm_mat) == "PC304_3_SemiBin_21", "R", "NR")
)
rownames(annotation_col) <- colnames(gmm_mat)

# Plot
pheatmap(
  gmm_mat,
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  annotation_col = annotation_col,
  main = "GMM abundance heatmap"
)

library(pheatmap)

# Define groups
annotation_col <- data.frame(
  Group = ifelse(colnames(gmm_mat) == "PC304_3_SemiBin_21", "R", "NR")
)
rownames(annotation_col) <- colnames(gmm_mat)

# Define colors with alpha
# Add alpha by appending two hex digits: 80 ~ 50% transparency
group_colors <- c("NR" = "#CC79A780", "R" = "#009E7380")  

# Plot heatmap
plot<-pheatmap(
  gmm_mat,
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_colors),
  main = "GMM abundance heatmap"
)

ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/GMM_heatmap.pdf", plot = plot, width = 6, height = 6)

