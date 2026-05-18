library(omixerRpm)
#Import NGLess outputs
eggnog_R <- read_tsv(
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_abundance_normalized_R.tsv",
  col_types = cols(.default = "c")
)
eggnog_NR <- read_tsv(
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_abundance_normalized_NR.tsv",
  col_types = cols(.default = "c")
)

library(dplyr)

# Add Group column to each table
eggnog_NR <- eggnog_NR %>% 
  mutate(Group = "NR")

eggnog_R <- eggnog_R %>% 
  mutate(Group = "R")

# Combine tables
eggnog_all <- bind_rows(eggnog_NR, eggnog_R)

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

# Aggregate gene counts per KEGG KO and Group
ko_group_summary <- eggnog_all %>%
  group_by(KEGG_ko, Group) %>%      # group by KO and Group
  summarise(total_rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") 

# Pivot to wide format: rows = KEGG_ko, columns = Group (R and NR)
ko_group_wide <- ko_group_summary %>%
  pivot_wider(
    names_from = Group,
    values_from = total_rel_abundance,
    values_fill = 0                 # fill missing combinations with 0
  )

# Optional: set rownames to KEGG_ko
ko_matrix <- ko_group_wide %>%
  column_to_rownames(var = "KEGG_ko")

# Check result
head(ko_matrix)

write.table(ko_matrix, file="/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KEGG_input_GMM.tsv", quote=FALSE, sep='\t', col.names = NA)

# Load the mapping database (v.1.09 created by me adding modules from 104 to 109)
listDB()
db <- loadDB("GMMs.v1.09")

# Run the module mapping on the loaded table.
mods <- rpm("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KEGG_input_GMM.tsv", minimum.coverage=0.3, annotation = 1, module.db = db)

# get the name of the first predicted module
getNames(db, mods@annotation[1,]) 

GMMs<-data.frame(mods@abundance) 
names<-as.data.frame(mods@annotation)
rownames(GMMs)<-names$Module
long_names<-as.data.frame(mods@db@module.names)
GMMs$names = long_names$V2[match(rownames(GMMs),rownames(long_names))]
colnames(GMMs) <- gsub("filtered$", "", colnames(GMMs))

#Write the file to a csv to save it. 
write.csv(GMMs, file = "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/GMMs.csv")