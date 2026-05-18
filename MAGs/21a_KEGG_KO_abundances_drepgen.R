library(tidyverse)
# Responders
# Directory with eggNOG annotation files
annot_dir <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/eggnog_mapper_derep_gen"

# Find all annotation files
annot_files <- list.files(
  annot_dir,
  pattern = "\\.emapper\\.annotations$",
  full.names = TRUE
)

# Function to read and clean one file
read_clean_eggnog <- function(file) {
  
  # MAG name from filename
  mag_id <- basename(file) %>%
    str_remove("\\.emapper\\.annotations$")
  
  # Read table (ignores ## lines automatically)
  df <- read_tsv(
    file,
    comment = "##",
    col_types = cols(.default = "c")
  )
  
  # Remove leading # from first column name if present
  colnames(df)[1] <- str_remove(colnames(df)[1], "^#")
  
  # Add MAG column
  df %>%
    mutate(MAG = mag_id) %>%
    relocate(MAG)
}

# Read all files and merge
clean_eggnog_R <- map_dfr(annot_files, read_clean_eggnog)

# Quick sanity check
print(dim(clean_eggnog_R))
glimpse(clean_eggnog_R)

# Write clean table
write_tsv(
  clean_eggnog_R,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/all_strains_eggnog_clean.tsv"
)


# Sanity check
stopifnot(all(c("MAG", "query", "KEGG_ko") %in% colnames(clean_eggnog_R)))

# ===============================
# 2. Clean KEGG KO column
# ===============================

ko_long <- clean_eggnog_R %>%
  filter(
    !is.na(KEGG_ko),
    KEGG_ko != "-",
    KEGG_ko != ""
  ) %>%
  # remove "ko:" prefix if present
  mutate(KEGG_ko = str_remove_all(KEGG_ko, "ko:")) %>%
  # split multiple KOs into long format
  separate_rows(KEGG_ko, sep = ",") %>%
  distinct(MAG, query, KEGG_ko)

# ===============================
# 3. KO abundance per MAG
# ===============================

ko_abundance <- ko_long %>%
  count(MAG, KEGG_ko, name = "gene_count")

# Save long-format table
write_tsv(
  ko_abundance,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_abundance_strains.tsv"
)

# ===============================
# 4. KO × MAG matrix (wide)
# ===============================

ko_matrix <- ko_abundance %>%
  pivot_wider(
    names_from = MAG,
    values_from = gene_count,
    values_fill = 0
  )

write_tsv(
  ko_matrix,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_abundance_matrix_strains.tsv"
)

# ===============================
# 5. Normalized KO abundance (optional)
# ===============================

# Total genes per MAG
mag_gene_counts <- clean_eggnog_R %>%
  count(MAG, name = "total_genes")

ko_abundance_norm <- ko_abundance %>%
  left_join(mag_gene_counts, by = "MAG") %>%
  mutate(
    rel_abundance = gene_count / total_genes
  )

write_tsv(
  ko_abundance_norm,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_abundance_normalized_strains.tsv"
)

# ===============================
# 6. Presence / absence table
# ===============================

ko_presence <- ko_matrix %>%
  mutate(across(-KEGG_ko, ~ ifelse(. > 0, 1, 0)))

write_tsv(
  ko_presence,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/KO_presence_absence_matrix_strains.tsv"
)




####################################### CAZy #######################################
library(tidyverse)


# Responders
eggnogR <- read_tsv(
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/all_strains_eggnog_clean.tsv",
  col_types = cols(.default = "c")
)

stopifnot("CAZy" %in% colnames(eggnogR))

# ===============================
# 2. Clean and expand CAZy column
# ===============================

cazy_longR <- eggnogR %>%
  filter(
    !is.na(CAZy),
    CAZy != "-",
    CAZy != ""
  ) %>%
  separate_rows(CAZy, sep = ",") %>%
  distinct(MAG, query, CAZy)

# ===============================
# 3. Count CAZy families per MAG
# ===============================

cazy_countsR <- cazy_longR %>%
  count(MAG, CAZy, name = "gene_count")

write_tsv(
  cazy_countsR,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/CAZy_family_counts_strains.tsv"
)

# ===============================
# 4. Wide CAZy matrix
# ===============================

cazy_matrixR <- cazy_countsR %>%
  pivot_wider(
    names_from = MAG,
    values_from = gene_count,
    values_fill = 0
  )

write_tsv(
  cazy_matrixR,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/CAZy_family_matrix_strains.tsv"
)

# ===============================
# 5. Normalized CAZy counts
# ===============================

mag_gene_countsR <- eggnogR %>%
  count(MAG, name = "total_genes")

cazy_counts_normR <- cazy_countsR %>%
  left_join(mag_gene_countsR, by = "MAG") %>%
  mutate(
    rel_abundance = gene_count / total_genes
  )

write_tsv(
  cazy_counts_normR,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/CAZy_family_counts_normalized_strains.tsv"
)

# ===============================
# 6. Group CAZy into major classes
# ===============================

cazy_class_countsR <- cazy_longR %>%
  mutate(
    CAZy_class = str_extract(CAZy, "^[A-Z]+")
  ) %>%
  count(MAG, CAZy_class, name = "gene_count")

write_tsv(
  cazy_class_countsR,
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/outputs_R/Tables/CAZy_class_counts_strains.tsv"
)

ggplot(cazy_class_countsR, aes(x = MAG, y = gene_count, fill = CAZy_class)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CAZy class composition per strain", y = "Gene count", x = "MAG")
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/CAZy_barplot.pdf", width = 6, height = 4)

library(pheatmap)

mat <- cazy_matrixR %>%
  column_to_rownames("CAZy") %>%
  as.matrix()

plot<-pheatmap(mat, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "CAZy family abundance heatmap")
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/CAZy_heatmap.pdf", plot = plot, width = 4, height = 5)

