# Load necessary libraries
library(tidyverse)
library(pheatmap)

# Set path to your Panaroo RTAB file
panaroo_file <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/panaroo/panaroo_output/gene_presence_absence.Rtab"

# Read the matrix
gene_matrix <- read.delim(panaroo_file, header = TRUE, row.names = 1, check.names = FALSE)

# Convert to numeric (0/1)
gene_matrix <- as.matrix(gene_matrix)
gene_matrix <- apply(gene_matrix, 2, as.numeric)
rownames(gene_matrix) <- rownames(read.delim(panaroo_file, header = TRUE)) # keep gene IDs
# Simple presence/absence heatmap
final_plot<-pheatmap(gene_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         color = c("white", "steelblue"),
         main = "S. copri Pangenome: Gene Presence/Absence")
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Panaroo_presence_absence_genes.pdf", plot = final_plot, width = 8, height = 4)

# Core genes: present in all MAGs
core_genes <- rownames(gene_matrix)[rowSums(gene_matrix) == ncol(gene_matrix)]
length(core_genes)  # number of core genes

# Accessory genes: present in some but not all MAGs
accessory_genes <- rownames(gene_matrix)[rowSums(gene_matrix) < ncol(gene_matrix)]
length(accessory_genes)  # number of accessory genes

library(UpSetR)

# Convert matrix to list of genes per MAG
gene_list <- apply(gene_matrix, 2, function(x) rownames(gene_matrix)[x==1])

# Upset plot
pdf("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Panaroo_venn_for7.pdf", width = 8, height = 4)
upset(fromList(gene_list), nsets = ncol(gene_matrix), order.by = "freq")
dev.off()

# Compute Jaccard distance between MAGs
dist_mat <- dist(t(gene_matrix), method = "binary")

# Hierarchical clustering
hc <- hclust(dist_mat, method = "average")

# Plot dendrogram
pdf("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Dendogram_by_gene_content_panaroos.pdf", width = 8, height = 6)
plot(hc, main = "Clustering of MAGs by gene content", xlab = "", sub = "")
dev.off()


# Plot core genes vs accessory
library(dplyr)
library(ggplot2)

# gene_matrix: your 614x7 presence/absence table (0/1)
num_genomes <- ncol(gene_matrix)

# Count in how many genomes each gene is present
gene_matrix$presence_count <- rowSums(gene_matrix, na.rm = TRUE)

# Classify genes
gene_matrix$category <- ifelse(
  gene_matrix$presence_count == num_genomes, "Core",
  "Accessory"
)

# Optional: if you want, you can define "soft core" (present in >=70% genomes)
threshold <- ceiling(0.7 * num_genomes)
gene_matrix$category_soft <- ifelse(
  gene_matrix$presence_count >= threshold, "Core",
  "Accessory"
)
gene_summary <- gene_matrix %>%
  group_by(category) %>%
  summarise(gene_count = n())

gene_summary_soft <- gene_matrix %>%
  group_by(category_soft) %>%
  summarise(gene_count = n())
print(gene_summary)
print(gene_summary_soft)

# Strict core
ggplot(gene_summary, aes(x = category, y = gene_count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Panaroo Pangenome: Core vs Accessory Genes",
       x = "Gene Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Core" = "#1b9e77", "Accessory" = "#d95f02"))

#Soft core
ggplot(gene_summary_soft, aes(x = category_soft, y = gene_count, fill = category_soft)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Panaroo Pangenome (Soft Core ≥70% MAGs)",
       x = "Gene Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Core" = "#1b9e77", "Accessory" = "#d95f02"))
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Panaroo_core_accessory_genes.pdf", width = 6, height = 4)