# Load libraries
library(ape)        # For reading trees
library(ggtree)     # For visualization
library(tidyverse)  # For data manipulation
library(dplyr)

tree_file <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/phylophlan3_2/RAxML_bestTree.all_mags_refined.tre"
tree <- read.tree(tree_file)

metadata_file <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/phylophlan3/MAG_metadata.csv"
metadata <- read_csv(metadata_file)
metadata <- metadata[c(-10),]

# Clean MAG_ID and prepare metadata
metadata <- metadata %>%
  mutate(MAG_ID = gsub(".fa$", "", MAG_ID)) 


checkm_file <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/checkm2_out_all/quality_report.tsv"
checkm <- read_tsv(checkm_file, show_col_types = FALSE)

# Inspect
str(checkm)

checkm_sel <- checkm %>%
  select(
    MAG_ID = Name,
    Completeness,
    Contamination,
    Genome_Size,
    GC_Content,
    Total_Contigs,
    Contig_N50
  )

metadata <- metadata %>%
  left_join(checkm_sel, by = "MAG_ID")



metadata <- metadata %>%  # remove .fa
  filter(MAG_ID %in% tree$tip.label) %>%        # keep only tree tips
  arrange(match(MAG_ID, tree$tip.label)) %>%   # order as in tree
  rename(label = MAG_ID)                        # <- THIS IS THE KEY

# Now plot
p <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label = label), size = 3, align = TRUE) +
  geom_tippoint(aes(color = Group_cutoff_4), size = 4) +
  scale_color_manual(values = c("R" = "#009E73", "NR" = "#CC79A7")) +
  theme_tree2() +
  ggtitle("PhyloPhlAn Tree Annotated with Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

p
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Phylophlan_MAGs_dendogram.pdf", plot = p, width = 10, height = 6)

pcircular <- ggtree(tree, layout = "circular") %<+% metadata +
  geom_tippoint(aes(color = Group_cutoff_4), size = 4) +
  scale_color_manual(values = c("R" = "#009E73", "NR" = "#CC79A7")) +
  theme_void() +
  ggtitle("Phylophlan Tree Annotated with Group") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )
pcircular
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Phylophlan_MAGs_circular.pdf", plot = p, width = 10, height = 6)
