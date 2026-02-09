# Load libraries
library(ape)        # For reading trees
library(ggtree)     # For visualization
library(tidyverse)  # For data manipulation
library(dplyr)

tree_file <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/strainphlan/strainphlan_out/RAxML_bestTree.t__SGB1626.StrainPhlAn4.tre"
tree <- read.tree(tree_file)

ps_clr_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_clr_corrected.Rds")
metadata<-data.frame(sample_data(ps_clr_corrected))
rm(ps_clr_corrected)


metadata <- metadata %>%  # remove .fa
  filter(SampleID %in% tree$tip.label) %>%        # keep only tree tips
  arrange(match(SampleID, tree$tip.label)) %>%   # order as in tree
  rename(label = SampleID)                        # <- THIS IS THE KEY

# Now plot
p <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label = label), size = 3, align = TRUE) +
  geom_tippoint(aes(color = Group_cutoff_4), size = 4) +
  scale_color_manual(values = c("R" = "#009E73", "NR" = "#CC79A7")) +
  theme_tree2() +
  ggtitle("StrainPhlAn4 Tree Annotated with Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

p

ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/StrainPhlAn4_dendogram.pdf", plot = p, width = 10, height = 6)


pcircular <- ggtree(tree, layout = "circular") %<+% metadata +
  geom_tippoint(aes(color = Group_cutoff_4), size = 4) +
  scale_color_manual(values = c("R" = "#009E73", "NR" = "#CC79A7")) +
  theme_void() +
  ggtitle("StrainPhlAn4 Tree Annotated with Group") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )

pcircular

ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/StrainPhlAn4_circular.pdf", plot = p, width = 10, height = 6)
