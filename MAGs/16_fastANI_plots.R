library(tidyverse)
library(reshape2)
library(patchwork)

metadata_file <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/phylophlan3/MAG_metadata.csv"
metadata <- read_csv(metadata_file)
metadata <- metadata[c(-10),]

# Clean MAG_ID and prepare metadata
metadata <- metadata %>%
  mutate(MAG_ID = gsub(".fa$", "", MAG_ID)) 
ani <- read.table(
  "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/MAGs/Prevotella_copri/RandNRonly/ANI_all_vs_all.txt",
  header = FALSE
)

colnames(ani) <- c("Genome1","Genome2","ANI","Fragments","Total")

ani <- ani %>%
  mutate(
    Genome1 = gsub(".fa","",Genome1),
    Genome2 = gsub(".fa","",Genome2)
  )
ani_rev <- ani %>%
  select(Genome1 = Genome2,
         Genome2 = Genome1,
         ANI)

ani_full <- bind_rows(
  ani %>% select(Genome1, Genome2, ANI),
  ani_rev
)

genomes <- unique(c(ani_full$Genome1, ani_full$Genome2))

self <- data.frame(
  Genome1 = genomes,
  Genome2 = genomes,
  ANI = 100
)

ani_full <- bind_rows(ani_full, self)
ani_matrix <- acast(ani_full, Genome1 ~ Genome2, value.var="ANI")

dist_matrix <- dist(100 - ani_matrix)

hc <- hclust(dist_matrix)

genome_order <- hc$labels[hc$order]
ani_full$Genome1 <- factor(ani_full$Genome1, levels=genome_order)
ani_full$Genome2 <- factor(ani_full$Genome2, levels=genome_order)

metadata <- metadata %>%
  mutate(MAG_ID = factor(MAG_ID, levels=genome_order))

library(tidyverse)
library(reshape2)
library(cowplot)

ani_sym <- ani_full %>%
  group_by(Genome1, Genome2) %>%
  summarise(ANI = max(ANI), .groups="drop")

#--- Prepare heatmap ---#
heatmap_plot <- ggplot(ani_sym, aes(Genome1, Genome2, fill=ANI)) +
  geom_tile(color="white") +
  scale_fill_gradientn(
    colours=c("white","yellow","orange","red"),
    limits=c(90,100),
    name="ANI (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

#--- Row annotation (left) with labels ---#
annotation_y <- ggplot(metadata, aes(x=1, y=MAG_ID, fill=Group_cutoff_4)) +
  geom_tile(alpha = 0.5) +
  scale_fill_manual(
    values=c(R="#009E73", NR="#CC79A7"),
    name="Group"
  ) +
  theme_void() +
  theme(
    legend.position="left",
    axis.text.y = element_text(hjust=1)
  )

#--- Column annotation (bottom) with labels ---#
annotation_x <- ggplot(metadata, aes(x=MAG_ID, y=1, fill=Group_cutoff_4)) +
  geom_tile(alpha = 0.5) +
  scale_fill_manual(
    values=c(R="#009E73", NR="#CC79A7"),
    guide = "none"   # legend only on left
  ) +
  theme_void() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text.y = element_blank()
  )

#--- Blank plot for top-left corner ---#
empty <- ggplot() + theme_void()

#--- Combine with cowplot::plot_grid ---#
# layout:
#   top row: empty | x annotation
#   bottom row: y annotation | heatmap
final_plot <- plot_grid(annotation_y, heatmap_plot, rel_widths=c(0.5,1), align = "h", axis = "lr", ncol = 2, rel_heights = c(1.9,1))
final_plot
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Fast_any_cont.pdf", plot = final_plot, width = 8, height = 4)
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Fast_any_cont.svg", plot = final_plot, width = 8, height = 4)


# Create categorical column based on ANI thresholds
ani_sym <- ani_sym %>%
  mutate(
    ANI_category = case_when(
      ANI > 95 ~ "Same strain (ANI>95%)",
      ANI >= 90 & ANI <= 95 ~ "Same species (95%<ANI>90%)",
      ANI < 90 ~ "Different species (90%>ANI)"
    )
  )

# Heatmap with categorical colors and legend
heatmap_plot <- ggplot(ani_sym, aes(Genome1, Genome2, fill = ANI_category)) +
  geom_tile(color="white") +
  scale_fill_manual(
    values = c(
      "Same strain (ANI>95%)" = "orange",
      "Same species (95%<ANI>90%)" = "lightblue",
      "Different species (90%>ANI)" = "lightgray"
    ),
    name = "ANI category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
final_plot <- plot_grid(annotation_y, heatmap_plot, rel_widths=c(0.5,1), align = "h", axis = "lr", ncol = 2, rel_heights = c(1.9,1))
final_plot

ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Fast_any_classes.pdf", plot = final_plot, width = 8, height = 4)
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Phylogeny/Fast_any_classes.svg", plot = final_plot, width = 8, height = 4)


