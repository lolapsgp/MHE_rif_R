library(tidyverse)

# Directory where all per-sample counts are stored
results_dir <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/rpoB_variant"
output_file <- file.path(results_dir, "_full_sequence_counts.tsv")

files <- list.files(
  path = results_dir, pattern = "_full_sequence_counts.txt$",
  recursive = TRUE,
  full.names = TRUE
)
cat("Found", length(files), "abundance files\n")
abundance_list <- lapply(files, function(f) {
  df <- read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  # Safety check
  if (!"Sample" %in% colnames(df)) {
    stop(paste("Missing 'Sample' column in", f))
  }
  
  return(df)
})

merged_abundance <- do.call(rbind, abundance_list)

# SORT 
merged_abundance <- merged_abundance[order(merged_abundance$Sample), ]


# WRITE OUTPUT
write.table(
  merged_abundance,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Merged table written to:", output_file, "\n")

merged_abundance <- merged_abundance %>%
  mutate(Sample = if_else(str_starts(Sample, "RS_") & !str_ends(Sample, "_ST"),
                          paste0(Sample, "_ST"),
                          Sample))

ps_met_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_met_corrected.Rds")

merged_abundance_filt <- merged_abundance

merged_abundance_filt$rpoB_rpk<-merged_abundance_filt$RPK
merged_abundance_filt$SampleID<-merged_abundance_filt$Sample
merged_abundance_filt$RPK<-NULL
merged_abundance_filt$Sample<-NULL

metadata<- data.frame(sample_data(ps_met_corrected)) 

data <- merge(merged_abundance_filt, metadata,
              by = "SampleID",
              all = TRUE)
data_plot <- data[, c("Seq_A_Count", "Seq_C_Count", "Group_cutoff_4", "Timepoint", "Study", "Volunteer")]


data_plot$Group_Timepoint <- paste(
  data_plot$Group_cutoff_4,
  data_plot$Timepoint,
  sep = "_"
)

library(ggplot2)
library(dplyr)
library(tidyr)

# Keep only samples with Seq_A_Count > 0 OR Seq_C_Count > 0
data_filtered <- data_plot %>% 
  filter(Seq_A_Count > 0 | Seq_C_Count > 0)

# Convert to long format for plotting
data_long <- data_filtered %>%
  pivot_longer(cols = c(Seq_A_Count, Seq_C_Count),
               names_to = "Allele",
               values_to = "Count") %>%
  mutate(Allele = recode(Allele, 
                         Seq_A_Count = "N mutation",
                         Seq_C_Count = "H reference"))

# Plot grouped bar chart per patient and timepoint
library(ggplot2)
library(ggtext)  # needed for element_markdown()

p <- ggplot(data_long, aes(x = Timepoint, y = Count, fill = Allele)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Volunteer, scales = "free_y") +
  scale_fill_manual(values = c("N mutation" = "#D55E00", "H reference" = "#999999")) +
  theme_minimal() +
  labs(title = "*rpoB* allele counts in *Segatella copri* (patients with either allele)",
       x = "Timepoint",
       y = "Read counts",
       fill = "Allele") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_markdown())  # <-- render markdown in title

# Save the plot
ggsave(filename = file.path("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Specific_motus_and_GMMs", 
                            "rpoB_allele_counts.pdf"),
       plot = p,
       width = 6, height = 6)