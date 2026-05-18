#!/usr/bin/env Rscript

# ---------------------------
# USER INPUT
# ---------------------------
results_dir <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/cfxA2_diamond"
output_file <- file.path(results_dir, "cfxA2_merged_abundance.tsv")

# ---------------------------
# FIND FILES RECURSIVELY
# ---------------------------
files <- list.files(
  path = results_dir,
  pattern = "_cfxA2_abundance\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(files), "abundance files\n")

# Print them so we can verify
print(files)

if (length(files) == 0) {
  stop("ERROR: No abundance files found. Check results_dir.")
}

# ---------------------------
# READ AND MERGE
# ---------------------------
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

ps_met_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_met_corrected.Rds")

merged_abundance$Sample <- sub("filtered$", "", merged_abundance$Sample)

merged_abundance_filt <- merged_abundance[
  merged_abundance$Sample %in% sample_names(ps_met_corrected),
]

merged_abundance_filt$cfxA2_rpk<-merged_abundance_filt$RPK
merged_abundance_filt$SampleID<-merged_abundance_filt$Sample
merged_abundance_filt$RPK<-NULL
merged_abundance_filt$Sample<-NULL

metadata<- data.frame(sample_data(ps_met_corrected)) 

data <- merge(merged_abundance_filt, metadata,
              by = "SampleID",
              all = TRUE)
data_plot <- data[, c("cfxA2_rpk", "Group_cutoff_4", "Timepoint", "Study", "Volunteer")]


data_plot$Group_Timepoint <- paste(
  data_plot$Group_cutoff_4,
  data_plot$Timepoint,
  sep = "_"
)


library(emmeans)
library(lmerTest)
model<-lmer(cfxA2_rpk ~ Timepoint * Group_cutoff_4 + (1 | Volunteer), data = data_plot)
anova(model)
emm <- emmeans(model, ~ Group_cutoff_4 )
pairs(emm)

# Reorder the x-axis
data_plot$Group_Timepoint <- factor(
  data_plot$Group_Timepoint,
  levels = c("NR_T0", "NR_T2", "R_T0", "R_T2")
)


ggplot(data_plot, aes(x=Group_Timepoint , y = cfxA2_rpk, fill = Group_cutoff_4)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = Study), alpha = 0.6,
              position = position_jitter(width = 0.2)) +
  geom_pwc(label = "{p.adj}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_fill_manual(values = c("R" = "#009E73", "NR" = "#CC79A7"),
                    name = "Group cutoff 4") +
  scale_color_manual(values = c("Spain" = "#E69F00", "UK" = "#0072B2"),
                     name = "Cohort") +
  labs(x = "Group", y = "rpk-transformed cfxA2 abundance") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_signif(
    comparisons = list(c("R_T0", "R_T2"), c("NR_T0", "NR_T2")),
    annotations = c("*", "."),
    y_position = c(max(data_plot$cfxA2_rpk) +1),
    tip_length = 0.03,
    textsize = 5
  )
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Specific_motus_and_GMMs/cfxA2_Group_Timepoint.pdf", height = 4, width = 4)



ggplot(data_plot, aes(x=Group_cutoff_4 , y = cfxA2_rpk, fill = Group_cutoff_4)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = Study), alpha = 0.6,
              position = position_jitter(width = 0.2)) +
  geom_pwc(label = "{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_fill_manual(values = c("R" =  "#009E73", "NR" = "#CC79A7"),
                    name = "Group_cutoff_4") +
  scale_color_manual(values = c("Spain" = "#E69F00", "UK" = "#0072B2"),
                     name = "Cohort") +
  labs(x = "Group", y = "rpk-transformed cfxA2 abundance") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Figures_merged/Specific_motus_and_GMMs/cfxA2_Group.pdf", height = 4, width = 4)


