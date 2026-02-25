####################### Use K0 raw abundance ##################################
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
metadata<-as.data.frame(sample_data(GMMs_corrected))

#Import NGLess outputs
KEGGs<-read.table(file = '/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs/KEGG_fpkm.profiles.txt',
                  header = TRUE, 
                  sep = "\t", 
                  row.names = 1, 
                  skip = 0, 
                  check.names = FALSE,
                  fill = TRUE,
                  quote = "")
# Clean "ko:" from rownames
rownames(KEGGs) <- sub("^ko:", "", rownames(KEGGs))
KEGGs<-KEGGs[-1,]
KEGGs_asnA<-KEGGs[c("K01914"),]
KEGGsUK<-read.table(file = '/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Patel_2021_09_10/KEGG_fpkm.profiles.txt',
                    header = TRUE, 
                    sep = "\t", 
                    row.names = 1, 
                    skip = 0, 
                    check.names = FALSE,
                    fill = TRUE,
                    quote = "")
# Clean "ko:" from rownames
rownames(KEGGsUK) <- sub("^ko:", "", rownames(KEGGsUK))
KEGGsUK<-KEGGsUK[-1,]
KEGGs_asnAUK<-KEGGsUK[c("K01914"),]

#merge
KEGGs_asnA_merged <- cbind(KEGGs_asnA, KEGGs_asnAUK)
cohort <- c(
  rep("Spain", ncol(KEGGs_asnA)),
  rep("UK", ncol(KEGGs_asnAUK))
)

colData <- data.frame(
  sample = colnames(KEGGs_asnA_merged),
  cohort = cohort
)

library(tidyr)
library(dplyr)

KEGGs_long <- KEGGs_asnA_merged %>%
  pivot_longer(cols = everything(),
               names_to = "sample",
               values_to = "asnA_count") %>%
  mutate(cohort = ifelse(sample %in% colnames(KEGGs_asnA),
                         "Spain", "UK"))
KEGGs_long$SampleID <- gsub("filtered", "", KEGGs_long$sample)
KEGGs_long$sample<-NULL
KEGGs_long_filtered <- KEGGs_long[match(metadata$SampleID, KEGGs_long$SampleID), ]



library(dplyr)

final_df <- metadata %>%
  left_join(KEGGs_long_filtered, by = c("SampleID"))

library(emmeans)
library(lmerTest)

final_df<-data.frame(final_df)
model<-lmer(asnA_count ~ Timepoint * Group_cutoff_4 + (1 | Volunteer), data = final_df)
anova(model)
emm <- emmeans(model, ~ Timepoint | Group_cutoff_4 )
pairs(emm)


final_df$Group_Timepoint <- paste(
  final_df$Group_cutoff_4,
  final_df$Timepoint,
  sep = "_"
)

# Reorder the x-axis
final_df$Group_Timepoint <- factor(
  final_df$Group_Timepoint,
  levels = c("NR_T0", "NR_T2", "R_T0", "R_T2")
)

ggplot(final_df, aes(x=Group_Timepoint , y = asnA_count, fill = Group_cutoff_4)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = Study), alpha = 0.6,
              position = position_jitter(width = 0.2))+
  scale_fill_manual(values = c("R" = "#009E73", "NR" = "#CC79A7"),
                    name = "Group cutoff 4") +
  scale_color_manual(values = c("Spain" = "#E69F00", "UK" = "#0072B2"),
                     name = "Cohort") +
  labs(x = "Group", y = "rpk-transformed asnA abundance") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))


####################### Use K0 corrected abundance ##################################
GMM_corrected<-readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected_asnA.Rds")
metadata<-as.data.frame(sample_data(GMMs_corrected))
otu<-data.frame(otu_table(GMMs_corrected))
otu<-otu[c("MF0150"),]
asnA<-data.frame(t(otu))
asnA$SampleID<-rownames(asnA)
final_df <- metadata %>%
  left_join(asnA, by = c("SampleID"))


final_df<-data.frame(final_df)
model<-lmer(asnA_count ~ Timepoint * Group_cutoff_4 + (1 | Volunteer), data = final_df)
anova(model)
emm <- emmeans(model, ~ Timepoint | Group_cutoff_4 )
pairs(emm)
emm <- emmeans(model, ~ Group_cutoff_4 | Timepoint  )
pairs(emm)



final_df$Group_Timepoint <- paste(
  final_df$Group_cutoff_4,
  final_df$Timepoint,
  sep = "_"
)

# Reorder the x-axis
final_df$Group_Timepoint <- factor(
  final_df$Group_Timepoint,
  levels = c("NR_T0", "NR_T2", "R_T0", "R_T2")
)

ggplot(final_df, aes(x=Group_Timepoint , y = asnA_count, fill = Group_cutoff_4)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = Study), alpha = 0.6,
              position = position_jitter(width = 0.2))+
  scale_fill_manual(values = c("R" = "#009E73", "NR" = "#CC79A7"),
                    name = "Group cutoff 4") +
  scale_color_manual(values = c("Spain" = "#E69F00", "UK" = "#0072B2"),
                     name = "Cohort") +
  labs(x = "Group", y = "rpk-transformed asnA abundance") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))