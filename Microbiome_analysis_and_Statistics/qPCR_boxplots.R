ps_clr_corrected <- readRDS("~/MHE_rif/outputs/merged/ps_clr_corrected.Rds")
ps_merged_raw <- readRDS("~/MHE_rif/outputs/merged/ps_merged_raw.Rds")

######################### qPCR results 16S #############################

library(readxl)

metadata <- as.data.frame(sample_data(ps_clr_corrected))
qPCR_results <- read_excel("~/Segatella_copri/Calculo_num_copias.xlsx", 
                                 sheet = "16S", range = "A1:H30")
qPCR_results<-data.frame(qPCR_results)

library(ggpubr)
library(ggplot2)
library(dplyr)
qPCR_results <- qPCR_results %>%
  filter(Volunteer != "CN")

p16S_Group<-ggplot(qPCR_results, aes(x = Group_cutoff_4, y = Copies_num_16S, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("p16S_Group.pdf", plot = p16S_Group, device = "pdf",   width = 4,
       height = 6)

eps <- 0.5  
p16S_sacaled_Group<-ggplot(qPCR_results, aes(x = Group_cutoff_4, y = Copies_num_16S + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_y_log10() +
  labs(y = expression(paste("Copies 16S ", italic("Segatella copri"), " (log10 scale)")))+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))


ggsave("p16S_sacaled_Group.pdf", plot = p16S_sacaled_Group, device = "pdf",   width = 4,
       height = 6)

qPCR_results <- qPCR_results %>%
  filter(Timepoint != "T11")

p16S_Timepoint<-ggplot(qPCR_results, aes(x = Timepoint, y = Copies_num_16S, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("p16S_Timepoint.pdf", plot = p16S_Timepoint, device = "pdf",   width = 4,
       height = 6)

eps <- 0.5  
p16S_sacaled_Timepoint<-ggplot(qPCR_results, aes(x = Timepoint, y = Copies_num_16S + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_y_log10() +
  labs(y = expression(paste("Copies 16S ", italic("Segatella copri"), " (log10 scale)")))+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))


ggsave("p16S_sacaled_Timepoint.pdf", plot = p16S_sacaled_Timepoint, device = "pdf",   width = 4,
       height = 6)





qPCR_results$Group_Timepoint <- factor(
  paste(qPCR_results$Group_cutoff_4, qPCR_results$Timepoint, sep = "_"),
  levels = c("R_T0", "R_T1", "R_T2", "NR_T0", "NR_T1", "NR_T2")
)
qPCR_results$Group_Timepoint <- factor(
  qPCR_results$Group_Timepoint,
  levels = c("R_T0", "R_T1", "R_T2", "NR_T0", "NR_T1", "NR_T2")
)

p16S_ordered_Timepoint<-ggplot(
  qPCR_results,
  aes(x = Group_Timepoint, y = Copies_num_16S + eps, fill = Group_cutoff_4)
) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90)
  ) +
  geom_pwc(
    label = "{p.adj}{p.adj.signif}",
    hide.ns = TRUE,
    method = "wilcox_test",
    p.adjust.method = "BH"
  ) +
  scale_y_log10() +
  labs(y = expression(paste("Copies primer 3 (", italic("Segatella copri"), " log10 scale)"))) +
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4") +
  geom_jitter(
    aes(group = Group_cutoff_4),
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  )

ggsave("p16S_ordered_Timepoint.pdf", plot = p16S_ordered_Timepoint, device = "pdf",   width = 4,
       height = 6)


qPCR_results <- qPCR_results %>%
  filter(Timepoint != "T2")

p16S_ordered_Timepoint_filt.pdf<-ggplot(
  qPCR_results,
  aes(x = Group_Timepoint, y = Copies_num_16S + eps, fill = Group_cutoff_4)
) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90)
  ) +
  geom_pwc(
    label = "{p.adj}{p.adj.signif}",
    hide.ns = TRUE,
    method = "wilcox_test",
    p.adjust.method = "BH"
  ) +
  scale_y_log10() +
  labs(y = expression(paste("Copies 16S (", italic("Segatella copri"), " log10 scale)"))) +
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4") +
  geom_jitter(
    aes(group = Group_cutoff_4),
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  )

ggsave("p16S_ordered_Timepoint_filt.pdf", plot = p16S_ordered_Timepoint_filt, device = "pdf",   width = 4,
       height = 6)

######################### qPCR results asnA #############################

library(readxl)

qPCR_results <- read_excel("~/Segatella_copri/Calculo_num_copias.xlsx", 
                           sheet = "asnA", range = "A1:H30")
qPCR_results<-data.frame(qPCR_results)

library(ggpubr)
library(ggplot2)
library(dplyr)
qPCR_results <- qPCR_results %>%
  filter(Volunteer != "CN")


plot1_Group_asna<-ggplot(qPCR_results, aes(x = Group_cutoff_4, y = Copies_num_asnA, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE, method = "wilcox_test", p.adjust.method = "BH")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("pasnA_Group.pdf", plot = plot1_Group_asna, device = "pdf",   width = 4,
       height = 6)


eps <- 0.5  
plot_sacaled_Group_asnA<-ggplot(qPCR_results, aes(x = Group_cutoff_4, y = Copies_num_asnA + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_y_log10() +
  labs(y = "Copies asnA (log10 scale)")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("pasnA_scaled_Group.pdf", plot = plot_sacaled_Group_asnA, device = "pdf",   width = 4,
       height = 6.5)


qPCR_results <- qPCR_results %>%
  filter(Timepoint != "T11")

plot1_Timepoint_asna<-ggplot(qPCR_results, aes(x = Timepoint, y = Copies_num_asnA + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  labs(y = "Copies asnA")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("pasnA_Timepoint.pdf", plot = plot1_Timepoint_asna, device = "pdf",   width = 4,
       height = 6.5)


eps <- 0.5  
plot_sacaled_Timepoint_asnA<-ggplot(qPCR_results, aes(x = Timepoint, y = Copies_num_asnA + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_y_log10() +
  labs(y = "Copies asnA (log10 scale)")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("pasnA_sacaled_Timepoint.pdf", plot = plot_sacaled_Timepoint_asnA, device = "pdf",   width = 4,
       height = 6)


qPCR_results$Group_Timepoint <- factor(
  paste(qPCR_results$Group_cutoff_4, qPCR_results$Timepoint, sep = "_"),
  levels = c("R_T0", "R_T1", "R_T2", "NR_T0", "NR_T1", "NR_T2")
)
qPCR_results$Group_Timepoint <- factor(
  qPCR_results$Group_Timepoint,
  levels = c("R_T0", "R_T1", "R_T2", "NR_T0", "NR_T1", "NR_T2")
)

pasnA_ordered_Timepoint.pdf<- ggplot(
  qPCR_results,
  aes(x = Group_Timepoint, y = Copies_num_asnA + eps, fill = Group_cutoff_4)
) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90)
  ) +
  geom_pwc(
    label = "{p.adj}{p.adj.signif}",
    hide.ns = TRUE,
    method = "wilcox_test",
    p.adjust.method = "BH"
  ) +
  scale_y_log10() +
  labs(y = expression(paste("Copies asnA (", italic("Segatella copri"), " log10 scale)"))) +
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4") +
  geom_jitter(
    aes(group = Group_cutoff_4),
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  )

ggsave("pasnA_ordered_Timepoint.pdf", plot = pasnA_ordered_Timepoint, device = "pdf",   width = 4,
       height = 6)


qPCR_results <- qPCR_results %>%
  filter(Timepoint != "T2")

pasnA_ordered_Timepoint_filt<-ggplot(
  qPCR_results,
  aes(x = Group_Timepoint, y = Copies_num_asnA + eps, fill = Group_cutoff_4)
) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90)
  ) +
  geom_pwc(
    label = "{p.adj}{p.adj.signif}",
    hide.ns = TRUE,
    method = "wilcox_test",
    p.adjust.method = "BH"
  ) +
  scale_y_log10() +
  labs(y = expression(paste("Copies num asnA (", italic("Segatella copri"), " log10 scale)"))) +
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4") +
  geom_jitter(
    aes(group = Group_cutoff_4),
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  )

ggsave("pasnA_ordered_Timepoint_filt.pdf", plot = pasnA_ordered_Timepoint_filt, device = "pdf",   width = 4,
       height = 6)


######################### qPCR results primer3 #############################
qPCR_results <- read_excel("~/Segatella_copri/Calculo_num_copias.xlsx", 
                           sheet = "primer3", range = "A1:H30")
qPCR_results<-data.frame(qPCR_results)

library(ggpubr)
library(ggplot2)
library(dplyr)
qPCR_results <- qPCR_results %>%
  filter(Volunteer != "CN")

p3_Group<-ggplot(qPCR_results, aes(x = Group_cutoff_4, y = Copies_num_primer3, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))

ggsave("p3_Group.pdf", plot = p3_Group, device = "pdf",   width = 4,
       height = 6)

eps <- 0.5  
p3_sacaled_Group<-ggplot(qPCR_results, aes(x = Group_cutoff_4, y = Copies_num_primer3 + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_y_log10() +
  labs(y = expression(paste("Copies primer 3 (", italic("Segatella copri"), " log10 scale)")))+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))


ggsave("p3_sacaled_Group.pdf", plot = p3_sacaled_Group, device = "pdf",   width = 4,
       height = 6)

qPCR_results <- qPCR_results %>%
  filter(Timepoint != "T11")

p3_Timepoint<-ggplot(qPCR_results, aes(x = Timepoint, y = Copies_num_primer3, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH")+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")

ggsave("p3_Timepoint.pdf", plot = p3_Timepoint, device = "pdf",   width = 4,
       height = 6)

eps <- 0.5  
p3_scaled_Timepoint<-ggplot(qPCR_results, aes(x = Timepoint, y = Copies_num_primer3 + eps, fill = Group_cutoff_4)) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE,
           method = "wilcox_test", p.adjust.method = "BH") +
  scale_y_log10() +
  labs(y = expression(paste("Copies primer 3 (", italic("Segatella copri"), " log10 scale)")))+
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4")+
  geom_jitter(aes(group = Group_cutoff_4), alpha =0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))


ggsave("p3_scaled_Timepoint.pdf", plot = p3_scaled_Timepoint, device = "pdf",   width = 4,
       height = 6)


qPCR_results$Group_Timepoint <- factor(
  paste(qPCR_results$Group_cutoff_4, qPCR_results$Timepoint, sep = "_"),
  levels = c("R_T0", "R_T1", "R_T2", "NR_T0", "NR_T1", "NR_T2")
)
qPCR_results$Group_Timepoint <- factor(
  qPCR_results$Group_Timepoint,
  levels = c("R_T0", "R_T1", "R_T2", "NR_T0", "NR_T1", "NR_T2")
)

ggplot(
  qPCR_results,
  aes(x = Group_Timepoint, y = Copies_num_primer3 + eps, fill = Group_cutoff_4)
) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90)
  ) +
  geom_pwc(
    label = "{p.adj}{p.adj.signif}",
    hide.ns = TRUE,
    method = "wilcox_test",
    p.adjust.method = "BH"
  ) +
  scale_y_log10() +
  labs(y = expression(paste("Copies primer 3 (", italic("Segatella copri"), " log10 scale)"))) +
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4") +
  geom_jitter(
    aes(group = Group_cutoff_4),
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  )

ggsave("p3_ordered_Timepoint.pdf", plot = p3_scaled_Timepoint, device = "pdf",   width = 4,
       height = 6)

qPCR_results <- qPCR_results %>%
  filter(Timepoint != "T2")

ggplot(
  qPCR_results,
  aes(x = Group_Timepoint, y = Copies_num_primer3 + eps, fill = Group_cutoff_4)
) +
  geom_boxplot(position = position_dodge(), alpha = 0.3, outliers = FALSE) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90)
  ) +
  geom_pwc(
    label = "{p.adj}{p.adj.signif}",
    hide.ns = TRUE,
    method = "wilcox_test",
    p.adjust.method = "BH"
  ) +
  scale_y_log10() +
  labs(y = expression(paste("Copies primer 3 (", italic("Segatella copri"), " log10 scale)"))) +
  scale_fill_manual(values = c("NR" = "#CC79A7", "R" = "#009E73"), name = "Group cutoff 4") +
  geom_jitter(
    aes(group = Group_cutoff_4),
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  )

ggsave("p3_ordered_Timepoint_filt.pdf", plot = p3_scaled_Timepoint, device = "pdf",   width = 4,
       height = 6)





