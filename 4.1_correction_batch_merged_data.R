library(vegan)
library(permute)
library(phyloseq)
library(microbiome)

####################################### UK and SPAIN ARSyN#########################################################
####################################### Taxonomy #########################################################
ps_merged_filtered <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_merged_filtered.Rds")
ps_merged_filtered<- subset_samples(ps_merged_filtered, Group_cutoff_4 %in% c("R", "NR"))
ps_merged_filtered<- subset_samples(ps_merged_filtered, Timepoint %in% c("T0", "T2"))

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_merged_long<- subset_samples(ps_merged_filtered, Timepoint %in% c("T0"))
ps_merged_long<- subset_samples(ps_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(ps_merged_long)
meta_tax1 <- meta(ps_merged_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_merged_long<- subset_samples(ps_merged_filtered, Timepoint %in% c("T2"))
ps_merged_long<- subset_samples(ps_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(ps_merged_long)
meta_tax2 <- meta(ps_merged_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")



#Remove batch effect taxonomy data and then apply clr transformation
#ARSyN
#Indicating group data
otu <- data.frame(otu_table(ps_merged_filtered))
metadata<-data.frame(sample_data(ps_merged_filtered))
rownames(metadata) <- gsub("-", ".", rownames(metadata))
metadata<-metadata[order(rownames(metadata)), ]
otu<-otu[,order(colnames(otu))]
all(rownames(metadata) == colnames(otu))

library(MultiBaC)

# transformation
# Batch correction
my_mbac <- createMbac (inputOmics = list(otu_S = otu[,metadata$Study == "Spain"],
                                         otu_UK = otu[,metadata$Study == "UK"]),
                       batchFactor = c("Spain", "UK"),
                       experimentalDesign = list("Spain" = metadata$Group_cutoff_4[metadata$Study == "Spain"],
                                                 "UK" = metadata$Group_cutoff_4[metadata$Study == "UK"]),
                       omicNames = "otu")
arsyn_1 <- ARSyNbac(my_mbac, modelName = "metadata", Variability = 0.95, 
                    batchEstimation = TRUE, Interaction = FALSE, beta=2)     

otu_corrected <- t(cbind(MultiAssayExperiment::assay(arsyn_1$CorrectedData$Spain),
                         MultiAssayExperiment::assay(arsyn_1$CorrectedData$UK)))
OTU = otu_table(t(otu_corrected), taxa_are_rows = TRUE)
TAX = tax_table(ps_merged_filtered)
SAMPLE <- sample_data(metadata)

ps_corrected <- phyloseq(OTU, TAX, SAMPLE)


saveRDS(ps_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_corrected.Rds")


# Euclidean distance
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_corrected_T0<- subset_samples(ps_corrected, Timepoint %in% c("T0"))
ps_corrected_T0<- subset_samples(ps_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_c0 <- abundances(ps_corrected_T0)
meta_tax_c0 <- meta(ps_corrected_T0)

# PERMANOVA Euclidean distance T0
adonis2(formula = t(otu_tax_c0) ~ Group_cutoff_4 + Study, 
        data = meta_tax_c0, permutations = 999, 
        method = "euclidean", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_corrected_T2<- subset_samples(ps_corrected, Timepoint %in% c("T2"))
ps_corrected_T2<- subset_samples(ps_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_c2 <- abundances(ps_corrected_T2)
meta_tax_c2 <- meta(ps_corrected_T2)

# PERMANOVA Euclidean distance T0
adonis2(formula = t(otu_tax_c2) ~ Group_cutoff_4 + Study, 
        data = meta_tax_c2, permutations = 999, 
        method = "euclidean", by = "margin")

# Euclidean distance because the data was clr-transformed
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_clr_corrected_T0<- subset_samples(ps_clr_corrected, Timepoint %in% c("T0"))
ps_clr_corrected_T0<- subset_samples(ps_clr_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_clr0 <- abundances(ps_clr_corrected_T0)
meta_tax_clr0 <- meta(ps_clr_corrected_T0)

# PERMANOVA euclidean distance T0
adonis2(formula = t(otu_tax_clr0) ~ Group_cutoff_4 + Study, 
        data = meta_tax_clr0, permutations = 999, 
        method = "euclidean", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_clr_corrected_T2<- subset_samples(ps_clr_corrected, Timepoint %in% c("T2"))
ps_clr_corrected_T2<- subset_samples(ps_clr_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_clr2 <- abundances(ps_clr_corrected_T2)
meta_tax_clr2 <- meta(ps_clr_corrected_T2)

# PERMANOVA euclidean distance T2
adonis2(formula = t(otu_tax_clr2) ~ Group_cutoff_4 + Study, 
        data = meta_tax_clr2, permutations = 999, 
        method = "euclidean", by = "margin")

dist = phyloseq::distance(ps_clr_corrected, method="euclidean")
ordination = ordinate(ps_clr_corrected, method="PCoA", distance=dist)
plot_ordination(ps_clr_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)



####################################### GMMs #########################################################
GMMs_merged<- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_merged.Rds")
GMMs_comp<- microbiome::transform(GMMs_merged, "compositional")
otu_gmm <- abundances(GMMs_comp)
meta_gmm <- meta(GMMs_comp)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_comp_long<- subset_samples(GMMs_comp, Timepoint %in% c("T0"))
GMMs_comp_long<- subset_samples(GMMs_comp_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(GMMs_comp_long)
meta_tax1 <- meta(GMMs_comp_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_comp_long<- subset_samples(GMMs_comp, Timepoint %in% c("T2"))
GMMs_comp_long<- subset_samples(GMMs_comp_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(GMMs_comp_long)
meta_tax2 <- meta(GMMs_comp_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")


# Remove batch effect indicating group data
set.seed(321)
otu_fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = otu_gmm,
  covariates = c("Group_cutoff_4", "Age", "Timepoint"),
  batch = "Study",  # Name of the column in meta_gmm
  data = meta_gmm,
  control = list(verbose = TRUE, diagnostic_plot = "adjust_batch_diagnostic.pdf")
)

GMMs_corrected<-GMMs_comp
otu_table(GMMs_corrected)<-otu_table(as.matrix(otu_fit_adjust_batch$feature_abd_adj), 
                                     taxa_are_rows = TRUE)
saveRDS(GMMs_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")


#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_corrected_long<- subset_samples(GMMs_corrected, Timepoint %in% c("T0"))
GMMs_corrected_long<- subset_samples(GMMs_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(GMMs_corrected_long)
meta_tax1 <- meta(GMMs_corrected_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_corrected_long<- subset_samples(GMMs_corrected, Timepoint %in% c("T2"))
GMMs_corrected_long<- subset_samples(GMMs_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(GMMs_corrected_long)
meta_tax2 <- meta(GMMs_corrected_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")

dist = phyloseq::distance(GMMs_corrected, method="bray")
ordination = ordinate(GMMs_corrected, method="PCoA", distance=dist)
plot_ordination(GMMs_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)



#PERMANOVA adjusting permutation scheme to see if we have differences between Timepoints
################## Confirm nesting is valid - no- so we use timepoints ############################################
table(meta_tax$Study, meta_tax$Volunteer)  # Should have only one Study per Volunteer

# Number of samples
n_samples <- nrow(meta_tax)
# Number of unique Volunteers
n_volunteers <- length(unique(meta_tax$Volunteer))

# Create the permutation scheme for later use in all comparisons
perm_scheme <- how(
  blocks = meta_tax$Study,              # Study is the outer block
  within = Within(type = "free"),       # Permute volunteers freely within each study
  plots = Plots(strata = meta_tax$Volunteer),
  nperm = 999
)
# Bray-Curtis distance
dist_mat <- vegdist(t(otu_tax), method = "bray")

# Run PERMANOVA with custom permutations
adonis2(
  formula = dist_mat ~ Group_cutoff_4 + Study + Timepoint,
  data = meta_tax,
  permutations = perm_scheme,
  by = "margin",
)



########################### Only UK and SPAIN MMuphin #########################################################
####################################### Taxonomy #########################################################
ps_merged_filtered <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_merged_filtered.Rds")
otu_tax <- abundances(ps_merged_filtered)
meta_tax <- meta(ps_merged_filtered)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_merged_long<- subset_samples(ps_merged_filtered, Timepoint %in% c("T0"))
ps_merged_long<- subset_samples(ps_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(ps_merged_long)
meta_tax1 <- meta(ps_merged_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_merged_long<- subset_samples(ps_merged_filtered, Timepoint %in% c("T2"))
ps_merged_long<- subset_samples(ps_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(ps_merged_long)
meta_tax2 <- meta(ps_merged_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")



#Remove batch effect taxonomy data and then apply clr transformation
#MMUPHin
#Indicating group data
set.seed(321)
otu_fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = otu_tax,
  covariates = c("Group_cutoff_4", "Age", "Timepoint"),
  batch = "Study",  # Name of the column in meta_tax
  data = meta_tax,
  control = list(verbose = TRUE, diagnostic_plot = "adjust_batch_diagnostic.pdf")
)

ps_corrected<-ps_merged_filtered
otu_table(ps_corrected)<-otu_table(as.matrix(otu_fit_adjust_batch$feature_abd_adj), 
                                   taxa_are_rows = TRUE)
saveRDS(ps_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_corrected.Rds")


# Bray-Curtis distance
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_corrected_T0<- subset_samples(ps_corrected, Timepoint %in% c("T0"))
ps_corrected_T0<- subset_samples(ps_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_c0 <- abundances(ps_corrected_T0)
meta_tax_c0 <- meta(ps_corrected_T0)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax_c0) ~ Group_cutoff_4 + Study, 
        data = meta_tax_c0, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_corrected_T2<- subset_samples(ps_corrected, Timepoint %in% c("T2"))
ps_corrected_T2<- subset_samples(ps_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_c2 <- abundances(ps_corrected_T2)
meta_tax_c2 <- meta(ps_corrected_T2)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax_c2) ~ Group_cutoff_4 + Study, 
        data = meta_tax_c2, permutations = 999, 
        method = "bray", by = "margin")

#Transform to clr
set.seed(156)
ps_clr_corrected<- microbiome::transform(ps_corrected, "clr")
saveRDS(ps_clr_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_clr_corrected.Rds")

# Euclidean distance because the data was clr-transformed
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_clr_corrected_T0<- subset_samples(ps_clr_corrected, Timepoint %in% c("T0"))
ps_clr_corrected_T0<- subset_samples(ps_clr_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_clr0 <- abundances(ps_clr_corrected_T0)
meta_tax_clr0 <- meta(ps_clr_corrected_T0)

# PERMANOVA euclidean distance T0
adonis2(formula = t(otu_tax_clr0) ~ Group_cutoff_4 + Study, 
        data = meta_tax_clr0, permutations = 999, 
        method = "euclidean", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_clr_corrected_T2<- subset_samples(ps_clr_corrected, Timepoint %in% c("T2"))
ps_clr_corrected_T2<- subset_samples(ps_clr_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_clr2 <- abundances(ps_clr_corrected_T2)
meta_tax_clr2 <- meta(ps_clr_corrected_T2)

# PERMANOVA euclidean distance T2
adonis2(formula = t(otu_tax_clr2) ~ Group_cutoff_4 + Study, 
        data = meta_tax_clr2, permutations = 999, 
        method = "euclidean", by = "margin")

dist = phyloseq::distance(ps_clr_corrected, method="euclidean")
ordination = ordinate(ps_clr_corrected, method="PCoA", distance=dist)
plot_ordination(ps_clr_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)



####################################### GMMs #########################################################
GMMs_merged<- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_merged.Rds")
GMMs_comp<- microbiome::transform(GMMs_merged, "compositional")
otu_gmm <- abundances(GMMs_comp)
meta_gmm <- meta(GMMs_comp)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_comp_long<- subset_samples(GMMs_comp, Timepoint %in% c("T0"))
GMMs_comp_long<- subset_samples(GMMs_comp_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(GMMs_comp_long)
meta_tax1 <- meta(GMMs_comp_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_comp_long<- subset_samples(GMMs_comp, Timepoint %in% c("T2"))
GMMs_comp_long<- subset_samples(GMMs_comp_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(GMMs_comp_long)
meta_tax2 <- meta(GMMs_comp_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")


# Remove batch effect indicating group data
set.seed(321)
otu_fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = otu_gmm,
  covariates = c("Group_cutoff_4", "Age", "Timepoint"),
  batch = "Study",  # Name of the column in meta_gmm
  data = meta_gmm,
  control = list(verbose = TRUE, diagnostic_plot = "adjust_batch_diagnostic.pdf")
)

GMMs_corrected<-GMMs_comp
otu_table(GMMs_corrected)<-otu_table(as.matrix(otu_fit_adjust_batch$feature_abd_adj), 
                                     taxa_are_rows = TRUE)
saveRDS(GMMs_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")


#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_corrected_long<- subset_samples(GMMs_corrected, Timepoint %in% c("T0"))
GMMs_corrected_long<- subset_samples(GMMs_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(GMMs_corrected_long)
meta_tax1 <- meta(GMMs_corrected_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_corrected_long<- subset_samples(GMMs_corrected, Timepoint %in% c("T2"))
GMMs_corrected_long<- subset_samples(GMMs_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(GMMs_corrected_long)
meta_tax2 <- meta(GMMs_corrected_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")

dist = phyloseq::distance(GMMs_corrected, method="bray")
ordination = ordinate(GMMs_corrected, method="PCoA", distance=dist)
plot_ordination(GMMs_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)



#PERMANOVA adjusting permutation scheme to see if we have differences between Timepoints
################## Confirm nesting is valid - no- so we use timepoints ############################################
table(meta_tax$Study, meta_tax$Volunteer)  # Should have only one Study per Volunteer

# Number of samples
n_samples <- nrow(meta_tax)
# Number of unique Volunteers
n_volunteers <- length(unique(meta_tax$Volunteer))

# Create the permutation scheme for later use in all comparisons
perm_scheme <- how(
  blocks = meta_tax$Study,              # Study is the outer block
  within = Within(type = "free"),       # Permute volunteers freely within each study
  plots = Plots(strata = meta_tax$Volunteer),
  nperm = 999
)
# Bray-Curtis distance
dist_mat <- vegdist(t(otu_tax), method = "bray")

# Run PERMANOVA with custom permutations
adonis2(
  formula = dist_mat ~ Group_cutoff_4 + Study + Timepoint,
  data = meta_tax,
  permutations = perm_scheme,
  by = "margin",
)
