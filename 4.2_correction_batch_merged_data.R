library(vegan)
library(permute)
library(phyloseq)
library(microbiome)

####################################### ARSyN UK and SPAIN #########################################################
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

#+225 because clr tranfromation need no negative data
OTU = otu_table(t(otu_corrected+225), taxa_are_rows = TRUE)
TAX = tax_table(ps_merged_filtered)
SAMPLE <- sample_data(metadata)
ps_corrected <- phyloseq(OTU, TAX, SAMPLE)
ps_clr_corrected<- microbiome::transform(ps_corrected, "clr")
saveRDS(ps_clr_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_clr_corrected.Rds")

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
GMMs_merged<- subset_samples(GMMs_merged, Timepoint %in% c("T0", "T2"))
GMMs_merged<- subset_samples(GMMs_merged, Group_cutoff_4 %in% c("R", "NR"))
otu_gmm <- abundances(GMMs_merged)
meta_gmm <- meta(GMMs_merged)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_merged_long<- subset_samples(GMMs_merged, Timepoint %in% c("T0"))
GMMs_merged_long<- subset_samples(GMMs_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(GMMs_merged_long)
meta_tax1 <- meta(GMMs_merged_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_merged_long<- subset_samples(GMMs_merged, Timepoint %in% c("T2"))
GMMs_merged_long<- subset_samples(GMMs_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(GMMs_merged_long)
meta_tax2 <- meta(GMMs_merged_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")


# Remove batch effect: ARSyN
#Indicating group data
otu <- data.frame(otu_table(GMMs_merged))
metadata<-data.frame(sample_data(GMMs_merged))
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
TAX = tax_table(GMMs_merged)
SAMPLE <- sample_data(metadata)

GMMs_corrected <- phyloseq(OTU, TAX, SAMPLE)

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

####################################### ARGs rgi #########################################################
ARGs_rgi <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ARGs_rgi.Rds")
ARGs_rgi<- subset_samples(ARGs_rgi, Group_cutoff_4 %in% c("R", "NR"))
ARGs_rgi<- subset_samples(ARGs_rgi, Timepoint %in% c("T0", "T2"))

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_long<- subset_samples(ARGs_rgi, Timepoint %in% c("T0"))
ARGs_rgi_long<- subset_samples(ARGs_rgi_long, Group_cutoff_4 %in% c("R", "NR"))
otu_arg_rgi1 <- abundances(ARGs_rgi_long)
meta_arg_rgi1 <- meta(ARGs_rgi_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_arg_rgi1) ~ Group_cutoff_4 + Study, 
        data = meta_arg_rgi1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_long<- subset_samples(ARGs_rgi, Timepoint %in% c("T2"))
ARGs_rgi_long<- subset_samples(ARGs_rgi_long, Group_cutoff_4 %in% c("R", "NR"))
otu_arg_rgi2 <- abundances(ARGs_rgi_long)
meta_arg_rgi2 <- meta(ARGs_rgi_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_arg_rgi2) ~ Group_cutoff_4 + Study, 
        data = meta_arg_rgi2, permutations = 999, 
        method = "bray", by = "margin")



#Remove batch effect arg_rgionomy data and then apply clr transformation
#ARSyN
#Indicating group data
otu <- data.frame(otu_table(ARGs_rgi))
metadata<-data.frame(sample_data(ARGs_rgi))
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
TAX = tax_table(ARGs_rgi)
SAMPLE <- sample_data(metadata)

ARGs_rgi_corrected <- phyloseq(OTU, TAX, SAMPLE)


saveRDS(ARGs_rgi_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ARGs_rgi_corrected.Rds")


# Euclidean distance
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_corrected_T0<- subset_samples(ARGs_rgi_corrected, Timepoint %in% c("T0"))
ARGs_rgi_corrected_T0<- subset_samples(ARGs_rgi_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_arg_rgi_c0 <- abundances(ARGs_rgi_corrected_T0)
meta_arg_rgi_c0 <- meta(ARGs_rgi_corrected_T0)

# PERMANOVA Euclidean distance T0
adonis2(formula = t(otu_arg_rgi_c0) ~ Group_cutoff_4 + Study, 
        data = meta_arg_rgi_c0, permutations = 999, 
        method = "euclidean", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_corrected_T2<- subset_samples(ARGs_rgi_corrected, Timepoint %in% c("T2"))
ARGs_rgi_corrected_T2<- subset_samples(ARGs_rgi_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_arg_rgi_c2 <- abundances(ARGs_rgi_corrected_T2)
meta_arg_rgi_c2 <- meta(ARGs_rgi_corrected_T2)

# PERMANOVA Euclidean distance T0
adonis2(formula = t(otu_arg_rgi_c2) ~ Group_cutoff_4 + Study, 
        data = meta_arg_rgi_c2, permutations = 999, 
        method = "euclidean", by = "margin")

dist = phyloseq::distance(ARGs_rgi_corrected, method="euclidean")
ordination = ordinate(ARGs_rgi_corrected, method="PCoA", distance=dist)
plot_ordination(ARGs_rgi_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)


####################################### ARGs ResFinder #########################################################
#All ARGs in the UK cohort only in NR group, so ARSyN is not applicable







########################### Not used ############################################
######################## MMuphin UK and SPAIN ###################################
######################### Taxonomy ###############################################
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
otu_gmm <- abundances(GMMs_merged)
meta_gmm <- meta(GMMs_merged)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_merged_long<- subset_samples(GMMs_merged, Timepoint %in% c("T0"))
GMMs_merged_long<- subset_samples(GMMs_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(GMMs_merged_long)
meta_tax1 <- meta(GMMs_merged_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
GMMs_merged_long<- subset_samples(GMMs_merged, Timepoint %in% c("T2"))
GMMs_merged_long<- subset_samples(GMMs_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(GMMs_merged_long)
meta_tax2 <- meta(GMMs_merged_long)

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

GMMs_corrected<-GMMs_merged
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

####################################### ARGs rgi #########################################################
ARGs_rgi<- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ARGs_rgi.Rds")
ARGs_rgi<- subset_samples(ARGs_rgi, Group_cutoff_4 %in% c("R", "NR"))
ARGs_rgi<- subset_samples(ARGs_rgi, Timepoint %in% c("T0", "T2"))
otu_ARG_rgi <- as.data.frame(abundances(ARGs_rgi))
#We remove decimals
otu_ARG_rgi <- otu_ARG_rgi %>%
  mutate(across(where(is.numeric), as.integer))
meta_ARG_rgi <- meta(ARGs_rgi)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_long<- subset_samples(ARGs_rgi, Timepoint %in% c("T0"))
ARGs_rgi_long<- subset_samples(ARGs_rgi_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_rgi1 <- abundances(ARGs_rgi_long)
meta_ARG_rgi1 <- meta(ARGs_rgi_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_ARG_rgi1) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_rgi1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_long<- subset_samples(ARGs_rgi, Timepoint %in% c("T2"))
ARGs_rgi_long<- subset_samples(ARGs_rgi_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_rgi2 <- abundances(ARGs_rgi_long)
meta_ARG_rgi2 <- meta(ARGs_rgi_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_ARG_rgi2) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_rgi2, permutations = 999, 
        method = "bray", by = "margin")


# Remove batch effect indicating group data
set.seed(321)
otu_fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = otu_ARG_rgi,
  covariates = c("Group_cutoff_4", "Age", "Timepoint"),
  batch = "Study",  # Name of the column in meta_ARG_rgi
  data = meta_ARG_rgi,
  control = list(verbose = TRUE, diagnostic_plot = "adjust_batch_diagnostic.pdf")
)

ARGs_rgi_corrected<-ARGs_rgi
otu_table(ARGs_rgi_corrected)<-otu_table(as.matrix(otu_fit_adjust_batch$feature_abd_adj), 
                                         taxa_are_rows = TRUE)
saveRDS(ARGs_rgi_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ARGs_rgi_corrected.Rds")


#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_corrected_long<- subset_samples(ARGs_rgi_corrected, Timepoint %in% c("T0"))
ARGs_rgi_corrected_long<- subset_samples(ARGs_rgi_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_rgi1 <- abundances(ARGs_rgi_corrected_long)
meta_ARG_rgi1 <- meta(ARGs_rgi_corrected_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_ARG_rgi1) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_rgi1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_rgi_corrected_long<- subset_samples(ARGs_rgi_corrected, Timepoint %in% c("T2"))
ARGs_rgi_corrected_long<- subset_samples(ARGs_rgi_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_rgi2 <- abundances(ARGs_rgi_corrected_long)
meta_ARG_rgi2 <- meta(ARGs_rgi_corrected_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_ARG_rgi2) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_rgi2, permutations = 999, 
        method = "euclidean", by = "margin")

dist = phyloseq::distance(ARGs_rgi_corrected, method="euclidean")
ordination = ordinate(ARGs_rgi_corrected, method="PCoA", distance=dist)
plot_ordination(ARGs_rgi_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)

####################################### ARGs ResFinder #########################################################
ARGs_ResFinder<- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ARGs_ResFinder.Rds")
ARGs_ResFinder<- subset_samples(ARGs_ResFinder, Group_cutoff_4 %in% c("R", "NR"))
ARGs_ResFinder<- subset_samples(ARGs_ResFinder, Timepoint %in% c("T0", "T2"))
otu_ARG_resf <- as.data.frame(abundances(ARGs_ResFinder))
#We remove decimals
otu_ARG_resf <- otu_ARG_resf %>%
  mutate(across(where(is.numeric), as.integer))
meta_ARG_resf <- meta(ARGs_ResFinder)

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_ResFinder_long<- subset_samples(ARGs_ResFinder, Timepoint %in% c("T0"))
ARGs_ResFinder_long<- subset_samples(ARGs_ResFinder_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_resf1 <- abundances(ARGs_ResFinder_long)
meta_ARG_resf1 <- meta(ARGs_ResFinder_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_ARG_resf1) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_resf1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_ResFinder_long<- subset_samples(ARGs_ResFinder, Timepoint %in% c("T2"))
ARGs_ResFinder_long<- subset_samples(ARGs_ResFinder_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_resf2 <- abundances(ARGs_ResFinder_long)
meta_ARG_resf2 <- meta(ARGs_ResFinder_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_ARG_resf2) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_resf2, permutations = 999, 
        method = "bray", by = "margin")


# Remove batch effect indicating group data
set.seed(321)
otu_fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = otu_ARG_resf,
  covariates = c("Group_cutoff_4", "Age", "Timepoint"),
  batch = "Study",  # Name of the column in meta_ARG_resf
  data = meta_ARG_resf,
  control = list(verbose = TRUE, diagnostic_plot = "adjust_batch_diagnostic.pdf")
)

ARGs_ResFinder_corrected<-ARGs_ResFinder
otu_table(ARGs_ResFinder_corrected)<-otu_table(as.matrix(otu_fit_adjust_batch$feature_abd_adj), 
                                               taxa_are_rows = TRUE)
saveRDS(ARGs_ResFinder_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ARGs_ResFinder_corrected.Rds")


#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_ResFinder_corrected_long<- subset_samples(ARGs_ResFinder_corrected, Timepoint %in% c("T0"))
ARGs_ResFinder_corrected_long<- subset_samples(ARGs_ResFinder_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_resf1 <- abundances(ARGs_ResFinder_corrected_long)
meta_ARG_resf1 <- meta(ARGs_ResFinder_corrected_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_ARG_resf1) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_resf1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ARGs_ResFinder_corrected_long<- subset_samples(ARGs_ResFinder_corrected, Timepoint %in% c("T2"))
ARGs_ResFinder_corrected_long<- subset_samples(ARGs_ResFinder_corrected_long, Group_cutoff_4 %in% c("R", "NR"))
otu_ARG_resf2 <- abundances(ARGs_ResFinder_corrected_long)
meta_ARG_resf2 <- meta(ARGs_ResFinder_corrected_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_ARG_resf2) ~ Group_cutoff_4 + Study, 
        data = meta_ARG_resf2, permutations = 999, 
        method = "euclidean", by = "margin")

dist = phyloseq::distance(ARGs_ResFinder_corrected, method="euclidean")
ordination = ordinate(ARGs_ResFinder_corrected, method="PCoA", distance=dist)
plot_ordination(ARGs_ResFinder_corrected, ordination, color="Study") + 
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





(ARGs_ResFinder, file = "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Patel_analysis/outputs/ARGs_ResFinder.Rds")
