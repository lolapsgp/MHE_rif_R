library(vegan)
library(permute)
library(phyloseq)
library(microbiome)

####################################### ARSyN UK and SPAIN #########################################################
####################################### Taxonomy #########################################################
met_merged_filtered <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_merged_filtered_met.Rds")
met_merged_filtered<- subset_samples(met_merged_filtered, Group_cutoff_4 %in% c("R", "NR"))
met_merged_filtered<- subset_samples(met_merged_filtered, Timepoint %in% c("T0", "T2"))

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_met_merged_long<- subset_samples(met_merged_filtered, Timepoint %in% c("T0"))
ps_met_merged_long<- subset_samples(ps_met_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax1 <- abundances(ps_met_merged_long)
meta_tax1 <- meta(ps_met_merged_long)

# PERMANOVA Bray-Curtis distance T0
adonis2(formula = t(otu_tax1) ~ Group_cutoff_4 + Study, 
        data = meta_tax1, permutations = 999, 
        method = "bray", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_met_merged_long<- subset_samples(met_merged_filtered, Timepoint %in% c("T2"))
ps_met_merged_long<- subset_samples(ps_met_merged_long, Group_cutoff_4 %in% c("R", "NR"))
otu_tax2 <- abundances(ps_met_merged_long)
meta_tax2 <- meta(ps_met_merged_long)

# PERMANOVA Bray-Curtis distance T2
adonis2(formula = t(otu_tax2) ~ Group_cutoff_4 + Study, 
        data = meta_tax2, permutations = 999, 
        method = "bray", by = "margin")



#Remove batch effect taxonomy data and then apply clr transformation
#ARSyN
#Indicating group data
otu <- data.frame(otu_table(met_merged_filtered))
metadata<-data.frame(sample_data(met_merged_filtered))
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
TAX = tax_table(met_merged_filtered)
SAMPLE <- sample_data(metadata)

ps_met_corrected <- phyloseq(OTU, TAX, SAMPLE)
saveRDS(ps_met_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_met_corrected.Rds")

#+225 because clr tranfromation need no negative data
OTU = otu_table(t(otu_corrected+225), taxa_are_rows = TRUE)
TAX = tax_table(met_merged_filtered)
SAMPLE <- sample_data(metadata)
ps_met_corrected <- phyloseq(OTU, TAX, SAMPLE)
ps_met_clr_corrected<- microbiome::transform(ps_met_corrected, "clr")
saveRDS(ps_met_clr_corrected, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_met_clr_corrected.Rds")

# Euclidean distance
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_met_corrected_T0<- subset_samples(ps_met_corrected, Timepoint %in% c("T0"))
ps_met_corrected_T0<- subset_samples(ps_met_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_c0 <- abundances(ps_met_corrected_T0)
meta_tax_c0 <- meta(ps_met_corrected_T0)

# PERMANOVA Euclidean distance T0
adonis2(formula = t(otu_tax_c0) ~ Group_cutoff_4 + Study, 
        data = meta_tax_c0, permutations = 999, 
        method = "euclidean", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_met_corrected_T2<- subset_samples(ps_met_corrected, Timepoint %in% c("T2"))
ps_met_corrected_T2<- subset_samples(ps_met_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_c2 <- abundances(ps_met_corrected_T2)
meta_tax_c2 <- meta(ps_met_corrected_T2)

# PERMANOVA Euclidean distance T0
adonis2(formula = t(otu_tax_c2) ~ Group_cutoff_4 + Study, 
        data = meta_tax_c2, permutations = 999, 
        method = "euclidean", by = "margin")

# Euclidean distance because the data was clr-transformed
#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_met_clr_corrected_T0<- subset_samples(ps_met_clr_corrected, Timepoint %in% c("T0"))
ps_met_clr_corrected_T0<- subset_samples(ps_met_clr_corrected_T0, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_clr0 <- abundances(ps_met_clr_corrected_T0)
meta_tax_clr0 <- meta(ps_met_clr_corrected_T0)

# PERMANOVA euclidean distance T0
adonis2(formula = t(otu_tax_clr0) ~ Group_cutoff_4 + Study, 
        data = meta_tax_clr0, permutations = 999, 
        method = "euclidean", by = "margin")

#PERMANOVA at T0 and at T2 to check differences by Study and Group
ps_met_clr_corrected_T2<- subset_samples(ps_met_clr_corrected, Timepoint %in% c("T2"))
ps_met_clr_corrected_T2<- subset_samples(ps_met_clr_corrected_T2, Group_cutoff_4 %in% c("R", "NR"))
otu_tax_clr2 <- abundances(ps_met_clr_corrected_T2)
meta_tax_clr2 <- meta(ps_met_clr_corrected_T2)

# PERMANOVA euclidean distance T2
adonis2(formula = t(otu_tax_clr2) ~ Group_cutoff_4 + Study, 
        data = meta_tax_clr2, permutations = 999, 
        method = "euclidean", by = "margin")

dist = phyloseq::distance(ps_met_clr_corrected, method="euclidean")
ordination = ordinate(ps_met_clr_corrected, method="PCoA", distance=dist)
plot_ordination(ps_met_clr_corrected, ordination, color="Study") + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  stat_ellipse(aes(group = Study), linetype = 3, level = 0.68) +
  geom_line(aes(group = Volunteer, color = Study)) +
  geom_text(aes(label = Timepoint), vjust = -1, size = 3)
