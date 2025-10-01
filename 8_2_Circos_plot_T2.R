library(metadeconfoundR)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(vegan)
library(dplyr)
library(tidyr)
library(stringr)
library(phyloseq)
library(circlize)
library(ComplexHeatmap)

library(here)
library(rtk)
library(MetBrewer)
library(umap)
library(DirichletMultinomial)
library(caret)
library(Hmisc)
library(plyr)

#####                                                                #####
##########                                                      ##########
###############              T2 metadeconf runs            ###############
##########                                                      ##########
#####                                                                #####

ps_clr_corrected <- readRDS("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/INCLIVA_2024/MHE_rif/outputs/merged/ps_clr_corrected.Rds")
ps_clr_corrected<- subset_samples(ps_clr_corrected, Group_cutoff_4 %in% c("R", "NR"))
ps_clr_corrected<- subset_samples(ps_clr_corrected, Timepoint %in% c("T2"))
ps_clr_corrected<- subset_samples(ps_clr_corrected, Longitudinal %in% c("Yes"))

ps_genus <- tax_glom(ps_clr_corrected, taxrank = 'Genus', NArm = FALSE)
Genera<-as.data.frame(t(otu_table(ps_genus)))
Genera<-Genera[, order(colnames(Genera))]
Genera<-Genera[order(rownames(Genera)),]

GMMs_corrected <- readRDS("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Group_cutoff_4 %in% c("R", "NR"))
GMMs_corrected<- subset_samples(GMMs_corrected, Timepoint %in% c("T2"))
GMMs_corrected<- subset_samples(GMMs_corrected, Longitudinal %in% c("Yes"))

GMMs<-as.data.frame(t(otu_table(GMMs_corrected)))
GMMs<-GMMs[, order(colnames(GMMs))]
GMMs<-GMMs[order(rownames(GMMs)),]

ARGs_rgi_corrected <- readRDS("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/INCLIVA_2024/MHE_rif/outputs/merged/ARGs_rgi_corrected.Rds")
ARGs_rgi_corrected<- subset_samples(ARGs_rgi_corrected, Group_cutoff_4 %in% c("R", "NR"))
ARGs_rgi_corrected<- subset_samples(ARGs_rgi_corrected, Timepoint %in% c("T2"))
ARGs_rgi_corrected<- subset_samples(ARGs_rgi_corrected, Longitudinal %in% c("Yes"))

ARGs_rgi<-as.data.frame(t(otu_table(ARGs_rgi_corrected)))
ARGs_rgi<-ARGs_rgi[, order(colnames(ARGs_rgi))]
ARGs_rgi<-ARGs_rgi[order(rownames(ARGs_rgi)),]

metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata<-metadata[,c("DmGenderSexID", "Group_cutoff_4", "Age", "PHES", "Ammonia",
                      "IL6", "Haemoglobin", "Leukocytes", "Absolute_Neutrophils",
                      "Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                      "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                      "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                      "Sodium", "Potassium", "Study")]
metadata$Albumin[metadata$Study == "Spain"] <- metadata$Albumin[metadata$Study == "Spain"] * 10
metadata$Neutrophils <- (metadata$Absolute_Neutrophils/metadata$Leukocytes)*100
metadata$Lymphocytes <- (metadata$Absolute_Lymphocytes/metadata$Leukocytes)*100
metadata$Monocytes <- (metadata$Absolute_Monocytes/metadata$Leukocytes)*100
metadata$Eosinophils <- (metadata$Absolute_Eosinophils/metadata$Leukocytes)*100

Imm<- subset(metadata,
             select = c("IL6", "Leukocytes", "Absolute_Neutrophils",
                        "Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                        "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils"))

Liver<- subset(metadata, select = c("Urea", "Creatinine", "Total_bilirubin",
                                    "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase"))

Hem<- subset(metadata, select = c("Haemoglobin", "INR", "Fibrinogen", "Ammonia",
                                  "Sodium", "Potassium"))
featureList<-list()
featureList$Liver<-Liver
featureList$Immunology<-Imm
featureList$Hematology<-Hem
featureList$Genera<-Genera
featureList$GMMs<-GMMs
featureList$ARGs<-ARGs_rgi

metadata_all <- metadata
metadata <- metadata[, c("Study", "DmGenderSexID","Group_cutoff_4", "Age", "PHES")]

metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata <- metadata %>% mutate(Study = ifelse(Study=="Spain",1,0))
metadata <- metadata %>% mutate(Group_cutoff_4 = ifelse(Group_cutoff_4=="R",1,0))
metadata <- metadata %>% mutate(Timepoint = ifelse(Timepoint=="T2",1,0))

meta_V1 <- metadata

#deconfListSafe <- deconfList
deconfList <- list()
for (i in seq_along(featureList)) {
  features <- featureList[[i]]
  features <- features[order(rownames(as.vector(features))), ]
  rmeovers <- c()
  for (j in colnames(features)) {
    if (is.numeric(features[, j])) {
      rmeovers <- c(rmeovers, j)
    }}
  print(rmeovers)
  #saveNames <- rownames(features)
  features <- as.data.frame(features[, rmeovers])
  #rownames(features) <- saveNames
  print(dim(features))
  if (nrow(features)<10 || ncol(features) < 2) {
    next
  }
  metamat <- meta_V1[rownames(meta_V1) %in% rownames(features), ]
  metamat <- metamat[order(rownames(metamat)), ]
  features <- features[rownames(features) %in% rownames(meta_V1), ]
  features <- features[order(rownames(features)), ]
  
  print(dim(metamat))
  print(dim(features))
  print(paste0("starting metadeconf for ", i))
  deconfList[[i]] <- MetaDeconfound(featureMat = features, 
                                    metaMat = metamat, 
                                    nnodes = 3,
                                    returnLong = T, 
                                    randomVar = c("Study"))
  names(deconfList)[i] <- names(featureList)[i]
}


# for (i in seq_along(deconfList)) {
#   print(i)
#   print(BuildHeatmap(deconfList[[i]]))
# }

# B1_deconf <- deconfList[[7]]
# B1_deconf <- MetaDeconfound(featureMat = features, 
#                metaMat = metamat, 
#                nnodes = 3,
#                logfile = here::here("intermediate/PoP_integrate_metadevonfs.log"), robustCutoff = 3)
#BuildHeatmap(deconfList[[7]])

pheno_dec_long <- NA
for (i in names(deconfList)) {
  interMed <- deconfList[[i]]
  interMed$feature <- paste0(i, "_", interMed$feature)
  interMed$metaVariable <- paste0("pheno_", interMed$metaVariable)
  if (is.null(nrow(pheno_dec_long))) {
    pheno_dec_long <- interMed
  } else {
    pheno_dec_long <- rbind(pheno_dec_long, interMed)
  }
  
  print(i)
  print(summary(as.factor(interMed$status)))
}

dim(pheno_dec_long)
pheno_dec_long <- pheno_dec_long[pheno_dec_long$status != "NS", ]
dim(pheno_dec_long)
pheno_dec_long <- pheno_dec_long[!is.infinite(pheno_dec_long$Ds), ]
dim(pheno_dec_long)
pheno_dec_long <- pheno_dec_long[!is.na(pheno_dec_long$feature), ]
dim(pheno_dec_long)

pheno_dec_long <- pheno_dec_long[, c(2,1,3:6)]
colnames(pheno_dec_long)[1] <- "metaVariable"
colnames(pheno_dec_long)[2] <- "feature"

# what ist max nrow
for (i in seq_along(featureList)) {
  print(dim(featureList[[i]]))
} 
# 17  8
# 17 10
# 17  6
# 17 88
# 17 143
# 11 122

set.seed(51)
randomBin13 <- sample(x = c(0,1), replace = T, prob = c(0.5, 0.5), size = 13)
randomBin13 <- rep(c(0,1), 7)[1:13]
randomBin12 <- randomBin13[1:12]

randomBin13 <- sample(x = c(0,1), replace = T, prob = c(0.5, 0.5), size = 363)

#outerList_safe <- outerList

# creat nested list of metadeconfoundR restults for all the possible combination of two datasets
outerList <- list()
for (i in seq_along(featureList)) {
  innerList <- list()
  print(paste("outerList", names(featureList)[i]))
  
  premetamat <- featureList[[i]]
  
  # Add random binary variable with minimum 5 samples per group
  N <- nrow(premetamat)
  min_per_group <- 5
  if (N >= 2 * min_per_group) {
    randomBin <- c(rep(0, min_per_group), rep(1, min_per_group),
                   sample(c(0,1), N - 2 * min_per_group, replace=TRUE))
    randomBin <- sample(randomBin)  # Shuffle to randomize order
    premetamat$randomBin <- randomBin
  } else {
    stop(paste("Not enough samples in dataset", names(featureList)[i], "to generate binary variable with minimum 5 per group"))
  }
  
  premetamat <- premetamat[, c(ncol(premetamat), 1:(ncol(premetamat)-1))]
  rmeovers <- c()
  
  for (j in colnames(premetamat)) {
    if (is.numeric(premetamat[, j])) {
      rmeovers <- c(rmeovers, j)
    }}
  #print(rmeovers)
  #saveNames <- rownames(features)
  premetamat <- as.data.frame(premetamat[, rmeovers])
  print(dim(premetamat))
  for (j in seq_along(featureList)) {
    inner_name <- names(featureList)[j]
    if (i <= j) {
      innerList[[j]] <- NA
      names(innerList)[j] <- inner_name
      print(paste("skipped", i, j))
      next
    }
    print(paste0("now doing ", i, " vs. ", j))
    features <- featureList[[j]]
    rmeovers <- c()
    for (k in colnames(features)) {
      if (is.numeric(features[, k])) {
        rmeovers <- c(rmeovers, k)
      }}
    #print(rmeovers)
    #saveNames <- rownames(features)
    features <- as.data.frame(features[, rmeovers])
    #rownames(features) <- saveNames
    #print(dim(features))
    if (nrow(features)<10 || ncol(features) < 2) {
      innerList[[j]] <- NA
      names(innerList)[j] <- inner_name
      print(paste("skipped", i, j))
      next
    }
    metamat <- premetamat[rownames(premetamat) %in% rownames(features),]
    features <- features[rownames(features) %in% rownames(metamat),]
    features <- features[order(rownames(features)),]
    metamat <- metamat[order(rownames(metamat)),]
    print(dim(features))
    print(dim(metamat))
    innerList[[j]] <- MetaDeconfound(featureMat = features, 
                                     metaMat = metamat, 
                                     nnodes = 3,
                                     deconfF = c("randomBin"),
                                     #logfile = here::here("intermediate/PoP_integrate_metadevonfs.log"), 
                                     returnLong = T)
    print(names(innerList))
    names(innerList)[j] <- inner_name
  }
  outerList[[i]] <- innerList
  names(outerList)[i] <- names(featureList)[i]
}
saveRDS(outerList, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/input/outerList_circos_plot_T2.rds")
outerList <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/input/outerList_circos_plot_T2.rds")
#VS_longSafe <- VS_long

# merge all individual metadeconfoundR outputs into one long-format
VS_long <- NA
for (i in names(outerList)) {
  innerList <- outerList[[i]]
  for (j in names(innerList)) {
    if (is.na(innerList[[j]])[1]) {
      next
    }
    print(paste(i,j))
    print(sum(!(innerList[[j]]$status %in% c("NA", "NS"))))
    print(sum((innerList[[j]]$Qs < 0.1), na.rm = T))
    innerList[[j]]$feature <- paste0(j, "_", innerList[[j]]$feature)
    innerList[[j]]$metaVariable <- paste0(i, "_", innerList[[j]]$metaVariable)
    if (is.na(VS_long[[1]])[1]) {
      VS_long <- innerList[[j]]
    } else {
      VS_long <- rbind(VS_long, innerList[[j]])
    }
  }
}

raw_VS_long <- VS_long
# remove all non-significant associations
# remove associations with the "random binary" variable, that was added to datasets 
# in order to fulfill metadeconfoundR requirment of 
# having a binary first metadata column
dim(VS_long)
VS_long <- VS_long[VS_long$status != "NS", ]
dim(VS_long)
rmRandBin <- grep(pattern = "randomBin", x = VS_long$metaVariable)
if (length(rmRandBin) > 0) {
  VS_long <- VS_long[-c(rmRandBin) , ]
}
dim(VS_long)


# add metadeconfoundR runs against the phenotipc data
VS_long <- rbind(VS_long, pheno_dec_long)
VS_long <- VS_long[!grepl("^C: ", VS_long$status), ]
data_all_1 <- VS_long %>%
  separate(feature, into = c("Group_1", "feature"), sep = "_", extra = "merge")
data_all <- data_all_1 %>%
  separate(metaVariable, into = c("Group_2", "metaVariable"), sep = "_", extra = "merge")
rm(data_all_1)
# export raw associations without removal of NS
rmRandBin <- grep(pattern = "randomBin", x = VS_long$metaVariable)
if (length(rmRandBin) > 0) {
  raw_VS_long <- raw_VS_long[-c(rmRandBin) , ]
}
raw_VS_long <- rbind(raw_VS_long, pheno_dec_long)
#write.table(raw_VS_long, here::here("output/allTestedAssociations.tsv"), sep = "\t", row.names = F)

#write.table(VS_long, here::here("output/allSignificantAssociations.tsv"), sep = "\t", row.names = F)


#####                                                                #####
##########                                                      ##########
###############                 break point                ###############
##########                                                      ##########
#####                                                                #####

#####                                                                #####
##########                                                      ##########
###############           start here for ninimal run       ###############
##########                                                      ##########
#####                                                                #####

# use these two commads to load all neccesary data to plot a new plot
#VS_long <- read.table("allSignificantAssociations.tsv", sep = "\t", header = T)
#load("allNiceNames.r")


Genera <- as.data.frame((tax_table(ps_genus))[,c("Genus")])
Genera$names<-Genera$Genus
Genera$Genus<-NULL
GMMs<- as.data.frame(tax_table(GMMs_corrected))
GMMs$names<-GMMs$V2
GMMs$V2<-NULL
ARGs_rgi<- as.data.frame(tax_table(ARGs_rgi_corrected))
ARGs_rgi$names<-ARGs_rgi$Best_Hit_ARO
ARGs_rgi<-ARGs_rgi[,c("names", "Best_Hit_ARO")]
ARGs_rgi$Best_Hit_ARO<-NULL
names<-colnames(metadata)
pheno<- as.data.frame(names)
Imm_names <- data.frame(
  names = c("IL6", "Leukocytes", "Absolute_Neutrophils",
            "Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
            "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils"))
rownames(Imm_names)<-colnames(Imm)

Liver_names <- data.frame(
  names = c("Urea", "Creatinine", "Total_bilirubin",
            "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase"))
rownames(Liver_names)<-colnames(Liver)

Hem_names <- data.frame(
  names = c("Haemoglobin", "INR", "Fibrinogen", "Ammonia",
            "Sodium", "Potassium"))
rownames(Hem_names)<-colnames(Hem)

allNiceNames <- list("Liver" = Liver_names,
                     "Immunology" = Imm_names,
                     "Hematology" = Hem_names,
                     "ARGs" = ARGs_rgi,
                     "Genera" = Genera, 
                     "GMMs" = GMMs, 
                     "pheno" = pheno)


# replace "machine readable" names with the pretty names supplied by the sub groups
nameMap <- dplyr::bind_rows(allNiceNames)
nameMap$old_names<-rownames(nameMap)
name_lookup <- setNames(nameMap$names, nameMap$old_names)

# Replace the features in data_all$feature using the lookup vector
data_all$feature <- ifelse(data_all$feature %in% names(name_lookup) & !is.na(name_lookup[data_all$feature]),
                           name_lookup[data_all$feature],
                           data_all$feature)
data_all$metaVariable <- ifelse(data_all$metaVariable %in% names(name_lookup) & !is.na(name_lookup[data_all$metaVariable]),
                                name_lookup[data_all$metaVariable],
                                data_all$metaVariable)

# assign a group and color to each feature
# this is needed for the circos plot
color_palette <- c("Liver" = "#999999", 
                   "Immunology" = "#E69F00", 
                   "Hematology" = "#56B4E9", 
                   "ARGs" = "#009E73", 
                   "Genera" = "#CC79A7", 
                   "GMMs" = "#0072B2",
                   "pheno" = "#BEBADA")

# Too many associations so that we cann see something, we filter by p-value
data_all <- data_all[data_all$Ps <= 0.01, ]

data_all<-data_all%>%
  mutate(metaVariable = case_when(
    metaVariable %in% c("ASV206") ~ "Lachnospiraceae_genusNA)",
    metaVariable %in% c("ASV312") ~ "Ruminococcaceae_genusNA",
    metaVariable %in% c("ASV198") ~ "Rhodospirillales_famNA_genusNA",
    metaVariable %in% c("ASV99") ~ "Clostridia UCG-014",
    metaVariable %in% c("GENERO") ~ "Gender",
    metaVariable %in% c("Type_Control") ~ "Subgroup: Control",
    metaVariable %in% c("Type_MCI") ~ "Subgroup: MCI",
    metaVariable %in% c("Type_MHE") ~ "Subgroup: MHE",
    metaVariable %in% c("Type_NMCI") ~ "Subgroup: NMCI",
    metaVariable %in% c("Type_NMHE") ~ "Subgroup: NMHE",
    TRUE ~ metaVariable))
data_all<-data_all%>%
  mutate(feature = case_when(
    feature %in% c("ASV390") ~ "Izemoplasmatales_famNA_genusNA",
    feature %in% c("ASV1297") ~ "Oscillospirales UCG-010",
    feature %in% c("ASV20") ~ "[Eubacterium] coprostanoligenes group",
    feature %in% c("ASV99") ~ "Clostridia UCG-014",
    feature %in% c("ASV206") ~ "Lachnospiraceae_genusNA",
    feature %in% c("ASV198") ~ "Rhodospirillales_famNA_genusNA",
    feature %in% c("ASV328") ~ "Barnesiellaceae_genusNA",
    TRUE ~ feature))

data_all <- data_all %>%
  mutate(colours = color_palette[Group_1])
color_palette1 <- data_all[, c("feature", "colours")]
group <- setNames(color_palette1$colours, color_palette1$feature)

data_all <- data_all %>%
  mutate(colours2 = color_palette[Group_2])
color_palette2 <- data_all[, c("metaVariable", "colours2")]
group2 <- setNames(color_palette2$colours2, color_palette2$metaVariable)
# Combine the names from both vectors
all_names <- c(names(group), names(group2))
unique_names <- unique(all_names)
# Combine the values from both vectors
all_values <- c(group, group2)
merged_group <- all_values[match(unique_names, all_names)]
names(merged_group) <- unique_names

# split data into positive and negative associations, minimum of 35% effect size
VS_long_pos <- data_all[data_all$Ds > 0.4, ]
VS_long_neg <- data_all[data_all$Ds < -0.4, ]
# Remove rows where Ds is NA
VS_long_pos <- VS_long_pos[!is.na(VS_long_pos$Ds), ]
VS_long_neg <- VS_long_neg[!is.na(VS_long_neg$Ds), ]
# create names for figure legend
new_names <- c("Liver function", 
               "Immunology", 
               "Hematology", 
               "Microbiome (ARGs)", 
               "Microbiome (genera)", 
               "Microbiome (GMMs)",
               "Phenotype"
)

color_palette <- c("Liver function" = "#999999", 
                   "Immunology" = "#E69F00", 
                   "Hematology" = "#56B4E9", 
                   "Microbiome (ARGs)" = "#009E73", 
                   "Microbiome (genera)" = "#CC79A7", 
                   "Microbiome (GMMs)" = "#0072B2",
                   "Phenotype" = "#BEBADA")

lgd_points <- Legend(at = new_names, 
                     type = "points", 
                     legend_gp = gpar(col = color_palette[new_names]), 
                     title_position = "topleft", 
                     title = "Feature Spaces")

# create 2 circos plots and save them as svg

svg(filename = "circos_plot_40_01_T2.svg", width = 10, height = 8)
par(mfrow=c(1,2))
circos.clear()
circos.par(start.degree = 60)
chordDiagram(VS_long_pos[, c(2,4,7)], annotationTrack = "grid",grid.col = merged_group, group = merged_group,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(VS_long_pos[, c(2,4,7)]))))*0.5))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.3)
}, bg.border = NA)

title(main="A", col="black", font=2)

circos.clear()
circos.par(start.degree = 60)
chordDiagram(VS_long_neg[, c(2,4,7)], annotationTrack = "grid", grid.col = merged_group, group = merged_group,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(VS_long_neg[, c(2,4,7)]))))*0.5))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.3)
}, bg.border = NA)
title(outer=T,main="B",col="black",font=2)
draw(lgd_points, just = c("left", "bottom"))
dev.off() 
#Remember to move the legend using Inkscape