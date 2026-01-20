library(dplyr)
library(stringr)

#Load taxonomy table
load(url("https://github.com/AlessioMilanese/motus_taxonomy/blob/master/data/motus_taxonomy_3.0.1.Rdata?raw=true"))
# Remove digits and space at the beginning of each rowname
rownames(motus3.0_taxonomy) <- sub("^\\d+\\s", "", rownames(motus3.0_taxonomy))
motus3.0_taxonomy$Kingdom <- sub("^\\d+\\s", "", motus3.0_taxonomy$Kingdom)
motus3.0_taxonomy$Phylum <- sub("^\\d+\\s", "", motus3.0_taxonomy$Phylum)
motus3.0_taxonomy$Class <- sub("^\\d+\\s", "", motus3.0_taxonomy$Class)
motus3.0_taxonomy$Order <- sub("^\\d+\\s", "", motus3.0_taxonomy$Order)
motus3.0_taxonomy$Family <- sub("^\\d+\\s", "", motus3.0_taxonomy$Family)
motus3.0_taxonomy$Genus <- sub("^\\d+\\s", "", motus3.0_taxonomy$Genus)
motus3.0_taxonomy$Species <- sub("^\\d+\\s", "", motus3.0_taxonomy$Species)
motus3.0_taxonomy$profiled <- sub("^\\d+\\s", "", motus3.0_taxonomy$profiled)

motus_taxonomy<-motus3.0_taxonomy %>%
  rownames_to_column("rownames")%>%
  mutate(mOTUs_ID = gsub("_v3_", "_v31_", mOTUs_ID)) 

#Import NGLess outputs
motus<-read.table(file = '/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Muto_2022_08_01/output_motus.tsv',
                  header = TRUE, 
                  sep = "\t", 
                  row.names = 1, 
                  skip = 0, 
                  check.names = FALSE,
                  fill = TRUE,
                  quote = "")

motus <- motus[, order(names(motus))]
motus <- motus %>%
  select(-ends_with(".1"))
motus <- na.omit(motus)

# Remove zeros
motus_ok <- motus[rowSums(motus != 0) > 0, ]

# Extract Taxa IDs and add as a new column
motus_df <- motus_ok %>%
  rownames_to_column(var = "rownames") %>%
  mutate(mOTUs_ID = str_extract(rownames, "(?<=\\[)(ref_mOTU_v31_\\d+|ext_mOTU_v31_\\d+|meta_mOTU_v31_\\d+)(?=\\])"))

# Check if mOTUs_ID in motus_df are contained in motus_taxonomy
motus_df$in_taxonomy <- motus_df$mOTUs_ID %in% motus_taxonomy$mOTUs_ID
summary(motus_df$in_taxonomy)

#Merge rows with in_taxonomy = FALSE with "unassigned"
unassigned_rows <- motus_df[motus_df$in_taxonomy == FALSE, ]
unassigned_rows <- unassigned_rows[,-c(1, 14, 15)]

#Sum the counts of the filtered columns for unassigned rows
unassigned_counts <- colSums(unassigned_rows)  # Exclude rownames and in_taxonomy columns

#Create a new row for unassigned
unassigned_row <- data.frame(
  rownames = "unassigned",
  t(unassigned_counts),  # Transpose to match the original data frame structure
  mOTUs_ID = "unassigned",
  in_taxonomy = TRUE,    # Set in_taxonomy to TRUE for the unassigned row
  stringsAsFactors = FALSE
)

#Keep TRUE and Combine the data frames
motus_df <- motus_df[motus_df$in_taxonomy == TRUE, ]
colnames(unassigned_row) <- gsub("\\.", "-", colnames(unassigned_row))
motus_df <- rbind(motus_df, unassigned_row)
summary(motus_df$in_taxonomy)


#Create phyloseq object
library(phyloseq)
#Tax table adjust
motus3.0_taxonomy_ok<-motus_taxonomy
rownames(motus3.0_taxonomy_ok)<-motus_taxonomy$mOTUs_ID
motus3.0_taxonomy_ok$mOTUs_ID<-NULL
motus3.0_taxonomy_ok$rownames<-NULL


#mOTU table adjust
motus_df_ok<-motus_df
rownames(motus_df_ok)<-motus_df$mOTUs_ID
motus_df_ok$mOTUs_ID<-NULL
motus_df_ok$rownames<-NULL
motus_df_ok$in_taxonomy<-NULL

# SAMPLE Data
# let's create the metadata data frame with SampleID
metadata <- data.frame(SampleID = colnames(motus_df_ok))

# We remove "_HE_Shotgun" from the end of each SampleID and extract the Volunteer and Timepoint information
metadata$Volunteer <- sub("-([0-9]+)_.*", "", metadata$SampleID)
metadata$Timepoint <- gsub("R\\d+-([0-9]+)_.*", "T\\1", metadata$SampleID)

metadata$Group <- ifelse(metadata$Volunteer %in% c("R05", "R10", "R13"), "Responder",
                         ifelse(metadata$Volunteer %in% c("R02", "R06", "R07"), "Non-Responder", NA))
rownames(metadata)<- metadata$SampleID

TAX = tax_table(as.matrix(motus3.0_taxonomy_ok))
OTU = otu_table(motus_df_ok, taxa_are_rows = TRUE)
SAMPLE = sample_data(metadata)

ps_raw_Muto<- phyloseq(OTU, TAX, SAMPLE)

#save complete 
saveRDS(motus3.0_taxonomy_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Muto_analysis/outputs/tax_table_ngless.Rds")
saveRDS(motus_df_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Muto_analysis/outputs/otu_table_ngless.Rds")  
saveRDS(ps_raw_Muto, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Muto_analysis/outputs/ps_raw_Muto.Rds") 
saveRDS(metadata, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Secuences_bibliography/Muto_analysis/outputs/Metadata_Muto.Rds") 