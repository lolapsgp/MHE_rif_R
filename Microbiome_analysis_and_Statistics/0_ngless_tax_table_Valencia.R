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
motus<-read.table(file = '/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/output_motus.tsv',
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
unassigned_rows <-unassigned_rows[, -c(1, 21, 22)]

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
colnames(motus_df_ok) <- gsub("filtered$", "", colnames(motus_df_ok))
motus_df_ok$mOTUs_ID<-NULL
motus_df_ok$rownames<-NULL
motus_df_ok$in_taxonomy<-NULL

#SAMPLE Data
metadata_rif <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/metadata_rif.Rds")
rownames(metadata_rif)<-metadata_rif$SampleID
rownames(metadata_rif) <- gsub("/", "_", rownames(metadata_rif))
metadata_rif$SampleID<-rownames(metadata_rif)
metadata_rif <- metadata_rif %>%
  mutate(Timepoint = case_when(
    grepl("PC239|PC251|PC286|PC287", SampleID) ~ "T9plus",
    grepl("_3$", SampleID) ~ "T6",
    grepl("_2$", SampleID) ~ "T3",
    TRUE ~ "T0"
  ))
metadata_rif <- metadata_rif %>%
  mutate(Longitudinal = case_when(
    grepl("PC239|PC251|PC286|PC287", SampleID) ~ "No",
    TRUE ~ "Yes"
  ))

TAX = tax_table(as.matrix(motus3.0_taxonomy_ok))
OTU = otu_table(motus_df_ok, taxa_are_rows = TRUE)
SAMPLE = sample_data(metadata_rif)

ps_raw_ngless<- phyloseq(OTU, TAX, SAMPLE)

#save complete 
saveRDS(motus3.0_taxonomy_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/ngless/tax_table_ngless.Rds")
saveRDS(motus_df_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/ngless/otu_table_ngless.Rds")  
saveRDS(ps_raw_ngless, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/ngless/ps_raw_ngless.Rds") 

