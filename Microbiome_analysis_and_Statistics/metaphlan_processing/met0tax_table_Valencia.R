library(tidyr)
library(dplyr)
library(microbiome)
metaphlan_out<-read.table(file = '/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/merged_abundance_table.txt',
                          header = TRUE, 
                          sep = "\t", 
                          row.names = 1, 
                          skip = 0, 
                          check.names = FALSE,
                          fill = TRUE,
                          quote = "")

metaphlan_out$Taxonomy<-rownames(metaphlan_out)
tax_df <- data.frame(Taxonomy = metaphlan_out$Taxonomy, stringsAsFactors = FALSE)
tax.clean <- tax_df %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "strain"), 
           sep = "\\|", fill = "right", extra = "drop")

# Remove the prefix from the relevant columns
tax.clean$Kingdom <- gsub("^[a-z]__", "", tax.clean$Kingdom)
tax.clean$Phylum <- gsub("^[a-z]__", "", tax.clean$Phylum)
tax.clean$Class   <- gsub("^[a-z]__", "", tax.clean$Class)
tax.clean$Order <- gsub("^[a-z]__", "", tax.clean$Order)
tax.clean$Family <- gsub("^[a-z]__", "", tax.clean$Family)
tax.clean$Genus <- gsub("^[a-z]__", "", tax.clean$Genus)
tax.clean$Species <- gsub("^[a-z]__", "", tax.clean$Species)
tax.clean$strain <- gsub("^[a-z]__", "", tax.clean$strain)

all<-cbind(metaphlan_out, tax.clean)
#Deletion of not strain taxons (higher taxons with abundance data duplicated)
all <- all[!is.na(all$strain), ]
strain_table <- all[,startsWith(colnames(all), "PC") ]
tax.clean <- tax.clean[!is.na(tax.clean$strain), ]


strain_table <- strain_table[rowSums(strain_table != 0) > 0, ]
colnames(strain_table) <- gsub("filtered$", "", colnames(strain_table))
rownames(strain_table) <- tax.clean$strain
rownames(tax.clean) <- tax.clean$strain

#SAMPLE Data
ps_clr_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/ps_clr_corrected.Rds")
ps_clr_corrected<-subset_samples(ps_clr_corrected, Study == "Spain")
metadata<-(sample_data(ps_clr_corrected))

TAX = tax_table(as.matrix(tax.clean))
OTU = otu_table(strain_table, taxa_are_rows = TRUE)
SAMPLE = sample_data(metadata)

#Sample data is missing SAMPLE = sample_data(metadata)

ps_raw_metaphlan<- phyloseq(OTU, TAX, SAMPLE)

#save complete saveRDS(ps_feces, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/outputs/R/ps_raw.Rds")  
saveRDS(tax.clean, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/metaphlan/tax_table_metaphlan.Rds")
saveRDS(metaphlan_out, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/metaphlan/otu_table_metaphlan.Rds") 
saveRDS(ps_raw_metaphlan, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/metaphlan/ps_raw_metaphlan.Rds") 
