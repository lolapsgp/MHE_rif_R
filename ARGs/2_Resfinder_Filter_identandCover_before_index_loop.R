library(tidyverse)
library(Biostrings)

# Set root directory for all samples
FASTA_ROOT <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Resfinder_results"

# Get all sample directories under FASTA_ROOT
sample_dirs <- list.dirs(FASTA_ROOT, full.names = TRUE, recursive = FALSE)

for (sample_dir in sample_dirs) {
  # Define result file paths for this sample
  tabfile <- file.path(sample_dir, "ResFinder_results_tab.txt")
  fastafile <- file.path(sample_dir, "ResFinder_Resistance_gene_seq_renamed.fsa")
  outfile <- file.path(sample_dir, "ResFinder_Resistance_gene_seq_filtered.fsa")
  
  # Skip if files are missing
  if (!file.exists(tabfile) | !file.exists(fastafile)) next
  
  # Read and filter table
  tab <- read.delim(tabfile, header = TRUE, sep = "\t")
  tab_filt <- tab %>% filter(Identity >= 80 & Coverage >= 80)
  
  # Create headers matching the FASTA format (>gene_accession_1, etc.)
  filt_ids <- tab_filt %>%
    mutate(hdr = paste0(Resistance.gene, "_", Accession.no.)) %>%
    pull(hdr)
  
  # Read FASTA and select those passing the filter
  fa <- readDNAStringSet(fastafile, format="fasta")
  # Escape regex metacharacters in filt_ids
  filt_ids_escaped <- str_replace_all(filt_ids, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")
  
  # Now safely use in regex
  sel <- str_detect(names(fa), paste(filt_ids_escaped, collapse = "|"))
  writeXStringSet(fa[sel], outfile, format="fasta")
}
