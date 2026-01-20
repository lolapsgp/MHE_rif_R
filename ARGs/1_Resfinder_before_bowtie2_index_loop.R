# Load necessary library
library(Biostrings)

# Set root directory for all samples
FASTA_ROOT <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Resfinder_results"

# Get all sample directories under FASTA_ROOT
sample_dirs <- list.dirs(FASTA_ROOT, full.names = TRUE, recursive = FALSE)

# Function to rename duplicate sequence names in a FASTA file
rename_duplicates_fasta <- function(input_fasta, output_fasta) {
  # Read DNA sequences
  fasta <- readDNAStringSet(input_fasta)
  
  # Extract original sequence names
  names_orig <- names(fasta)
  
  # Create a named vector and count duplicates
  counts <- integer(length(names_orig))
  names(counts) <- names_orig
  
  # Track counts for each name
  count_tracker <- list()
  
  # New names vector
  new_names <- character(length(names_orig))
  
  for (i in seq_along(names_orig)) {
    name_i <- names_orig[i]
    # Increment count for each occurrence
    if (is.null(count_tracker[[name_i]])) {
      count_tracker[[name_i]] <- 1
    } else {
      count_tracker[[name_i]] <- count_tracker[[name_i]] + 1
    }
    # Append count suffix
    suffix <- paste0("_", count_tracker[[name_i]])
    new_names[i] <- paste0(name_i, suffix)
  }
  
  # Set new unique names back
  names(fasta) <- new_names
  
  # Write the updated FASTA file
  writeXStringSet(fasta, filepath = output_fasta)
  
  message(paste0("Processed file saved as: ", output_fasta))
}

for (sample in sample_names) {
  input_fasta <- file.path(sample_dirs, "ResFinder_Resistance_gene_seq.fsa")
  output_fasta <- file.path(sample_dirs, "ResFinder_Resistance_gene_seq_renamed.fsa")
  rename_duplicates_fasta(input_fasta, output_fasta)
}