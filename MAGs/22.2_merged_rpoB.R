#!/usr/bin/env Rscript

# ---------------------------
# USER INPUT
# ---------------------------
results_dir <- "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/rpoB_in_MAGs"
output_file <- file.path(results_dir, "rpoB_all_hits_MAGs.tsv")

# ---------------------------
# FIND FILES
# ---------------------------
presence_files <- list.files(
  path = results_dir,
  pattern = "_rpoB_presence\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# ---------------------------
# PROCESS FUNCTION
# ---------------------------
process_mag <- function(pres_file) {
  
  mag_name <- basename(dirname(pres_file))
  
  # -------- Presence --------
  presence_text <- readLines(pres_file, warn = FALSE)
  
  if (length(presence_text) == 0) {
    presence <- NA
  } else if (grepl("NOT detected", presence_text)) {
    presence <- 0
  } else if (grepl("detected", presence_text)) {
    presence <- 1
  } else {
    presence <- NA
  }
  
  # -------- Hits file --------
  hits_file <- file.path(dirname(pres_file),
                         paste0(mag_name, "_rpoB_hits.tsv"))
  
  # If missing or empty → return NA row
  if (!file.exists(hits_file) || file.size(hits_file) == 0) {
    
    return(data.frame(
      MAG = mag_name,
      qseqid = NA,
      sseqid = NA,
      pident = NA,
      length = NA,
      mismatch = NA,
      gapopen = NA,
      qstart = NA,
      qend = NA,
      sstart = NA,
      send = NA,
      evalue = NA,
      bitscore = NA,
      presence = presence,
      n_hits = 0
    ))
  }
  
  # -------- Read hits --------
  hits <- read.table(hits_file, header = FALSE, sep = "\t")
  
  colnames(hits) <- c(
    "qseqid","sseqid","pident","length","mismatch",
    "gapopen","qstart","qend","sstart","send",
    "evalue","bitscore"
  )
  
  n_hits <- nrow(hits)
  
  # Add metadata columns
  hits$MAG <- mag_name
  hits$presence <- presence
  hits$n_hits <- n_hits
  
  return(hits)
}

# ---------------------------
# APPLY
# ---------------------------
all_hits_list <- lapply(presence_files, process_mag)

all_hits <- do.call(rbind, all_hits_list)

# Reorder columns nicely
all_hits <- all_hits[, c(
  "MAG",
  "qseqid","sseqid","pident","length","mismatch",
  "gapopen","qstart","qend","sstart","send",
  "evalue","bitscore",
  "presence","n_hits"
)]

# ---------------------------
# WRITE OUTPUT
# ---------------------------
write.table(
  all_hits,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Merged hit table written to:", output_file, "\n")