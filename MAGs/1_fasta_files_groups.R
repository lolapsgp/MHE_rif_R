# Assuming your data frame is named metadata2
# Filter SampleIDs where Group_cutoff_4 is "R"
R_samples <- metadata2$SampleID[metadata2$Group_cutoff_4 == "R"]

# Filter SampleIDs where Group_cutoff_4 is "NR"
NR_samples <- metadata2$SampleID[metadata2$Group_cutoff_4 == "NR"]

# Write the SampleIDs to "R.txt"
writeLines(R_samples, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Responders.txt")

# Write the SampleIDs to "NR.txt"
writeLines(NR_samples, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/NonR.txt")


responders <- readLines("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Responders.txt")
nonr <- readLines("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/NonR.txt")
src_dirs <- c("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/UK", "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/Spain")

dir.create("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/Responders", showWarnings=FALSE)
dir.create("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/NonR", showWarnings=FALSE)

copy_fasta <- function(samples, target_dir) {
  for (sample in samples) {
    for (src_dir in src_dirs) {
      file_path <- file.path(src_dir, paste0(sample, ".fcontigs.fasta"))
      if (file.exists(file_path)) {
        file.copy(file_path, target_dir)
      }
    }
  }
}

copy_fasta(responders, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/Responders")
copy_fasta(nonr, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/assembled_filtered_contigs/NonR")
