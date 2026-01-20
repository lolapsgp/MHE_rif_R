# Load necessary library
library(fs)
setwd("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/01.RawData/")
# Function to get the current working directory
get_cwd <- function() {
  return(getwd())
}

# Function to list directories in the current working directory
list_directories <- function() {
  return(dir_ls(path = ".", type = "directory"))
}

# Function to list .fq.gz files in a directory
list_fq_files <- function(dir_path) {
  return(dir_ls(path = dir_path, glob = "*.fq.gz"))
}

create_txt_file <- function(output_file = "output.txt") {
  # Get current working directory
  cwd <- get_cwd()
  
  # Get list of directories
  directories <- list_directories()
  
  # Open file connection
  con <- file(output_file, "w")
  
  # Write basedir
  writeLines(paste0("basedir: ", cwd), con)
  
  # Write samples header
  writeLines("samples:", con)
  
  # Write directory names followed by their .fq.gz file names
  for (dir in directories) {
    dir_name <- basename(dir)
    writeLines(paste0("   ", dir_name, ":"), con)
    
    fq_files <- list_fq_files(dir)
    if (length(fq_files) > 0) {
      # Group files by their prefix (everything before _1.fq.gz or _2.fq.gz)
      file_groups <- split(fq_files, sub("_[12]\\.fq\\.gz$", "", basename(fq_files)))
      
      for (group in file_groups) {
        if (length(group) == 2) {
          writeLines("   - paired:", con)
          for (file in sort(basename(group))) {
            writeLines(paste0("      - ", file), con)
          }
        }
      }
    } else {
      writeLines("    No .fq.gz files found", con)
    }
  }
  
  # Close file connection
  close(con)
  
  cat("File", output_file, "has been created successfully.\n")
}

# Call the main function
create_txt_file("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/output.txt")
create_txt_file("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/output.yaml")

