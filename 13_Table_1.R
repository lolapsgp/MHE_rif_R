library(dplyr)
library(phyloseq)
library(tidyr)
library(purrr)
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Timepoint %in% c("T0"))
metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata_ok<-metadata[,c("DmGenderSexID", "Group_cutoff_4", "Age", "PHES", "Ammonia",
                                 "IL6", "Haemoglobin", "Leukocytes", "Lymphocytes","Neutrophils", "Monocytes","Absolute_Neutrophils",
                                 "Eosinophils","Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                                 "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                                 "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                                 "Sodium", "Potassium")]
rm(GMMs_corrected)
# Function to calculate statistics
calc_stats <- function(x) {
  return(c(median(x, na.rm = TRUE), 
           Q1 = quantile(x, 0.25, na.rm = TRUE), 
           Q3 = quantile(x, 0.75, na.rm = TRUE)))
}

# Function to calculate Wilcoxon p-value adjusted by BH
calc_wilcoxon_bh <- function(data, variable) {
  control_data <- data %>% filter(Group_cutoff_4 == "R") %>% pull(!!sym(variable))
  p_values <- data %>%
    filter(Group_cutoff_4 != "R") %>%
    group_by(Group_cutoff_4) %>%
    summarize(p_value = wilcox.test(!!sym(variable), control_data, na.rm = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  return(p_values)
}

# Function to convert p-values to asterisks
p_value_to_asterisk <- function(p) {
  if (p > 0.05) {
    return("")
  } else if (p <= 0.05 && p > 0.01) {
    return("*")
  } else if (p <= 0.01 && p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

# Calculate the statistics for each variable grouped by Group_cutoff_4
stats_table <- metadata_ok %>%
  group_by(Group_cutoff_4) %>%
  summarize(across(everything(), list(stats = ~ calc_stats(.)))) %>%
  pivot_longer(cols = -Group_cutoff_4, names_to = c("Variable", "Statistic"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = c(Group_cutoff_4, Statistic), values_from = Value, names_sep = "_")

# Get list of variables to iterate over
variables <- names(metadata_ok)[!names(metadata_ok) %in% "Group_cutoff_4"]

# Initialize an empty list to store final statistics
final_stats <- list()

# Loop through each variable to calculate stats and p-values
for (variable in variables) {
  # Calculate stats
  group_stats <- metadata_ok %>%
    group_by(Group_cutoff_4) %>%
    summarize(stats = list(calc_stats(!!sym(variable)))) %>%
    ungroup()
  
  # Calculate Wilcoxon p-values
  wilcoxon_results <- calc_wilcoxon_bh(metadata_ok, variable)
  
  # Integrate p-values into the stats lists
  group_stats <- group_stats %>%
    mutate(stats = map2(stats, Group_cutoff_4, ~ {
      if (.y == "R") {
        .x  # Return only stats for R
      } else {
        # Add the general adjusted p-value
        adj_p_value <- wilcoxon_results$adj_p_value[wilcoxon_results$Group_cutoff_4 == .y]
        stats_with_p_values <- c(.x, p_value_adj = round(adj_p_value, 4), p_value_to_asterisk(adj_p_value))
        
        stats_with_p_values
      }
    }))
  
  # Append the results to final_stats
  final_stats[[variable]] <- group_stats
}

# Combine all results into a final table
final_table <- bind_rows(final_stats, .id = "Variable")

# Optionally, reshape the final_table for better viewing
final_table_long <- final_table %>%
  pivot_wider(names_from = Group_cutoff_4, values_from = stats, names_sep = "_")

# Display the final table
View(final_table_long)

#Save as character
library(tidyverse)
table_ok <- final_table_long %>%
  mutate(
    R = map_chr(R, ~ paste0(.[1], " (", .[2], ", ", .[3], ")")),
    NR = map_chr(NR, ~ paste0(.[1], " (", .[2], ", ", .[3], "), p_value_adj ", .[4], .[5])),
  ) %>%
  select(Variable, R, NR) %>% 
  data.frame()

table_ok <- final_table_long %>%
  mutate(
    R = map_chr(R, ~ paste0(
      round(.x[1], 2), " (", round(.x[2], 2), ", ", round(.x[3], 2), ")"
    )),
    NR = map_chr(NR, ~ paste0(
      round(as.numeric(.x[1]), 2), " (", round(as.numeric(.x[2]), 2), ", ", round(as.numeric(.x[3]), 2), ") ",
      .x[5]  # Keep p-value (already a string)
    ))
  ) %>%
  select(Variable, R, NR) %>%
  data.frame()
writexl::write_xlsx(table_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Tables/Table_1_t0.xlsx")



#For Timepoint 2
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Timepoint %in% c("T2"))
metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata_ok<-metadata[,c("DmGenderSexID", "Group_cutoff_4", "Age", "PHES", "Ammonia",
                         "IL6", "Haemoglobin", "Leukocytes", "Lymphocytes","Neutrophils", "Monocytes","Absolute_Neutrophils",
                         "Eosinophils","Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                         "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                         "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                         "Sodium", "Potassium")]
rm(GMMs_corrected)
# Function to calculate statistics
calc_stats <- function(x) {
  return(c(median(x, na.rm = TRUE), 
           Q1 = quantile(x, 0.25, na.rm = TRUE), 
           Q3 = quantile(x, 0.75, na.rm = TRUE)))
}

# Function to calculate Wilcoxon p-value adjusted by BH
calc_wilcoxon_bh <- function(data, variable) {
  control_data <- data %>% filter(Group_cutoff_4 == "R") %>% pull(!!sym(variable))
  p_values <- data %>%
    filter(Group_cutoff_4 != "R") %>%
    group_by(Group_cutoff_4) %>%
    summarize(p_value = wilcox.test(!!sym(variable), control_data, na.rm = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  return(p_values)
}

# Function to convert p-values to asterisks
p_value_to_asterisk <- function(p) {
  if (p > 0.05) {
    return("")
  } else if (p <= 0.05 && p > 0.01) {
    return("*")
  } else if (p <= 0.01 && p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

# Calculate the statistics for each variable grouped by Group_cutoff_4
stats_table <- metadata_ok %>%
  group_by(Group_cutoff_4) %>%
  summarize(across(everything(), list(stats = ~ calc_stats(.)))) %>%
  pivot_longer(cols = -Group_cutoff_4, names_to = c("Variable", "Statistic"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = c(Group_cutoff_4, Statistic), values_from = Value, names_sep = "_")

# Get list of variables to iterate over
variables <- names(metadata_ok)[!names(metadata_ok) %in% "Group_cutoff_4"]

# Initialize an empty list to store final statistics
final_stats <- list()

# Loop through each variable to calculate stats and p-values
for (variable in variables) {
  # Calculate stats
  group_stats <- metadata_ok %>%
    group_by(Group_cutoff_4) %>%
    summarize(stats = list(calc_stats(!!sym(variable)))) %>%
    ungroup()
  
  # Calculate Wilcoxon p-values
  wilcoxon_results <- calc_wilcoxon_bh(metadata_ok, variable)
  
  # Integrate p-values into the stats lists
  group_stats <- group_stats %>%
    mutate(stats = map2(stats, Group_cutoff_4, ~ {
      if (.y == "R") {
        .x  # Return only stats for R
      } else {
        # Add the general adjusted p-value
        adj_p_value <- wilcoxon_results$adj_p_value[wilcoxon_results$Group_cutoff_4 == .y]
        stats_with_p_values <- c(.x, p_value_adj = round(adj_p_value, 4), p_value_to_asterisk(adj_p_value))
        
        stats_with_p_values
      }
    }))
  
  # Append the results to final_stats
  final_stats[[variable]] <- group_stats
}

# Combine all results into a final table
final_table <- bind_rows(final_stats, .id = "Variable")

# Optionally, reshape the final_table for better viewing
final_table_long <- final_table %>%
  pivot_wider(names_from = Group_cutoff_4, values_from = stats, names_sep = "_")

# Display the final table
View(final_table_long)

#Save as character
library(tidyverse)
table_ok <- final_table_long %>%
  mutate(
    R = map_chr(R, ~ paste0(.[1], " (", .[2], ", ", .[3], ")")),
    NR = map_chr(NR, ~ paste0(.[1], " (", .[2], ", ", .[3], "), p_value_adj ", .[4], .[5])),
  ) %>%
  select(Variable, R, NR) %>% 
  data.frame()

table_ok <- final_table_long %>%
  mutate(
    R = map_chr(R, ~ paste0(
      round(.x[1], 2), " (", round(.x[2], 2), ", ", round(.x[3], 2), ")"
    )),
    NR = map_chr(NR, ~ paste0(
      round(as.numeric(.x[1]), 2), " (", round(as.numeric(.x[2]), 2), ", ", round(as.numeric(.x[3]), 2), ") ",
      .x[5]  # Keep p-value (already a string)
    ))
  ) %>%
  select(Variable, R, NR) %>%
  data.frame()
writexl::write_xlsx(table_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Tables/Table_1_t2.xlsx")


####################### Study ############################
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Timepoint %in% c("T0"))
metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata_ok<-metadata[,c("DmGenderSexID", "Study", "Age", "PHES", "Ammonia",
                         "IL6", "Haemoglobin", "Leukocytes", "Lymphocytes","Neutrophils", "Monocytes","Absolute_Neutrophils",
                         "Eosinophils","Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                         "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                         "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                         "Sodium", "Potassium")]
rm(GMMs_corrected)
# Function to calculate statistics
calc_stats <- function(x) {
  return(c(median(x, na.rm = TRUE), 
           Q1 = quantile(x, 0.25, na.rm = TRUE), 
           Q3 = quantile(x, 0.75, na.rm = TRUE)))
}

# Function to calculate Wilcoxon p-value adjusted by BH
calc_wilcoxon_bh <- function(data, variable) {
  control_data <- data %>% filter(Study == "Spain") %>% pull(!!sym(variable))
  p_values <- data %>%
    filter(Study != "Spain") %>%
    group_by(Study) %>%
    summarize(p_value = wilcox.test(!!sym(variable), control_data, na.rm = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  return(p_values)
}

# Function to convert p-values to asterisks
p_value_to_asterisk <- function(p) {
  if (p > 0.05) {
    return("")
  } else if (p <= 0.05 && p > 0.01) {
    return("*")
  } else if (p <= 0.01 && p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

# Calculate the statistics for each variable grouped by Study
stats_table <- metadata_ok %>%
  group_by(Study) %>%
  summarize(across(everything(), list(stats = ~ calc_stats(.)))) %>%
  pivot_longer(cols = -Study, names_to = c("Variable", "Statistic"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = c(Study, Statistic), values_from = Value, names_sep = "_")

# Get list of variables to iterate over
variables <- names(metadata_ok)[!names(metadata_ok) %in% "Study"]

# Initialize an empty list to store final statistics
final_stats <- list()

# Loop through each variable to calculate stats and p-values
for (variable in variables) {
  # Calculate stats
  group_stats <- metadata_ok %>%
    group_by(Study) %>%
    summarize(stats = list(calc_stats(!!sym(variable)))) %>%
    ungroup()
  
  # Calculate Wilcoxon p-values
  wilcoxon_results <- calc_wilcoxon_bh(metadata_ok, variable)
  
  # Integrate p-values into the stats lists
  group_stats <- group_stats %>%
    mutate(stats = map2(stats, Study, ~ {
      if (.y == "Spain") {
        .x  # Return only stats for R
      } else {
        # Add the general adjusted p-value
        adj_p_value <- wilcoxon_results$adj_p_value[wilcoxon_results$Study == .y]
        stats_with_p_values <- c(.x, p_value_adj = round(adj_p_value, 4), p_value_to_asterisk(adj_p_value))
        
        stats_with_p_values
      }
    }))
  
  # Append the results to final_stats
  final_stats[[variable]] <- group_stats
}

# Combine all results into a final table
final_table <- bind_rows(final_stats, .id = "Variable")

# Optionally, reshape the final_table for better viewing
final_table_long <- final_table %>%
  pivot_wider(names_from = Study, values_from = stats, names_sep = "_")

# Display the final table
View(final_table_long)

#Save as character
library(tidyverse)
table_ok <- final_table_long %>%
  mutate(
    Spain = map_chr(Spain, ~ paste0(.[1], " (", .[2], ", ", .[3], ")")),
    UK = map_chr(UK, ~ paste0(.[1], " (", .[2], ", ", .[3], "), p_value_adj ", .[4], .[5])),
  ) %>%
  select(Variable, Spain, UK) %>% 
  data.frame()

table_ok <- final_table_long %>%
  mutate(
    Spain = map_chr(Spain, ~ paste0(
      round(.x[1], 2), " (", round(.x[2], 2), ", ", round(.x[3], 2), ")"
    )),
    UK = map_chr(UK, ~ paste0(
      round(as.numeric(.x[1]), 2), " (", round(as.numeric(.x[2]), 2), ", ", round(as.numeric(.x[3]), 2), ") ",
      .x[5]  # Keep p-value (already a string)
    ))
  ) %>%
  select(Variable, Spain, UK) %>%
  data.frame()
writexl::write_xlsx(table_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Tables/Table_Study_t0.xlsx")



#For Timepoint 2
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Timepoint %in% c("T2"))
metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata_ok<-metadata[,c("DmGenderSexID", "Study", "Age", "PHES", "Ammonia",
                         "IL6", "Haemoglobin", "Leukocytes", "Lymphocytes","Neutrophils", "Monocytes","Absolute_Neutrophils",
                         "Eosinophils","Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                         "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                         "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                         "Sodium", "Potassium")]
rm(GMMs_corrected)
# Function to calculate statistics
calc_stats <- function(x) {
  return(c(median(x, na.rm = TRUE), 
           Q1 = quantile(x, 0.25, na.rm = TRUE), 
           Q3 = quantile(x, 0.75, na.rm = TRUE)))
}

# Function to calculate Wilcoxon p-value adjusted by BH
calc_wilcoxon_bh <- function(data, variable) {
  control_data <- data %>% filter(Study == "Spain") %>% pull(!!sym(variable))
  p_values <- data %>%
    filter(Study != "Spain") %>%
    group_by(Study) %>%
    summarize(p_value = wilcox.test(!!sym(variable), control_data, na.rm = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  return(p_values)
}

# Function to convert p-values to asterisks
p_value_to_asterisk <- function(p) {
  if (p > 0.05) {
    return("")
  } else if (p <= 0.05 && p > 0.01) {
    return("*")
  } else if (p <= 0.01 && p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

# Calculate the statistics for each variable grouped by Study
stats_table <- metadata_ok %>%
  group_by(Study) %>%
  summarize(across(everything(), list(stats = ~ calc_stats(.)))) %>%
  pivot_longer(cols = -Study, names_to = c("Variable", "Statistic"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = c(Study, Statistic), values_from = Value, names_sep = "_")

# Get list of variables to iterate over
variables <- names(metadata_ok)[!names(metadata_ok) %in% "Study"]

# Initialize an empty list to store final statistics
final_stats <- list()

# Loop through each variable to calculate stats and p-values
for (variable in variables) {
  # Calculate stats
  group_stats <- metadata_ok %>%
    group_by(Study) %>%
    summarize(stats = list(calc_stats(!!sym(variable)))) %>%
    ungroup()
  
  # Calculate Wilcoxon p-values
  wilcoxon_results <- calc_wilcoxon_bh(metadata_ok, variable)
  
  # Integrate p-values into the stats lists
  group_stats <- group_stats %>%
    mutate(stats = map2(stats, Study, ~ {
      if (.y == "Spain") {
        .x  # Return only stats for R
      } else {
        # Add the general adjusted p-value
        adj_p_value <- wilcoxon_results$adj_p_value[wilcoxon_results$Study == .y]
        stats_with_p_values <- c(.x, p_value_adj = round(adj_p_value, 4), p_value_to_asterisk(adj_p_value))
        
        stats_with_p_values
      }
    }))
  
  # Append the results to final_stats
  final_stats[[variable]] <- group_stats
}

# Combine all results into a final table
final_table <- bind_rows(final_stats, .id = "Variable")

# Optionally, reshape the final_table for better viewing
final_table_long <- final_table %>%
  pivot_wider(names_from = Study, values_from = stats, names_sep = "_")

# Display the final table
View(final_table_long)

#Save as character
library(tidyverse)
table_ok <- final_table_long %>%
  mutate(
    Spain = map_chr(Spain, ~ paste0(.[1], " (", .[2], ", ", .[3], ")")),
    UK = map_chr(UK, ~ paste0(.[1], " (", .[2], ", ", .[3], "), p_value_adj ", .[4], .[5])),
  ) %>%
  select(Variable, Spain, UK) %>% 
  data.frame()

table_ok <- final_table_long %>%
  mutate(
    Spain = map_chr(Spain, ~ paste0(
      round(.x[1], 2), " (", round(.x[2], 2), ", ", round(.x[3], 2), ")"
    )),
    UK = map_chr(UK, ~ paste0(
      round(as.numeric(.x[1]), 2), " (", round(as.numeric(.x[2]), 2), ", ", round(as.numeric(.x[3]), 2), ") ",
      .x[5]  # Keep p-value (already a string)
    ))
  ) %>%
  select(Variable, Spain, UK) %>%
  data.frame()
writexl::write_xlsx(table_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Tables/Table_Study_t2.xlsx")

########################### Timepoint #################################

library(dplyr)
library(phyloseq)
library(tidyr)
library(purrr)
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Group_cutoff_4 %in% c("R"))
metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata_ok<-metadata[,c("Timepoint", "Age", "PHES", "Ammonia",
                         "IL6", "Haemoglobin", "Leukocytes", "Lymphocytes","Neutrophils", "Monocytes","Absolute_Neutrophils",
                         "Eosinophils","Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                         "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                         "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                         "Sodium", "Potassium")]
rm(GMMs_corrected)
# Function to calculate statistics
calc_stats <- function(x) {
  return(c(median(x, na.rm = TRUE), 
           Q1 = quantile(x, 0.25, na.rm = TRUE), 
           Q3 = quantile(x, 0.75, na.rm = TRUE)))
}

# Function to calculate Wilcoxon p-value adjusted by BH
calc_wilcoxon_bh <- function(data, variable) {
  control_data <- data %>% filter(Timepoint == "T0") %>% pull(!!sym(variable))
  p_values <- data %>%
    filter(Timepoint != "T0") %>%
    group_by(Timepoint) %>%
    summarize(p_value = wilcox.test(!!sym(variable), control_data, na.rm = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  return(p_values)
}

# Function to convert p-values to asterisks
p_value_to_asterisk <- function(p) {
  if (p > 0.05) {
    return("")
  } else if (p <= 0.05 && p > 0.01) {
    return("*")
  } else if (p <= 0.01 && p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

# Calculate the statistics for each variable grouped by Timepoint
stats_table <- metadata_ok %>%
  group_by(Timepoint) %>%
  summarize(across(everything(), list(stats = ~ calc_stats(.)))) %>%
  pivot_longer(cols = -Timepoint, names_to = c("Variable", "Statistic"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = c(Timepoint, Statistic), values_from = Value, names_sep = "_")

# Get list of variables to iterate over
variables <- names(metadata_ok)[!names(metadata_ok) %in% "Timepoint"]

# Initialize an empty list to store final statistics
final_stats <- list()

# Loop through each variable to calculate stats and p-values
for (variable in variables) {
  # Calculate stats
  group_stats <- metadata_ok %>%
    group_by(Timepoint) %>%
    summarize(stats = list(calc_stats(!!sym(variable)))) %>%
    ungroup()
  
  # Calculate Wilcoxon p-values
  wilcoxon_results <- calc_wilcoxon_bh(metadata_ok, variable)
  
  # Integrate p-values into the stats lists
  group_stats <- group_stats %>%
    mutate(stats = map2(stats, Timepoint, ~ {
      if (.y == "T0") {
        .x  # Return only stats for T0
      } else {
        # Add the general adjusted p-value
        adj_p_value <- wilcoxon_results$adj_p_value[wilcoxon_results$Timepoint == .y]
        stats_with_p_values <- c(.x, p_value_adj = round(adj_p_value, 4), p_value_to_asterisk(adj_p_value))
        
        stats_with_p_values
      }
    }))
  
  # Append the results to final_stats
  final_stats[[variable]] <- group_stats
}

# Combine all results into a final table
final_table <- bind_rows(final_stats, .id = "Variable")

# Optionally, reshape the final_table for better viewing
final_table_long <- final_table %>%
  pivot_wider(names_from = Timepoint, values_from = stats, names_sep = "_")

# Display the final table
View(final_table_long)

#Save as character
library(tidyverse)
table_ok <- final_table_long %>%
  mutate(
    T0 = map_chr(T0, ~ paste0(.[1], " (", .[2], ", ", .[3], ")")),
    T2 = map_chr(T2, ~ paste0(.[1], " (", .[2], ", ", .[3], "), p_value_adj ", .[4], .[5])),
  ) %>%
  select(Variable, T0, T2) %>% 
  data.frame()

table_ok <- final_table_long %>%
  mutate(
    T0 = map_chr(T0, ~ paste0(
      round(.x[1], 2), " (", round(.x[2], 2), ", ", round(.x[3], 2), ")"
    )),
    T2 = map_chr(T2, ~ paste0(
      round(as.numeric(.x[1]), 2), " (", round(as.numeric(.x[2]), 2), ", ", round(as.numeric(.x[3]), 2), ") ",
      .x[5]  # Keep p-value (already a string)
    ))
  ) %>%
  select(Variable, T0, T2) %>%
  data.frame()
writexl::write_xlsx(table_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Tables/Table_Timepoint_R.xlsx")



#For NR
library(dplyr)
library(phyloseq)
library(tidyr)
library(purrr)
GMMs_corrected <- readRDS("/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/GMMs_corrected.Rds")
GMMs_corrected<- subset_samples(GMMs_corrected, Group_cutoff_4 %in% c("NR"))
metadata<-data.frame(sample_data(GMMs_corrected))
rownames(metadata)<-metadata$SampleID
metadata$DmGenderSexID[metadata$DmGenderSexID == 2] <- 0
metadata_ok<-metadata[,c("Timepoint", "Age", "PHES", "Ammonia",
                         "IL6", "Haemoglobin", "Leukocytes", "Lymphocytes","Neutrophils", "Monocytes","Absolute_Neutrophils",
                         "Eosinophils","Absolute_Lymphocytes", "Absolute_Monocytes", "Absolute_Eosinophils",
                         "INR", "Fibrinogen", "Urea", "Creatinine", "Total_bilirubin",
                         "Albumin", "AST", "ALT", "GGT", "Alkaline_phosphatase",
                         "Sodium", "Potassium")]
rm(GMMs_corrected)
# Function to calculate statistics
calc_stats <- function(x) {
  return(c(median(x, na.rm = TRUE), 
           Q1 = quantile(x, 0.25, na.rm = TRUE), 
           Q3 = quantile(x, 0.75, na.rm = TRUE)))
}

# Function to calculate Wilcoxon p-value adjusted by BH
calc_wilcoxon_bh <- function(data, variable) {
  control_data <- data %>% filter(Timepoint == "T0") %>% pull(!!sym(variable))
  p_values <- data %>%
    filter(Timepoint != "T0") %>%
    group_by(Timepoint) %>%
    summarize(p_value = wilcox.test(!!sym(variable), control_data, na.rm = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  return(p_values)
}

# Function to convert p-values to asterisks
p_value_to_asterisk <- function(p) {
  if (p > 0.05) {
    return("")
  } else if (p <= 0.05 && p > 0.01) {
    return("*")
  } else if (p <= 0.01 && p > 0.001) {
    return("**")
  } else {
    return("***")
  }
}

# Calculate the statistics for each variable grouped by Timepoint
stats_table <- metadata_ok %>%
  group_by(Timepoint) %>%
  summarize(across(everything(), list(stats = ~ calc_stats(.)))) %>%
  pivot_longer(cols = -Timepoint, names_to = c("Variable", "Statistic"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = c(Timepoint, Statistic), values_from = Value, names_sep = "_")

# Get list of variables to iterate over
variables <- names(metadata_ok)[!names(metadata_ok) %in% "Timepoint"]

# Initialize an empty list to store final statistics
final_stats <- list()

# Loop through each variable to calculate stats and p-values
for (variable in variables) {
  # Calculate stats
  group_stats <- metadata_ok %>%
    group_by(Timepoint) %>%
    summarize(stats = list(calc_stats(!!sym(variable)))) %>%
    ungroup()
  
  # Calculate Wilcoxon p-values
  wilcoxon_results <- calc_wilcoxon_bh(metadata_ok, variable)
  
  # Integrate p-values into the stats lists
  group_stats <- group_stats %>%
    mutate(stats = map2(stats, Timepoint, ~ {
      if (.y == "T0") {
        .x  # Return only stats for T0
      } else {
        # Add the general adjusted p-value
        adj_p_value <- wilcoxon_results$adj_p_value[wilcoxon_results$Timepoint == .y]
        stats_with_p_values <- c(.x, p_value_adj = round(adj_p_value, 4), p_value_to_asterisk(adj_p_value))
        
        stats_with_p_values
      }
    }))
  
  # Append the results to final_stats
  final_stats[[variable]] <- group_stats
}

# Combine all results into a final table
final_table <- bind_rows(final_stats, .id = "Variable")

# Optionally, reshape the final_table for better viewing
final_table_long <- final_table %>%
  pivot_wider(names_from = Timepoint, values_from = stats, names_sep = "_")

# Display the final table
View(final_table_long)

#Save as character
library(tidyverse)
table_ok <- final_table_long %>%
  mutate(
    T0 = map_chr(T0, ~ paste0(.[1], " (", .[2], ", ", .[3], ")")),
    T2 = map_chr(T2, ~ paste0(.[1], " (", .[2], ", ", .[3], "), p_value_adj ", .[4], .[5])),
  ) %>%
  select(Variable, T0, T2) %>% 
  data.frame()

table_ok <- final_table_long %>%
  mutate(
    T0 = map_chr(T0, ~ paste0(
      round(.x[1], 2), " (", round(.x[2], 2), ", ", round(.x[3], 2), ")"
    )),
    T2 = map_chr(T2, ~ paste0(
      round(as.numeric(.x[1]), 2), " (", round(as.numeric(.x[2]), 2), ", ", round(as.numeric(.x[3]), 2), ") ",
      .x[5]  # Keep p-value (already a string)
    ))
  ) %>%
  select(Variable, T0, T2) %>%
  data.frame()
writexl::write_xlsx(table_ok, "/fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/MHE_rif/outputs/merged/Tables/Table_Timepoint_NR.xlsx")


