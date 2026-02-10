# ============================================================================
# STEP 2: PREPROCESSING - CPTAC COLORECTAL CANCER PROTEOMICS (R Version)
# ============================================================================
# 
# Objective: Clean sample metadata, prepare expression matrix, document decisions
# Scope: Data cleaning only (NO deep analysis, filtering, or visualization)
# Output: 3 clean CSV files + summary statistics for report writing
#
# ============================================================================

library(readxl)
library(dplyr)
library(stringr)

# SETUP
# ============================================================================
base_dir <- "~/Downloads/bm "  # Note: space after "bm"
sample_file <- file.path(base_dir, "BMS5021_CPTAC_sample_information.xlsx")
data_file <- file.path(base_dir, "BMS5021_CPTAC_data.xlsx")

cat(paste(rep("=", 80), collapse=""), "\n")
cat("STEP 2: PREPROCESSING - CPTAC DATA\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# 1. LOAD RAW DATA
# ============================================================================
cat("1. Loading raw data...\n")
sample_df <- read_excel(sample_file, sheet = "Sheet1", col_names = FALSE)
data_df <- read_excel(data_file, sheet = "Results")

cat(sprintf("   Sample info table: %d rows × %d cols\n", nrow(sample_df), ncol(sample_df)))
cat(sprintf("   Expression matrix: %d proteins × %d cols\n\n", nrow(data_df), ncol(data_df)))

# 2. CLEAN SAMPLE INFORMATION TABLE
# ============================================================================
cat("2. Cleaning sample information...\n")

# Extract TMT labels (row 2)
tmt_labels <- as.character(sample_df[2, 3:ncol(sample_df)])

# Build mapping (start from row 3)
mapping <- tibble()

for (i in 3:(nrow(sample_df))) {
  set_id <- sample_df[[i, 2]]
  
  if (is.na(set_id)) next
  
  for (j in seq_along(tmt_labels)) {
    raw_val <- sample_df[[i, j + 2]]
    
    if (is.na(raw_val)) {
      tissue_type <- NA
      patient_id <- NA
    } else {
      # Clean format artifacts (_x000D_, \n)
      clean <- str_remove_all(as.character(raw_val), "_x000D_|\\\\n") %>% 
        str_trim()
      
      # Identify tissue type
      if (str_detect(clean, "ColonRef|Internal")) {
        tissue_type <- "Internal_Standard"
        patient_id <- "ColonRef"
      } else {
        patient_id <- str_extract(clean, "^[A-Z0-9]+")
        tissue_type <- if_else(
          str_detect(clean, "Normal|Solid"), 
          "Normal", 
          if_else(str_detect(clean, "Tumor|Tumour|Primary"), "Tumor", NA_character_)
        )
      }
    }
    
    col_name <- sprintf("%s_%s", set_id, tmt_labels[j])
    mapping <- bind_rows(mapping, tibble(
      SampleSet = set_id,
      TMT_Channel = tmt_labels[j],
      Patient_ID = patient_id,
      Tissue_Type = tissue_type,
      Column_Name = col_name
    ))
  }
}

# Remove internal standards (TMT_131)
mapping_clean <- mapping %>% filter(Tissue_Type != "Internal_Standard")

n_normal <- sum(mapping_clean$Tissue_Type == "Normal")
n_tumor <- sum(mapping_clean$Tissue_Type == "Tumor")
n_standards <- sum(mapping$Tissue_Type == "Internal_Standard")

cat(sprintf("   Normal: %d | Tumor: %d | Standards (removed): %d\n", 
            n_normal, n_tumor, n_standards))
cat(sprintf("   Cleaned metadata: %d samples\n\n", nrow(mapping_clean)))

# 3. PREPARE EXPRESSION MATRIX
# ============================================================================
cat("3. Preparing expression matrix...\n")

# Select only the TMT columns from mapping_clean
expr_cols <- intersect(mapping_clean$Column_Name, names(data_df))
id_cols <- c("Protein Group", "Accession")

expr_clean <- data_df %>% select(all_of(c(id_cols, expr_cols)))

n_proteins <- nrow(expr_clean)
n_samples <- length(expr_cols)

cat(sprintf("   Proteins: %d\n", n_proteins))
cat(sprintf("   Samples: %d (54 = 60 total - 6 internal standards)\n\n", n_samples))

# 4. ASSESS ZERO VALUES (basic check)
# ============================================================================
cat("4. Assessing zero values...\n")

expr_vals <- as.matrix(expr_clean[, -c(1:2)])
total_data_points <- prod(dim(expr_vals))
total_zeros <- sum(expr_vals == 0, na.rm = TRUE)
zero_pct <- 100 * total_zeros / total_data_points

proteins_all_zeros <- sum(rowSums(expr_vals == 0, na.rm = TRUE) == ncol(expr_vals))

cat(sprintf("   Total data points: %d\n", total_data_points))
cat(sprintf("   Zero values: %d (%.1f%%)\n", total_zeros, zero_pct))
cat(sprintf("   Proteins with all zeros: %d\n\n", proteins_all_zeros))

# 5. IDENTIFY PAIRED SAMPLES (for reference)
# ============================================================================
cat("5. Identifying paired samples...\n")

paired_count <- 0
for (patient in unique(mapping_clean$Patient_ID)) {
  patient_data <- mapping_clean %>% filter(Patient_ID == patient)
  tissue_types <- n_distinct(patient_data$Tissue_Type)
  
  if (tissue_types > 1) {
    paired_count <- paired_count + 1
  }
}

cat(sprintf("   Paired patients (both Normal & Tumor): %d\n\n", paired_count))

# 6. SAVE CLEANED FILES
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Saving cleaned files...\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# File 1: Sample metadata
write.csv(mapping_clean, 
          file.path(base_dir, "CPTAC_Step2_sample_metadata.csv"), 
          row.names = FALSE)
cat("✓ CPTAC_Step2_sample_metadata.csv\n")

# File 2: Expression matrix (proteins + 54 samples)
write.csv(expr_clean, 
          file.path(base_dir, "CPTAC_Step2_expression_matrix_cleaned.csv"), 
          row.names = FALSE)
cat("✓ CPTAC_Step2_expression_matrix_cleaned.csv\n")

# File 3: Column annotation
col_annot <- mapping_clean %>% 
  select(Column_Name, SampleSet, TMT_Channel, Patient_ID, Tissue_Type)
write.csv(col_annot, 
          file.path(base_dir, "CPTAC_Step2_column_annotation.csv"), 
          row.names = FALSE)
cat("✓ CPTAC_Step2_column_annotation.csv\n\n")

# 7. SUMMARY FOR REPORT
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("PREPROCESSING SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

summary_df <- tibble(
  Metric = c(
    "Raw proteins",
    "Clean proteins",
    "Total samples",
    "Normal samples",
    "Tumor samples",
    "Paired patients",
    "Zero values",
    "Zero percentage (%)"
  ),
  Value = c(
    nrow(data_df),
    n_proteins,
    n_samples,
    n_normal,
    n_tumor,
    paired_count,
    total_zeros,
    round(zero_pct, 2)
  )
)

print(summary_df)

cat("\n✅ STEP 2 COMPLETE\n")
cat(sprintf("Output location: %s\n\n", base_dir))
```

