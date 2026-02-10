# ============================================================================
# STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS - CPTAC COLORECTAL CANCER
# ============================================================================
#
# ANALYSIS DECISIONS & JUSTIFICATION:
#
# Method Selection:
# • Unpaired Welch t-test: For all 54 samples (Normal n=21 vs Tumor n=33)
#   Rationale: Handles unequal sample sizes and variances
# • Paired t-test: For 9 paired patients (within-subject comparison)
#   Rationale: Reduces inter-individual variation, increases statistical power
#
# Data Transformation:
# • log2-transformation: log2(expression + 1)
#   Rationale: Normalize wide expression range, stabilize variance for t-tests
#
# Multiple Testing Correction:
# • Benjamini-Hochberg (BH) FDR control
#   Rationale: Standard in omics analyses, balances Type I/II error
#
# Significance Thresholds:
# • |log2FC| > 1 (≥2-fold change)
# • FDR < 0.05 (multiple testing corrected)
#   Rationale: Conservative thresholds for biological relevance
#
# ============================================================================

library(dplyr)
library(ggplot2)

# SETUP
# ============================================================================
base_dir   <- "~/Downloads/bm "
expr_file  <- file.path(base_dir, "CPTAC_Step2_expression_matrix_cleaned.csv")
meta_file  <- file.path(base_dir, "CPTAC_Step2_sample_metadata.csv")

cat(paste(rep("=", 80), collapse=""), "\n")
cat("STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# 1. LOAD & PREPARE DATA
# ============================================================================
cat("1. Loading and preparing data...\n")

expr <- read.csv(expr_file, check.names = FALSE)
meta <- read.csv(meta_file)

cat(sprintf("   Expression matrix: %d proteins × %d cols\n", 
            nrow(expr), ncol(expr)))
cat(sprintf("   Metadata: %d samples\n", nrow(meta)))

# 获取列名（第1-2列是蛋白注释，其余是样本）
protein_col <- "Protein Group"  # 注意有空格！
accession_col <- "Accession"
sample_cols_raw <- colnames(expr)[-c(1,2)]

# 只保留在 metadata 中出现的样本
sample_cols <- intersect(sample_cols_raw, meta$Column_Name)

# 准备表达矩阵
expr_matrix <- expr[, c(protein_col, accession_col, sample_cols)]
expr_values <- expr_matrix[, sample_cols]
meta_use <- meta %>% filter(Column_Name %in% sample_cols)

# Apply the log2 transformation
expr_log2 <- log2(expr_values + 1)

normal_samples <- meta_use$Column_Name[meta_use$Tissue_Type == "Normal"]
tumor_samples  <- meta_use$Column_Name[meta_use$Tissue_Type == "Tumor"]

cat(sprintf("   Samples: Normal=%d, Tumor=%d\n", 
            length(normal_samples), length(tumor_samples)))
cat(sprintf("   Data transformation: log2(expression + 1)\n"))
cat(sprintf("✓ Ready for statistical analysis\n\n"))

# 2. UNPAIRED DIFFERENTIAL EXPRESSION (Welch t-test on all 54 samples)
# ============================================================================
cat("2. UNPAIRED ANALYSIS (Welch t-test, all 54 samples)\n")
cat(paste(rep("-", 80), collapse=""), "\n")
cat(sprintf("Testing: Tumor (n=%d) vs Normal (n=%d)\n", 
            length(tumor_samples), length(normal_samples)))
cat("Method: Welch t-test (t = (mean1 - mean2) / SE difference)\n\n")

results_unpaired <- data.frame(
  Protein   = expr_matrix[[protein_col]],
  Accession = expr_matrix[[accession_col]],
  mean_normal = NA_real_,
  mean_tumor  = NA_real_,
  log2FC      = NA_real_,
  pvalue      = NA_real_
)

for (i in seq_len(nrow(expr_log2))) {
  normal_vals <- as.numeric(expr_log2[i, normal_samples])
  tumor_vals  <- as.numeric(expr_log2[i, tumor_samples])
  
  results_unpaired$mean_normal[i] <- mean(normal_vals, na.rm = TRUE)
  results_unpaired$mean_tumor[i]  <- mean(tumor_vals, na.rm = TRUE)
  results_unpaired$log2FC[i] <- results_unpaired$mean_tumor[i] - results_unpaired$mean_normal[i]
      
      # Welch t-test
      tt <- t.test(tumor_vals, normal_vals, var.equal = FALSE)
      results_unpaired$pvalue[i] <- tt$p.value
    }
    
    # FDR correction
    results_unpaired$FDR <- p.adjust(results_unpaired$pvalue, method = "BH")
    
    # 统计显著蛋白
    sig_unpaired <- results_unpaired %>%
      filter(FDR < 0.05 & abs(log2FC) > 1)
    
    n_up <- sum(sig_unpaired$log2FC > 0)
    n_down <- sum(sig_unpaired$log2FC < 0)
    
    cat(sprintf("Total proteins tested: %d\n", nrow(results_unpaired)))
    cat(sprintf("Significant proteins (FDR<0.05, |log2FC|>1): %d\n", nrow(sig_unpaired)))
    cat(sprintf("  ↑ Upregulated in tumor: %d\n", n_up))
    cat(sprintf("  ↓ Downregulated in tumor: %d\n\n", n_down))
    
    # Save results
    write.csv(results_unpaired %>% arrange(FDR),
              file.path(base_dir, "CPTAC_Step4_DE_unpaired_all.csv"), 
              row.names = FALSE)
    write.csv(sig_unpaired %>% arrange(FDR),
              file.path(base_dir, "CPTAC_Step4_DE_unpaired_significant.csv"), 
              row.names = FALSE)
    
    cat("✓ Results saved:\n")
    cat("  - CPTAC_Step4_DE_unpaired_all.csv\n")
    cat("  - CPTAC_Step4_DE_unpaired_significant.csv\n\n")
    
    # 3. PAIRED DIFFERENTIAL EXPRESSION (9 paired patients)
    # ============================================================================
    cat("3. PAIRED ANALYSIS (Paired t-test, 9 paired patients)\n")
    cat(paste(rep("-", 80), collapse=""), "\n")
    cat("Method: Paired t-test (within-subject comparison)\n")
    cat("Rationale: Reduces inter-individual variation\n\n")
    
    # Build a pairing table: For each patient, select one Normal and one Tumor.
    pair_table <- meta_use %>%
      filter(Tissue_Type %in% c("Normal", "Tumor"),
             !is.na(Patient_ID)) %>%
      group_by(Patient_ID) %>%
      filter(n_distinct(Tissue_Type) == 2) %>%
      summarise(
        Normal = Column_Name[Tissue_Type == "Normal"][1],
        Tumor  = Column_Name[Tissue_Type == "Tumor"][1],
        .groups = "drop"
      )
    
    cat(sprintf("Paired patients identified: %d\n\n", nrow(pair_table)))
    
    # Paired t-test
    if (nrow(pair_table) > 0) {
      results_paired <- data.frame(
        Protein   = expr_matrix[["Protein Group"]],  # ✅ 修复：用双括号
        Accession = expr_matrix[["Accession"]],      # ✅ 修复：用双括号
        log2FC_paired = NA_real_,
        pvalue_paired = NA_real_
      )
      
      for (i in seq_len(nrow(expr_log2))) {
        vals_normal <- sapply(pair_table$Normal, function(col) expr_log2[i, col])
        vals_tumor  <- sapply(pair_table$Tumor,  function(col) expr_log2[i, col])
        
        # Paired t-test
        tt <- t.test(vals_tumor, vals_normal, paired = TRUE)
        
        results_paired$log2FC_paired[i] <- mean(vals_tumor - vals_normal)
        results_paired$pvalue_paired[i] <- tt$p.value
      }
      
      results_paired$FDR_paired <- p.adjust(results_paired$pvalue_paired, method = "BH")
      
      sig_paired <- results_paired %>%
        filter(FDR_paired < 0.05 & abs(log2FC_paired) > 1)
      
      n_up_paired <- sum(sig_paired$log2FC_paired > 0)
      n_down_paired <- sum(sig_paired$log2FC_paired < 0)
      
      cat(sprintf("Significant proteins (FDR<0.05, |log2FC|>1): %d\n", nrow(sig_paired)))
      cat(sprintf("  ↑ Upregulated in tumor: %d\n", n_up_paired))
      cat(sprintf("  ↓ Downregulated in tumor: %d\n\n", n_down_paired))
      
      # Save paired results
      write.csv(results_paired %>% arrange(FDR_paired),
                file.path(base_dir, "CPTAC_Step4_DE_paired_all.csv"), 
                row.names = FALSE)
      write.csv(sig_paired %>% arrange(FDR_paired),
                file.path(base_dir, "CPTAC_Step4_DE_paired_significant.csv"), 
                row.names = FALSE)
      
      cat("✓ Results saved:\n")
      cat("  - CPTAC_Step4_DE_paired_all.csv\n")
      cat("  - CPTAC_Step4_DE_paired_significant.csv\n\n")
      
      # Find proteins significant in BOTH analyses
      sig_both <- results_unpaired %>%
        left_join(results_paired, by = "Accession") %>%
        filter(FDR < 0.05 & abs(log2FC) > 1 &
                 FDR_paired < 0.05 & abs(log2FC_paired) > 1)
      
      cat(sprintf("Proteins significant in BOTH analyses: %d (highest confidence)\n\n", 
                  nrow(sig_both)))
    } else {
      cat("⚠ No paired patients found – paired analysis skipped.\n\n")
    }
    # 4. TOP DIFFERENTIALLY EXPRESSED PROTEINS
    # ============================================================================
    cat("4. TOP DIFFERENTIALLY EXPRESSED PROTEINS\n")
    cat(paste(rep("-", 80), collapse=""), "\n\n")
    
    top_up <- sig_unpaired %>%
      arrange(desc(log2FC)) %>%
      head(5)
    
    top_down <- sig_unpaired %>%
      arrange(log2FC) %>%
      head(5)
    
    cat("Top 5 Upregulated in Tumor:\n")
    print(top_up %>% select(Protein, log2FC, FDR))
    
    cat("\nTop 5 Downregulated in Tumor:\n")
    print(top_down %>% select(Protein, log2FC, FDR))
    cat("\n")
    
    # 5. CREATE VISUALIZATIONS
    # ============================================================================
    cat("5. Creating publication-ready visualizations...\n\n")
    
    output_dir <- base_dir
    
    # ---- FIGURE 1: Volcano Plot ----
    cat("Generating Volcano Plot...\n")
    
    volcano_df <- results_unpaired %>%
      mutate(
        negLog10P = -log10(pvalue),
        sig = ifelse(FDR < 0.05 & abs(log2FC) > 1, "Significant", "Not significant")
      )
    
    p_volcano <- ggplot(volcano_df, aes(x = log2FC, y = negLog10P, color = sig)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = c("Significant" = "#e74c3c", "Not significant" = "#95a5a6")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", size = 0.8) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", size = 0.8) +
      annotate("text", x = -6, y = max(volcano_df$negLog10P)*0.95, 
               label = paste0("Down: ", sum(sig_unpaired$log2FC < 0)), 
               fontface = "bold", size = 4, color = "#3498db") +
      annotate("text", x = 6, y = max(volcano_df$negLog10P)*0.95, 
               label = paste0("Up: ", sum(sig_unpaired$log2FC > 0)), 
               fontface = "bold", size = 4, color = "#e74c3c") +
      theme_minimal() +
      theme(legend.position = "bottom",
            panel.grid.major = element_line(color = "gray90")) +
      labs(
        title = "Volcano Plot: Tumor vs Normal (Unpaired Welch t-test)",
        x = "log2(Fold Change)",
        y = "-log10(p-value)",
        color = "Significance"
      )
    
    ggsave(file.path(output_dir, "CPTAC_Step4_Fig1_volcano_plot.png"),
           p_volcano, width = 9, height = 7, dpi = 300)
    
    cat("   ✓ Volcano Plot saved\n")
    
    # ---- FIGURE 2: MA Plot ----
    cat("Generating MA Plot...\n")
    
    ma_df <- results_unpaired %>%
      mutate(
        A = (log2(mean_normal + 1) + log2(mean_tumor + 1)) / 2,
        M = log2FC,
        sig = ifelse(FDR < 0.05 & abs(log2FC) > 1, "Significant", "Not significant")
      )
    
    p_ma <- ggplot(ma_df, aes(x = A, y = M, color = sig)) +
      geom_point(alpha = 0.5, size = 1.5) +
      scale_color_manual(values = c("Significant" = "#e74c3c", "Not significant" = "#95a5a6")) +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) +
      theme_minimal() +
      theme(panel.grid.major = element_line(color = "gray90")) +
      labs(
        title = "MA Plot: Average Expression vs Fold Change",
        x = "Average log2 Expression (A)",
        y = "log2(Fold Change) (M)",
        color = "Significance"
      )
    
    ggsave(file.path(output_dir, "CPTAC_Step4_Fig2_ma_plot.png"),
           p_ma, width = 8, height = 6, dpi = 300)
    
    cat("   ✓ MA Plot saved\n")
    
    # ---- FIGURE 3: log2FC Distribution ----
    cat("Generating log2FC Distribution...\n")
    
    p_hist <- ggplot(results_unpaired, aes(x = log2FC, fill = ifelse(abs(log2FC) > 1 & FDR < 0.05, "Significant", "Not significant"))) +
      geom_histogram(bins = 60, alpha = 0.8) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", size = 1) +
      scale_fill_manual(values = c("Significant" = "#e74c3c", "Not significant" = "#95a5a6")) +
      theme_minimal() +
      theme(panel.grid.major = element_line(color = "gray90")) +
      labs(
        title = "Distribution of log2(Fold Changes)",
        x = "log2(Fold Change)",
        y = "Frequency",
        fill = "Significance"
      )
    
    ggsave(file.path(output_dir, "CPTAC_Step4_Fig3_fc_distribution.png"),
           p_hist, width = 8, height = 5, dpi = 300)
    
    cat("   ✓ FC Distribution saved\n\n")
    
    # 6. SUMMARY TABLE
    # ============================================================================
    cat("6. ANALYSIS SUMMARY\n")
    cat(paste(rep("-", 80), collapse=""), "\n\n")
    
    summary_step4 <- tibble(
      Metric = c(
        "Total proteins analyzed",
        "Data transformation",
        "Statistical test (unpaired)",
        "Significant proteins (unpaired)",
        "  - Upregulated",
        "  - Downregulated",
        "Normal samples",
        "Tumor samples",
        "Paired patients",
        "Multiple testing correction",
        "FDR threshold",
        "log2FC threshold"
      ),
      Value = c(
        nrow(results_unpaired),
        "log2(expression + 1)",
        "Welch t-test",
        nrow(sig_unpaired),
        n_up,
        n_down,
        length(normal_samples),
        length(tumor_samples),
        ifelse(nrow(pair_table) > 0, nrow(pair_table), "N/A"),
        "Benjamini-Hochberg (BH)",
        "0.05",
        ">1"
      )
    )
    
    print(summary_step4)
    
    write.csv(summary_step4,
              file.path(base_dir, "CPTAC_Step4_summary.csv"),
              row.names = FALSE)
    
    cat("\n✅ STEP 4 COMPLETE\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("Output files:\n")
    cat("  ✓ CPTAC_Step4_DE_unpaired_all.csv\n")
    cat("  ✓ CPTAC_Step4_DE_unpaired_significant.csv\n")
    if (nrow(pair_table) > 0) {
      cat("  ✓ CPTAC_Step4_DE_paired_all.csv\n")
      cat("  ✓ CPTAC_Step4_DE_paired_significant.csv\n")
    }
    cat("  ✓ CPTAC_Step4_Fig1_volcano_plot.png\n")
    cat("  ✓ CPTAC_Step4_Fig2_ma_plot.png\n")
    cat("  ✓ CPTAC_Step4_Fig3_fc_distribution.png\n")
    cat("  ✓ CPTAC_Step4_summary.csv\n\n")
    cat(sprintf("Location: %s\n\n", base_dir))
    
    