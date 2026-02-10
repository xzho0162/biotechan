# CPTAC Colorectal Cancer Proteomics Analysis

Differential protein expression analysis of colorectal cancer (CRC) using TMT-labelled quantitative proteomics data from the [Clinical Proteomic Tumor Analysis Consortium (CPTAC)](https://proteomics.cancer.gov/programs/cptac).

## Overview

This project implements a complete proteomics data analysis pipeline in R, covering data preprocessing, exploratory visualization, and differential expression testing between matched normal and tumour colorectal tissue.

**Key findings:**
- 1,656 differentially expressed proteins identified (unpaired Welch t-test, FDR < 0.05, |log2FC| > 1)
- 1,504 upregulated and 152 downregulated in tumour tissue
- 1,066 proteins validated across both unpaired and paired analyses (highest confidence)

## Dataset

- **Source:** CPTAC TMT-labelled proteomics (colorectal cancer cohort)
- **Samples:** 54 experimental samples (21 Normal, 33 Tumour) across 6 TMT batches
- **Proteins:** ~10,294 protein groups
- **Design:** 9 patients with matched normal/tumour tissue enabling paired analysis

## Pipeline

```
01_preprocessing.R       Data cleaning, metadata extraction, expression matrix alignment
02_exploration.R         Descriptive statistics, distribution plots, batch effect assessment
03_differential_expression.R   Welch t-test, paired t-test, volcano/MA plots
```

### 1. Preprocessing
- Parse TMT channel-to-sample mapping from raw Excel files
- Clean formatting artifacts (`_x000D_`, line breaks) from sample annotations
- Remove internal standards (TMT_131 / ColonRef)
- Validate column alignment between metadata and expression matrix

### 2. Exploration & Visualization
- Per-sample median intensity profiles (log scale)
- Normal vs Tumour global expression comparison
- Cross-batch consistency assessment (CV across 6 TMT sets)
- Sample composition summary

### 3. Differential Expression
- **Transformation:** log2(expression + 1)
- **Unpaired test:** Welch t-test across all 54 samples
- **Paired test:** Paired t-test for 9 matched patients
- **Correction:** Benjamini-Hochberg FDR
- **Thresholds:** FDR < 0.05, |log2FC| > 1
- Cross-validation of results between both statistical approaches

## Project Structure

```
R/
  01_preprocessing.R
  02_exploration.R
  03_differential_expression.R
data/
  raw/                  # Raw CPTAC Excel files (not included)
  processed/            # Cleaned CSV outputs
output/
  figures/              # Publication-ready plots
  results/              # DE results tables
```

## Requirements

- R >= 4.0
- Packages: `readxl`, `dplyr`, `stringr`, `ggplot2`, `tidyr`

## Usage

```r
# Run the pipeline in order:
source("R/01_preprocessing.R")
source("R/02_exploration.R")
source("R/03_differential_expression.R")
```

Place raw CPTAC Excel files in `data/raw/` before running. Output files and figures will be generated in `output/`.

## License

MIT
