# CPTAC Colorectal Cancer Proteomics Analysis

End-to-end differential protein expression analysis of colorectal cancer (CRC) using TMT-labelled quantitative proteomics data from the [Clinical Proteomic Tumor Analysis Consortium (CPTAC)](https://proteomics.cancer.gov/programs/cptac).

## Quick Project Pitch (HR-Friendly)

This is an **independently completed bioinformatics project** where I designed and implemented a full analysis workflow in R, from raw data cleaning to statistical testing and result visualization.

- **What I built:** A 3-stage reproducible analysis pipeline (preprocessing, exploration, differential expression)
- **Core methods:** Welch t-test + paired t-test, BH-FDR multiple testing correction
- **Outcome:** 1,656 significant proteins identified; 1,066 cross-validated in both unpaired and paired analyses
- **Skills demonstrated:** data wrangling, statistical analysis, visualization, research-style reporting, reproducible scripting

## Overview

This project implements a complete proteomics data analysis pipeline in R, covering data preprocessing, exploratory visualization, and differential expression testing between matched normal and tumour colorectal tissue.

## Key Findings

- **1,656** differentially expressed proteins identified (unpaired Welch t-test, FDR < 0.05, |log2FC| > 1)
- **1,504** upregulated and **152** downregulated in tumour tissue
- **1,066** proteins validated across both unpaired and paired analyses (highest confidence)

## My Role

I independently completed the full workflow, including:

1. Parsing and cleaning messy sample annotation from raw Excel files
2. Building metadata-expression alignment checks
3. Designing both unpaired and paired differential expression analyses
4. Applying multiple-testing control and significance thresholds
5. Producing publication-style figures and summary tables

## Dataset

- **Source:** CPTAC TMT-labelled proteomics (colorectal cancer cohort)
- **Samples:** 54 experimental samples (21 Normal, 33 Tumour) across 6 TMT batches
- **Proteins:** ~10,294 protein groups
- **Design:** 9 patients with matched normal/tumour tissue enabling paired analysis

## Pipeline

```text
01_preprocessing.R              Data cleaning, metadata extraction, expression matrix alignment
02_exploration.R                Descriptive statistics, distribution plots, batch effect assessment
03_differential_expression.R    Welch t-test, paired t-test, volcano/MA plots
```

### 1) Preprocessing
- Parse TMT channel-to-sample mapping from raw Excel files
- Clean formatting artifacts (`_x000D_`, line breaks) from sample annotations
- Remove internal standards (TMT_131 / ColonRef)
- Validate column alignment between metadata and expression matrix

### 2) Exploration & Visualization
- Per-sample median intensity profiles (log scale)
- Normal vs Tumour global expression comparison
- Cross-batch consistency assessment (CV across 6 TMT sets)
- Sample composition summary

### 3) Differential Expression
- **Transformation:** log2(expression + 1)
- **Unpaired test:** Welch t-test across all 54 samples
- **Paired test:** Paired t-test for 9 matched patients
- **Correction:** Benjamini-Hochberg FDR
- **Thresholds:** FDR < 0.05, |log2FC| > 1
- Cross-validation of results between both statistical approaches

## Project Structure

```text
R/
  01_preprocessing.R
  02_exploration.R
  03_differential_expression.R
data/
  raw/                  # Raw CPTAC Excel files (not included)
  processed/            # Cleaned CSV outputs
output/
  figures/              # Publication-ready plots
  results/              # Differential expression result tables
```

## Requirements

- R >= 4.0
- Packages: `readxl`, `dplyr`, `stringr`, `ggplot2`, `tidyr`

Install dependencies:

```r
install.packages(c("readxl", "dplyr", "stringr", "ggplot2", "tidyr"))
```

## Usage

```r
# Run the pipeline in order
source("R/01_preprocessing.R")
source("R/02_exploration.R")
source("R/03_differential_expression.R")
```

Place raw CPTAC Excel files in `data/raw/` before running. Output files and figures will be generated under `data/processed/` and `output/`.

## Notes for Recruiters / Non-Technical Reviewers

If you are reviewing this project from a hiring perspective, the most relevant points are:

- It demonstrates end-to-end ownership of a real-world analytical task
- It combines domain data processing with statistical rigor and clear reporting
- It includes traceable logic from raw input to final biological insights

## License

MIT
