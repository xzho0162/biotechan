# CLAUDE.md

## Project Overview

**biotechan** is a CPTAC colorectal cancer proteomics analysis project. It implements a complete TMT-labelled quantitative proteomics pipeline in R: preprocessing, exploratory visualization, and differential expression analysis between matched normal and tumour tissue.

## Project Structure

```
biotechan/
  R/
    01_preprocessing.R           # Data cleaning, metadata extraction, matrix alignment
    02_exploration.R             # Descriptive stats, distribution plots, batch assessment
    03_differential_expression.R # Welch t-test, paired t-test, volcano/MA plots
  data/
    raw/                         # Raw CPTAC Excel files (not in git)
    processed/                   # Cleaned CSV outputs (not in git)
  output/
    figures/                     # Publication-ready plots (not in git)
    results/                     # DE results tables (not in git)
  .gitignore
  README.md
  CLAUDE.md
```

## Tech Stack

- **Language**: R (>= 4.0)
- **Packages**: readxl, dplyr, stringr, ggplot2, tidyr

## Development Guidelines

### R Conventions

- Scripts are numbered in execution order (`01_`, `02_`, `03_`)
- Use `check.names = FALSE` when reading CSVs with non-standard column names
- Use `[[ ]]` for column access when names contain spaces (e.g., `expr[["Protein Group"]]`)
- All paths are relative — `data/raw`, `data/processed`, `output/`

### Data Conventions

- Raw data files go in `data/raw/` and are never modified
- Processed outputs go in `data/processed/`
- Figures and result CSVs go in `output/`
- None of these are committed to git (see `.gitignore`)

### Git Conventions

- Clear, descriptive commit messages
- Do not commit data files, Excel files, or generated outputs
- Do not commit secrets or credentials

## Running the Pipeline

```r
source("R/01_preprocessing.R")
source("R/02_exploration.R")
source("R/03_differential_expression.R")
```

Place raw CPTAC Excel files in `data/raw/` before running.

## Key Analysis Parameters

- **Transformation**: log2(expression + 1)
- **Unpaired test**: Welch t-test (Normal n=21 vs Tumor n=33)
- **Paired test**: Paired t-test (9 matched patients)
- **FDR correction**: Benjamini-Hochberg
- **Significance**: FDR < 0.05, |log2FC| > 1
