#!/usr/bin/env Rscript

# 02-data-processing.R
# Dataset combination and comprehensive data quality analysis
# NOTE: Comprehensive summary printed at the end of script execution
#
# This script:
# 1. Sources datasets from 01-data-preprocessing.R (with backup redundancy)
# 2. Combines MMC1 (clinical) and MMC2 (gene expression) datasets
# 3. Performs comprehensive column completeness analysis
# 4. Stores results in R objects (no file exports)

# Load required libraries
suppressMessages(suppressWarnings({
  if (!require(tidyverse, quietly = TRUE)) {
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
  }
  if (!require(readxl, quietly = TRUE)) {
    install.packages("readxl", repos = "http://cran.us.r-project.org")
  }
  if (!require(here, quietly = TRUE)) {
    install.packages("here", repos = "http://cran.us.r-project.org")
  }

  library(tidyverse)
  library(readxl)
  library(readr)
  library(here)
}))

# ===============================================================================
# DATA SOURCE CONFIGURATION
# ===============================================================================

# Configure preferred data source for development
# Options: "auto" (full fallback chain), "preprocessing", "local", "github", "sciencedirect"
PREFERRED_DATA_SOURCE <- "auto" # Change this to control data loading

# ===============================================================================
# DATA LOADING VIA SOURCE SCRIPT (WITH BACKUP REDUNDANCY)
# ===============================================================================

# Initialise variables
clinical_data <- NULL
gene_expression_data <- NULL
data_source_type <- ""

# Data loading based on preference
if (PREFERRED_DATA_SOURCE == "local") {
  # Skip to local loading directly
  data_dir <- here("data")
  mmc1_file <- file.path(data_dir, "mmc1-ClinicalInformation.csv")
  mmc2_file <- file.path(data_dir, "mmc2-GeneExpression.csv")

  if (file.exists(mmc1_file) && file.exists(mmc2_file)) {
    clinical_data <- read_csv(mmc1_file, show_col_types = FALSE)
    gene_expression_data <- read_csv(mmc2_file, show_col_types = FALSE)

    if (
      nrow(clinical_data) == 1001 &&
        ncol(clinical_data) == 35 &&
        nrow(gene_expression_data) == 775 &&
        ncol(gene_expression_data) == 21
    ) {
      data_source_type <- "local"
      cat("✓ Using local CSV files\n")
    } else {
      stop("Local files have invalid dimensions")
    }
  } else {
    stop(
      "Local CSV files not found. Run 01-data-preprocessing.R first or change PREFERRED_DATA_SOURCE to 'auto'"
    )
  }
} else {
  # Use original fallback chain for other options

  # Attempt 1: Source from preprocessing script
  tryCatch(
    {
      suppressMessages(suppressWarnings({
        source(here("supplementary-materials", "01-data-preprocessing.R"))
      }))

      # Check if data was loaded successfully from preprocessing script
      data_dir <- here("data")
      mmc1_file <- file.path(data_dir, "mmc1-ClinicalInformation.csv")
      mmc2_file <- file.path(data_dir, "mmc2-GeneExpression.csv")

      if (file.exists(mmc1_file) && file.exists(mmc2_file)) {
        clinical_data <- read_csv(mmc1_file, show_col_types = FALSE)
        gene_expression_data <- read_csv(mmc2_file, show_col_types = FALSE)

        # Validate dimensions
        if (
          nrow(clinical_data) == 1001 &&
            ncol(clinical_data) == 35 &&
            nrow(gene_expression_data) == 775 &&
            ncol(gene_expression_data) == 21
        ) {
          data_source_type <- "preprocessing_script"
        } else {
          stop("Invalid dimensions from preprocessing script")
        }
      } else {
        stop("Preprocessing script did not create expected files")
      }
    },
    error = function(e) {
      # Backup redundancy: Direct data loading if preprocessing script fails
      suppressMessages(suppressWarnings({
        # MMC1 Clinical Data Loading
        direct_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc1.xlsx"
        repo_url <- "https://raw.githubusercontent.com/joelyu/HDS-W1/refs/heads/main/data/mmc1-ClinicalInformation.csv"

        clinical_data <<- NULL

        # Try ScienceDirect → GitHub → Local
        tryCatch(
          {
            temp_file <- tempfile(fileext = ".xlsx")
            download.file(direct_url, temp_file, mode = "wb", quiet = TRUE)

            if (file.exists(temp_file) && file.size(temp_file) > 0) {
              clinical_data <<- read_excel(
                temp_file,
                sheet = "Clinical Information",
                skip = 3
              )

              if (nrow(clinical_data) == 1001 && ncol(clinical_data) == 35) {
                data_source_type <<- "backup_original"
              } else {
                stop("Wrong dimensions")
              }
            }
          },
          error = function(e2) {
            tryCatch(
              {
                clinical_data <<- read_csv(repo_url, show_col_types = FALSE)
                if (nrow(clinical_data) == 1001 && ncol(clinical_data) == 35) {
                  data_source_type <<- "backup_github"
                } else {
                  stop("Wrong dimensions")
                }
              },
              error = function(e3) {
                local_path <- here("data", "mmc1-ClinicalInformation.csv")
                if (file.exists(local_path)) {
                  clinical_data <<- read_csv(local_path, show_col_types = FALSE)
                  if (
                    nrow(clinical_data) == 1001 && ncol(clinical_data) == 35
                  ) {
                    data_source_type <<- "backup_local"
                  }
                }
              }
            )
          }
        )

        # MMC2 Gene Expression Data Loading
        direct_url_mmc2 <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc2.xlsx"
        repo_url_mmc2 <- "https://raw.githubusercontent.com/joelyu/HDS-W1/refs/heads/main/data/mmc2-GeneExpression.csv"

        gene_expression_data <<- NULL

        tryCatch(
          {
            temp_file_mmc2 <- tempfile(fileext = ".xlsx")
            download.file(
              direct_url_mmc2,
              temp_file_mmc2,
              mode = "wb",
              quiet = TRUE
            )

            if (file.exists(temp_file_mmc2) && file.size(temp_file_mmc2) > 0) {
              gene_expression_data <<- read_excel(
                temp_file_mmc2,
                sheet = "Gene Expression",
                skip = 3
              )

              if (
                nrow(gene_expression_data) != 775 ||
                  ncol(gene_expression_data) != 21
              ) {
                stop("Wrong dimensions")
              }
            }
          },
          error = function(e2) {
            tryCatch(
              {
                gene_expression_data <<- read_csv(
                  repo_url_mmc2,
                  show_col_types = FALSE
                )
                if (
                  nrow(gene_expression_data) != 775 ||
                    ncol(gene_expression_data) != 21
                ) {
                  stop("Wrong dimensions")
                }
              },
              error = function(e3) {
                local_path_mmc2 <- here("data", "mmc2-GeneExpression.csv")
                if (file.exists(local_path_mmc2)) {
                  gene_expression_data <<- read_csv(
                    local_path_mmc2,
                    show_col_types = FALSE
                  )
                }
              }
            )
          }
        )

        # Final check
        if (is.null(clinical_data) || is.null(gene_expression_data)) {
          clinical_data <<- data.frame()
          gene_expression_data <<- data.frame()
          data_source_type <<- "all-failed"
        }
      }))
    }
  )
} # End of data source preference else block

# ===============================================================================
# INDIVIDUAL DATASET COMPLETENESS ANALYSIS
# ===============================================================================

# Clinical data completeness
clinical_completeness <- clinical_data %>%
  summarise_all(~ sum(!is.na(.))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "complete_cases"
  ) %>%
  mutate(
    total_cases = nrow(clinical_data),
    completeness_pct = round(100 * complete_cases / total_cases, 1),
    missing_cases = total_cases - complete_cases,
    dataset = "MMC1_Clinical"
  ) %>%
  arrange(desc(completeness_pct), variable)

# Gene expression data completeness
expression_completeness <- gene_expression_data %>%
  summarise_all(~ sum(!is.na(.))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "complete_cases"
  ) %>%
  mutate(
    total_cases = nrow(gene_expression_data),
    completeness_pct = round(100 * complete_cases / total_cases, 1),
    missing_cases = total_cases - complete_cases,
    dataset = "MMC2_Expression"
  ) %>%
  arrange(desc(completeness_pct), variable)

# ===============================================================================
# DATASET COMBINATION (REPLICATING RMD APPROACH)
# ===============================================================================

# Select key variables from clinical dataset (MMC1)
clinical_selected <- clinical_data %>%
  select(
    sample_id = `Sample  ID`,
    gender = Gender,
    response_therapy = `Response to initial therapy`,
    age_at_diagnosis = `age at diagnosis`,
    MYC = `log2 MYC expr`,
    BCL2 = `log2 BCL2 expr`,
    log2_bcl6_mmc1 = `log2 BCL6 expr`
  )

# Select gene expression data from MMC2 (exclude Samples and metadata columns)
expression_selected <- gene_expression_data %>%
  select(-Samples, -`RNA-seq core set of samples`) %>%
  rename(log2_bcl6_mmc2 = BCL6) %>%
  bind_cols(sample_id = gene_expression_data$Samples, .)

# Inner join to keep only patients with both clinical and expression data
combined_data <- clinical_selected %>%
  inner_join(expression_selected, by = "sample_id") %>%
  mutate(
    # Create binary therapy response variable
    complete_response = case_when(
      response_therapy == "Complete response" ~ "Complete response",
      response_therapy %in%
        c("Partial response", "No response") ~ "Incomplete response",
      TRUE ~ NA_character_
    ),
    complete_response = factor(
      complete_response,
      levels = c("Complete response", "Incomplete response")
    ),
    # Create age groups based on NCCN-IPI classification
    age_group_nccn = case_when(
      age_at_diagnosis <= 40 ~ "≤40 years",
      age_at_diagnosis <= 60 ~ "41-60 years",
      age_at_diagnosis <= 75 ~ "61-75 years",
      age_at_diagnosis > 75 ~ ">75 years",
      TRUE ~ NA_character_
    ),
    age_group_nccn = factor(
      age_group_nccn,
      levels = c("≤40 years", "41-60 years", "61-75 years", ">75 years")
    )
  ) %>%
  # Reorder columns for logical flow
  relocate(complete_response, .after = response_therapy) %>%
  relocate(age_group_nccn, .after = age_at_diagnosis)

# Validate BCL6 consistency
bcl6_validation <- combined_data %>%
  filter(!is.na(log2_bcl6_mmc1), !is.na(log2_bcl6_mmc2)) %>%
  mutate(bcl6_difference = abs(log2_bcl6_mmc1 - log2_bcl6_mmc2))

# Remove duplicate BCL6 column
combined_data <- combined_data %>%
  select(-log2_bcl6_mmc2) %>%
  rename(BCL6 = log2_bcl6_mmc1)

# ===============================================================================
# COMBINED DATASET COMPLETENESS ANALYSIS
# ===============================================================================

combined_completeness <- combined_data %>%
  summarise_all(~ sum(!is.na(.))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "complete_cases"
  ) %>%
  mutate(
    total_cases = nrow(combined_data),
    completeness_pct = round(100 * complete_cases / total_cases, 1),
    missing_cases = total_cases - complete_cases,
    dataset = "Combined",
    variable_type = case_when(
      variable == "sample_id" ~ "Identifier",
      variable %in%
        c(
          "gender",
          "response_therapy",
          "age_group_nccn"
        ) ~ "Clinical_Categorical",
      variable %in% c("age_at_diagnosis") ~ "Clinical_Continuous",
      variable %in% c("MYC", "BCL2", "BCL6") ~ "Primary_Genes",
      TRUE ~ "Expression_Genes"
    )
  ) %>%
  arrange(variable_type, desc(completeness_pct), variable)

# ===============================================================================
# SUMMARY STATISTICS
# ===============================================================================

completeness_by_type <- combined_completeness %>%
  group_by(variable_type) %>%
  summarise(
    num_variables = n(),
    mean_completeness = round(mean(completeness_pct), 1),
    min_completeness = min(completeness_pct),
    max_completeness = max(completeness_pct),
    variables_100_pct = sum(completeness_pct == 100),
    variables_95_plus = sum(completeness_pct >= 95),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_completeness))

key_vars_completeness <- combined_completeness %>%
  filter(
    variable %in%
      c(
        "sample_id",
        "gender",
        "age_at_diagnosis",
        "age_group_nccn",
        "response_therapy",
        "MYC",
        "BCL2",
        "BCL6"
      )
  ) %>%
  select(variable, complete_cases, total_cases, completeness_pct, missing_cases)

# Create comprehensive completeness report
full_completeness_report <- bind_rows(
  clinical_completeness,
  expression_completeness,
  combined_completeness
) %>%
  arrange(dataset, desc(completeness_pct), variable)

# Store data source information
data_sources <- list(
  mmc1_source = data_source_type,
  mmc2_source = data_source_type # Both datasets loaded from same source
)

# ===============================================================================
# STORE OBJECTS FOR POTENTIAL SOURCING
# ===============================================================================

# All key objects are now available in R environment:
# - clinical_data, gene_expression_data, combined_data
# - clinical_completeness, expression_completeness, combined_completeness
# - completeness_by_type, key_vars_completeness, full_completeness_report
# - bcl6_validation, data_sources
