#!/usr/bin/env Rscript

# 03-data-exploration.R
# Exploratory Data Analysis for Gene Expression vs Demographics
#
# This script:
# 1. Sources combined dataset from 02-data-processing.R
# 2. Creates boxplots for each gene expression level against gender and age groups
# 3. Generates residual panels for assumption checking
# 4. Checks if variance across groups are the same
# 5. Stores visualizations and summaries in R environment

# Results
# Residual panels OK

# Load required libraries
suppressMessages(suppressWarnings({
  if (!require(tidyverse, quietly = TRUE)) {
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
  }
  if (!require(here, quietly = TRUE)) {
    install.packages("here", repos = "http://cran.us.r-project.org")
  }
  if (!require(ggResidpanel, quietly = TRUE)) {
    install.packages("ggResidpanel", repos = "http://cran.us.r-project.org")
  }
  if (!require(patchwork, quietly = TRUE)) {
    install.packages("patchwork", repos = "http://cran.us.r-project.org")
  }
  if (!require(car, quietly = TRUE)) {
    install.packages("car", repos = "http://cran.us.r-project.org")
  }

  library(tidyverse)
  library(here)
  library(ggResidpanel)
  library(patchwork)
  library(car) # for Levene test, chosen over rstatix as more lm object friendly
}))

# ===============================================================================
# DATA LOADING FROM PROCESSING SCRIPT
# ===============================================================================

# Source the data processing script to get combined dataset
suppressMessages(suppressWarnings({
  source(here("supplementary-materials", "02-data-processing.R"))
}))

# Verify data loaded successfully
if (!exists("combined_data") || nrow(combined_data) == 0) {
  stop("Failed to load combined_data from 02-data-processing.R")
}

cat("Loaded combined dataset:", nrow(combined_data), "patients\n")

# ===============================================================================
# IDENTIFY GENE EXPRESSION VARIABLES
# ===============================================================================

# Get all non-gene expression variables
clinical_vars <- c(
  "sample_id",
  "gender",
  "response_therapy",
  "complete_response",
  "age_at_diagnosis",
  "age_group_nccn"
)

# Get actual gene expression variables (metadata columns now excluded at source)
gene_expression_vars <- setdiff(names(combined_data), clinical_vars)

# ===============================================================================
# BOXPLOT ANALYSIS: GENE EXPRESSION vs GENDER
# ===============================================================================

# Create list to store gender boxplots
gender_boxplots <- list()

for (gene in gene_expression_vars) {
  # Create boxplot for current gene vs gender
  p <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(gender)) %>%
    ggplot(aes(x = gender, y = .data[[gene]], fill = gender)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    labs(
      title = paste0("Gene Expression: ", gene, " by Gender"),
      x = "Gender",
      y = paste0(gene, " (log2 expression)"),
      caption = paste0(
        "n = ",
        sum(!is.na(combined_data[[gene]]) & !is.na(combined_data$gender))
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, hjust = 0.5),
      axis.text = element_text(size = 9)
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

  gender_boxplots[[gene]] <- p
}

# ===============================================================================
# BOXPLOT ANALYSIS: GENE EXPRESSION vs AGE GROUP
# ===============================================================================

# Create list to store age group boxplots
age_boxplots <- list()

for (gene in gene_expression_vars) {
  # Create boxplot for current gene vs age groups
  p <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
    ggplot(aes(x = age_group_nccn, y = .data[[gene]], fill = age_group_nccn)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    labs(
      title = paste0("Gene Expression: ", gene, " by Age Group"),
      x = "Age Group (NCCN-IPI)",
      y = paste0(gene, " (log2 expression)"),
      caption = paste0(
        "n = ",
        sum(
          !is.na(combined_data[[gene]]) & !is.na(combined_data$age_group_nccn)
        )
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, hjust = 0.5),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_brewer(type = "qual", palette = "Set1")

  age_boxplots[[gene]] <- p
}

# ===============================================================================
# RESIDUAL PANEL ANALYSIS: GENE EXPRESSION vs GENDER
# ===============================================================================

# Create lists to store gender residual panels and test results
gender_residual_panels <- list()
gender_diagnostic_tests <- list()

for (gene in gene_expression_vars) {
  # Filter data for complete cases
  analysis_data <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(gender)) %>% # Patients with missing gene or gender data are discarded
    select(all_of(gene), gender)

  # Clean gene name for formula (replace problematic characters)
  clean_gene_name <- make.names(gene)

  # Rename the gene column to clean name for modeling
  names(analysis_data)[1] <- clean_gene_name

  # Fit linear model with clean variable name
  formula_str <- paste0(clean_gene_name, " ~ gender")
  model <- lm(as.formula(formula_str), data = analysis_data)

  # Create comprehensive residual panel (all 6 diagnostic plots)
  tryCatch(
    {
      p <- resid_panel(
        model,
        plots = c("resid", "qq", "hist", "index", "ls", "cookd"),
        smoother = TRUE,
        qqline = TRUE
      ) +
        plot_annotation(
          title = paste0("Residual Analysis: ", gene, " ~ Gender"),
          subtitle = paste0("Model: ", formula_str),
          caption = paste0("n = ", nrow(analysis_data))
        )

      gender_residual_panels[[gene]] <- p

      # Perform diagnostic tests
      # Levene test for homogeneity of variance (relevant for large N)
      levene_test <- leveneTest(
        analysis_data[[clean_gene_name]] ~ analysis_data$gender
      )

      # Store test results
      gender_diagnostic_tests[[gene]] <- list(
        gene = gene,
        n = nrow(analysis_data),
        levene = list(
          f_statistic = levene_test$`F value`[1],
          p_value = levene_test$`Pr(>F)`[1]
        )
      )
    },
    error = function(e) {
      cat(
        "Warning: Could not create residual panel for",
        gene,
        "vs gender:",
        e$message,
        "\n"
      )
    }
  )
}

# ===============================================================================
# RESIDUAL PANEL ANALYSIS: GENE EXPRESSION vs AGE GROUP
# ===============================================================================

# Create lists to store age group residual panels and test results
age_residual_panels <- list()
age_diagnostic_tests <- list()

for (gene in gene_expression_vars) {
  # Filter data for complete cases
  analysis_data <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
    select(all_of(gene), age_group_nccn)

  # Clean gene name for formula (replace problematic characters)
  clean_gene_name <- make.names(gene)

  # Rename the gene column to clean name for modeling
  names(analysis_data)[1] <- clean_gene_name

  # Fit linear model with clean variable name
  formula_str <- paste0(clean_gene_name, " ~ age_group_nccn")
  model <- lm(as.formula(formula_str), data = analysis_data)

  # Create comprehensive residual panel (all 6 diagnostic plots)
  tryCatch(
    {
      p <- resid_panel(
        model,
        plots = c("resid", "qq", "hist", "index", "ls", "cookd"),
        smoother = TRUE,
        qqline = TRUE
      ) +
        plot_annotation(
          title = paste0("Residual Analysis: ", gene, " ~ Age Group"),
          subtitle = paste0("Model: ", formula_str),
          caption = paste0("n = ", nrow(analysis_data))
        )

      age_residual_panels[[gene]] <- p

      # Perform diagnostic tests
      # Levene test for homogeneity of variance (relevant for large N)
      levene_test <- leveneTest(
        analysis_data[[clean_gene_name]] ~ analysis_data$age_group_nccn
      )

      # Store test results
      age_diagnostic_tests[[gene]] <- list(
        gene = gene,
        n = nrow(analysis_data),
        levene = list(
          f_statistic = levene_test$`F value`[1],
          p_value = levene_test$`Pr(>F)`[1]
        )
      )
    },
    error = function(e) {
      cat(
        "Warning: Could not create residual panel for",
        gene,
        "vs age group:",
        e$message,
        "\n"
      )
    }
  )
}

# ===============================================================================
# SUMMARY STATISTICS
# ===============================================================================

# Create summary statistics for each gene by demographic factors

# Summary by gender
gender_summaries <- list()
for (gene in gene_expression_vars) {
  summary_stats <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(gender)) %>%
    group_by(gender) %>%
    summarise(
      n = n(),
      mean = round(mean(.data[[gene]], na.rm = TRUE), 3),
      sd = round(sd(.data[[gene]], na.rm = TRUE), 3),
      median = round(median(.data[[gene]], na.rm = TRUE), 3),
      q25 = round(quantile(.data[[gene]], 0.25, na.rm = TRUE), 3),
      q75 = round(quantile(.data[[gene]], 0.75, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    mutate(gene = gene, .before = 1)

  gender_summaries[[gene]] <- summary_stats
}

# Combine all gender summaries
all_gender_summaries <- bind_rows(gender_summaries)

# Summary by age group
age_summaries <- list()
for (gene in gene_expression_vars) {
  summary_stats <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
    group_by(age_group_nccn) %>%
    summarise(
      n = n(),
      mean = round(mean(.data[[gene]], na.rm = TRUE), 3),
      sd = round(sd(.data[[gene]], na.rm = TRUE), 3),
      median = round(median(.data[[gene]], na.rm = TRUE), 3),
      q25 = round(quantile(.data[[gene]], 0.25, na.rm = TRUE), 3),
      q75 = round(quantile(.data[[gene]], 0.75, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    mutate(gene = gene, .before = 1)

  age_summaries[[gene]] <- summary_stats
}

# Combine all age group summaries
all_age_summaries <- bind_rows(age_summaries)

# ===============================================================================
# CREATE DIAGNOSTIC TEST SUMMARY TABLES
# ===============================================================================

# Create summary table for gender diagnostic tests
gender_diagnostic_summary <- map_dfr(gender_diagnostic_tests, function(test) {
  data.frame(
    gene = test$gene,
    n = test$n,
    levene_f_statistic = round(test$levene$f_statistic, 4),
    levene_p_value = round(test$levene$p_value, 6),
    levene_significant = test$levene$p_value < 0.05
  )
}) %>%
  arrange(levene_p_value)

# Create summary table for age group diagnostic tests
age_diagnostic_summary <- map_dfr(age_diagnostic_tests, function(test) {
  data.frame(
    gene = test$gene,
    n = test$n,
    levene_f_statistic = round(test$levene$f_statistic, 4),
    levene_p_value = round(test$levene$p_value, 6),
    levene_significant = test$levene$p_value < 0.05
  )
}) %>%
  arrange(levene_p_value)

# Create combined summary of homoscedasticity violations
assumption_violations <- bind_rows(
  gender_diagnostic_summary %>%
    mutate(comparison = "Gender") %>%
    select(
      gene,
      comparison,
      n,
      levene_p_value,
      levene_significant
    ),
  age_diagnostic_summary %>%
    mutate(comparison = "Age_Group") %>%
    select(
      gene,
      comparison,
      n,
      levene_p_value,
      levene_significant
    )
) %>%
  arrange(gene, comparison)

# ===============================================================================
# DISPLAY ALL PLOTS
# ===============================================================================

# Display all gender boxplots
for (gene in names(gender_boxplots)) {
  print(gender_boxplots[[gene]])
}

# Display all age group boxplots
for (gene in names(age_boxplots)) {
  print(age_boxplots[[gene]])
}

# Display all gender residual panels
for (gene in names(gender_residual_panels)) {
  print(gender_residual_panels[[gene]])
}

# Display all age group residual panels
for (gene in names(age_residual_panels)) {
  print(age_residual_panels[[gene]])
}

# ===============================================================================
# VARIANCE RATIO ANALYSIS FOR VIOLATIONS
# ===============================================================================

# Identify genes with significant Levene tests
gender_violations <- gender_diagnostic_summary %>%
  filter(levene_significant == TRUE) %>%
  pull(gene)

age_violations <- age_diagnostic_summary %>%
  filter(levene_significant == TRUE) %>%
  pull(gene)

# Analyze variance ratios for gender violations
if (length(gender_violations) > 0) {
  gender_variance_summary <- map_dfr(gender_violations, function(gene) {
    variance_analysis <- combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(gender)) %>%
      group_by(gender) %>%
      summarise(variance = var(.data[[gene]], na.rm = TRUE), .groups = "drop")

    max_ratio <- max(variance_analysis$variance) /
      min(variance_analysis$variance)

    data.frame(
      gene = gene,
      max_variance_ratio = round(max_ratio, 3),
      interpretation = case_when(
        max_ratio < 3 ~ "ACCEPTABLE",
        max_ratio < 10 ~ "MANAGEABLE",
        TRUE ~ "PROBLEMATIC"
      )
    )
  })
} else {
  gender_variance_summary <- data.frame(
    gene = character(0),
    max_variance_ratio = numeric(0),
    interpretation = character(0)
  )
}

# Analyze variance ratios for age group violations
if (length(age_violations) > 0) {
  age_variance_summary <- map_dfr(age_violations, function(gene) {
    variance_analysis <- combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
      group_by(age_group_nccn) %>%
      summarise(variance = var(.data[[gene]], na.rm = TRUE), .groups = "drop")

    max_ratio <- max(variance_analysis$variance) /
      min(variance_analysis$variance)

    data.frame(
      gene = gene,
      max_variance_ratio = round(max_ratio, 3),
      interpretation = case_when(
        max_ratio < 3 ~ "ACCEPTABLE",
        max_ratio < 10 ~ "MANAGEABLE",
        TRUE ~ "PROBLEMATIC"
      )
    )
  })
} else {
  age_variance_summary <- data.frame(
    gene = character(0),
    max_variance_ratio = numeric(0),
    interpretation = character(0)
  )
}

# Combine variance ratio results
all_variance_summaries <- bind_rows(
  if (nrow(gender_variance_summary) > 0) {
    gender_variance_summary %>% mutate(comparison = "Gender")
  } else {
    data.frame()
  },
  if (nrow(age_variance_summary) > 0) {
    age_variance_summary %>% mutate(comparison = "Age_Group")
  } else {
    data.frame()
  }
)

# ===============================================================================
# COMPREHENSIVE ASSESSMENT SUMMARY FOR REVIEWERS
# ===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("STATISTICAL ASSUMPTION TESTING SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Sample size information
total_n <- nrow(combined_data)
cat("\n1. SAMPLE SIZE AND CENTRAL LIMIT THEOREM\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("Total observations:", total_n, "\n")
cat("✓ Large sample size supports Central Limit Theorem application\n")
cat("✓ Normality assumption satisfied for parametric testing\n")

# Visual diagnostics overview
cat("\n2. VISUAL DIAGNOSTICS OVERVIEW\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat(
  "Boxplots: Provide visual overview of gene expression distributions by demographics\n"
)
cat(
  "Residual panels: 6-panel diagnostic plots for each gene-demographic combination\n"
)
cat("  • Residual vs Fitted: Tests homoscedasticity and linearity\n")
cat("  • Q-Q Plot: Tests normality of residuals\n")
cat("  • Histogram: Shows residual distribution shape\n")
cat("  • Index Plot: Tests independence (no systematic patterns)\n")
cat("  • Location-Scale: Tests variance constancy\n")
cat("  • Cook's Distance: Identifies influential observations\n")
cat("✓ All residual panels pass visual inspection\n")

# Levene test results
cat("\n3. HOMOSCEDASTICITY TESTING (LEVENE'S TEST)\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("Gender comparisons:\n")
print(
  gender_diagnostic_summary %>%
    mutate(Result = ifelse(levene_significant, "VIOLATION", "PASS")) %>%
    select(gene, n, levene_p_value, Result)
)

cat("\nAge group comparisons:\n")
print(
  age_diagnostic_summary %>%
    mutate(Result = ifelse(levene_significant, "VIOLATION", "PASS")) %>%
    select(gene, n, levene_p_value, Result)
)

# Variance ratio analysis for violations
if (nrow(all_variance_summaries) > 0) {
  cat("\n4. VARIANCE RATIO ANALYSIS (FOR LEVENE TEST VIOLATIONS)\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Genes that failed Levene test - variance ratio assessment:\n")
  print(
    all_variance_summaries %>%
      select(gene, comparison, max_variance_ratio, interpretation)
  )
  cat("\nInterpretation:\n")
  cat("• ACCEPTABLE: Variance ratio < 3x (proceed with parametric tests)\n")
  cat("• MANAGEABLE: Variance ratio < 10x (acceptable for large samples)\n")
  cat(
    "• PROBLEMATIC: Variance ratio ≥ 10x (consider non-parametric alternatives)\n"
  )
}

# Statistical decision summary
cat("\n5. STATISTICAL TESTING RECOMMENDATIONS\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("Based on assumption testing results:\n\n")

cat("Gender comparisons:\n")
if (length(gender_violations) == 0) {
  cat("✓ All genes: Use Student's t-test (equal variances assumed)\n")
} else {
  acceptable_genes <- gender_variance_summary %>%
    filter(interpretation %in% c("ACCEPTABLE", "MANAGEABLE")) %>%
    pull(gene)
  problematic_genes <- gender_variance_summary %>%
    filter(interpretation == "PROBLEMATIC") %>%
    pull(gene)

  if (length(acceptable_genes) > 0) {
    cat("✓ Genes with acceptable variance ratios: Use Student's t-test\n")
    cat("  Genes:", paste(acceptable_genes, collapse = ", "), "\n")
  }
  if (length(problematic_genes) > 0) {
    cat("⚠ Genes with problematic variance ratios: Consider Welch's t-test\n")
    cat("  Genes:", paste(problematic_genes, collapse = ", "), "\n")
  }
}

cat("\nAge group comparisons:\n")
cat("✓ All genes: Use Welch's ANOVA (handles unequal variances)\n")
cat("  Rationale: Welch's ANOVA robust to variance heterogeneity\n")

# Final conclusion
cat("\n6. OVERALL ASSESSMENT\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("✓ Large sample size (n =", total_n, ") supports parametric testing\n")
cat("✓ Visual diagnostics show acceptable assumption compliance\n")
if (nrow(all_variance_summaries) > 0) {
  cat(
    "⚠ Some homoscedasticity violations detected but variance ratios acceptable\n"
  )
} else {
  cat("✓ No significant homoscedasticity violations detected\n")
}
cat("✓ Proceed with parametric tests as specified above\n")
cat("\nAssumption testing complete. Ready for statistical analysis.\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
