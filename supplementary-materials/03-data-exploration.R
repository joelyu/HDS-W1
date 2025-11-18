#!/usr/bin/env Rscript

# 03-data-exploration.R
# Exploratory Data Analysis for Gene Expression vs Demographics
#
# This script:
# 1. Sources combined dataset from 02-data-processing.R
# 2. Creates boxplots for each gene expression level against gender and age groups
# 3. Generates residual panels for assumption checking
# 4. Stores visualizations and summaries in R environment

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

cat("Found", length(gene_expression_vars), "gene expression variables:\n")
cat(paste(gene_expression_vars, collapse = ", "), "\n\n")

# ===============================================================================
# BOXPLOT ANALYSIS: GENE EXPRESSION vs GENDER
# ===============================================================================

cat("Creating boxplots for gene expression vs gender...\n")

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

cat("Creating boxplots for gene expression vs age groups...\n")

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

cat("Creating residual panels for gene expression vs gender...\n")

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

      cat(
        "✓",
        gene,
        "- n =",
        nrow(analysis_data),
        "| Levene p =",
        round(levene_test$`Pr(>F)`[1], 4),
        "\n"
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

cat("Creating residual panels for gene expression vs age groups...\n")

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

      cat(
        "✓",
        gene,
        "- n =",
        nrow(analysis_data),
        "| Levene p =",
        round(levene_test$`Pr(>F)`[1], 4),
        "\n"
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
cat("Generating summary statistics...\n")

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

cat("\nDiagnostic Test Summary (Large N Analysis):\n")
cat("- Genes with homoscedasticity violations (Levene p < 0.05):\n")
variance_violations <- assumption_violations %>%
  filter(levene_significant) %>%
  group_by(gene) %>%
  summarise(violations = paste(comparison, collapse = ", "), .groups = "drop")
if (nrow(variance_violations) > 0) {
  for (i in 1:nrow(variance_violations)) {
    cat(
      "  ",
      variance_violations$gene[i],
      "(",
      variance_violations$violations[i],
      ")\n"
    )
  }
} else {
  cat("   None\n")
}
cat(
  "Note: Normality testing omitted for large N analysis (Central Limit Theorem applies)\n"
)

# ===============================================================================
# STORE RESULTS FOR POTENTIAL SOURCING
# ===============================================================================

cat("Analysis complete! Available objects:\n")
cat("VISUALIZATION OBJECTS:\n")
cat(
  "- gender_boxplots: List of",
  length(gender_boxplots),
  "boxplots (gene vs gender)\n"
)
cat(
  "- age_boxplots: List of",
  length(age_boxplots),
  "boxplots (gene vs age group)\n"
)
cat(
  "- gender_residual_panels: List of",
  length(gender_residual_panels),
  "residual panels (6-panel diagnostics: gene vs gender)\n"
)
cat(
  "- age_residual_panels: List of",
  length(age_residual_panels),
  "residual panels (6-panel diagnostics: gene vs age group)\n"
)
cat("\nDIAGNOSTIC TEST OBJECTS:\n")
cat(
  "- gender_diagnostic_tests: List of",
  length(gender_diagnostic_tests),
  "detailed test results (gene vs gender)\n"
)
cat(
  "- age_diagnostic_tests: List of",
  length(age_diagnostic_tests),
  "detailed test results (gene vs age group)\n"
)
cat(
  "- gender_diagnostic_summary: Data frame with Levene test results (gender)\n"
)
cat(
  "- age_diagnostic_summary: Data frame with Levene test results (age groups)\n"
)
cat(
  "- assumption_violations: Combined summary table of homoscedasticity violations\n"
)
cat("\nDESCRIPTIVE STATISTICS:\n")
cat("- all_gender_summaries: Summary statistics by gender\n")
cat("- all_age_summaries: Summary statistics by age group\n")
cat("- gene_expression_vars: Vector of gene expression variable names\n\n")

cat("KEY TABLES FOR ANALYSIS:\n")
cat(
  "View(gender_diagnostic_summary)  # Levene test results for gender models\n"
)
cat(
  "View(age_diagnostic_summary)     # Levene test results for age group models\n"
)
cat("View(assumption_violations)      # Homoscedasticity violations summary\n")

# ===============================================================================
# DISPLAY ALL PLOTS
# ===============================================================================

cat("Displaying all plots...\n\n")

# Display all gender boxplots
cat("=== GENDER BOXPLOTS ===\n")
for (gene in names(gender_boxplots)) {
  cat("Displaying:", gene, "vs Gender\n")
  print(gender_boxplots[[gene]])
  cat("\n")
}

# Display all age group boxplots
cat("=== AGE GROUP BOXPLOTS ===\n")
for (gene in names(age_boxplots)) {
  cat("Displaying:", gene, "vs Age Groups\n")
  print(age_boxplots[[gene]])
  cat("\n")
}

# Display all gender residual panels
cat("=== GENDER RESIDUAL PANELS ===\n")
for (gene in names(gender_residual_panels)) {
  cat("Displaying residual analysis:", gene, "~ Gender\n")
  print(gender_residual_panels[[gene]])
  cat("\n")
}

# Display all age group residual panels
cat("=== AGE GROUP RESIDUAL PANELS ===\n")
for (gene in names(age_residual_panels)) {
  cat("Displaying residual analysis:", gene, "~ Age Groups\n")
  print(age_residual_panels[[gene]])
  cat("\n")
}

cat("All plots displayed!\n")
cat("Summary data available in: all_gender_summaries, all_age_summaries\n")

# ===============================================================================
# FOLLOW-UP ANALYSIS: VARIANCE RATIO INVESTIGATION
# ===============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("INVESTIGATING HOMOSCEDASTICITY VIOLATIONS\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Identify genes with significant Levene tests
gender_violations <- gender_diagnostic_summary %>%
  filter(levene_significant == TRUE) %>%
  pull(gene)

age_violations <- age_diagnostic_summary %>%
  filter(levene_significant == TRUE) %>%
  pull(gene)

cat("Genes with significant Levene test for GENDER:\n")
if (length(gender_violations) > 0) {
  cat(paste(gender_violations, collapse = ", "), "\n")
} else {
  cat("None\n")
}

cat("\nGenes with significant Levene test for AGE GROUPS:\n")
if (length(age_violations) > 0) {
  cat(paste(age_violations, collapse = ", "), "\n")
} else {
  cat("None\n")
}

# ===============================================================================
# VARIANCE RATIO ANALYSIS FOR GENDER VIOLATIONS
# ===============================================================================

if (length(gender_violations) > 0) {
  cat(
    "Analyzing variance ratios for",
    length(gender_violations),
    "genes with gender violations...\n"
  )

  # Create detailed variance analysis table
  gender_variance_details <- map_dfr(gender_violations, function(gene) {
    combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(gender)) %>%
      group_by(gene = !!gene, gender) %>%
      summarise(
        n = n(),
        variance = var(.data[[gene]], na.rm = TRUE),
        sd = sd(.data[[gene]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variance_ratio = variance / min(variance))
  })

  # Create summary table
  gender_variance_summary <- map_dfr(gender_violations, function(gene) {
    variance_analysis <- combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(gender)) %>%
      group_by(gender) %>%
      summarise(
        variance = var(.data[[gene]], na.rm = TRUE),
        .groups = "drop"
      )

    max_ratio <- max(variance_analysis$variance) /
      min(variance_analysis$variance)

    interpretation <- case_when(
      max_ratio < 3 ~ "ACCEPTABLE",
      max_ratio < 10 ~ "MANAGEABLE",
      TRUE ~ "PROBLEMATIC"
    )

    data.frame(
      gene = gene,
      max_variance_ratio = round(max_ratio, 3),
      interpretation = interpretation,
      recommendation = case_when(
        max_ratio < 3 ~ "Continue with parametric tests",
        max_ratio < 10 ~ "Acceptable for large samples",
        TRUE ~ "Consider non-parametric alternatives"
      )
    )
  })

  gender_variance_results <- list(
    summary = gender_variance_summary,
    details = gender_variance_details
  )
} else {
  cat("✓ No significant gender-related homoscedasticity violations\n")
  gender_variance_summary <- data.frame(
    gene = character(0),
    max_variance_ratio = numeric(0),
    interpretation = character(0),
    recommendation = character(0)
  )
  gender_variance_details <- data.frame()
  gender_variance_results <- list(
    summary = gender_variance_summary,
    details = gender_variance_details
  )
}

# ===============================================================================
# VARIANCE RATIO ANALYSIS FOR AGE GROUP VIOLATIONS
# ===============================================================================

if (length(age_violations) > 0) {
  cat(
    "Analyzing variance ratios for",
    length(age_violations),
    "genes with age group violations...\n"
  )

  # Create detailed variance analysis table
  age_variance_details <- map_dfr(age_violations, function(gene) {
    combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
      group_by(gene = !!gene, age_group_nccn) %>%
      summarise(
        n = n(),
        variance = var(.data[[gene]], na.rm = TRUE),
        sd = sd(.data[[gene]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variance_ratio = variance / min(variance))
  })

  # Create summary table
  age_variance_summary <- map_dfr(age_violations, function(gene) {
    variance_analysis <- combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
      group_by(age_group_nccn) %>%
      summarise(
        variance = var(.data[[gene]], na.rm = TRUE),
        .groups = "drop"
      )

    max_ratio <- max(variance_analysis$variance) /
      min(variance_analysis$variance)

    interpretation <- case_when(
      max_ratio < 3 ~ "ACCEPTABLE",
      max_ratio < 10 ~ "MANAGEABLE",
      TRUE ~ "PROBLEMATIC"
    )

    data.frame(
      gene = gene,
      max_variance_ratio = round(max_ratio, 3),
      interpretation = interpretation,
      recommendation = case_when(
        max_ratio < 3 ~ "Continue with parametric tests",
        max_ratio < 10 ~ "Acceptable for large samples",
        TRUE ~ "Consider non-parametric alternatives"
      )
    )
  })

  age_variance_results <- list(
    summary = age_variance_summary,
    details = age_variance_details
  )
} else {
  cat("✓ No significant age-related homoscedasticity violations\n")
  age_variance_summary <- data.frame(
    gene = character(0),
    max_variance_ratio = numeric(0),
    interpretation = character(0),
    recommendation = character(0)
  )
  age_variance_details <- data.frame()
  age_variance_results <- list(
    summary = age_variance_summary,
    details = age_variance_details
  )
}

# ===============================================================================
# SUMMARY RECOMMENDATIONS
# ===============================================================================

# ===============================================================================
# COMBINED VARIANCE RATIO SUMMARY
# ===============================================================================

# Combine all variance ratio results
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
) %>%
  arrange(desc(max_variance_ratio))

# Create interpretation summary
variance_interpretation_summary <- if (nrow(all_variance_summaries) > 0) {
  all_variance_summaries %>%
    count(interpretation, name = "count") %>%
    mutate(percentage = round(100 * count / sum(count), 1))
} else {
  data.frame(
    interpretation = character(0),
    count = numeric(0),
    percentage = numeric(0)
  )
}

cat("\n=== VARIANCE RATIO INVESTIGATION COMPLETE ===\n")
cat(
  "Total violations found:",
  length(c(gender_violations, age_violations)),
  "\n"
)

if (nrow(all_variance_summaries) > 0) {
  cat("Breakdown by interpretation:\n")
  for (i in 1:nrow(variance_interpretation_summary)) {
    cat(
      "- ",
      variance_interpretation_summary$interpretation[i],
      ": ",
      variance_interpretation_summary$count[i],
      " genes (",
      variance_interpretation_summary$percentage[i],
      "%)\n",
      sep = ""
    )
  }
} else {
  cat(
    "✓ No homoscedasticity violations detected. Parametric tests appropriate.\n"
  )
}

cat("\nData frames available in workspace:\n")
cat("- gender_variance_summary: Summary of gender variance ratios\n")
cat("- gender_variance_details: Detailed gender variance analysis\n")
cat("- age_variance_summary: Summary of age group variance ratios\n")
cat("- age_variance_details: Detailed age group variance analysis\n")
cat("- all_variance_summaries: Combined summary table\n")
cat("- variance_interpretation_summary: Count by interpretation category\n")
cat("\nUse View() to examine tables: View(all_variance_summaries)\n")
