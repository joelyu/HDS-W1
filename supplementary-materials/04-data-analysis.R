#!/usr/bin/env Rscript

# 04-data-analysis.R
# Statistical Analysis for Gene Expression vs Demographics
#
# This script:
# 1. Sources exploratory analysis results from 03-data-exploration.R
# 2. Implements Welch ANOVA for age group comparisons (handles unequal variances)
# 3. Implements standard t-tests for gender comparisons (Levene test showed statistically significant unequal variance for certain genes, but variance ratio and group sizes were still appropriate for standard t-test)
# 4. Performs Games-Howell post-hoc tests for significant ANOVA results
# 5. Stores comprehensive statistical results

# Load required libraries
suppressMessages(suppressWarnings({
  if (!require(tidyverse, quietly = TRUE)) {
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
  }
  if (!require(here, quietly = TRUE)) {
    install.packages("here", repos = "http://cran.us.r-project.org")
  }
  if (!require(rstatix, quietly = TRUE)) {
    install.packages("rstatix", repos = "http://cran.us.r-project.org")
  }

  library(tidyverse)
  library(here)
  library(rstatix)
}))

# ===============================================================================
# DATA LOADING FROM EXPLORATION SCRIPT
# ===============================================================================

cat("Loading data and exploratory analysis results...\n")
source(here("supplementary-materials", "03-data-exploration.R"))

# Verify data loaded successfully
if (!exists("combined_data") || nrow(combined_data) == 0) {
  stop("Failed to load combined_data from 03-data-exploration.R")
}

if (!exists("gene_expression_vars") || length(gene_expression_vars) == 0) {
  stop("Failed to load gene_expression_vars from 03-data-exploration.R")
}

cat("✓ Loaded combined dataset:", nrow(combined_data), "patients\n")
cat("✓ Loaded", length(gene_expression_vars), "gene expression variables\n\n")

# ===============================================================================
# STATISTICAL ANALYSIS STRATEGY
# ===============================================================================

cat("=== STATISTICAL ANALYSIS STRATEGY ===\n")
cat("Based on Levene test results and variance ratio analysis:\n")
cat("- GENDER COMPARISONS: Standard t-tests (equal variances assumed)\n")
cat("- AGE GROUP COMPARISONS: Welch ANOVA (unequal variances handled)\n")
cat("- POST-HOC: Games-Howell for significant Welch ANOVA results\n\n")

# ===============================================================================
# GENDER COMPARISONS: STANDARD T-TESTS
# ===============================================================================

cat("=== PERFORMING GENDER COMPARISONS (T-TESTS) ===\n")

# Initialize results storage
gender_ttest_results <- list()
gender_ttest_summary <- list()

for (gene in gene_expression_vars) {
  cat("Analyzing:", gene, "\n")

  # Prepare data for t-test
  analysis_data <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(gender)) %>%
    select(all_of(gene), gender)

  if (nrow(analysis_data) < 10) {
    cat("  Skipping - insufficient data (n =", nrow(analysis_data), ")\n")
    next
  }

  # Perform t-test
  tryCatch(
    {
      # Two-sample t-test assuming equal variances
      ttest_result <- t.test(
        analysis_data[[gene]] ~ analysis_data$gender,
        var.equal = TRUE # Based on variance ratio analysis
      )

      # Calculate effect size (Cohen's d)
      group_stats <- analysis_data %>%
        group_by(gender) %>%
        summarise(
          n = n(),
          mean = mean(.data[[gene]], na.rm = TRUE),
          sd = sd(.data[[gene]], na.rm = TRUE),
          .groups = "drop"
        )

      pooled_sd <- sqrt(
        ((group_stats$n[1] - 1) *
          group_stats$sd[1]^2 +
          (group_stats$n[2] - 1) * group_stats$sd[2]^2) /
          (sum(group_stats$n) - 2)
      )

      cohens_d <- abs(diff(group_stats$mean)) / pooled_sd

      # Store detailed results
      gender_ttest_results[[gene]] <- list(
        gene = gene,
        n_total = nrow(analysis_data),
        group_stats = group_stats,
        ttest = ttest_result,
        cohens_d = cohens_d
      )

      # Store summary results
      gender_ttest_summary[[gene]] <- data.frame(
        gene = gene,
        n_total = nrow(analysis_data),
        n_female = group_stats$n[group_stats$gender == "F"],
        n_male = group_stats$n[group_stats$gender == "M"],
        mean_female = round(group_stats$mean[group_stats$gender == "F"], 3),
        mean_male = round(group_stats$mean[group_stats$gender == "M"], 3),
        mean_difference = round(abs(diff(group_stats$mean)), 3),
        t_statistic = round(ttest_result$statistic, 3),
        df = ttest_result$parameter,
        p_value = ttest_result$p.value,
        p_adjusted = NA, # Will be adjusted later
        significant = ttest_result$p.value < 0.05,
        cohens_d = round(cohens_d, 3),
        effect_size = case_when(
          cohens_d < 0.2 ~ "Small",
          cohens_d < 0.5 ~ "Small-Medium",
          cohens_d < 0.8 ~ "Medium",
          TRUE ~ "Large"
        )
      )

      cat(
        "  ✓ t =",
        round(ttest_result$statistic, 3),
        "| p =",
        round(ttest_result$p.value, 4),
        "| d =",
        round(cohens_d, 3),
        "\n"
      )
    },
    error = function(e) {
      cat("  Warning: Could not perform t-test for", gene, ":", e$message, "\n")
    }
  )
}

# Combine t-test summaries and adjust p-values
gender_ttest_table <- bind_rows(gender_ttest_summary) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  mutate(significant_adjusted = p_adjusted < 0.05) %>%
  arrange(p_value)

cat("Gender analysis complete:", nrow(gender_ttest_table), "genes analyzed\n\n")

# ===============================================================================
# AGE GROUP COMPARISONS: WELCH ANOVA
# ===============================================================================

cat("=== PERFORMING AGE GROUP COMPARISONS (WELCH ANOVA) ===\n")

# Initialize results storage
age_anova_results <- list()
age_anova_summary <- list()

for (gene in gene_expression_vars) {
  cat("Analyzing:", gene, "\n")

  # Prepare data for ANOVA
  analysis_data <- combined_data %>%
    filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
    select(all_of(gene), age_group_nccn)

  if (nrow(analysis_data) < 20) {
    cat("  Skipping - insufficient data (n =", nrow(analysis_data), ")\n")
    next
  }

  # Perform Welch ANOVA
  tryCatch(
    {
      # Welch ANOVA (handles unequal variances)
      anova_result <- oneway.test(
        analysis_data[[gene]] ~ analysis_data$age_group_nccn,
        var.equal = FALSE # Welch correction for unequal variances
      )

      # Calculate group statistics
      group_stats <- analysis_data %>%
        group_by(age_group_nccn) %>%
        summarise(
          n = n(),
          mean = mean(.data[[gene]], na.rm = TRUE),
          sd = sd(.data[[gene]], na.rm = TRUE),
          .groups = "drop"
        )

      # Calculate eta-squared (effect size) - approximate for Welch ANOVA
      grand_mean <- mean(analysis_data[[gene]], na.rm = TRUE)
      ss_between <- sum(group_stats$n * (group_stats$mean - grand_mean)^2)
      ss_total <- sum((analysis_data[[gene]] - grand_mean)^2, na.rm = TRUE)
      eta_squared <- ss_between / ss_total

      # Store detailed results
      age_anova_results[[gene]] <- list(
        gene = gene,
        n_total = nrow(analysis_data),
        group_stats = group_stats,
        anova = anova_result,
        eta_squared = eta_squared
      )

      # Store summary results
      age_anova_summary[[gene]] <- data.frame(
        gene = gene,
        n_total = nrow(analysis_data),
        f_statistic = round(anova_result$statistic, 3),
        df_num = anova_result$parameter[1],
        df_den = round(anova_result$parameter[2], 1),
        p_value = anova_result$p.value,
        p_adjusted = NA, # Will be adjusted later
        significant = anova_result$p.value < 0.05,
        eta_squared = round(eta_squared, 3),
        effect_size = case_when(
          eta_squared < 0.01 ~ "Small",
          eta_squared < 0.06 ~ "Small-Medium",
          eta_squared < 0.14 ~ "Medium",
          TRUE ~ "Large"
        )
      )

      cat(
        "  ✓ F =",
        round(anova_result$statistic, 3),
        "| p =",
        round(anova_result$p.value, 4),
        "| η² =",
        round(eta_squared, 3),
        "\n"
      )
    },
    error = function(e) {
      cat("  Warning: Could not perform ANOVA for", gene, ":", e$message, "\n")
    }
  )
}

# Combine ANOVA summaries and adjust p-values
age_anova_table <- bind_rows(age_anova_summary) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  mutate(significant_adjusted = p_adjusted < 0.05) %>%
  arrange(p_value)

cat("Age group analysis complete:", nrow(age_anova_table), "genes analyzed\n\n")

# ===============================================================================
# POST-HOC ANALYSIS: GAMES-HOWELL FOR SIGNIFICANT WELCH ANOVA RESULTS
# ===============================================================================

cat("=== PERFORMING POST-HOC ANALYSIS (GAMES-HOWELL) ===\n")

# Identify genes with significant ANOVA results
significant_anova_genes <- age_anova_table %>%
  filter(significant_adjusted == TRUE) %>%
  pull(gene)

if (length(significant_anova_genes) > 0) {
  cat(
    "Found",
    length(significant_anova_genes),
    "genes with significant age group differences\n"
  )

  # Initialize post-hoc results storage
  games_howell_results <- list()
  games_howell_summary <- list()

  for (gene in significant_anova_genes) {
    cat("Post-hoc analysis for:", gene, "\n")

    # Prepare data
    analysis_data <- combined_data %>%
      filter(!is.na(.data[[gene]]), !is.na(age_group_nccn)) %>%
      select(all_of(gene), age_group_nccn)

    tryCatch(
      {
        # Games-Howell test for unequal variances (appropriate for Welch ANOVA)
        games_howell_result <- analysis_data %>%
          games_howell_test(formula(paste(gene, "~ age_group_nccn")))

        # Convert to consistent format
        games_howell_tidy <- games_howell_result %>%
          select(group1, group2, estimate = estimate, p.adj, p.adj.signif)

        # Store results
        games_howell_results[[gene]] <- list(
          gene = gene,
          comparisons = games_howell_result,
          tidy_results = games_howell_tidy
        )

        # Create summary table
        games_howell_summary[[gene]] <- games_howell_tidy %>%
          mutate(
            gene = gene,
            comparison = paste(group1, "vs", group2),
            mean_difference = round(estimate, 3),
            p_value = p.adj,
            significant = p.adj < 0.05,
            .before = 1
          ) %>%
          select(
            gene,
            comparison,
            mean_difference,
            p_value,
            significant
          )

        # Display significant comparisons
        significant_comparisons <- sum(
          games_howell_tidy$p.adj < 0.05,
          na.rm = TRUE
        )
        cat(
          "  ✓",
          significant_comparisons,
          "significant pairwise comparisons found\n"
        )
      },
      error = function(e) {
        cat(
          "  Warning: Could not perform post-hoc analysis for",
          gene,
          ":",
          e$message,
          "\n"
        )
      }
    )
  }

  # Combine Games-Howell summaries
  if (length(games_howell_summary) > 0) {
    games_howell_table <- bind_rows(games_howell_summary) %>%
      arrange(gene, p_value)
  } else {
    games_howell_table <- data.frame()
  }
} else {
  cat(
    "No genes with significant age group differences - skipping post-hoc analysis\n"
  )
  games_howell_results <- list()
  games_howell_table <- data.frame()
}

# ===============================================================================
# RESULTS SUMMARY
# ===============================================================================

cat("\n=== STATISTICAL ANALYSIS COMPLETE ===\n")

# Gender analysis summary
significant_gender <- sum(gender_ttest_table$significant_adjusted, na.rm = TRUE)
cat("GENDER COMPARISONS (t-tests):\n")
cat("- Total genes analyzed:", nrow(gender_ttest_table), "\n")
cat("- Significant after FDR correction:", significant_gender, "\n")

# Age group analysis summary
significant_age <- sum(age_anova_table$significant_adjusted, na.rm = TRUE)
cat("\nAGE GROUP COMPARISONS (Welch ANOVA):\n")
cat("- Total genes analyzed:", nrow(age_anova_table), "\n")
cat("- Significant after FDR correction:", significant_age, "\n")

# Post-hoc summary
if (length(games_howell_results) > 0) {
  total_comparisons <- sum(sapply(games_howell_summary, nrow))
  significant_comparisons <- sum(games_howell_table$significant, na.rm = TRUE)
  cat("\nPOST-HOC ANALYSIS (Games-Howell):\n")
  cat("- Genes with significant ANOVA:", length(games_howell_results), "\n")
  cat("- Total pairwise comparisons:", total_comparisons, "\n")
  cat("- Significant pairwise comparisons:", significant_comparisons, "\n")
}

cat("\nDATA FRAMES AVAILABLE:\n")
cat("- gender_ttest_table: Complete t-test results for gender comparisons\n")
cat(
  "- age_anova_table: Complete Welch ANOVA results for age group comparisons\n"
)
if (nrow(games_howell_table) > 0) {
  cat(
    "- games_howell_table: Post-hoc pairwise comparisons for significant ANOVA results\n"
  )
}
cat("\nDETAILED RESULTS AVAILABLE:\n")
cat("- gender_ttest_results: Detailed t-test objects and statistics\n")
cat("- age_anova_results: Detailed ANOVA objects and statistics\n")
if (length(games_howell_results) > 0) {
  cat("- games_howell_results: Detailed post-hoc analysis objects\n")
}

cat("\nUse View() to examine results: View(gender_ttest_table)\n")
cat("Statistical analysis complete!\n")
