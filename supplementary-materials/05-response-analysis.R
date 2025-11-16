#!/usr/bin/env Rscript

# 05-response-analysis.R
# Chi-squared Analysis of Therapy Response vs Demographics
#
# This script:
# 1. Sources combined data from 02-data-processing.R
# 2. Performs Chi-squared tests for therapy response vs gender and age groups
# 3. Checks expected frequency assumptions
# 4. Creates contingency tables with counts and percentages
# 5. Provides summary tables of associations

# Load required libraries
suppressMessages(suppressWarnings({
  if (!require(tidyverse, quietly = TRUE)) {
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
  }
  if (!require(here, quietly = TRUE)) {
    install.packages("here", repos = "http://cran.us.r-project.org")
  }
  if (!require(broom, quietly = TRUE)) {
    install.packages("broom", repos = "http://cran.us.r-project.org")
  }

  library(tidyverse)
  library(here)
  library(broom)
}))

# ===============================================================================
# DATA LOADING FROM PROCESSING SCRIPT
# ===============================================================================

cat("Loading combined data with therapy response variables...\n")
source(here("supplementary-materials", "02-data-processing.R"))

# Verify data loaded successfully with new complete_response column
if (!exists("combined_data") || nrow(combined_data) == 0) {
  stop("Failed to load combined_data from 02-data-processing.R")
}

if (!"complete_response" %in% names(combined_data)) {
  stop("complete_response column not found. Please check 02-data-processing.R")
}

cat("✓ Loaded combined dataset:", nrow(combined_data), "patients\n")
cat(
  "✓ Found complete_response column with levels:",
  levels(combined_data$complete_response),
  "\n\n"
)

# Flag to control output display (set to TRUE to show analysis outputs)
show_outputs <- FALSE

# ===============================================================================
# DATA EXPLORATION: THERAPY RESPONSE DISTRIBUTION
# ===============================================================================

# Overall therapy response distribution
overall_response <- combined_data %>%
  filter(!is.na(complete_response)) %>%
  count(complete_response) %>%
  mutate(
    percentage = round(100 * n / sum(n), 1),
    label = paste0(complete_response, ": ", n, " (", percentage, "%)")
  )

# Store total sample size for analysis
n_total_response <- sum(overall_response$n)

# ===============================================================================
# CHI-SQUARED ANALYSIS: THERAPY RESPONSE VS GENDER
# ===============================================================================

# Create contingency table for therapy response vs gender
gender_response_table <- combined_data %>%
  filter(!is.na(complete_response), !is.na(gender)) %>%
  count(gender, complete_response) %>%
  pivot_wider(names_from = complete_response, values_from = n, values_fill = 0)

# Convert to matrix for chi-squared test
gender_matrix <- as.matrix(gender_response_table[, -1])
rownames(gender_matrix) <- gender_response_table$gender

# Check expected frequencies
gender_expected <- chisq.test(gender_matrix)$expected

# Check if expected frequencies meet assumptions (all cells >= 5)
min_expected_gender <- min(gender_expected)

# Perform Chi-squared test
gender_chisq_result <- chisq.test(gender_matrix)

# Calculate Cramér's V for effect size
cramers_v_gender <- sqrt(
  gender_chisq_result$statistic /
    (sum(gender_matrix) * (min(dim(gender_matrix)) - 1))
)

# Interpret effect size
effect_interpretation_gender <- case_when(
  cramers_v_gender < 0.1 ~ "Negligible",
  cramers_v_gender < 0.3 ~ "Small",
  cramers_v_gender < 0.5 ~ "Medium",
  TRUE ~ "Large"
)

# Create percentage table for clinical interpretation
gender_percentage_table <- combined_data %>%
  filter(!is.na(complete_response), !is.na(gender)) %>%
  count(gender, complete_response) %>%
  group_by(gender) %>%
  mutate(
    total = sum(n),
    percentage = round(100 * n / total, 1)
  ) %>%
  ungroup() %>%
  select(gender, complete_response, n, percentage) %>%
  pivot_wider(
    names_from = complete_response,
    values_from = c(n, percentage),
    names_sep = "_"
  )

# ===============================================================================
# CHI-SQUARED ANALYSIS: THERAPY RESPONSE VS AGE GROUP
# ===============================================================================

# Create contingency table for therapy response vs age group
age_response_table <- combined_data %>%
  filter(!is.na(complete_response), !is.na(age_group_nccn)) %>%
  count(age_group_nccn, complete_response) %>%
  pivot_wider(names_from = complete_response, values_from = n, values_fill = 0)

# Convert to matrix for chi-squared test
age_matrix <- as.matrix(age_response_table[, -1])
rownames(age_matrix) <- age_response_table$age_group_nccn

# Check expected frequencies
age_expected <- chisq.test(age_matrix)$expected

# Check if expected frequencies meet assumptions
min_expected_age <- min(age_expected)

# Perform Chi-squared test
age_chisq_result <- chisq.test(age_matrix)

# Calculate Cramér's V for effect size
cramers_v_age <- sqrt(
  age_chisq_result$statistic /
    (sum(age_matrix) * (min(dim(age_matrix)) - 1))
)

# Interpret effect size
effect_interpretation_age <- case_when(
  cramers_v_age < 0.1 ~ "Negligible",
  cramers_v_age < 0.3 ~ "Small",
  cramers_v_age < 0.5 ~ "Medium",
  TRUE ~ "Large"
)

# Create percentage table for clinical interpretation
age_percentage_table <- combined_data %>%
  filter(!is.na(complete_response), !is.na(age_group_nccn)) %>%
  count(age_group_nccn, complete_response) %>%
  group_by(age_group_nccn) %>%
  mutate(
    total = sum(n),
    percentage = round(100 * n / total, 1)
  ) %>%
  ungroup() %>%
  select(age_group_nccn, complete_response, n, percentage) %>%
  pivot_wider(
    names_from = complete_response,
    values_from = c(n, percentage),
    names_sep = "_"
  )

# ===============================================================================
# VISUALIZATION: THERAPY RESPONSE BY DEMOGRAPHICS
# ===============================================================================

# Gender vs therapy response visualization
gender_response_plot <- combined_data %>%
  filter(!is.na(complete_response), !is.na(gender)) %>%
  ggplot(aes(x = gender, fill = complete_response)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "Complete response" = "#2E86AB",
      "Incomplete response" = "#E63946"
    ),
    name = "Therapy Response"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Therapy Response Distribution by Gender",
    subtitle = paste0(
      "Chi-squared test: p = ",
      round(gender_chisq_result$p.value, 4)
    ),
    x = "Gender",
    y = "Proportion of Patients"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Age group vs therapy response visualization
age_response_plot <- combined_data %>%
  filter(!is.na(complete_response), !is.na(age_group_nccn)) %>%
  ggplot(aes(x = age_group_nccn, fill = complete_response)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "Complete response" = "#2E86AB",
      "Incomplete response" = "#E63946"
    ),
    name = "Therapy Response"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Therapy Response Distribution by Age Group",
    subtitle = paste0(
      "Chi-squared test: p = ",
      round(age_chisq_result$p.value, 4)
    ),
    x = "Age Group (NCCN-IPI Classification)",
    y = "Proportion of Patients"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display plots (controlled by show_outputs flag)
if (show_outputs) {
  print(gender_response_plot)
  print(age_response_plot)
}

# ===============================================================================
# RESULTS SUMMARY TABLES
# ===============================================================================

# Create comprehensive results summary
chi_squared_summary <- data.frame(
  demographic = c("Gender", "Age Group"),
  chi_squared = round(
    c(gender_chisq_result$statistic, age_chisq_result$statistic),
    3
  ),
  df = c(gender_chisq_result$parameter, age_chisq_result$parameter),
  p_value = round(c(gender_chisq_result$p.value, age_chisq_result$p.value), 4),
  significant = c(
    gender_chisq_result$p.value < 0.05,
    age_chisq_result$p.value < 0.05
  ),
  cramers_v = round(c(cramers_v_gender, cramers_v_age), 3),
  effect_size = c(effect_interpretation_gender, effect_interpretation_age),
  min_expected_freq = round(c(min_expected_gender, min_expected_age), 2),
  assumptions_met = c(min_expected_gender >= 5, min_expected_age >= 5)
)

# Store results for potential further analysis
therapy_response_results <- list(
  gender = list(
    contingency_table = gender_response_table,
    percentage_table = gender_percentage_table,
    chi_squared = gender_chisq_result,
    expected_frequencies = gender_expected,
    cramers_v = cramers_v_gender,
    plot = gender_response_plot
  ),
  age_group = list(
    contingency_table = age_response_table,
    percentage_table = age_percentage_table,
    chi_squared = age_chisq_result,
    expected_frequencies = age_expected,
    cramers_v = cramers_v_age,
    plot = age_response_plot
  ),
  summary = chi_squared_summary
)

# ===============================================================================
# CLINICAL INTERPRETATION (Prepare for conditional output)
# ===============================================================================

# Gender interpretation logic (without output)
gender_interpretation <- if (gender_chisq_result$p.value < 0.05) {
  # Find which gender has higher complete response rate
  gender_rates <- gender_percentage_table %>%
    select(gender, ends_with("Complete response")) %>%
    rename(complete_rate = `percentage_Complete response`)

  best_gender <- gender_rates$gender[which.max(gender_rates$complete_rate)]
  best_rate <- max(gender_rates$complete_rate)

  list(
    significant = TRUE,
    message = paste0(
      "GENDER: Statistically significant association with therapy response\n- Higher complete response rate in ",
      best_gender,
      " patients (",
      best_rate,
      "%)"
    )
  )
} else {
  list(
    significant = FALSE,
    message = "GENDER: No statistically significant association with therapy response\n- Treatment response appears similar between males and females"
  )
}

# Age group interpretation logic (without output)
age_interpretation <- if (age_chisq_result$p.value < 0.05) {
  # Find which age group has highest complete response rate
  age_rates <- age_percentage_table %>%
    select(age_group_nccn, ends_with("Complete response")) %>%
    rename(complete_rate = `percentage_Complete response`)

  best_age_group <- age_rates$age_group_nccn[which.max(age_rates$complete_rate)]
  best_age_rate <- max(age_rates$complete_rate)
  worst_age_group <- age_rates$age_group_nccn[which.min(
    age_rates$complete_rate
  )]
  worst_age_rate <- min(age_rates$complete_rate)

  list(
    significant = TRUE,
    message = paste0(
      "AGE GROUP: Statistically significant association with therapy response\n",
      "- Highest complete response rate: ",
      best_age_group,
      " (",
      best_age_rate,
      "%)\n",
      "- Lowest complete response rate: ",
      worst_age_group,
      " (",
      worst_age_rate,
      "%)"
    )
  )
} else {
  list(
    significant = FALSE,
    message = "AGE GROUP: No statistically significant association with therapy response\n- Treatment response appears similar across age groups"
  )
}

# ===============================================================================
# CONDITIONAL OUTPUT SECTION
# ===============================================================================

# Check for errors or warnings that need user attention
has_errors <- FALSE
error_messages <- character(0)

# Check for assumption violations
if (min_expected_gender < 5) {
  has_errors <- TRUE
  error_messages <- c(
    error_messages,
    paste0(
      "WARNING: Gender analysis has expected frequencies < 5 (minimum: ",
      round(min_expected_gender, 2),
      "). Consider Fisher's exact test."
    )
  )
}

if (min_expected_age < 5) {
  has_errors <- TRUE
  error_messages <- c(
    error_messages,
    paste0(
      "WARNING: Age group analysis has expected frequencies < 5 (minimum: ",
      round(min_expected_age, 2),
      "). Consider Fisher's exact test or combine categories."
    )
  )
}

# Display output conditionally
if (has_errors || show_outputs) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("THERAPY RESPONSE ANALYSIS RESULTS\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  if (has_errors) {
    cat("\n=== WARNINGS/ERRORS ===\n")
    for (msg in error_messages) {
      cat(msg, "\n")
    }
  }

  if (show_outputs) {
    cat("\n=== THERAPY RESPONSE DISTRIBUTION ===\n")
    cat("Overall therapy response distribution:\n")
    for (i in 1:nrow(overall_response)) {
      cat("-", overall_response$label[i], "\n")
    }
    cat("\nTotal patients with therapy response data:", n_total_response, "\n")

    cat("\n=== THERAPY RESPONSE vs GENDER ===\n")
    cat("Contingency table (observed frequencies):\n")
    print(gender_response_table)
    cat("\nExpected frequencies:\n")
    print(round(gender_expected, 2))
    cat("\nExpected frequency check:\n")
    cat("- Minimum expected frequency:", round(min_expected_gender, 2), "\n")
    if (min_expected_gender >= 5) {
      cat("- ✓ All expected frequencies ≥ 5: Chi-squared test appropriate\n")
    } else {
      cat("- ⚠ Some expected frequencies < 5: Consider Fisher's exact test\n")
    }
    cat("\nChi-squared test results:\n")
    cat(
      "- Chi-squared statistic:",
      round(gender_chisq_result$statistic, 3),
      "\n"
    )
    cat("- Degrees of freedom:", gender_chisq_result$parameter, "\n")
    cat("- p-value:", round(gender_chisq_result$p.value, 4), "\n")
    cat("- Cramér's V (effect size):", round(cramers_v_gender, 3), "\n")
    cat("- Effect size interpretation:", effect_interpretation_gender, "\n")
    cat("\nResponse rates by gender:\n")
    print(gender_percentage_table)

    cat("\n=== THERAPY RESPONSE vs AGE GROUP ===\n")
    cat("Contingency table (observed frequencies):\n")
    print(age_response_table)
    cat("\nExpected frequencies:\n")
    print(round(age_expected, 2))
    cat("\nExpected frequency check:\n")
    cat("- Minimum expected frequency:", round(min_expected_age, 2), "\n")
    if (min_expected_age >= 5) {
      cat("- ✓ All expected frequencies ≥ 5: Chi-squared test appropriate\n")
    } else {
      cat(
        "- ⚠ Some expected frequencies < 5: Consider Fisher's exact test or combine categories\n"
      )
    }
    cat("\nChi-squared test results:\n")
    cat("- Chi-squared statistic:", round(age_chisq_result$statistic, 3), "\n")
    cat("- Degrees of freedom:", age_chisq_result$parameter, "\n")
    cat("- p-value:", round(age_chisq_result$p.value, 4), "\n")
    cat("- Cramér's V (effect size):", round(cramers_v_age, 3), "\n")
    cat("- Effect size interpretation:", effect_interpretation_age, "\n")
    cat("\nResponse rates by age group:\n")
    print(age_percentage_table)

    cat("\n=== SUMMARY RESULTS ===\n")
    print(chi_squared_summary)

    cat("\n=== CLINICAL INTERPRETATION ===\n")
    cat(gender_interpretation$message, "\n\n")
    cat(age_interpretation$message, "\n")

    cat("\n=== VISUALIZATIONS ===\n")
    print(gender_response_plot)
    print(age_response_plot)
  }

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("DATA OBJECTS AVAILABLE FOR FURTHER ANALYSIS:\n")
  cat("- therapy_response_results: Complete analysis results including plots\n")
  cat("- chi_squared_summary: Summary table of all Chi-squared tests\n")
  cat("- gender_response_table: Gender contingency table\n")
  cat("- age_response_table: Age group contingency table\n")
  cat("- gender_percentage_table: Gender response rates\n")
  cat("- age_percentage_table: Age group response rates\n")
  cat("- overall_response: Overall therapy response distribution\n")
  cat("\nUse View() to examine tables: View(chi_squared_summary)\n")
  cat("Set show_outputs <- TRUE to display full analysis results\n")
  cat("Therapy response analysis complete!\n")
}
