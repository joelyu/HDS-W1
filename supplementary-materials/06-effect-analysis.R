#!/usr/bin/env Rscript

# 06-effect-analysis-clean.R
# Logistic Regression Analysis: Gene Expression × Demographics Predicting Therapy Response
#
# This script:
# 1. Sources combined data with therapy response from 02-data-processing.R
# 2. Fits logistic regression models with interaction terms
# 3. Tests whether demographic factors modify gene expression effects on therapy response
# 4. Provides odds ratios and clinical interpretation
# 5. Identifies genes with significant demographic interactions

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
  if (!require(ResourceSelection, quietly = TRUE)) {
    install.packages("ResourceSelection", repos = "http://cran.us.r-project.org")
  }
  if (!require(performance, quietly = TRUE)) {
    install.packages("performance", repos = "http://cran.us.r-project.org")
  }

  library(tidyverse)
  library(here)
  library(broom)
  library(ResourceSelection)
  library(performance)
}))

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

# Function to fit interaction models for a demographic variable
fit_interaction_models <- function(data, gene_vars, demographic_var, demographic_name, silent = TRUE) {
  results <- list()

  for (gene in gene_vars) {
    tryCatch({
      # Filter data for this gene
      gene_data <- data %>%
        filter(!is.na(.data[[gene]]), !is.na(.data[[demographic_var]]))

      if (nrow(gene_data) < 10) next

      # Fit nested models
      null_model <- glm(response_binary ~ 1,
                       data = gene_data, family = binomial)

      gene_model <- glm(reformulate(gene, "response_binary"),
                       data = gene_data, family = binomial)

      main_effects_model <- glm(reformulate(c(gene, demographic_var), "response_binary"),
                               data = gene_data, family = binomial)

      interaction_model <- glm(reformulate(c(gene, demographic_var,
                                           paste(gene, demographic_var, sep = ":")),
                                         "response_binary"),
                              data = gene_data, family = binomial)

      # Likelihood ratio tests
      lrt_gene_vs_null <- anova(null_model, gene_model, test = "Chisq")
      lrt_main_vs_gene <- anova(gene_model, main_effects_model, test = "Chisq")
      lrt_interaction_vs_main <- anova(main_effects_model, interaction_model, test = "Chisq")

      # Extract p-values
      gene_p <- lrt_gene_vs_null$`Pr(>Chi)`[2]
      demographic_p <- lrt_main_vs_gene$`Pr(>Chi)`[2]
      interaction_p <- lrt_interaction_vs_main$`Pr(>Chi)`[2]

      # Store results
      results[[gene]] <- list(
        null_model = null_model,
        gene_model = gene_model,
        main_effects_model = main_effects_model,
        interaction_model = interaction_model,
        gene_p = gene_p,
        demographic_p = demographic_p,
        interaction_p = interaction_p,
        n_patients = nrow(gene_data)
      )

      if (!silent) {
        cat("Analyzing:", gene, "×", demographic_name, "interaction\n")
        cat("  ✓ Gene p =", format(gene_p, digits = 2),
            "| Demo p =", format(demographic_p, digits = 4),
            "| Interaction p =", format(interaction_p, digits = 4), "\n")
      }

    }, error = function(e) {
      if (!silent) cat("Error analyzing", gene, ":", e$message, "\n")
    })
  }

  return(results)
}

# Function to extract model fit statistics
extract_model_fit <- function(model, model_name, gene, demographic) {
  if (is.null(model)) return(NULL)

  data.frame(
    gene = gene,
    demographic = demographic,
    model_type = model_name,
    deviance = round(model$deviance, 3),
    df_residual = model$df.residual,
    chi_sq_p = round(pchisq(model$deviance, model$df.residual, lower.tail = FALSE), 4),
    aic = round(model$aic, 1),
    fit_quality = ifelse(
      pchisq(model$deviance, model$df.residual, lower.tail = FALSE) >= 0.05,
      "Good", "Poor"
    ),
    stringsAsFactors = FALSE
  )
}

# Function to extract Hosmer-Lemeshow test results
extract_hl_test <- function(model, model_name, gene, demographic) {
  if (is.null(model)) return(NULL)

  tryCatch({
    observed_outcomes <- model$y
    fitted_probs <- fitted(model)
    hl_test <- hoslem.test(observed_outcomes, fitted_probs, g = 10)

    data.frame(
      gene = gene,
      demographic = demographic,
      model_type = model_name,
      hl_statistic = round(hl_test$statistic, 3),
      hl_df = hl_test$parameter,
      hl_p_value = round(hl_test$p.value, 4),
      hl_fit_quality = ifelse(hl_test$p.value >= 0.05, "Good", "Poor"),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    return(NULL)
  })
}

# Function to extract AIC scores
extract_aic_scores <- function(gene, models, demographic) {
  data.frame(
    gene = gene,
    demographic = demographic,
    null_aic = round(models$null_model$aic, 2),
    gene_aic = round(models$gene_model$aic, 2),
    main_effects_aic = round(models$main_effects_model$aic, 2),
    interaction_aic = round(models$interaction_model$aic, 2),
    stringsAsFactors = FALSE
  )
}

# ===============================================================================
# DATA LOADING AND PREPARATION
# ===============================================================================

suppressMessages({
  source(here("supplementary-materials", "02-data-processing.R"))
})

# Verify data loaded successfully
if (!exists("combined_data") || nrow(combined_data) == 0) {
  stop("Failed to load combined_data from 02-data-processing.R")
}

if (!"complete_response" %in% names(combined_data)) {
  stop("complete_response column not found. Please check 02-data-processing.R")
}

# Prepare data for logistic regression
analysis_data <- combined_data %>%
  filter(!is.na(complete_response)) %>%
  mutate(response_binary = ifelse(complete_response == "Complete response", 1, 0))

# Identify gene expression variables
non_gene_cols <- c(
  "sample_id", "gender", "response_therapy", "complete_response",
  "age_at_diagnosis", "os_years", "censored", "age_group_nccn", "response_binary"
)
gene_expression_vars <- setdiff(names(analysis_data), non_gene_cols)

# Verify binary coding
response_counts <- table(analysis_data$response_binary, analysis_data$complete_response)

# ===============================================================================
# MAIN ANALYSIS: GENE × DEMOGRAPHIC INTERACTIONS
# ===============================================================================

# Analyze gene × gender interactions
gender_interaction_results <- fit_interaction_models(
  analysis_data, gene_expression_vars, "gender", "Gender", silent = TRUE
)

# Analyze gene × age group interactions
age_interaction_results <- fit_interaction_models(
  analysis_data, gene_expression_vars, "age_group_nccn", "Age", silent = TRUE
)

# ===============================================================================
# MULTIPLE TESTING CORRECTION AND SIGNIFICANCE IDENTIFICATION
# ===============================================================================

# Gender interactions: FDR correction
gender_interaction_p_values <- sapply(gender_interaction_results, function(x) x$interaction_p)
gender_interaction_p_values <- gender_interaction_p_values[!is.na(gender_interaction_p_values)]
gender_fdr_adj <- p.adjust(gender_interaction_p_values, method = "fdr")
sig_gender_genes <- names(gender_fdr_adj)[gender_fdr_adj < 0.05]

# Age interactions: FDR correction
age_interaction_p_values <- sapply(age_interaction_results, function(x) x$interaction_p)
age_interaction_p_values <- age_interaction_p_values[!is.na(age_interaction_p_values)]
age_fdr_adj <- p.adjust(age_interaction_p_values, method = "fdr")
sig_age_genes <- names(age_fdr_adj)[age_fdr_adj < 0.05]

# ===============================================================================
# COMPREHENSIVE LRT P-VALUES TABLE WITH FDR CORRECTION
# ===============================================================================

# Create comprehensive table with all LRT p-values for age analysis
age_lrt_comprehensive <- data.frame(
  Gene = names(age_interaction_results),
  # Raw LRT p-values
  Gene_vs_Null_P = sapply(age_interaction_results, function(x) round(x$gene_p, 6)),
  Age_vs_Gene_P = sapply(age_interaction_results, function(x) round(x$demographic_p, 6)),
  Interaction_vs_Main_P = sapply(age_interaction_results, function(x) round(x$interaction_p, 6)),
  N_Patients = sapply(age_interaction_results, function(x) x$n_patients),
  stringsAsFactors = FALSE
)

# Add FDR corrections for each test
age_lrt_comprehensive <- age_lrt_comprehensive %>%
  mutate(
    # FDR-adjusted p-values
    Gene_vs_Null_FDR = round(p.adjust(Gene_vs_Null_P, method = "fdr"), 6),
    Age_vs_Gene_FDR = round(p.adjust(Age_vs_Gene_P, method = "fdr"), 6),
    Interaction_vs_Main_FDR = round(p.adjust(Interaction_vs_Main_P, method = "fdr"), 6),
    # Significance indicators
    Gene_Significant = Gene_vs_Null_FDR < 0.05,
    Age_Significant = Age_vs_Gene_FDR < 0.05,
    Interaction_Significant = Interaction_vs_Main_FDR < 0.05
  ) %>%
  arrange(Interaction_vs_Main_P)

# ===============================================================================
# CREATE SUMMARY TABLES
# ===============================================================================

# Gender interaction table
gender_interaction_table <- data.frame(
  Gene = names(gender_interaction_results),
  Gene_P_Value = sapply(gender_interaction_results, function(x) round(x$gene_p, 4)),
  Gender_P_Value = sapply(gender_interaction_results, function(x) round(x$demographic_p, 4)),
  Interaction_P_Value = sapply(gender_interaction_results, function(x) round(x$interaction_p, 4)),
  Interaction_FDR = round(gender_fdr_adj[names(gender_interaction_results)], 4),
  Significant = names(gender_interaction_results) %in% sig_gender_genes,
  N_Patients = sapply(gender_interaction_results, function(x) x$n_patients),
  stringsAsFactors = FALSE
)

# Age interaction table
age_interaction_table <- data.frame(
  Gene = names(age_interaction_results),
  Gene_P_Value = sapply(age_interaction_results, function(x) round(x$gene_p, 4)),
  Age_P_Value = sapply(age_interaction_results, function(x) round(x$demographic_p, 4)),
  Interaction_P_Value = sapply(age_interaction_results, function(x) round(x$interaction_p, 4)),
  Interaction_FDR = round(age_fdr_adj[names(age_interaction_results)], 4),
  Significant = names(age_interaction_results) %in% sig_age_genes,
  N_Patients = sapply(age_interaction_results, function(x) x$n_patients),
  stringsAsFactors = FALSE
)

# Combined effects summary
all_gene_effects <- data.frame(
  Gene = gene_expression_vars,
  Gender_Gene_P = gender_interaction_table$Gene_P_Value[match(gene_expression_vars, gender_interaction_table$Gene)],
  Gender_Interaction_P = gender_interaction_table$Interaction_P_Value[match(gene_expression_vars, gender_interaction_table$Gene)],
  Gender_Interaction_FDR = gender_interaction_table$Interaction_FDR[match(gene_expression_vars, gender_interaction_table$Gene)],
  Age_Gene_P = age_interaction_table$Gene_P_Value[match(gene_expression_vars, age_interaction_table$Gene)],
  Age_Interaction_P = age_interaction_table$Interaction_P_Value[match(gene_expression_vars, age_interaction_table$Gene)],
  Age_Interaction_FDR = age_interaction_table$Interaction_FDR[match(gene_expression_vars, age_interaction_table$Gene)],
  Significant_Gender_Interaction = gene_expression_vars %in% sig_gender_genes,
  Significant_Age_Interaction = gene_expression_vars %in% sig_age_genes,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Significant_Age_Interaction), desc(Significant_Gender_Interaction), Age_Gene_P)

# ===============================================================================
# DETAILED ANALYSIS FOR SIGNIFICANT INTERACTIONS
# ===============================================================================

# Extract detailed results for significant genes
gender_detailed_results <- list()
if (length(sig_gender_genes) > 0) {
  for (gene in sig_gender_genes) {
    model_summary <- tidy(gender_interaction_results[[gene]]$interaction_model, conf.int = TRUE)
    model_summary$OR <- exp(model_summary$estimate)
    model_summary$OR_lower <- exp(model_summary$conf.low)
    model_summary$OR_upper <- exp(model_summary$conf.high)
    gender_detailed_results[[gene]] <- model_summary
  }
}

age_detailed_results <- list()
if (length(sig_age_genes) > 0) {
  for (gene in sig_age_genes) {
    model_summary <- tidy(age_interaction_results[[gene]]$interaction_model, conf.int = TRUE)
    model_summary$OR <- exp(model_summary$estimate)
    model_summary$OR_lower <- exp(model_summary$conf.low)
    model_summary$OR_upper <- exp(model_summary$conf.high)
    age_detailed_results[[gene]] <- model_summary
  }
}

# ===============================================================================
# MODEL FIT ASSESSMENT
# ===============================================================================

# Extract model fit statistics for all models
model_fit_results <- list()

for (gene in names(gender_interaction_results)) {
  gender_models <- gender_interaction_results[[gene]]
  model_fit_results <- append(model_fit_results, list(
    extract_model_fit(gender_models$null_model, "Null", gene, "Gender"),
    extract_model_fit(gender_models$gene_model, "Gene", gene, "Gender"),
    extract_model_fit(gender_models$main_effects_model, "Main Effects", gene, "Gender"),
    extract_model_fit(gender_models$interaction_model, "Interaction", gene, "Gender")
  ))
}

for (gene in names(age_interaction_results)) {
  age_models <- age_interaction_results[[gene]]
  model_fit_results <- append(model_fit_results, list(
    extract_model_fit(age_models$null_model, "Null", gene, "Age"),
    extract_model_fit(age_models$gene_model, "Gene", gene, "Age"),
    extract_model_fit(age_models$main_effects_model, "Main Effects", gene, "Age"),
    extract_model_fit(age_models$interaction_model, "Interaction", gene, "Age")
  ))
}

model_fit_table <- do.call(rbind, model_fit_results[!sapply(model_fit_results, is.null)])

# ===============================================================================
# HOSMER-LEMESHOW GOODNESS-OF-FIT TESTS
# ===============================================================================

# Extract Hosmer-Lemeshow test results
hl_results <- list()

for (gene in names(gender_interaction_results)) {
  gender_models <- gender_interaction_results[[gene]]
  hl_results <- append(hl_results, list(
    extract_hl_test(gender_models$null_model, "Null", gene, "Gender"),
    extract_hl_test(gender_models$gene_model, "Gene", gene, "Gender"),
    extract_hl_test(gender_models$main_effects_model, "Main Effects", gene, "Gender"),
    extract_hl_test(gender_models$interaction_model, "Interaction", gene, "Gender")
  ))
}

for (gene in names(age_interaction_results)) {
  age_models <- age_interaction_results[[gene]]
  hl_results <- append(hl_results, list(
    extract_hl_test(age_models$null_model, "Null", gene, "Age"),
    extract_hl_test(age_models$gene_model, "Gene", gene, "Age"),
    extract_hl_test(age_models$main_effects_model, "Main Effects", gene, "Age"),
    extract_hl_test(age_models$interaction_model, "Interaction", gene, "Age")
  ))
}

hosmer_lemeshow_table <- do.call(rbind, hl_results[!sapply(hl_results, is.null)])

# ===============================================================================
# AIC COMPARISON FOR BOTH DEMOGRAPHICS
# ===============================================================================

# Gender AIC comparison
gender_aic_results <- list()
for (gene in names(gender_interaction_results)) {
  gender_aic_results[[gene]] <- extract_aic_scores(gene, gender_interaction_results[[gene]], "Gender")
}
gender_aic_table <- do.call(rbind, gender_aic_results)

# Age AIC comparison
age_aic_results <- list()
for (gene in names(age_interaction_results)) {
  age_aic_results[[gene]] <- extract_aic_scores(gene, age_interaction_results[[gene]], "Age")
}
age_aic_table <- do.call(rbind, age_aic_results)

# Add delta AIC calculations for both tables
calculate_delta_aic <- function(aic_table) {
  aic_table %>%
    mutate(
      best_aic = pmin(null_aic, gene_aic, main_effects_aic, interaction_aic),
      delta_null = round(null_aic - best_aic, 2),
      delta_gene = round(gene_aic - best_aic, 2),
      delta_main_effects = round(main_effects_aic - best_aic, 2),
      delta_interaction = round(interaction_aic - best_aic, 2),
      best_model = case_when(
        null_aic == best_aic ~ "Null",
        gene_aic == best_aic ~ "Gene",
        main_effects_aic == best_aic ~ "Main Effects",
        interaction_aic == best_aic ~ "Interaction"
      ),
      evidence_strength = case_when(
        pmin(delta_null, delta_gene, delta_main_effects, delta_interaction) == 0 &
        pmax(delta_null, delta_gene, delta_main_effects, delta_interaction) <= 2 ~ "Weak evidence",
        pmax(delta_null, delta_gene, delta_main_effects, delta_interaction) <= 4 ~ "Moderate evidence",
        TRUE ~ "Strong evidence"
      )
    )
}

gender_aic_table <- calculate_delta_aic(gender_aic_table)
age_aic_table <- calculate_delta_aic(age_aic_table)

# ===============================================================================
# COOK'S DISTANCE ANALYSIS FOR BEST MODELS
# ===============================================================================

# Identify best models for Cook's distance analysis
myc_gene_model <- age_interaction_results[["MYC"]]$gene_model
ccnd2_interaction_model <- age_interaction_results[["CCND2"]]$interaction_model

# Cook's distance analysis for MYC
myc_cooks <- cooks.distance(myc_gene_model)
myc_data <- analysis_data %>%
  filter(!is.na(MYC), !is.na(age_group_nccn)) %>%
  mutate(
    observation = row_number(),
    cooks_distance = myc_cooks,
    influential = cooks_distance > 4 / nrow(.)
  )

# Cook's distance analysis for CCND2
ccnd2_cooks <- cooks.distance(ccnd2_interaction_model)
ccnd2_data <- analysis_data %>%
  filter(!is.na(CCND2), !is.na(age_group_nccn)) %>%
  mutate(
    observation = row_number(),
    cooks_distance = ccnd2_cooks,
    influential = cooks_distance > 4 / nrow(.)
  )

# Create Cook's distance plots
myc_cooks_plot <- ggplot(myc_data, aes(x = observation, y = cooks_distance)) +
  geom_point(aes(color = influential), alpha = 0.7) +
  geom_hline(yintercept = 4 / nrow(myc_data), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"), name = "Influential") +
  labs(
    title = "Cook's Distance: MYC Gene Model (Age Analysis)",
    subtitle = paste0("Threshold: ", round(4 / nrow(myc_data), 4),
                     " | Influential points: ", sum(myc_data$influential)),
    x = "Observation", y = "Cook's Distance"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ccnd2_cooks_plot <- ggplot(ccnd2_data, aes(x = observation, y = cooks_distance)) +
  geom_point(aes(color = influential), alpha = 0.7) +
  geom_hline(yintercept = 4 / nrow(ccnd2_data), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"), name = "Influential") +
  labs(
    title = "Cook's Distance: CCND2 Interaction Model (Age Analysis)",
    subtitle = paste0("Threshold: ", round(4 / nrow(ccnd2_data), 4),
                     " | Influential points: ", sum(ccnd2_data$influential)),
    x = "Observation", y = "Cook's Distance"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# ===============================================================================
# PERFORMANCE PACKAGE DIAGNOSTICS
# ===============================================================================

# Initialize performance package diagnostics
myc_outliers <- NULL
ccnd2_outliers <- NULL
myc_collinearity <- NULL
ccnd2_collinearity <- NULL

tryCatch({
  myc_outliers <- check_outliers(myc_gene_model, threshold = list('cook' = 0.1))
  ccnd2_outliers <- check_outliers(ccnd2_interaction_model, threshold = list('cook' = 0.1))
  myc_collinearity <- check_collinearity(myc_gene_model)
  ccnd2_collinearity <- check_collinearity(ccnd2_interaction_model)
}, error = function(e) {
  # Performance package diagnostics failed silently
})

# ===============================================================================
# RESULTS SUMMARY AND CONSOLE OUTPUT
# ===============================================================================

cat("=== GENE EXPRESSION × DEMOGRAPHIC INTERACTION ANALYSIS COMPLETE ===\n\n")

cat("DATA SUMMARY:\n")
cat("- Loaded dataset:", nrow(combined_data), "patients\n")
cat("- Analysis dataset:", nrow(analysis_data), "patients with complete response data\n")
cat("- Gene expression variables analyzed:", length(gene_expression_vars), "\n")
cat("- Complete response rate:", round(100 * sum(analysis_data$response_binary == 1) / nrow(analysis_data), 1), "%\n\n")

cat("RESPONSE VARIABLE VERIFICATION:\n")
print(response_counts)
cat("✓ 1 = Complete response,", sum(analysis_data$response_binary == 1), "patients\n")
cat("✓ 0 = Incomplete response,", sum(analysis_data$response_binary == 0), "patients\n\n")

cat("INTERACTION ANALYSIS RESULTS:\n")
cat("- Gender interactions analyzed:", length(gender_interaction_results), "genes\n")
cat("- Significant gender interactions (FDR < 0.05):", length(sig_gender_genes), "genes\n")
if (length(sig_gender_genes) > 0) {
  cat("  Significant genes:", paste(sig_gender_genes, collapse = ", "), "\n")
}

cat("- Age interactions analyzed:", length(age_interaction_results), "genes\n")
cat("- Significant age interactions (FDR < 0.05):", length(sig_age_genes), "genes\n")
if (length(sig_age_genes) > 0) {
  cat("  Significant genes:", paste(sig_age_genes, collapse = ", "), "\n")
}
cat("\n")

cat("MODEL FIT ASSESSMENT:\n")
cat("- Total models fitted:", nrow(model_fit_table), "\n")
cat("- Models with poor fit (chi-squared p < 0.05):", sum(model_fit_table$chi_sq_p < 0.05, na.rm = TRUE), "\n")
cat("- Models with good fit (chi-squared p >= 0.05):", sum(model_fit_table$chi_sq_p >= 0.05, na.rm = TRUE), "\n")

if (!is.null(hosmer_lemeshow_table) && nrow(hosmer_lemeshow_table) > 0) {
  cat("- Hosmer-Lemeshow poor fit (p < 0.05):", sum(hosmer_lemeshow_table$hl_p_value < 0.05, na.rm = TRUE), "\n")
  cat("- Hosmer-Lemeshow good fit (p >= 0.05):", sum(hosmer_lemeshow_table$hl_p_value >= 0.05, na.rm = TRUE), "\n")
}
cat("\n")

cat("COOK'S DISTANCE ANALYSIS:\n")
cat("- MYC model influential observations:", sum(myc_data$influential), "out of", nrow(myc_data), "\n")
cat("- CCND2 model influential observations:", sum(ccnd2_data$influential), "out of", nrow(ccnd2_data), "\n\n")

cat("PERFORMANCE PACKAGE DIAGNOSTICS:\n")
if (!is.null(myc_outliers)) {
  cat("- MYC model outliers: Performance package analysis completed\n")
} else {
  cat("- MYC model outliers: Analysis failed or no outliers detected\n")
}

if (!is.null(myc_collinearity)) {
  cat("- MYC model collinearity: No multicollinearity detected\n")
} else {
  cat("- MYC model collinearity: No issues or analysis failed\n")
}

if (!is.null(ccnd2_collinearity)) {
  cat("- CCND2 model collinearity: Multicollinearity detected (see ccnd2_collinearity)\n")
} else {
  cat("- CCND2 model collinearity: Analysis failed\n")
}
cat("\n")

cat("=== DATA OBJECTS AVAILABLE IN R ENVIRONMENT ===\n\n")

cat("MAIN RESULTS:\n")
cat("- gender_interaction_table: Gender × gene interaction results with FDR correction\n")
cat("- age_interaction_table: Age × gene interaction results with FDR correction\n")
cat("- age_lrt_comprehensive: Comprehensive LRT p-values (raw + FDR) for all age models\n")
cat("- all_gene_effects: Combined summary of all gene effects and interactions\n")
cat("- sig_gender_genes: Vector of genes with significant gender interactions\n")
cat("- sig_age_genes: Vector of genes with significant age interactions\n\n")

cat("DETAILED MODEL RESULTS:\n")
cat("- gender_interaction_results: List of all fitted models for gender analysis\n")
cat("- age_interaction_results: List of all fitted models for age analysis\n")
cat("- gender_detailed_results: Detailed coefficients for significant gender interactions\n")
cat("- age_detailed_results: Detailed coefficients for significant age interactions\n\n")

cat("MODEL DIAGNOSTICS:\n")
cat("- model_fit_table: Chi-squared goodness-of-fit tests for all models\n")
if (!is.null(hosmer_lemeshow_table)) {
  cat("- hosmer_lemeshow_table: Hosmer-Lemeshow goodness-of-fit tests\n")
}
cat("- gender_aic_table: AIC comparison and model selection for gender analysis\n")
cat("- age_aic_table: AIC comparison and model selection for age analysis\n\n")

cat("OUTLIER ANALYSIS:\n")
cat("- myc_data: MYC model data with Cook's distances\n")
cat("- ccnd2_data: CCND2 model data with Cook's distances\n")
cat("- myc_cooks_plot: Cook's distance plot for MYC model\n")
cat("- ccnd2_cooks_plot: Cook's distance plot for CCND2 model\n")
if (!is.null(myc_outliers)) cat("- myc_outliers: Performance package outlier detection for MYC\n")
if (!is.null(ccnd2_outliers)) cat("- ccnd2_outliers: Performance package outlier detection for CCND2\n\n")

cat("MULTICOLLINEARITY ANALYSIS:\n")
if (!is.null(myc_collinearity)) {
  cat("- myc_collinearity: VIF analysis for MYC model (NULL = no issues)\n")
} else {
  cat("- myc_collinearity: VIF analysis for MYC model (NULL)\n")
}
if (!is.null(ccnd2_collinearity)) {
  cat("- ccnd2_collinearity: VIF analysis for CCND2 model\n")
} else {
  cat("- ccnd2_collinearity: VIF analysis for CCND2 model (NULL)\n")
}
cat("\n")

cat("BEST MODEL OBJECTS:\n")
cat("- myc_gene_model: Best performing model for MYC (gene main effect)\n")
cat("- ccnd2_interaction_model: Best performing model for CCND2 (interaction model)\n\n")

cat("=== ANALYSIS COMPLETE ===\n")
cat("Use the data objects listed above for further analysis, plotting, or reporting.\n")
cat("Key findings are summarized in gender_interaction_table and age_interaction_table.\n")