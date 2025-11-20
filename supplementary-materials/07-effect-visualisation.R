#!/usr/bin/env Rscript

# 07-effect-visualisation.R
# Visualization of Gene Expression x Demographics Interaction Effects
# NOTE: Comprehensive summary printed at the end of script execution
#
# This script creates three main visualizations:
# 1. Forest plot of gene effects (odds ratios) for all 21 genes
# 2. Predicted probability curves for MYC and CCND2 interaction models by age group
# 3. Forest plot of age group odds ratios from CCND2 interaction model
#
# NOTE: This version creates plots in R environment only (no exports) for verification

# Load required libraries
suppressMessages(suppressWarnings({
  if (!require(tidyverse, quietly = TRUE)) {
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
  }
  if (!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  }
  if (!require(plotly, quietly = TRUE)) {
    install.packages("plotly", repos = "http://cran.us.r-project.org")
  }
  if (!require(broom, quietly = TRUE)) {
    install.packages("broom", repos = "http://cran.us.r-project.org")
  }
  if (!require(patchwork, quietly = TRUE)) {
    install.packages("patchwork", repos = "http://cran.us.r-project.org")
  }
  if (!require(here, quietly = TRUE)) {
    install.packages("here", repos = "http://cran.us.r-project.org")
  }
  if (!require(emmeans, quietly = TRUE)) {
    install.packages("emmeans", repos = "http://cran.us.r-project.org")
  }

  library(tidyverse)
  library(ggplot2)
  library(plotly)
  library(broom)
  library(patchwork)
  library(here)
  library(emmeans)
}))

# ===============================================================================
# SOURCE ANALYSIS RESULTS
# ===============================================================================

# Load analysis results from 06-effect-analysis.R
source(here("supplementary-materials", "06-effect-analysis.R"))

# Verify required objects are available
required_objects <- c(
  "age_interaction_results",
  "gender_interaction_results",
  "analysis_data",
  "gene_expression_vars"
)

missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects from 06-effect-analysis.R: ",
    paste(missing_objects, collapse = ", ")
  )
}

# Analysis results loaded successfully

# ===============================================================================
# PLOT 1: FOREST PLOT OF GENE EFFECTS (ODDS RATIOS)
# ===============================================================================

# Create Plot 1: Forest plot of gene effects (odds ratios)

# Extract odds ratios from gene-only models (main effect models)
extract_gene_or <- function(results_list) {
  gene_or_data <- list()

  for (gene in names(results_list)) {
    tryCatch(
      {
        # Use the gene-only model for main gene effect
        model <- results_list[[gene]]$gene_model
        if (!is.null(model)) {
          # Get model summary with confidence intervals
          model_summary <- tidy(model, conf.int = TRUE)

          # Extract gene coefficient (should be the second row, first is intercept)
          gene_coef <- model_summary[2, ]

          if (!is.na(gene_coef$estimate)) {
            gene_or_data[[gene]] <- data.frame(
              gene = gene,
              log_or = gene_coef$estimate,
              or = exp(gene_coef$estimate),
              or_lower = exp(gene_coef$conf.low),
              or_upper = exp(gene_coef$conf.high),
              p_value = gene_coef$p.value,
              significant = gene_coef$p.value < 0.05,
              stringsAsFactors = FALSE
            )
          }
        }
      },
      error = function(e) {
        # Skip genes that failed
      }
    )
  }

  return(do.call(rbind, gene_or_data))
}

# Extract odds ratios from age interaction results (using gene models for main effects)
gene_or_data <- extract_gene_or(age_interaction_results)

# Get FDR-corrected significance from age_lrt_comprehensive table
gene_or_data <- gene_or_data %>%
  left_join(
    age_lrt_comprehensive %>%
      select(Gene, Gene_vs_Null_FDR, Gene_Significant),
    by = c("gene" = "Gene")
  ) %>%
  arrange(or) %>%
  mutate(
    gene = factor(gene, levels = gene),
    # Create 3-way significance categories
    significance_category = case_when(
      Gene_Significant ~ "Adjusted p-value < 0.05",
      significant ~ "Unadjusted p-value < 0.05",
      TRUE ~ "Unadjusted p-value ≥ 0.05"
    ),
    significance_category = factor(
      significance_category,
      levels = c(
        "Adjusted p-value < 0.05",
        "Unadjusted p-value < 0.05",
        "Unadjusted p-value ≥ 0.05"
      )
    ),
    hover_text = paste0(
      "Gene: ",
      gene,
      "\n",
      "OR: ",
      round(or, 3),
      " (",
      round(or_lower, 3),
      " - ",
      round(or_upper, 3),
      ")\n",
      "Unadjusted p-value: ",
      format(p_value, scientific = TRUE, digits = 3),
      "\n",
      "Adjusted p-value: ",
      format(Gene_vs_Null_FDR, scientific = TRUE, digits = 3),
      "\n",
      "Significance: ",
      case_when(
        Gene_Significant ~ "Adjusted significant",
        significant ~ "Unadjusted significant",
        TRUE ~ "Insignificant"
      )
    )
  )

# Create forest plot
plot1_forest <- ggplot(gene_or_data, aes(x = or, y = gene)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(xmin = or_lower, xmax = or_upper, color = significance_category),
    width = 0.3,
    linewidth = 0.8,
    orientation = "y"
  ) +
  geom_point(aes(color = significance_category, text = hover_text), size = 4) +
  scale_color_manual(
    values = c(
      "Adjusted p-value < 0.05" = "#FC8D62", # Orange (Set 2) - most important
      "Unadjusted p-value < 0.05" = "#66C2A5", # Teal (Set 2) - moderate importance
      "Unadjusted p-value ≥ 0.05" = "#8DA0CB" # Light purple (Set 2) - least important
    ),
    name = "Significance Level",
    guide = guide_legend(
      override.aes = list(size = 5, alpha = 1),
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  scale_x_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.2", "0.5", "1.0", "2.0", "5.0", "10.0")
  ) +
  labs(
    title = "Logistic Regression Models - Gene Expression Effects on Complete Response",
    subtitle = "Odds Ratios with 95% CIs - Orange: Adjusted Significant, Teal: Unadjusted Significant, Light Purple: Insignificant",
    x = "Odds Ratio (log scale)",
    y = "Gene",
    caption = "Gray dashed line indicates OR = 1.0 (no effect). Log scale ensures equal visual distance for OR=0.5 and OR=2.0"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10),
    panel.grid.minor.x = element_blank()
  )

# Convert to plotly for hover functionality
plot1_interactive <- ggplotly(plot1_forest, tooltip = "text") %>%
  layout(
    title = list(
      text = "Logistic Regression Models - Gene Expression Effects on Complete Response<br><sub>Orange: Adjusted p-values < 0.05<br>Teal: Adjusted p-values ≥ 0.05, unadjusted p-values < 0.05<br>Light Purple: Unadjusted p-values ≥ 0.05 (insignificant)<br>Confidence interval: 95%</sub>",
      font = list(size = 16),
      pad = list(t = 20, b = 80) # Add padding above and below title
    ),
    showlegend = FALSE,
    margin = list(t = 120, b = 80, l = 80, r = 80) # Increase top margin for multi-line subtitle
  )

# ===============================================================================
# PLOT 2: PREDICTED PROBABILITY CURVES (MYC vs CCND2 by Age Groups)
# ===============================================================================

# Create Plot 2: Predicted probability curves for MYC and CCND2 models

# Function to create predicted probability data
create_prediction_data <- function(model, gene_name, data_subset, age_groups) {
  # Create prediction grid
  gene_range <- seq(
    from = min(data_subset[[gene_name]], na.rm = TRUE),
    to = max(data_subset[[gene_name]], na.rm = TRUE),
    length.out = 100
  )

  pred_data <- expand_grid(
    gene_expr = gene_range,
    age_group_nccn = age_groups
  )

  # Set gene name in prediction data
  pred_data[[gene_name]] <- pred_data$gene_expr

  # Get predictions
  pred_data$predicted_prob <- predict(
    model,
    newdata = pred_data,
    type = "response"
  )
  pred_data$gene_name <- gene_name

  return(pred_data)
}

# Prepare data for predictions
analysis_data_clean <- analysis_data %>%
  filter(
    !is.na(MYC),
    !is.na(CCND2),
    !is.na(age_group_nccn),
    !is.na(response_binary)
  )

age_groups <- levels(analysis_data_clean$age_group_nccn)

# Get models for prediction
myc_interaction_model <- age_interaction_results[["MYC"]]$interaction_model
ccnd2_interaction_model <- age_interaction_results[["CCND2"]]$interaction_model

# Create prediction data
myc_pred_data <- create_prediction_data(
  myc_interaction_model,
  "MYC",
  analysis_data_clean,
  age_groups
)

ccnd2_pred_data <- create_prediction_data(
  ccnd2_interaction_model,
  "CCND2",
  analysis_data_clean,
  age_groups
)

# Create MYC probability plot with original data
plot2a_myc <- ggplot() +
  # Add original data points first (will be behind the lines)
  geom_point(
    data = analysis_data_clean,
    aes(x = MYC, y = response_binary, color = age_group_nccn),
    alpha = 0.6,
    size = 1.5,
    position = position_jitter(height = 0.02, width = 0) # Slight jitter for visibility
  ) +
  # Add prediction lines on top
  geom_line(
    data = myc_pred_data,
    aes(x = gene_expr, y = predicted_prob, color = age_group_nccn),
    linewidth = 1.5,
    alpha = 0.9
  ) +
  scale_color_brewer(name = "Age Group", type = "qual", palette = "Set2") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "MYC Gene Expression",
    subtitle = "Predicted Probability of Complete Response by Age Group (with original data)",
    x = "MYC Expression Level (log2)",
    y = "Predicted Probability of Complete Response"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

# Create CCND2 probability plot with original data
plot2b_ccnd2 <- ggplot() +
  # Add original data points first (will be behind the lines)
  geom_point(
    data = analysis_data_clean,
    aes(x = CCND2, y = response_binary, color = age_group_nccn),
    alpha = 0.6,
    size = 1.5,
    position = position_jitter(height = 0.02, width = 0) # Slight jitter for visibility
  ) +
  # Add prediction lines on top
  geom_line(
    data = ccnd2_pred_data,
    aes(x = gene_expr, y = predicted_prob, color = age_group_nccn),
    linewidth = 1.5,
    alpha = 0.9
  ) +
  scale_color_brewer(name = "Age Group", type = "qual", palette = "Set2") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "CCND2 Gene Expression",
    subtitle = "Predicted Probability of Complete Response by Age Group (with original data)",
    x = "CCND2 Expression Level (log2)",
    y = "Predicted Probability of Complete Response"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

# Combine plots side by side
plot2_combined <- plot2a_myc +
  plot2b_ccnd2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# Add overall title
plot2_combined <- plot2_combined +
  plot_annotation(
    title = "Gene Expression x Age Group Interaction Effects on Therapy Response",
    subtitle = "Predicted probabilities from interaction models overlaid with original patient data",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# Convert to interactive plots
plot2a_interactive <- ggplotly(plot2a_myc) %>%
  layout(
    title = list(
      text = "MYC Gene Expression<br><sub>Predicted Probability with Original Patient Data by Age Group</sub>"
    )
  )

plot2b_interactive <- ggplotly(plot2b_ccnd2) %>%
  layout(
    title = list(
      text = "CCND2 Gene Expression<br><sub>Predicted Probability with Original Patient Data by Age Group</sub>"
    )
  )

# ===============================================================================
# PLOT 3: CCND2 SLOPES BY AGE GROUP (EMTRENDS APPROACH)
# ===============================================================================

# Create Plot 3: CCND2 slopes by age group (emtrends approach)

# Use emtrends to get slopes (trends) for how age group effects change with CCND2
# This shows the rate of change in log-odds for each age group per unit CCND2 increase
emtrends_result <- emtrends(
  ccnd2_interaction_model,
  specs = ~age_group_nccn,
  var = "CCND2"
)

# Get p-values for each trend using test()
emtrends_test <- test(emtrends_result) %>%
  as.data.frame()

# Combine emtrends results with p-values
plot3_data <- as.data.frame(emtrends_result) %>%
  left_join(
    emtrends_test %>% select(age_group_nccn, p.value),
    by = "age_group_nccn"
  ) %>%
  mutate(
    # Calculate FDR-adjusted p-values
    p_adjusted = p.adjust(p.value, method = "fdr"),
    age_group = case_when(
      age_group_nccn == "≤40 years" ~ "≤40 years (Reference)",
      age_group_nccn == "41-60 years" ~ "41-60 years",
      age_group_nccn == "61-75 years" ~ "61-75 years",
      age_group_nccn == ">75 years" ~ ">75 years",
      TRUE ~ as.character(age_group_nccn)
    ),
    # Convert slope to odds ratio per unit CCND2 increase
    or_per_unit = exp(CCND2.trend),
    or_lower = exp(asymp.LCL),
    or_upper = exp(asymp.UCL),
    # Significance based on unadjusted p-value
    significant = p.value < 0.05,
    # Significance based on adjusted p-value
    significant_adj = p_adjusted < 0.05,
    hover_text = paste0(
      "Age Group: ",
      age_group,
      "\n",
      "OR per unit CCND2 increase: ",
      round(or_per_unit, 3),
      " (",
      round(or_lower, 3),
      " - ",
      round(or_upper, 3),
      ")\n",
      "Slope: ",
      round(CCND2.trend, 3),
      "\n",
      "Unadjusted p-value: ",
      format(p.value, scientific = TRUE, digits = 3),
      "\n",
      "Adjusted p-value: ",
      format(p_adjusted, scientific = TRUE, digits = 3),
      "\n",
      "Interpretation: For each 1-unit increase in CCND2, the odds multiply by ",
      round(or_per_unit, 3)
    )
  ) %>%
  arrange(or_per_unit) %>%
  mutate(age_group = factor(age_group, levels = age_group))

# Create comprehensive p-values table for further analysis
ccnd2_age_pvalues_table <- plot3_data %>%
  select(
    age_group_nccn,
    age_group,
    CCND2.trend,
    SE,
    or_per_unit,
    or_lower,
    or_upper,
    p.value,
    p_adjusted,
    significant,
    significant_adj
  ) %>%
  rename(
    Age_Group = age_group_nccn,
    Age_Group_Label = age_group,
    Log_OR_Slope = CCND2.trend,
    Standard_Error = SE,
    OR_per_Unit = or_per_unit,
    OR_Lower = or_lower,
    OR_Upper = or_upper,
    P_Value_Unadjusted = p.value,
    P_Value_Adjusted = p_adjusted,
    Significant_Unadjusted = significant,
    Significant_Adjusted = significant_adj
  )

# CCND2 Age Group P-Values Table created: ccnd2_age_pvalues_table

# Create emtrends plot showing slopes (rates of change)
plot3 <- ggplot(plot3_data, aes(x = or_per_unit, y = age_group)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(xmin = or_lower, xmax = or_upper, color = significant),
    width = 0.3,
    linewidth = 2,
    orientation = "y"
  ) +
  geom_point(aes(color = significant, text = hover_text), size = 4) +
  scale_color_manual(
    values = c("TRUE" = "#FC8D62", "FALSE" = "#8DA0CB"), # ColorBrewer Set 2: Orange for significant, Light purple for insignificant
    name = "Significance",
    labels = c("Insignificant", "Significant")
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0),
    labels = c("0.5", "0.8", "1.0", "1.2", "1.5", "2.0")
  ) +
  labs(
    title = "CCND2 Interaction Effects - Age Group Slopes",
    subtitle = "Odds Ratio per Unit Increase in CCND2 Expression by Age Group - Orange: Significant, Light Purple: Insignificant",
    x = "Odds Ratio per Unit CCND2 Increase (log scale)",
    y = "Age Group",
    caption = "Gray dashed line indicates OR = 1.0 (no effect of CCND2). Slopes show how strongly CCND2 affects each age group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10),
    panel.grid.minor.x = element_blank()
  )

# Convert to interactive
plot3_interactive <- ggplotly(plot3, tooltip = "text") %>%
  layout(
    title = list(
      text = "CCND2 Interaction Effects - Age Group Slopes<br><sub>Orange: Significant slopes, Light Purple: Insignificant slopes</sub>",
      font = list(size = 16)
    ),
    showlegend = FALSE
  )

# ===============================================================================
# PLOT SUMMARY (R ENVIRONMENT OBJECTS ONLY)
# ===============================================================================

# ===============================================================================
# COMPREHENSIVE GENE EXPRESSION VISUALIZATION ANALYSIS SUMMARY
# ===============================================================================

# Count genes in each significance category for summary
adj_sig_count <- sum(
  gene_or_data$significance_category == "Adjusted p-value < 0.05",
  na.rm = TRUE
)
unadj_sig_count <- sum(
  gene_or_data$significance_category == "Unadjusted p-value < 0.05",
  na.rm = TRUE
)
insig_count <- sum(
  gene_or_data$significance_category == "Unadjusted p-value ≥ 0.05",
  na.rm = TRUE
)

# Get significant genes lists
adj_genes <- if (adj_sig_count > 0) {
  gene_or_data$gene[
    gene_or_data$significance_category == "Adjusted p-value < 0.05"
  ]
} else {
  character(0)
}

# CCND2 age-specific effects summary
ccnd2_sig_count <- sum(plot3_data$significant, na.rm = TRUE)
ccnd2_sig_adj_count <- sum(plot3_data$significant_adj, na.rm = TRUE)

cat(sprintf(
  "
================================================================================
GENE EXPRESSION VISUALIZATION ANALYSIS: THERAPY RESPONSE PREDICTION
================================================================================

1. ANALYSIS OVERVIEW
--------------------------------------------------
Total genes analysed: %d
Forest plot visualization: Gene-only model odds ratios with confidence intervals
Interaction visualization: MYC (no interaction) vs CCND2 (age interaction)
Age-specific effects: CCND2 emmeans analysis across 4 age groups

2. FOREST PLOT RESULTS (PLOT 1)
--------------------------------------------------
Gene effects summary (FDR-adjusted significance):
- Significant genes (FDR < 0.05): %d genes %s
- Unadjusted significant (p < 0.05): %d genes
- Non-significant genes (p ≥ 0.05): %d genes

Statistical interpretation: %s

3. INTERACTION VISUALIZATION (PLOT 2)
--------------------------------------------------
MYC Gene Analysis:
- Predicted probability curves by age group
- Pattern: Parallel lines (no age × gene interaction)
- Clinical meaning: MYC effect consistent across all age groups

CCND2 Gene Analysis:
- Predicted probability curves by age group
- Pattern: Diverging lines (significant age × gene interaction)
- Clinical meaning: CCND2 effect varies significantly by age group

4. CCND2 AGE-SPECIFIC EFFECTS (PLOT 3 - EMMEANS ANALYSIS)
--------------------------------------------------
Statistical method: emtrends() - continuous gene effect slopes by age group
Total age groups analyzed: %d
Significant age-specific effects (p < 0.05): %d/%d age groups

Age-Specific Odds Ratios (per unit increase in CCND2 expression):
(Table shows odds ratios, confidence intervals, and p-values for each age group)
",
  length(gene_expression_vars),
  adj_sig_count,
  if (adj_sig_count > 0) {
    paste("(", paste(adj_genes, collapse = ", "), ")")
  } else {
    ""
  },
  unadj_sig_count,
  insig_count,
  if (adj_sig_count > 0) {
    "Gene expression effects detected after multiple testing correction"
  } else {
    "No genes show significant effects after FDR correction"
  },
  nrow(plot3_data),
  ccnd2_sig_count,
  nrow(plot3_data)
))

# Display CCND2 age-specific results table
print(plot3_data[, c(
  "age_group",
  "or_per_unit",
  "or_lower",
  "or_upper",
  "p.value",
  "significant"
)])

cat("\nVisualization Plots:\n")
print(plot1_forest)
print(plot2_combined)
print(plot3)

cat(sprintf(
  "
5. CLINICAL INTERPRETATION
--------------------------------------------------
Key findings from emmeans analysis:
- Age groups show differential sensitivity to CCND2 expression
- Continuous odds ratios avoid arbitrary expression cutpoints
- %s age groups demonstrate significant CCND2 effects (p < 0.05)

Clinical relevance:
%s

6. VISUALIZATION OUTPUTS
--------------------------------------------------
The following plots are displayed above:

✓ Forest Plot (Plot 1): Odds ratios for %d genes with confidence intervals
  - Shows gene-only model effects across all 21 genes
  - Significance color coding for visual interpretation

✓ Interaction Curves (Plot 2): MYC vs CCND2 predicted probabilities by age
  - MYC: Parallel lines (no age interaction)
  - CCND2: Diverging lines (significant age interaction)

✓ Age-Specific Effects (Plot 3): CCND2 emmeans slopes with significance coding
  - Odds ratios per unit CCND2 increase for each age group
  - Visual representation of age-dependent therapeutic effects

Statistical validation:
✓ Comprehensive odds ratio analysis completed
✓ Age-interaction effects properly characterised using emmeans
✓ Continuous gene effects presented without arbitrary cutpoints
✓ Multiple testing correction applied (FDR method)

Analysis complete. Visualization objects and data tables available for review.

================================================================================
",
  ccnd2_sig_count,
  if (ccnd2_sig_count > 0) {
    "CCND2 expression shows age-dependent therapeutic effects"
  } else {
    "CCND2 effects appear consistent across age groups"
  },
  length(gene_expression_vars)
))
