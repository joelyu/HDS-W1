#!/usr/bin/env Rscript

# Data pre-processing script for Reddy's DLBCL Dataset
# Converts Excel supplementary tables to CSV format for analysis
#
# Source: Reddy et al. (2017) Cell, Volume 171, Issue 2
# DOI: 10.1016/j.cell.2017.09.027

# Load required libraries
if (!require(readxl, quietly = TRUE)) {
  install.packages("readxl", repos = "http://cran.us.r-project.org")
  library(readxl)
}

if (!require(here, quietly = TRUE)) {
  install.packages("here", repos = "http://cran.us.r-project.org")
  library(here)
}

# Create data directory if it doesn't exist
data_dir <- here("data")
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
  cat("Created data directory:", data_dir, "\n")
} else {
  cat("Found existing data directory:", data_dir, "\n")
}

# Define source URLs and temporary file paths
mmc1_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc1.xlsx"
mmc2_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc2.xlsx"

mmc1_local <- tempfile(fileext = ".xlsx")
mmc2_local <- tempfile(fileext = ".xlsx")

# Helper function to convert sheet name to PascalCase
to_pascal_case <- function(text) {
  # Split by spaces and convert each word to title case
  words <- strsplit(text, " ")[[1]]
  words <- tools::toTitleCase(tolower(words))
  # Remove spaces and concatenate
  paste0(words, collapse = "")
}

# Function to process Excel supplementary table to CSV
process_supplementary_table <- function(
  url,
  sheet_name,
  expected_rows,
  expected_cols,
  description
) {
  tryCatch(
    {
      # Download to temporary file
      temp_file <- tempfile(fileext = ".xlsx")
      cat("Downloading", description, "...\n")
      download.file(url, temp_file, mode = "wb", quiet = FALSE)

      # Verify download
      if (!file.exists(temp_file) || file.size(temp_file) == 0) {
        stop("Download failed or file is empty")
      }
      cat("Downloaded", description, "successfully", "\n")

      # Read Excel file, skip first 3 rows
      cat("Reading Excel sheet:", sheet_name, "\n")
      data <- read_excel(temp_file, sheet = sheet_name, skip = 3)

      # Check dimensions
      actual_rows <- nrow(data)
      actual_cols <- ncol(data)

      cat("Loaded data:", actual_rows, "rows x", actual_cols, "columns\n")

      # Verify expected dimensions
      if (actual_rows != expected_rows || actual_cols != expected_cols) {
        warning(paste0(
          "Dimension failed, expected: ",
          expected_rows,
          "x",
          expected_cols,
          ", Downloaded: ",
          actual_rows,
          "x",
          actual_cols
        ))
      } else {
        cat("Dimensions verified for", description, "\n")
      }

      # Generate output filename from URL and sheet name
      # Extract mmc number from URL (e.g., "mmc1", "mmc2")
      mmc_match <- regmatches(url, regexpr("mmc\\d+", url))
      sheet_pascal <- to_pascal_case(sheet_name)
      output_filename <- paste0(mmc_match, "-", sheet_pascal, ".csv")

      # Export to CSV
      output_path <- here("data", output_filename)
      write.csv(data, output_path, row.names = FALSE)

      cat("Exported to:", output_path, "\n")

      return(TRUE)
    },
    error = function(e) {
      cat("Error processing", description, ":", e$message, "\n")
      cat(
        "Please download manually from:",
        url,
        "and follow the pre-processing steps in the README",
        "\n\n"
      )
      return(FALSE)
    }
  )
}

# Process both tables using the unified function
cat("Downloading tables...\n\n")

# Process MMC1 - Clinical Information (1001 rows x 35 columns expected)
mmc1_success <- process_supplementary_table(
  url = mmc1_url,
  sheet_name = "Clinical Information",
  expected_rows = 1001,
  expected_cols = 35,
  description = "MMC1 (Clinical Information)"
)

# Process MMC2 - Gene Expression (775 rows x 21 columns expected)
mmc2_success <- process_supplementary_table(
  url = mmc2_url,
  sheet_name = "Gene Expression",
  expected_rows = 775,
  expected_cols = 21,
  description = "MMC2 (Gene Expression)"
)

# Temporary files are automatically cleaned up by the OS

# Summary
cat("\nDownload and pre-processing complete!\n")
cat("-------------------------------------\n")
cat("File are ready in", data_dir, ":\n")

output_files <- list.files(data_dir, pattern = "*.csv", full.names = TRUE)
for (file in output_files) {
  if (file.exists(file)) {
    size_mb <- round(file.size(file) / 1024^2, 2)
    cat("  -", basename(file), "(", size_mb, "MB )\n")
  }
}
