# =============================================================================
# MASTER SCRIPT: Run All Analyses
# Financial Econometrics Coursework 2025-2026
# =============================================================================
# This script runs the complete analysis pipeline:
# 1. Part A: Monte Carlo Power Comparison
# 2. Part B: Data Download and QTE Analysis
# 3. Generate Report
# =============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
cat("║     FINANCIAL ECONOMETRICS COURSEWORK 2025-2026                           ║\n")
cat("║     Complete Analysis Pipeline                                            ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Set working directory (modify if needed)
# setwd("/path/to/your/project")

# Record start time
start_time <- Sys.time()

# =============================================================================
# PART A: Monte Carlo Power Comparison
# =============================================================================
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("PART A: Monte Carlo Power Comparison\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("This will take approximately 30-60 minutes for R=1000, B=499\n")
cat("For quick testing, modify N, R, B in Part_A_Monte_Carlo.R\n\n")

run_part_a <- readline(prompt = "Run Part A simulation? (y/n): ")

if (tolower(run_part_a) == "y") {
  cat("\n>>> Starting Part A...\n\n")
  source("Part_A_Monte_Carlo.R")
  cat("\n>>> Part A completed!\n")
} else {
  cat(">>> Skipping Part A\n")
}

# =============================================================================
# PART B: QTE Analysis
# =============================================================================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("PART B: Quantile Transfer Entropy Analysis\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

run_part_b <- readline(prompt = "Run Part B analysis? (y/n): ")

if (tolower(run_part_b) == "y") {
  # First, prepare data
  cat("\n>>> Step 1: Preparing data...\n\n")
  source("Part_B_Data_Download.R")

  # Then run QTE analysis
  cat("\n>>> Step 2: Running QTE analysis...\n\n")
  source("Part_B_QTE.R")
  cat("\n>>> Part B completed!\n")
} else {
  cat(">>> Skipping Part B\n")
}

# =============================================================================
# GENERATE REPORT
# =============================================================================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("GENERATE PDF REPORT\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

generate_report <- readline(prompt = "Generate PDF report? (y/n): ")

if (tolower(generate_report) == "y") {
  if (require("rmarkdown")) {
    cat("\n>>> Rendering PDF report...\n")
    tryCatch({
      rmarkdown::render("Coursework_Report.Rmd", output_format = "pdf_document")
      cat(">>> Report generated: Coursework_Report.pdf\n")
    }, error = function(e) {
      cat(">>> Error generating PDF. Make sure you have LaTeX installed.\n")
      cat(">>> Try: tinytex::install_tinytex()\n")
      cat(">>> Error message:", e$message, "\n")
    })
  } else {
    cat(">>> Please install rmarkdown: install.packages('rmarkdown')\n")
  }
} else {
  cat(">>> Skipping report generation\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total time: %.2f minutes\n", as.numeric(duration)))
cat("\nOutput files:\n")
cat("  - Part_A_Results.csv       : Monte Carlo simulation results\n")
cat("  - Part_B_QTE_Results.csv   : QTE analysis results\n")
cat("  - Coursework_Report.pdf    : Complete report (if generated)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
