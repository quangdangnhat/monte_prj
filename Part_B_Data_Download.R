# =============================================================================
# PART B: DATA DOWNLOAD SCRIPT
# =============================================================================
# This script downloads real data for QTE analysis:
# 1. Economic Policy Uncertainty (EPU) Index from policyuncertainty.com
# 2. S&P 500 data from Yahoo Finance
# =============================================================================

rm(list = ls())

# Install and load required packages
packages <- c("quantmod", "zoo", "xts", "dplyr", "readr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("=============================================================================\n")
cat("DATA DOWNLOAD FOR QTE ANALYSIS\n")
cat("=============================================================================\n\n")

# =============================================================================
# 1. DOWNLOAD S&P 500 DATA
# =============================================================================
cat(">>> Downloading S&P 500 data from Yahoo Finance...\n")

tryCatch({
  # Download S&P 500 (^GSPC)
  getSymbols("^GSPC", src = "yahoo", from = "2000-01-01", to = Sys.Date())

  # Calculate monthly returns
  sp500_monthly <- to.monthly(GSPC)
  sp500_returns <- monthlyReturn(Cl(sp500_monthly), type = "log") * 100  # Percentage

  # Convert to data frame
  sp500_df <- data.frame(
    date = index(sp500_returns),
    returns = as.numeric(sp500_returns)
  )
  colnames(sp500_df) <- c("date", "returns")

  cat(sprintf("   Downloaded %d months of S&P 500 data\n", nrow(sp500_df)))
  cat(sprintf("   Date range: %s to %s\n",
              min(sp500_df$date), max(sp500_df$date)))

  # Save to CSV
  write.csv(sp500_df, "sp500_data.csv", row.names = FALSE)
  cat("   Saved to 'sp500_data.csv'\n\n")

}, error = function(e) {
  cat("   Error downloading S&P 500 data. Creating simulated data instead.\n")
  cat("   Error message:", e$message, "\n\n")

  # Create simulated S&P 500 returns
  n <- 300  # 25 years of monthly data
  dates <- seq(as.Date("2000-01-01"), by = "month", length.out = n)
  returns <- rnorm(n, mean = 0.5, sd = 4)  # Typical stock market returns

  sp500_df <- data.frame(date = dates, returns = returns)
  write.csv(sp500_df, "sp500_data.csv", row.names = FALSE)
  cat("   Simulated data saved to 'sp500_data.csv'\n\n")
})

# =============================================================================
# 2. EPU DATA INSTRUCTIONS
# =============================================================================
cat(">>> Economic Policy Uncertainty (EPU) Data\n")
cat("   \n")
cat("   EPU data must be downloaded manually from:\n")
cat("   https://www.policyuncertainty.com/us_monthly.html\n")
cat("   \n")
cat("   Steps:\n")
cat("   1. Go to the website above\n")
cat("   2. Download the US Monthly EPU Index\n")
cat("   3. Save as 'epu_raw.csv' in the working directory\n")
cat("   4. Run the processing code below\n\n")

# Check if EPU file exists
if (file.exists("epu_raw.csv")) {
  cat(">>> Processing EPU data...\n")

  epu_raw <- read.csv("epu_raw.csv")
  print(head(epu_raw))

  # Process EPU data (adjust column names based on actual file)
  # Typical format: Year, Month, EPU_Index
  if ("Year" %in% colnames(epu_raw) && "Month" %in% colnames(epu_raw)) {
    epu_df <- data.frame(
      date = as.Date(paste(epu_raw$Year, epu_raw$Month, "01", sep = "-")),
      uncertainty = epu_raw[, 3]  # Assuming EPU is in 3rd column
    )
  } else {
    # Try to detect date column
    epu_df <- epu_raw
    colnames(epu_df) <- c("date", "uncertainty")
    epu_df$date <- as.Date(epu_df$date)
  }

  write.csv(epu_df, "epu_data.csv", row.names = FALSE)
  cat("   Processed EPU data saved to 'epu_data.csv'\n\n")

} else {
  cat(">>> EPU file not found. Creating simulated uncertainty data...\n")

  # Create simulated EPU-like data
  if (exists("sp500_df")) {
    n <- nrow(sp500_df)
    dates <- sp500_df$date
  } else {
    n <- 300
    dates <- seq(as.Date("2000-01-01"), by = "month", length.out = n)
  }

  # Simulate EPU-like process (AR(1) with occasional spikes)
  set.seed(456)
  uncertainty <- numeric(n)
  uncertainty[1] <- 100
  for (t in 2:n) {
    shock <- ifelse(runif(1) < 0.05, rnorm(1, 0, 50), rnorm(1, 0, 10))
    uncertainty[t] <- 0.9 * uncertainty[t-1] + shock
    uncertainty[t] <- max(uncertainty[t], 20)  # Floor at 20
  }

  epu_df <- data.frame(date = dates, uncertainty = uncertainty)
  write.csv(epu_df, "epu_data.csv", row.names = FALSE)
  cat("   Simulated uncertainty data saved to 'epu_data.csv'\n\n")
}

# =============================================================================
# 3. MERGE DATA
# =============================================================================
cat(">>> Merging datasets...\n")

sp500_df <- read.csv("sp500_data.csv")
epu_df <- read.csv("epu_data.csv")

sp500_df$date <- as.Date(sp500_df$date)
epu_df$date <- as.Date(epu_df$date)

# Merge by date
merged_data <- merge(sp500_df, epu_df, by = "date")
merged_data <- merged_data[complete.cases(merged_data), ]

cat(sprintf("   Merged dataset: %d observations\n", nrow(merged_data)))
cat(sprintf("   Date range: %s to %s\n",
            min(merged_data$date), max(merged_data$date)))

# Calculate log changes for uncertainty (to make stationary)
merged_data$uncertainty_change <- c(NA, diff(log(merged_data$uncertainty + 1)) * 100)
merged_data <- merged_data[-1, ]  # Remove first NA row

# Save final dataset
final_data <- data.frame(
  date = merged_data$date,
  returns = merged_data$returns,
  uncertainty = merged_data$uncertainty_change
)

write.csv(final_data, "qte_analysis_data.csv", row.names = FALSE)
cat("   Final dataset saved to 'qte_analysis_data.csv'\n")

# =============================================================================
# 4. DATA SUMMARY
# =============================================================================
cat("\n")
cat("=============================================================================\n")
cat("DATA SUMMARY\n")
cat("=============================================================================\n")
cat(sprintf("Variable         | Mean      | Std Dev   | Min       | Max\n"))
cat("-----------------------------------------------------------------------------\n")
cat(sprintf("Returns          | %9.4f | %9.4f | %9.4f | %9.4f\n",
            mean(final_data$returns), sd(final_data$returns),
            min(final_data$returns), max(final_data$returns)))
cat(sprintf("Uncertainty      | %9.4f | %9.4f | %9.4f | %9.4f\n",
            mean(final_data$uncertainty), sd(final_data$uncertainty),
            min(final_data$uncertainty), max(final_data$uncertainty)))
cat("=============================================================================\n")

cat("\nData preparation complete. Now run 'Part_B_QTE.R' for analysis.\n")
