# =============================================================================
# FINANCIAL ECONOMETRICS COURSEWORK 2025-2026
# PART B: Quantile Transfer Entropy (QTE) Analysis
# =============================================================================
# This script analyzes the information flow from Economic Uncertainty (X)
# to Financial Returns (Y) using the QTE framework with Stationary Block Bootstrap
# =============================================================================

# Clear environment and set seed
rm(list = ls())
set.seed(123)

# Load required packages
packages <- c("quantreg", "tseries", "zoo", "xts", "readr", "dplyr", "ggplot2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("=============================================================================\n")
cat("PART B: Quantile Transfer Entropy (QTE) Analysis\n")
cat("=============================================================================\n")

# =============================================================================
# PARAMETERS
# =============================================================================
B <- 499                          # Bootstrap replications
tau_values <- c(0.10, 0.50, 0.90) # Quantiles to test
k_values <- c(1, 2, 3, 5)         # Lag values to consider
alpha <- 0.05                     # Significance level

cat(sprintf("Bootstrap Replications (B): %d\n", B))
cat(sprintf("Quantiles (tau): %s\n", paste(tau_values, collapse = ", ")))
cat(sprintf("Lags (k): %s\n", paste(k_values, collapse = ", ")))
cat("=============================================================================\n\n")

# =============================================================================
# SECTION 1: DATA LOADING AND PREPROCESSING
# =============================================================================
# Instructions for data:
# 1. Download EPU data from: https://www.policyuncertainty.com/
#    (Use US Monthly EPU Index or your preferred country)
# 2. Download S&P 500 or other financial index data
# 3. Save as CSV files in the working directory
# =============================================================================

# Function to simulate data if real data not available
# This generates data with known properties for testing
generate_simulated_data <- function(n = 2000, rho_x = 0.3, rho_y = 0.2, causality = 0.15) {
  cat(">>> Generating simulated data for demonstration...\n")
  cat("    (Replace with real EPU and financial data for actual analysis)\n\n")

  # Generate AR(1) process for uncertainty (X)
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (t in 2:n) {
    x[t] <- rho_x * x[t-1] + rnorm(1)
  }

  # Generate returns (Y) with potential causality from X
  y <- numeric(n)
  y[1] <- rnorm(1)
  for (t in 2:n) {
    # Y depends on its own lag and potentially on X's lag
    y[t] <- rho_y * y[t-1] + causality * x[t-1] + rnorm(1, sd = 1 + 0.3 * abs(x[t-1]))
  }

  # Make stationary (difference if needed)
  # Here we use percentage changes
  return(data.frame(
    date = 1:n,
    uncertainty = x,  # Already stationary (simulated)
    returns = y       # Already stationary (simulated)
  ))
}

# Try to load real data, otherwise use simulated
load_data <- function() {
  # Check for existing data files
  epu_file <- "epu_data.csv"
  sp500_file <- "sp500_data.csv"

  if (file.exists(epu_file) && file.exists(sp500_file)) {
    cat(">>> Loading real data from CSV files...\n")
    epu <- read.csv(epu_file)
    sp500 <- read.csv(sp500_file)
    # Process and merge data (adjust column names as needed)
    # This is a placeholder - adjust based on actual data format
    data <- merge(epu, sp500, by = "date")
    return(data)
  } else {
    cat(">>> Data files not found. Using simulated data.\n")
    cat("    To use real data, place 'epu_data.csv' and 'sp500_data.csv'\n")
    cat("    in the working directory.\n\n")
    return(generate_simulated_data(n = 2000))
  }
}

# Load data
data <- load_data()

# Extract series
X <- data$uncertainty  # Economic Uncertainty Index (e.g., EPU)
Y <- data$returns      # Financial Returns

n <- length(Y)
cat(sprintf("Sample size: %d observations\n", n))

# =============================================================================
# SECTION 2: STATIONARITY TESTS
# =============================================================================
cat("\n>>> Checking stationarity...\n")

# ADF test for X
adf_x <- adf.test(X, alternative = "stationary")
cat(sprintf("ADF test for Uncertainty (X): p-value = %.4f %s\n",
            adf_x$p.value, ifelse(adf_x$p.value < 0.05, "(Stationary)", "(Non-stationary)")))

# ADF test for Y
adf_y <- adf.test(Y, alternative = "stationary")
cat(sprintf("ADF test for Returns (Y): p-value = %.4f %s\n",
            adf_y$p.value, ifelse(adf_y$p.value < 0.05, "(Stationary)", "(Non-stationary)")))

# If non-stationary, take differences
if (adf_x$p.value >= 0.05) {
  cat(">>> Differencing X to achieve stationarity...\n")
  X <- diff(X)
  Y <- Y[-1]  # Align lengths
  n <- length(Y)
}

if (adf_y$p.value >= 0.05) {
  cat(">>> Differencing Y to achieve stationarity...\n")
  Y <- diff(Y)
  X <- X[-1]
  n <- length(Y)
}

# =============================================================================
# SECTION 3: QUANTILE TRANSFER ENTROPY (QTE) CALCULATION
# =============================================================================
# QTE = log(sum(rho_tau(u1))) - log(sum(rho_tau(u2)))
# where:
#   u1 = residuals from restricted model: Q_tau(Y_t | Y_{t-k})
#   u2 = residuals from unrestricted model: Q_tau(Y_t | Y_{t-k}, X_{t-k})
#   rho_tau(u) = u * (tau - I(u < 0)) is the check function
# =============================================================================

# Check function for quantile regression
check_function <- function(u, tau) {
  return(u * (tau - (u < 0)))
}

# Calculate QTE for given tau and k
calculate_qte <- function(Y, X, tau, k) {
  n <- length(Y)

  # Create lagged variables
  Y_t <- Y[(k+1):n]           # Y_t (response)
  Y_lag <- Y[1:(n-k)]         # Y_{t-k}
  X_lag <- X[1:(n-k)]         # X_{t-k}

  # Restricted model: Q_tau(Y_t | Y_{t-k})
  model_restricted <- rq(Y_t ~ Y_lag, tau = tau)
  u1 <- residuals(model_restricted)

  # Unrestricted model: Q_tau(Y_t | Y_{t-k}, X_{t-k})
  model_unrestricted <- rq(Y_t ~ Y_lag + X_lag, tau = tau)
  u2 <- residuals(model_unrestricted)

  # Calculate total quantile loss
  loss_restricted <- sum(check_function(u1, tau))
  loss_unrestricted <- sum(check_function(u2, tau))

  # QTE = log ratio of losses
  # Add small constant to avoid log(0)
  qte <- log(loss_restricted + 1e-10) - log(loss_unrestricted + 1e-10)

  return(qte)
}

# =============================================================================
# SECTION 4: STATIONARY BLOCK BOOTSTRAP
# =============================================================================
# The Stationary Bootstrap uses random block lengths from geometric distribution
# This preserves temporal dependence while imposing the null hypothesis
# =============================================================================

# Function to generate stationary bootstrap indices
stationary_bootstrap_indices <- function(n, expected_block_length) {
  # Probability of starting new block
  p <- 1 / expected_block_length

  indices <- integer(n)
  i <- 1

  while (i <= n) {
    # Random starting point
    start <- sample(1:n, 1)

    # Generate block (geometric distribution for length)
    j <- start
    while (i <= n) {
      indices[i] <- ((j - 1) %% n) + 1  # Wrap around
      i <- i + 1
      j <- j + 1

      # Probability of ending block
      if (runif(1) < p) break
    }
  }

  return(indices)
}

# Bootstrap test for QTE
bootstrap_qte_test <- function(Y, X, tau, k, B, expected_block_length = NULL) {
  n <- length(Y)

  # Default block length: n^(1/3) (common choice)
  if (is.null(expected_block_length)) {
    expected_block_length <- ceiling(n^(1/3))
  }

  # Calculate observed QTE
  qte_observed <- calculate_qte(Y, X, tau, k)

  # Bootstrap distribution under null hypothesis
  # Null: X does not cause Y
  # Impose null by resampling X and Y independently
  qte_bootstrap <- numeric(B)

  for (b in 1:B) {
    # Generate TWO independent bootstrap index samples
    # This breaks any dependence between X and Y (imposes null)
    indices_X <- stationary_bootstrap_indices(n, expected_block_length)
    indices_Y <- stationary_bootstrap_indices(n, expected_block_length)

    X_star <- X[indices_X]
    Y_star <- Y[indices_Y]

    # Calculate QTE on bootstrap sample
    qte_bootstrap[b] <- calculate_qte(Y_star, X_star, tau, k)
  }

  # Centered bootstrap statistic
  # Z_b = QTE_bootstrap - QTE_observed
  Z_b <- qte_bootstrap - qte_observed

  # P-value: proportion of Z_b >= QTE_observed
  p_value <- mean(Z_b >= qte_observed)

  return(list(
    qte = qte_observed,
    p_value = p_value,
    bootstrap_dist = qte_bootstrap,
    block_length = expected_block_length
  ))
}

# =============================================================================
# SECTION 5: RUN QTE ANALYSIS
# =============================================================================
cat("\n")
cat("=============================================================================\n")
cat("RUNNING QTE ANALYSIS\n")
cat("=============================================================================\n")

# Storage for results
results <- data.frame(
  tau = numeric(),
  k = numeric(),
  QTE = numeric(),
  p_value = numeric(),
  significant = character(),
  stringsAsFactors = FALSE
)

# Calculate expected block length
block_length <- ceiling(n^(1/3))
cat(sprintf("Expected block length for Stationary Bootstrap: %d\n", block_length))
cat("(Based on n^(1/3) rule of thumb)\n\n")

# Run analysis for each combination of tau and k
for (tau in tau_values) {
  cat(sprintf(">>> Processing tau = %.2f ...\n", tau))

  for (k in k_values) {
    cat(sprintf("    Lag k = %d: ", k))

    # Run bootstrap test
    test_result <- bootstrap_qte_test(Y, X, tau, k, B, block_length)

    # Determine significance
    sig <- ifelse(test_result$p_value < alpha, "***", "")

    cat(sprintf("QTE = %.4f, p-value = %.4f %s\n",
                test_result$qte, test_result$p_value, sig))

    # Store results
    results <- rbind(results, data.frame(
      tau = tau,
      k = k,
      QTE = test_result$qte,
      p_value = test_result$p_value,
      significant = sig,
      stringsAsFactors = FALSE
    ))
  }
  cat("\n")
}

# =============================================================================
# SECTION 6: RESULTS PRESENTATION
# =============================================================================
cat("=============================================================================\n")
cat("QTE RESULTS TABLE\n")
cat("=============================================================================\n")

# Reshape results for better presentation
cat("\nTable 1: QTE Estimates\n")
cat("-----------------------------------------------------------------------------\n")
cat(sprintf("%-10s", "Lag (k)"))
for (tau in tau_values) {
  cat(sprintf("| %-15s", paste0("tau = ", tau)))
}
cat("\n")
cat("-----------------------------------------------------------------------------\n")

for (k in k_values) {
  cat(sprintf("%-10d", k))
  for (tau in tau_values) {
    row <- results[results$tau == tau & results$k == k, ]
    qte_str <- sprintf("%.4f%s", row$QTE, row$significant)
    cat(sprintf("| %-15s", qte_str))
  }
  cat("\n")
}
cat("-----------------------------------------------------------------------------\n")
cat("*** indicates significance at 5% level\n\n")

cat("\nTable 2: Bootstrap P-values\n")
cat("-----------------------------------------------------------------------------\n")
cat(sprintf("%-10s", "Lag (k)"))
for (tau in tau_values) {
  cat(sprintf("| %-15s", paste0("tau = ", tau)))
}
cat("\n")
cat("-----------------------------------------------------------------------------\n")

for (k in k_values) {
  cat(sprintf("%-10d", k))
  for (tau in tau_values) {
    row <- results[results$tau == tau & results$k == k, ]
    cat(sprintf("| %-15.4f", row$p_value))
  }
  cat("\n")
}
cat("-----------------------------------------------------------------------------\n")

# =============================================================================
# SECTION 7: INTERPRETATION AND DISCUSSION
# =============================================================================
cat("\n")
cat("=============================================================================\n")
cat("INTERPRETATION\n")
cat("=============================================================================\n")

# Find significant results
sig_results <- results[results$significant == "***", ]

if (nrow(sig_results) > 0) {
  cat("Significant information transfer detected:\n\n")

  # Analyze by quantile
  for (tau in tau_values) {
    tau_sig <- sig_results[sig_results$tau == tau, ]
    if (nrow(tau_sig) > 0) {
      if (tau == 0.10) {
        cat(sprintf("- At tau = 0.10 (lower tail/losses): %d significant lag(s)\n",
                    nrow(tau_sig)))
        cat("  => Uncertainty helps predict EXTREME LOSSES\n")
      } else if (tau == 0.50) {
        cat(sprintf("- At tau = 0.50 (median): %d significant lag(s)\n",
                    nrow(tau_sig)))
        cat("  => Uncertainty helps predict AVERAGE returns\n")
      } else if (tau == 0.90) {
        cat(sprintf("- At tau = 0.90 (upper tail/gains): %d significant lag(s)\n",
                    nrow(tau_sig)))
        cat("  => Uncertainty helps predict EXTREME GAINS\n")
      }
    }
  }

  # Find strongest effect
  max_qte_row <- results[which.max(results$QTE), ]
  cat(sprintf("\nStrongest information transfer: tau = %.2f, k = %d (QTE = %.4f)\n",
              max_qte_row$tau, max_qte_row$k, max_qte_row$QTE))
} else {
  cat("No significant information transfer detected at the 5% level.\n")
  cat("This suggests that economic uncertainty does not significantly\n")
  cat("predict financial returns at any quantile tested.\n")
}

cat("\n")
cat("DISCUSSION:\n")
cat("- QTE > 0 indicates that past uncertainty provides predictive information\n")
cat("  for the current return's quantile beyond what past returns provide.\n")
cat("- Significance at lower quantiles (tau = 0.10) suggests uncertainty\n")
cat("  is particularly informative for predicting extreme losses.\n")
cat("- Significance at upper quantiles (tau = 0.90) suggests uncertainty\n")
cat("  helps predict extreme gains.\n")
cat("- The lag (k) with the strongest effect indicates the time horizon\n")
cat("  over which uncertainty affects returns.\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================
write.csv(results, "Part_B_QTE_Results.csv", row.names = FALSE)
cat("\nResults saved to 'Part_B_QTE_Results.csv'\n")

# =============================================================================
# END OF PART B
# =============================================================================
cat("\n=============================================================================\n")
cat("END OF QTE ANALYSIS\n")
cat("=============================================================================\n")
