# =============================================================================
# FINANCIAL ECONOMETRICS COURSEWORK 2025-2026
# PART A: Monte Carlo Power Comparison for Heteroscedasticity Tests
# =============================================================================
# This script compares the finite-sample power of:
#   1. Sup Goldfeld-Quandt (Sup-GQ) Test
#   2. Combined P-Value Test (Fisher's G statistic)
#   3. Breusch-Pagan Test (existing test)
#   4. White Test (existing test)
# =============================================================================

# Clear environment and set seed for reproducibility
rm(list = ls())
set.seed(123)

# Load required packages
if (!require("lmtest")) install.packages("lmtest", repos = "http://cran.r-project.org")
library(lmtest)

# =============================================================================
# SIMULATION PARAMETERS (as specified in coursework)
# =============================================================================
N <- 100            # Sample size
R <- 1000           # Monte Carlo replications
B <- 499            # Bootstrap replications
Delta <- c(0, 1, 3) # Heteroscedasticity levels (0 = size check)
trim <- 0.15        # Trimming percentage (15% on both ends)
alpha <- 0.05       # Significance level

cat("=============================================================================\n")
cat("PART A: Monte Carlo Power Comparison\n")
cat("=============================================================================\n")
cat(sprintf("Sample Size (N): %d\n", N))
cat(sprintf("Monte Carlo Replications (R): %d\n", R))
cat(sprintf("Bootstrap Replications (B): %d\n", B))
cat(sprintf("Trimming: %.0f%% on both ends\n", trim * 100))
cat(sprintf("Significance Level: %.0f%%\n", alpha * 100))
cat("=============================================================================\n\n")

# =============================================================================
# FUNCTION: Generate Data with Heteroscedasticity Break at tau = 0.5
# =============================================================================
# Model: y_t = 1 + x_t + epsilon_t
# where x_t ~ N(0,1)
# and epsilon_t ~ N(0, sigma_t^2)
# with sigma_t^2 = 1 + delta * I(t/N > 0.5)
# =============================================================================
generate_data <- function(n, delta) {
  x <- rnorm(n, mean = 0, sd = 1)

  # Variance structure: break at midpoint (tau = 0.5)
  # First half: variance = 1
  # Second half: variance = 1 + delta
  var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))

  # Generate errors with heteroscedastic variance
  epsilon <- rnorm(n, mean = 0, sd = sqrt(var_structure))

  # Generate response variable
  y <- 1 + x + epsilon

  return(data.frame(y = y, x = x))
}

# =============================================================================
# FUNCTION: Calculate Sup-GQ and Combined G Statistics
# =============================================================================
# Sup-GQ: Maximum F-statistic across all split points
# G (Fisher): -2 * sum(log(p_values)) - combines all p-values
# =============================================================================
calc_statistics <- function(y, x, trim = 0.15) {
  n <- length(y)

  # Define grid of split points (excluding trimmed portions)
  tau_min <- floor(n * trim)
  tau_max <- ceiling(n * (1 - trim))
  tau_grid <- tau_min:tau_max

  M <- length(tau_grid)  # Number of split points
  f_stats <- numeric(M)
  p_vals <- numeric(M)

  idx <- 1
  for (tau in tau_grid) {
    # Split data at tau
    y1 <- y[1:tau]
    x1 <- x[1:tau]
    y2 <- y[(tau + 1):n]
    x2 <- x[(tau + 1):n]

    # Fit separate regressions and get RSS
    rss1 <- sum(resid(lm(y1 ~ x1))^2)
    rss2 <- sum(resid(lm(y2 ~ x2))^2)

    # Degrees of freedom (n - k, where k = 2 for intercept + slope)
    df1 <- length(y1) - 2
    df2 <- length(y2) - 2

    # Goldfeld-Quandt F-statistic: ratio of variances
    # F = (RSS2/df2) / (RSS1/df1)
    F_GQ <- (rss2 / df2) / (rss1 / df1)

    # P-value from F-distribution (one-tailed, upper tail)
    p_val <- pf(F_GQ, df2, df1, lower.tail = FALSE)

    f_stats[idx] <- F_GQ
    p_vals[idx] <- p_val
    idx <- idx + 1
  }

  # Sup-GQ: Maximum F-statistic
  sup_gq <- max(f_stats)

  # Fisher's G statistic: -2 * sum(log(p_values))
  # Clip p-values to avoid log(0)
  p_vals_clipped <- pmax(p_vals, 1e-10)
  g_stat <- -2 * sum(log(p_vals_clipped))

  return(list(sup = sup_gq, g = g_stat, f_stats = f_stats, p_vals = p_vals))
}

# =============================================================================
# FUNCTION: Breusch-Pagan Test
# =============================================================================
# Tests H0: homoscedasticity vs Ha: variance is a function of regressors
# =============================================================================
breusch_pagan_test <- function(y, x) {
  model <- lm(y ~ x)
  bp_test <- bptest(model)
  return(bp_test$p.value)
}

# =============================================================================
# FUNCTION: White Test
# =============================================================================
# More general test that includes squared terms and cross-products
# =============================================================================
white_test <- function(y, x) {
  model <- lm(y ~ x)
  resid_sq <- resid(model)^2

  # Auxiliary regression: e^2 ~ x + x^2
  x_sq <- x^2
  aux_model <- lm(resid_sq ~ x + x_sq)

  # Test statistic: n * R^2 ~ Chi-squared(2)
  n <- length(y)
  r_squared <- summary(aux_model)$r.squared
  white_stat <- n * r_squared
  p_value <- pchisq(white_stat, df = 2, lower.tail = FALSE)

  return(p_value)
}

# =============================================================================
# WILD BOOTSTRAP PROCEDURE
# =============================================================================
# Uses Rademacher distribution: v ~ {-1, +1} with equal probability
# Bootstrap sample: y* = y_hat + e_hat * v
# This preserves heteroscedasticity structure under the null
# =============================================================================
wild_bootstrap <- function(y_obs, x_obs, B, stats_obs) {
  n <- length(y_obs)

  # Fit null model (homoscedastic)
  null_model <- lm(y_obs ~ x_obs)
  y_hat <- fitted(null_model)
  e_hat <- resid(null_model)

  # Storage for bootstrap statistics
  boot_sup_vals <- numeric(B)
  boot_g_vals <- numeric(B)

  for (b in 1:B) {
    # Rademacher weights
    v <- sample(c(-1, 1), n, replace = TRUE)

    # Bootstrap sample
    y_star <- y_hat + e_hat * v

    # Calculate test statistics on bootstrap sample
    stats_star <- calc_statistics(y_star, x_obs, trim = 0.15)
    boot_sup_vals[b] <- stats_star$sup
    boot_g_vals[b] <- stats_star$g
  }

  # Bootstrap p-values
  pval_sup <- mean(boot_sup_vals >= stats_obs$sup)
  pval_g <- mean(boot_g_vals >= stats_obs$g)

  return(list(pval_sup = pval_sup, pval_g = pval_g))
}

# =============================================================================
# MAIN MONTE CARLO SIMULATION
# =============================================================================
results_table <- data.frame(
  Delta = numeric(),
  Power_SupGQ = numeric(),
  Power_G = numeric(),
  Power_BP = numeric(),
  Power_White = numeric()
)

for (delta in Delta) {
  cat(sprintf("\n>>> Running simulation for Delta = %d ...\n", delta))

  # Counters for rejections
  reject_sup <- 0
  reject_g <- 0
  reject_bp <- 0
  reject_white <- 0

  # Progress tracking
  pb <- txtProgressBar(min = 0, max = R, style = 3)

  for (r in 1:R) {
    # Generate data with heteroscedasticity
    data <- generate_data(N, delta)
    y_obs <- data$y
    x_obs <- data$x

    # Calculate observed test statistics (Sup-GQ and G)
    stats_obs <- calc_statistics(y_obs, x_obs, trim = 0.15)

    # Wild Bootstrap for Sup-GQ and G
    boot_results <- wild_bootstrap(y_obs, x_obs, B, stats_obs)

    # Breusch-Pagan test (asymptotic)
    pval_bp <- breusch_pagan_test(y_obs, x_obs)

    # White test (asymptotic)
    pval_white <- white_test(y_obs, x_obs)

    # Count rejections at alpha = 0.05
    if (boot_results$pval_sup < alpha) reject_sup <- reject_sup + 1
    if (boot_results$pval_g < alpha) reject_g <- reject_g + 1
    if (pval_bp < alpha) reject_bp <- reject_bp + 1
    if (pval_white < alpha) reject_white <- reject_white + 1

    setTxtProgressBar(pb, r)
  }
  close(pb)

  # Calculate empirical power (rejection rate)
  power_sup <- reject_sup / R
  power_g <- reject_g / R
  power_bp <- reject_bp / R
  power_white <- reject_white / R

  cat(sprintf("   Sup-GQ Power: %.3f\n", power_sup))
  cat(sprintf("   G (Fisher) Power: %.3f\n", power_g))
  cat(sprintf("   Breusch-Pagan Power: %.3f\n", power_bp))
  cat(sprintf("   White Test Power: %.3f\n", power_white))

  # Store results
  results_table <- rbind(results_table, data.frame(
    Delta = delta,
    Power_SupGQ = power_sup,
    Power_G = power_g,
    Power_BP = power_bp,
    Power_White = power_white
  ))
}

# =============================================================================
# RESULTS PRESENTATION
# =============================================================================
cat("\n\n")
cat("=============================================================================\n")
cat("RESULTS: Empirical Size (delta=0) and Power (delta=1,3) at 5% Level\n")
cat("=============================================================================\n")
print(results_table)
cat("=============================================================================\n")

# Create formatted table
cat("\n\nFormatted Results Table:\n")
cat("-----------------------------------------------------------------------------\n")
cat(sprintf("%-8s | %-12s | %-12s | %-12s | %-12s\n",
            "Delta", "Sup-GQ", "G (Fisher)", "Breusch-Pagan", "White"))
cat("-----------------------------------------------------------------------------\n")
for (i in 1:nrow(results_table)) {
  row_label <- ifelse(results_table$Delta[i] == 0, "0 (Size)",
                      as.character(results_table$Delta[i]))
  cat(sprintf("%-8s | %-12.3f | %-12.3f | %-12.3f | %-12.3f\n",
              row_label,
              results_table$Power_SupGQ[i],
              results_table$Power_G[i],
              results_table$Power_BP[i],
              results_table$Power_White[i]))
}
cat("-----------------------------------------------------------------------------\n")

# =============================================================================
# INTERPRETATION NOTES
# =============================================================================
cat("\n\nINTERPRETATION:\n")
cat("- Delta = 0: Size check (should be close to nominal 5%)\n")
cat("- Delta = 1: Moderate heteroscedasticity (variance doubles after break)\n")
cat("- Delta = 3: Strong heteroscedasticity (variance quadruples after break)\n")
cat("\n")
cat("- Sup-GQ: Focuses on worst-case (maximum) split point\n")
cat("- G (Fisher): Combines information from ALL split points\n")
cat("- Breusch-Pagan: Tests if variance depends on regressors\n")
cat("- White: More general test (includes squared terms)\n")
cat("\n")
cat("Note: Sup-GQ and G use Wild Bootstrap; BP and White use asymptotic distribution\n")

# Save results to file
write.csv(results_table, "Part_A_Results.csv", row.names = FALSE)
cat("\nResults saved to 'Part_A_Results.csv'\n")

# =============================================================================
# END OF PART A
# =============================================================================
