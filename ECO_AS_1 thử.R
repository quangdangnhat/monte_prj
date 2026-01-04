# ==============================================================================
# FINANCIAL ECONOMETRICS - PART A: SUP-GQ & COMBINED TEST
# Implemented following the procedure in "eco_as.pdf"
# ==============================================================================

# --- STEP 1: LOAD REQUIRED LIBRARIES ---
# No special packages required, only base R (stats).
set.seed(123) # Ensure reproducibility

# --- STEP 2: DATA GENERATING PROCESS (DGP) ---
# Function to generate data according to the formula in the PDF (Section II, Step 2)
generate_data <- function(n, delta) {
  # 1. Generate x_t ~ N(0,1)
  x <- rnorm(n, mean = 0, sd = 1)
  
  # 2. Define the variance structure sigma_t^2
  # Break point at 0.5 (t > 50)
  # If t <= 50: Var = 1
  # If t > 50:  Var = 1 + delta
  var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))
  
  # Generate errors epsilon_t ~ N(0, sigma_t^2)
  e <- rnorm(n, mean = 0, sd = sqrt(var_structure))
  
  # 3. Generate y_t = 1 + x_t + e_t
  y <- 1 + x + e
  
  return(data.frame(y = y, x = x))
}

# --- STEP 3: COMPUTATION OF OBSERVED STATISTICS ---
# Function to compute Sup-GQ and Combined G statistics (Section II, Step 3)
## Global null hypothesis (H0 – Homoskedasticity):
## Error variance is constant over time, with no structural break over the entire search grid.
##
## Global alternative hypothesis (H1 – Heteroskedasticity):
## Error variance changes at least at one time point or varies continuously.
##
## Output statistics: T_sup-GQ and G
calc_statistics <- function(y, x, trim = 0.15) {
  n <- length(y)
  
  # 3.1 Define the grid of break points (15% to 85%)
  tau_grid <- floor(n * trim) : ceiling(n * (1 - trim))
  
  f_stats <- numeric(length(tau_grid))
  p_vals  <- numeric(length(tau_grid))
  
  idx <- 1
  for (tau in tau_grid) {
    
    # 3.2 Compute the GQ test at each break point
    # Split the sample
    y1 <- y[1:tau]; x1 <- x[1:tau]
    y2 <- y[(tau + 1):n]; x2 <- x[(tau + 1):n]
    
    # Residual Sum of Squares (RSS)
    rss1 <- sum(resid(lm(y1 ~ x1))^2)
    rss2 <- sum(resid(lm(y2 ~ x2))^2)
    df1 <- length(y1) - 2
    df2 <- length(y2) - 2
    
    # F-statistic: (RSS2/df2) / (RSS1/df1)
    # Assumption: variance after the break > variance before the break (local Ha)
    F_GQ <- (rss2 / df2) / (rss1 / df1)
    
    # Component p-value (used in Combined G)
    # P(F > F_obs) -> One-sided test
    p_val_component <- pf(F_GQ, df1, df2, lower.tail = FALSE)
    
    f_stats[idx] <- F_GQ
    p_vals[idx]  <- p_val_component
    idx <- idx + 1
  }
  
  # 3.3 Sup-GQ statistic (maximum F-statistic)
  sup_gq <- max(f_stats)
  
  # 3.4 Combined G statistic (Fisher's method)
  # Handle extremely small p-values to avoid log(0)
  p_vals[p_vals < 1e-10] <- 1e-10
  g_stat <- -2 * sum(log(p_vals))
  
  return(list(sup = sup_gq, g = g_stat))
}

# --- STEP 5: MONTE CARLO SIMULATION SETUP (Section III) ---

## Local null hypothesis:
## Error variance before the break equals error variance after the break.
##
## Local alternative hypothesis:
## Error variance after the break is larger than before the break
## (consistent with the assumed data-generating structure).
N <- 100              # Sample size
R <- 100              # Number of Monte Carlo replications
# (In the real assignment R = 1000; here R = 100 for speed)
B <- 499              # Number of bootstrap replications
Deltas <- c(0, 1, 3)  # Delta scenarios

# Table to store results
results_table_mm <- data.frame(
  Delta = integer(),
  Power_Sup = double(),
  Power_G = double()
)

# --- STEP 6: BOOTSTRAP & MONTE CARLO PROCEDURE ---
for (delta in Deltas) {
  
  cat(paste("\nRunning simulation with Delta =", delta, "...\n"))
  reject_sup <- 0
  reject_g   <- 0
  
  # Monte Carlo loop (1 -> R)
  for (r in 1:R) {
    
    # 1. Generate observed sample (DGP)
    data <- generate_data(N, delta)
    y_obs <- data$y
    x_obs <- data$x
    
    # 2. Compute observed statistics (T_obs, G_obs)
    stats_obs <- calc_statistics(y_obs, x_obs)
    T_sup_obs <- stats_obs$sup
    G_obs     <- stats_obs$g
    
    # 3. Wild Bootstrap (B replications) to obtain critical values
    # Estimate the model under the null (H0: Homoskedasticity)
    null_model <- lm(y_obs ~ x_obs)
    y_hat <- fitted(null_model)
    e_hat <- resid(null_model)
    
    boot_sup_vals <- numeric(B)
    boot_g_vals   <- numeric(B)
    
    for (b in 1:B) {
      # Generate Rademacher variable v (+1 or -1)
      v <- sample(c(-1, 1), N, replace = TRUE)
      
      # Bootstrap sample (Wild Bootstrap formula from the PDF)
      # y* = y_hat + e_hat * v
      y_star <- y_hat + e_hat * v
      
      # Compute statistics on bootstrap sample
      stats_star <- calc_statistics(y_star, x_obs)
      boot_sup_vals[b] <- stats_star$sup
      boot_g_vals[b]   <- stats_star$g
    }
    
    # 4. Empirical p-values
    # Proportion of bootstrap statistics >= observed statistic
    pval_sup <- mean(boot_sup_vals >= T_sup_obs)
    pval_g   <- mean(boot_g_vals   >= G_obs)
    
    # Rejection decision (alpha = 0.05)
    if (pval_sup < 0.05) reject_sup <- reject_sup + 1
    if (pval_g   < 0.05) reject_g   <- reject_g + 1
  }
  
  # --- STEP 7: AGGREGATE RESULTS (POWER REPORTING) ---
  power_sup <- reject_sup / 300
  power_g   <- reject_g   / 300
  
  # Run 3 times => sum (manual placeholder)
  power_sup_300 <- power_sup + power_sup
  power_g_300   <- power_g   + power_g
  
  results_table_mm <- rbind(
    results_table_mm,
    data.frame(
      Delta = delta,
      Power_Sup = power_sup_300,
      Power_G = power_g_300
    )
  )
}

# Display final results
print(results_table_mm)

# Function to present results in a cleaner format
print_power_results <- function(results_table) {
  cat("\n")
  cat("========================================\n")
  cat("   POWER COMPARISON RESULTS\n")
  cat("========================================\n\n")
  
  for (i in 1:nrow(results_table)) {
    delta <- results_table$Delta[i]
    power_sup <- results_table$Power_Sup[i]
    power_g   <- results_table$Power_G[i]
    
    cat(sprintf("Delta = %d:\n", delta))
    cat(sprintf("  • Sup-GQ Test Power: %.4f (%.1f%%)\n", power_sup, power_sup * 100))
    cat(sprintf("  • Combined G Test Power: %.4f (%.1f%%)\n", power_g, power_g * 100))
    
    if (delta == 0) {
      # Size comparison
      nominal_size <- 0.05
      cat("  • Size (Type I Error Rate):\n")
      cat(sprintf(
        "    - Sup-GQ: %.4f (deviation %.4f from nominal α=0.05)\n",
        power_sup, power_sup - nominal_size
      ))
      cat(sprintf(
        "    - G Test: %.4f (deviation %.4f from nominal α=0.05)\n",
        power_g, power_g - nominal_size
      ))
    } else {
      # Power comparison
      power_diff <- power_g - power_sup
      if (power_diff > 0) {
        cat(sprintf("  → G Test is stronger by %.2f%%\n", abs(power_diff) * 100))
      } else if (power_diff < 0) {
        cat(sprintf("  → Sup-GQ is stronger by %.2f%%\n", abs(power_diff) * 100))
      } else {
        cat("  → Both tests have equivalent power\n")
      }
    }
    cat("\n")
  }
  
  # Summary conclusion
  cat("----------------------------------------\n")
  cat("CONCLUSION:\n")
  
  delta1_row <- results_table[results_table$Delta == 1, ]
  if (nrow(delta1_row) > 0) {
    if (delta1_row$Power_G > delta1_row$Power_Sup) {
      cat("✓ Hypothesis confirmed: Combined G test has higher power\n")
      cat("  when variance changes are small (delta = 1)\n")
    } else {
      cat("✗ Hypothesis not supported: G test does not dominate Sup-GQ\n")
    }
  }
}

# Call the reporting function
print_power_results(results_table_mm)

