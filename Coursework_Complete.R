# ╔═══════════════════════════════════════════════════════════════════════════════╗
# ║                                                                               ║
# ║         FINANCIAL ECONOMETRICS COURSEWORK 2025-2026                           ║
# ║         Complete Analysis - All in One File                                   ║
# ║                                                                               ║
# ╚═══════════════════════════════════════════════════════════════════════════════╝
#
# This file contains:
#   Part A: Monte Carlo Power Comparison (Sup-GQ vs Fisher's G)
#   Part B: Quantile Transfer Entropy (QTE) Analysis
#
# Author: [Your Name]
# Date: [Date]
#
# ═══════════════════════════════════════════════════════════════════════════════

# ╔═══════════════════════════════════════════════════════════════════════════════╗
# ║  [0] CONFIGURATION - CHỈNH SỬA Ở ĐÂY                                          ║
# ╚═══════════════════════════════════════════════════════════════════════════════╝

# Chọn phần nào sẽ chạy (TRUE = chạy, FALSE = bỏ qua)
RUN_PART_A <- TRUE
RUN_PART_B <- TRUE

# Simulation parameters (theo yêu cầu đề bài)
# Thời gian chạy ước tính: 5-10 giờ
R_SIMS <- 1000    # Monte Carlo replications
B_BOOT <- 499     # Bootstrap replications

# ╔═══════════════════════════════════════════════════════════════════════════════╗
# ║  [1] SETUP & PACKAGES                                                         ║
# ╚═══════════════════════════════════════════════════════════════════════════════╝

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
cat("║     FINANCIAL ECONOMETRICS COURSEWORK 2025-2026                           ║\n")
cat("║     Complete Analysis                                                     ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Clear environment (giữ lại config)
config_vars <- c("RUN_PART_A", "RUN_PART_B", "R_SIMS", "B_BOOT")
rm(list = setdiff(ls(), config_vars))

# Set seed for reproducibility
set.seed(123)

# Install and load packages
cat(">>> [1] Loading packages...\n")

packages <- c("quantreg", "tseries", "zoo", "xts", "dplyr", "quantmod", "parallel")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

cat("    All packages loaded successfully!\n\n")

# Display configuration
cat(">>> Configuration:\n")
cat(sprintf("    RUN_PART_A: %s\n", RUN_PART_A))
cat(sprintf("    RUN_PART_B: %s\n", RUN_PART_B))
cat(sprintf("    R (Monte Carlo reps): %d\n", R_SIMS))
cat(sprintf("    B (Bootstrap reps): %d\n", B_BOOT))
cat("\n")

# ╔═══════════════════════════════════════════════════════════════════════════════╗
# ║                                                                               ║
# ║                              PART A                                           ║
# ║              MONTE CARLO POWER COMPARISON                                     ║
# ║                                                                               ║
# ╚═══════════════════════════════════════════════════════════════════════════════╝

if (RUN_PART_A) {

  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("                         PART A: MONTE CARLO                               \n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [A.1] PARAMETERS
  # ─────────────────────────────────────────────────────────────────────────────

  N <- 100              # Sample size
  R <- R_SIMS           # Monte Carlo replications
  B <- B_BOOT           # Bootstrap replications
  Delta <- c(2)   # CHỈ CHẠY DELTA=2 (vì delta=0 và delta=1 đã có kết quả)
  trim <- 0.15          # Trimming percentage (15% on both ends)
  alpha <- 0.05         # Significance level

  cat(">>> [A.1] Parameters:\n")
  cat(sprintf("    Sample Size (N): %d\n", N))
  cat(sprintf("    Monte Carlo Replications (R): %d\n", R))
  cat(sprintf("    Bootstrap Replications (B): %d\n", B))
  cat(sprintf("    Trimming: %.0f%% on both ends\n", trim * 100))
  cat(sprintf("    Delta values: %s\n", paste(Delta, collapse = ", ")))
  cat("\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [A.2] FUNCTIONS
  # ─────────────────────────────────────────────────────────────────────────────

  cat(">>> [A.2] Defining functions...\n")

  # Function: Generate data with heteroscedasticity break at tau = 0.5
  # Model: y_t = 1 + x_t + epsilon_t
  # where epsilon_t ~ N(0, sigma_t^2) and sigma_t^2 = 1 + delta * I(t/N > 0.5)
  generate_data <- function(n, delta) {
    x <- rnorm(n, mean = 0, sd = 1)
    var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))
    epsilon <- rnorm(n, mean = 0, sd = sqrt(var_structure))
    y <- 1 + x + epsilon
    return(data.frame(y = y, x = x))
  }

  # Function: Calculate Sup-GQ and Fisher's G statistics (OPTIMIZED)
  # Sup-GQ: Maximum F-statistic across all split points
  # G (Fisher): -2 * sum(log(p_values))
  # Using direct matrix algebra instead of lm() for speed
  calc_statistics <- function(y, x, trim = 0.15) {
    n <- length(y)
    tau_min <- floor(n * trim)
    tau_max <- ceiling(n * (1 - trim))
    tau_grid <- tau_min:tau_max
    M <- length(tau_grid)

    f_stats <- numeric(M)
    p_vals <- numeric(M)

    # Pre-compute design matrix
    X_full <- cbind(1, x)

    idx <- 1
    for (tau in tau_grid) {
      # Split indices
      idx1 <- 1:tau
      idx2 <- (tau + 1):n

      # First subsample - direct OLS
      X1 <- X_full[idx1, , drop = FALSE]
      y1 <- y[idx1]
      beta1 <- solve(t(X1) %*% X1) %*% t(X1) %*% y1
      rss1 <- sum((y1 - X1 %*% beta1)^2)
      df1 <- tau - 2

      # Second subsample - direct OLS
      X2 <- X_full[idx2, , drop = FALSE]
      y2 <- y[idx2]
      beta2 <- solve(t(X2) %*% X2) %*% t(X2) %*% y2
      rss2 <- sum((y2 - X2 %*% beta2)^2)
      df2 <- n - tau - 2

      # F-statistic
      F_GQ <- (rss2 / df2) / (rss1 / df1)
      p_val <- pf(F_GQ, df2, df1, lower.tail = FALSE)

      f_stats[idx] <- F_GQ
      p_vals[idx] <- p_val
      idx <- idx + 1
    }

    sup_gq <- max(f_stats)
    p_vals_clipped <- pmax(p_vals, 1e-10)
    g_stat <- -2 * sum(log(p_vals_clipped))

    return(list(sup = sup_gq, g = g_stat))
  }

  # Function: Wild Bootstrap
  # Uses Rademacher distribution: v ~ {-1, +1} with equal probability
  wild_bootstrap <- function(y_obs, x_obs, B, stats_obs) {
    n <- length(y_obs)
    null_model <- lm(y_obs ~ x_obs)
    y_hat <- fitted(null_model)
    e_hat <- resid(null_model)

    boot_sup_vals <- numeric(B)
    boot_g_vals <- numeric(B)

    for (b in 1:B) {
      v <- sample(c(-1, 1), n, replace = TRUE)
      y_star <- y_hat + e_hat * v
      stats_star <- calc_statistics(y_star, x_obs, trim = 0.15)
      boot_sup_vals[b] <- stats_star$sup
      boot_g_vals[b] <- stats_star$g
    }

    pval_sup <- mean(boot_sup_vals >= stats_obs$sup)
    pval_g <- mean(boot_g_vals >= stats_obs$g)

    return(list(pval_sup = pval_sup, pval_g = pval_g))
  }

  cat("    Functions defined.\n\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [A.3] MONTE CARLO SIMULATION (PARALLEL)
  # ─────────────────────────────────────────────────────────────────────────────

  cat(">>> [A.3] Running Monte Carlo simulation (PARALLEL)...\n")

  # Detect number of cores (giới hạn để tránh quá tải)
  n_cores <- min(detectCores() - 2, 6)  # Dùng tối đa 6 cores, để lại 2 cores cho hệ thống
  if (n_cores < 1) n_cores <- 1
  cat(sprintf("    Using %d CPU cores for parallel processing\n", n_cores))

  # Function for single Monte Carlo replication
  single_mc_rep <- function(r, N, delta, B, alpha) {
    # Generate data
    data <- generate_data(N, delta)
    y_obs <- data$y
    x_obs <- data$x

    # Calculate observed statistics
    stats_obs <- calc_statistics(y_obs, x_obs, trim = 0.15)

    # Wild Bootstrap for Sup-GQ and G
    boot_results <- wild_bootstrap(y_obs, x_obs, B, stats_obs)

    # Return rejection indicators
    return(c(
      reject_sup = as.integer(boot_results$pval_sup < alpha),
      reject_g = as.integer(boot_results$pval_g < alpha)
    ))
  }

  results_A <- data.frame(
    Delta = numeric(),
    Power_SupGQ = numeric(),
    Power_G = numeric()
  )

  for (delta in Delta) {
    cat(sprintf("\n    Delta = %d: ", delta))
    start_time <- Sys.time()

    # Create cluster for parallel processing
    cl <- makeCluster(n_cores)

    # Export functions and variables to cluster
    clusterExport(cl, c("generate_data", "calc_statistics", "wild_bootstrap",
                        "single_mc_rep", "N", "delta", "B", "alpha"), envir = environment())

    # Set seed for reproducibility on each worker
    clusterSetRNGStream(cl, 123 + which(Delta == delta))

    # Run parallel Monte Carlo
    mc_results <- parLapply(cl, 1:R, function(r) {
      single_mc_rep(r, N, delta, B, alpha)
    })

    # Stop cluster
    stopCluster(cl)

    # Aggregate results
    reject_sup <- sum(sapply(mc_results, function(x) x["reject_sup"]))
    reject_g <- sum(sapply(mc_results, function(x) x["reject_g"]))

    # Calculate power
    results_A <- rbind(results_A, data.frame(
      Delta = delta,
      Power_SupGQ = reject_sup / R,
      Power_G = reject_g / R
    ))

    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
    cat(sprintf(" Done in %.1f mins (Sup=%.3f, G=%.3f)\n",
                elapsed, reject_sup/R, reject_g/R))
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # [A.4] PART A RESULTS
  # ─────────────────────────────────────────────────────────────────────────────

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("                    PART A RESULTS                                         \n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("\n")
  cat("───────────────────────────────────────────\n")
  cat(sprintf("%8s %15s %12s\n", "Delta", "Power (SupGQ)", "Power (G)"))
  cat("───────────────────────────────────────────\n")

  for (i in 1:nrow(results_A)) {
    cat(sprintf("%8d %15.3f %12.3f\n",
                results_A$Delta[i],
                results_A$Power_SupGQ[i],
                results_A$Power_G[i]))
  }
  cat("───────────────────────────────────────────\n")

  cat("\nInterpretation:\n")
  cat("• Delta = 0: Size check (should be ≈ 0.05 = nominal level)\n")
  cat("• Delta = 1: Variance doubles after midpoint\n")
  cat("• Delta = 2: Variance triples after midpoint\n")
  cat("• Higher power = better at detecting heteroscedasticity\n")
  cat("\n")

} else {
  cat(">>> PART A: Skipped (RUN_PART_A = FALSE)\n\n")
}

# ╔═══════════════════════════════════════════════════════════════════════════════╗
# ║                                                                               ║
# ║                              PART B                                           ║
# ║              QUANTILE TRANSFER ENTROPY (QTE)                                  ║
# ║                                                                               ║
# ╚═══════════════════════════════════════════════════════════════════════════════╝

if (RUN_PART_B) {

  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("                         PART B: QTE ANALYSIS                              \n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [B.1] PARAMETERS
  # ─────────────────────────────────────────────────────────────────────────────

  B_qte <- B_BOOT                     # Bootstrap replications
  tau_values <- c(0.10, 0.50, 0.90)   # Quantiles to test
  k_values <- c(1, 2, 3, 5)           # Lag values
  alpha_qte <- 0.05                   # Significance level

  cat(">>> [B.1] Parameters:\n")
  cat(sprintf("    Bootstrap Replications (B): %d\n", B_qte))
  cat(sprintf("    Quantiles (tau): %s\n", paste(tau_values, collapse = ", ")))
  cat(sprintf("    Lags (k): %s\n", paste(k_values, collapse = ", ")))
  cat("\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [B.2] DATA PREPARATION
  # ─────────────────────────────────────────────────────────────────────────────

  cat(">>> [B.2] Preparing data...\n")

  # Generate simulated data (replace with real data if available)
  # This simulates EPU-like uncertainty and stock returns
  generate_qte_data <- function(n = 2000, causality = 0.15) {
    # AR(1) process for uncertainty (X)
    x <- numeric(n)
    x[1] <- rnorm(1)
    for (t in 2:n) {
      x[t] <- 0.3 * x[t-1] + rnorm(1)
    }

    # Returns (Y) with potential causality from X
    y <- numeric(n)
    y[1] <- rnorm(1)
    for (t in 2:n) {
      y[t] <- 0.2 * y[t-1] + causality * x[t-1] + rnorm(1, sd = 1 + 0.3*abs(x[t-1]))
    }

    return(data.frame(uncertainty = x, returns = y))
  }

  # Generate data
  set.seed(456)
  qte_data <- generate_qte_data(n = 2000, causality = 0.15)
  X <- qte_data$uncertainty
  Y <- qte_data$returns
  n_qte <- length(Y)

  cat(sprintf("    Sample size: %d observations\n", n_qte))
  cat("    (Using simulated data - replace with real EPU/S&P500 if needed)\n\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [B.3] QTE FUNCTIONS
  # ─────────────────────────────────────────────────────────────────────────────

  cat(">>> [B.3] Defining QTE functions...\n")

  # Check function for quantile regression
  check_function <- function(u, tau) {
    return(u * (tau - (u < 0)))
  }

  # Calculate QTE for given tau and k
  # QTE = log(loss_restricted) - log(loss_unrestricted)
  calculate_qte <- function(Y, X, tau, k) {
    n <- length(Y)

    Y_t <- Y[(k+1):n]
    Y_lag <- Y[1:(n-k)]
    X_lag <- X[1:(n-k)]

    # Restricted model: Q_tau(Y_t | Y_{t-k})
    model_r <- rq(Y_t ~ Y_lag, tau = tau)
    u1 <- residuals(model_r)

    # Unrestricted model: Q_tau(Y_t | Y_{t-k}, X_{t-k})
    model_u <- rq(Y_t ~ Y_lag + X_lag, tau = tau)
    u2 <- residuals(model_u)

    loss_r <- sum(check_function(u1, tau))
    loss_u <- sum(check_function(u2, tau))

    qte <- log(loss_r + 1e-10) - log(loss_u + 1e-10)
    return(qte)
  }

  # Stationary Bootstrap indices
  stationary_bootstrap_indices <- function(n, block_length) {
    p <- 1 / block_length
    indices <- integer(n)
    i <- 1

    while (i <= n) {
      start <- sample(1:n, 1)
      j <- start
      while (i <= n) {
        indices[i] <- ((j - 1) %% n) + 1
        i <- i + 1
        j <- j + 1
        if (runif(1) < p) break
      }
    }
    return(indices)
  }

  # Bootstrap test for QTE
  bootstrap_qte_test <- function(Y, X, tau, k, B) {
    n <- length(Y)
    block_length <- ceiling(n^(1/3))

    qte_obs <- calculate_qte(Y, X, tau, k)
    qte_boot <- numeric(B)

    for (b in 1:B) {
      # Independent resampling to impose null (X does not cause Y)
      idx_X <- stationary_bootstrap_indices(n, block_length)
      idx_Y <- stationary_bootstrap_indices(n, block_length)

      X_star <- X[idx_X]
      Y_star <- Y[idx_Y]

      qte_boot[b] <- calculate_qte(Y_star, X_star, tau, k)
    }

    # Centered statistic and p-value
    Z_b <- qte_boot - qte_obs
    p_value <- mean(Z_b >= qte_obs)

    return(list(qte = qte_obs, p_value = p_value))
  }

  cat("    Functions defined.\n\n")

  # ─────────────────────────────────────────────────────────────────────────────
  # [B.4] QTE ANALYSIS
  # ─────────────────────────────────────────────────────────────────────────────

  cat(">>> [B.4] Running QTE analysis...\n")

  results_B <- data.frame(
    tau = numeric(),
    k = numeric(),
    QTE = numeric(),
    p_value = numeric(),
    significant = character(),
    stringsAsFactors = FALSE
  )

  for (tau in tau_values) {
    cat(sprintf("\n    tau = %.2f: ", tau))

    for (k in k_values) {
      result <- bootstrap_qte_test(Y, X, tau, k, B_qte)
      sig <- ifelse(result$p_value < alpha_qte, "***", "")

      results_B <- rbind(results_B, data.frame(
        tau = tau,
        k = k,
        QTE = result$qte,
        p_value = result$p_value,
        significant = sig,
        stringsAsFactors = FALSE
      ))

      cat(".")
    }
    cat(" Done")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # [B.5] PART B RESULTS
  # ─────────────────────────────────────────────────────────────────────────────

  cat("\n\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("                    PART B RESULTS                                         \n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("\n")

  # Table 2: QTE Estimates
  cat("Table 2: QTE Estimates (*** = significant at 5%)\n")
  cat("─────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-8s | %-15s | %-15s | %-15s\n", "Lag (k)", "tau = 0.10", "tau = 0.50", "tau = 0.90"))
  cat("─────────────────────────────────────────────────────────────────────────\n")

  for (k in k_values) {
    cat(sprintf("%-8d |", k))
    for (tau in tau_values) {
      row <- results_B[results_B$tau == tau & results_B$k == k, ]
      cat(sprintf(" %6.4f%-3s      |", row$QTE, row$significant))
    }
    cat("\n")
  }
  cat("─────────────────────────────────────────────────────────────────────────\n")

  # Table 3: P-values
  cat("\nTable 3: Bootstrap P-values\n")
  cat("─────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-8s | %-15s | %-15s | %-15s\n", "Lag (k)", "tau = 0.10", "tau = 0.50", "tau = 0.90"))
  cat("─────────────────────────────────────────────────────────────────────────\n")

  for (k in k_values) {
    cat(sprintf("%-8d |", k))
    for (tau in tau_values) {
      row <- results_B[results_B$tau == tau & results_B$k == k, ]
      cat(sprintf(" %-15.4f |", row$p_value))
    }
    cat("\n")
  }
  cat("─────────────────────────────────────────────────────────────────────────\n")

  cat("\nInterpretation:\n")
  cat("• tau = 0.10: Lower tail (extreme losses)\n")
  cat("• tau = 0.50: Median (average returns)\n")
  cat("• tau = 0.90: Upper tail (extreme gains)\n")
  cat("• QTE > 0 with p < 0.05: Significant information transfer from X to Y\n")
  cat("\n")

} else {
  cat(">>> PART B: Skipped (RUN_PART_B = FALSE)\n\n")
}

# ╔═══════════════════════════════════════════════════════════════════════════════╗
# ║  FINAL SUMMARY                                                                ║
# ╚═══════════════════════════════════════════════════════════════════════════════╝

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                         ANALYSIS COMPLETE                                 ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

if (RUN_PART_A) {
  cat("PART A Summary:\n")
  cat("───────────────\n")
  print(results_A)
  cat("\n")
}

if (RUN_PART_B) {
  cat("PART B Summary:\n")
  cat("───────────────\n")
  sig_results <- results_B[results_B$significant == "***", ]
  if (nrow(sig_results) > 0) {
    cat(sprintf("Found %d significant QTE result(s):\n", nrow(sig_results)))
    for (i in 1:nrow(sig_results)) {
      cat(sprintf("  - tau=%.2f, k=%d: QTE=%.4f (p=%.4f)\n",
                  sig_results$tau[i], sig_results$k[i],
                  sig_results$QTE[i], sig_results$p_value[i]))
    }
  } else {
    cat("No significant QTE results at 5% level.\n")
  }
  cat("\n")
}

cat("═══════════════════════════════════════════════════════════════════════════\n")

# ═══════════════════════════════════════════════════════════════════════════════
# END OF FILE
# ═══════════════════════════════════════════════════════════════════════════════
