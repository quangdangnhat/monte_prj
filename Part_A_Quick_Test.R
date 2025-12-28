# =============================================================================
# QUICK TEST VERSION - Part A Monte Carlo
# Use smaller R and B for quick testing (full version takes ~1 hour)
# =============================================================================

rm(list = ls())
set.seed(123)

if (!require("lmtest")) install.packages("lmtest", repos = "http://cran.r-project.org")
library(lmtest)

# QUICK TEST PARAMETERS
N <- 100
R <- 50     # Reduced from 1000
B <- 99     # Reduced from 499
Delta <- c(0, 1, 3)
trim <- 0.15
alpha <- 0.05

cat("=============================================================================\n")
cat("QUICK TEST: Part A Monte Carlo (R=50, B=99)\n")
cat("=============================================================================\n")

# Functions (same as full version)
generate_data <- function(n, delta) {
  x <- rnorm(n, mean = 0, sd = 1)
  var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))
  epsilon <- rnorm(n, mean = 0, sd = sqrt(var_structure))
  y <- 1 + x + epsilon
  return(data.frame(y = y, x = x))
}

calc_statistics <- function(y, x, trim = 0.15) {
  n <- length(y)
  tau_min <- floor(n * trim)
  tau_max <- ceiling(n * (1 - trim))
  tau_grid <- tau_min:tau_max
  M <- length(tau_grid)
  f_stats <- numeric(M)
  p_vals <- numeric(M)

  idx <- 1
  for (tau in tau_grid) {
    y1 <- y[1:tau]; x1 <- x[1:tau]
    y2 <- y[(tau + 1):n]; x2 <- x[(tau + 1):n]
    rss1 <- sum(resid(lm(y1 ~ x1))^2)
    rss2 <- sum(resid(lm(y2 ~ x2))^2)
    df1 <- length(y1) - 2
    df2 <- length(y2) - 2
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

breusch_pagan_test <- function(y, x) {
  model <- lm(y ~ x)
  bp_test <- bptest(model)
  return(bp_test$p.value)
}

white_test <- function(y, x) {
  model <- lm(y ~ x)
  resid_sq <- resid(model)^2
  x_sq <- x^2
  aux_model <- lm(resid_sq ~ x + x_sq)
  n <- length(y)
  r_squared <- summary(aux_model)$r.squared
  white_stat <- n * r_squared
  p_value <- pchisq(white_stat, df = 2, lower.tail = FALSE)
  return(p_value)
}

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
    stats_star <- calc_statistics(y_star, x_obs)
    boot_sup_vals[b] <- stats_star$sup
    boot_g_vals[b] <- stats_star$g
  }

  pval_sup <- mean(boot_sup_vals >= stats_obs$sup)
  pval_g <- mean(boot_g_vals >= stats_obs$g)
  return(list(pval_sup = pval_sup, pval_g = pval_g))
}

# Main simulation
results_table <- data.frame(
  Delta = numeric(), Power_SupGQ = numeric(), Power_G = numeric(),
  Power_BP = numeric(), Power_White = numeric()
)

for (delta in Delta) {
  cat(sprintf("\n>>> Delta = %d ...\n", delta))
  reject_sup <- reject_g <- reject_bp <- reject_white <- 0

  for (r in 1:R) {
    data <- generate_data(N, delta)
    y_obs <- data$y; x_obs <- data$x
    stats_obs <- calc_statistics(y_obs, x_obs)
    boot_results <- wild_bootstrap(y_obs, x_obs, B, stats_obs)
    pval_bp <- breusch_pagan_test(y_obs, x_obs)
    pval_white <- white_test(y_obs, x_obs)

    if (boot_results$pval_sup < alpha) reject_sup <- reject_sup + 1
    if (boot_results$pval_g < alpha) reject_g <- reject_g + 1
    if (pval_bp < alpha) reject_bp <- reject_bp + 1
    if (pval_white < alpha) reject_white <- reject_white + 1

    if (r %% 10 == 0) cat(sprintf("   Completed %d/%d\n", r, R))
  }

  results_table <- rbind(results_table, data.frame(
    Delta = delta,
    Power_SupGQ = reject_sup / R,
    Power_G = reject_g / R,
    Power_BP = reject_bp / R,
    Power_White = reject_white / R
  ))
}

cat("\n=============================================================================\n")
cat("QUICK TEST RESULTS (not for submission - use full version)\n")
cat("=============================================================================\n")
print(results_table)
cat("=============================================================================\n")
cat("\nNote: For actual submission, run Part_A_Monte_Carlo.R with R=1000, B=499\n")
