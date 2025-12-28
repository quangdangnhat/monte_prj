set.seed(123)
N <- 100            
R <- 1000         
B <- 499
Delta<-c(0,1,3)
generate_data <- function(n, delta) {
  x <- rnorm(n, mean = 0, sd = 1)
  var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))
  e <- rnorm(n, mean = 0, sd = sqrt(var_structure))
  y <- 1 + x + e
  
  return(data.frame(y=y, x=x))
}



calc_statistics <- function(y, x, trim=0.15) {
  n <- length(y)
  tau_grid <- floor(n*trim) : ceiling(n*(1-trim))
  
  f_stats <- numeric(length(tau_grid))
  p_vals  <- numeric(length(tau_grid))
  
  idx <- 1
  for (tau in tau_grid) {
    y1 <- y[1:tau]; x1 <- x[1:tau]
    y2 <- y[(tau+1):n]; x2 <- x[(tau+1):n]
    rss1 <- sum(resid(lm(y1 ~ x1))^2)
    rss2 <- sum(resid(lm(y2 ~ x2))^2)
    df1 <- length(y1) - 2
    df2 <- length(y2) - 2
    F_GQ <- (rss2 / df2) / (rss1 / df1)
    p_val_component <- pf(F_GQ, df1, df2, lower.tail = FALSE)
    
    f_stats[idx] <- F_GQ
    p_vals[idx]  <- p_val_component
    idx <- idx + 1
  }
  sup_gq <- max(f_stats)
  p_vals[p_vals < 1e-10] <- 1e-10
  g_stat <- -2 * sum(log(p_vals))
  
  return(list(sup = sup_gq, g = g_stat))
}



results_table_mm <- data.frame(Delta=integer(), Power_Sup=double(), Power_G=double())

for (delta in Delta) {
  
  cat(paste("\nĐang chạy mô phỏng với Delta =", delta, "...\n"))
  reject_sup <- 0
  reject_g   <- 0
  for (r in 1:R) {
    data <- generate_data(N, delta)
    y_obs <- data$y
    x_obs <- data$x
    stats_obs <- calc_statistics(y_obs, x_obs)
    T_sup_obs <- stats_obs$sup
    G_obs     <- stats_obs$g
    null_model <- lm(y_obs ~ x_obs)
    y_hat <- fitted(null_model)
    e_hat <- resid(null_model)
    
    boot_sup_vals <- numeric(B)
    boot_g_vals   <- numeric(B)
    
    for (b in 1:B) {
      v <- sample(c(-1, 1), N, replace = TRUE)
      y_star <- y_hat + e_hat * v
      stats_star <- calc_statistics(y_star, x_obs)
      boot_sup_vals[b] <- stats_star$sup
      boot_g_vals[b]   <- stats_star$g
    }
    pval_sup <- mean(boot_sup_vals >= T_sup_obs)
    pval_g   <- mean(boot_g_vals   >= G_obs)
    if (pval_sup < 0.05) reject_sup <- reject_sup + 1
    if (pval_g   < 0.05) reject_g   <- reject_g + 1
  }
  power_sup <- reject_sup / R
  power_g   <- reject_g / R
  
  results_table_mm <- rbind(results_table_mm, data.frame(
    Delta = delta,
    Power_Sup = power_sup,
    
    Power_G = power_g
  ))
}
print(results_table_mm)

