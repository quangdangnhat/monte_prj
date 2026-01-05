# Giải Thích Chi Tiết: Coursework_Complete.R

## Dành cho Newbie - Financial Econometrics 2025-2026

---

## Mục Lục

1. [Tổng Quan](#1-tổng-quan)
2. [Cấu Hình](#2-cấu-hình-configuration)
3. [Part A: Monte Carlo Power Comparison](#3-part-a-monte-carlo-power-comparison)
4. [Part B: Quantile Transfer Entropy](#4-part-b-quantile-transfer-entropy-qte)
5. [Hướng Dẫn Chạy](#5-hướng-dẫn-chạy)
6. [Đọc Hiểu Kết Quả](#6-đọc-hiểu-kết-quả)

---

## 1. Tổng Quan

### File này làm gì?

| Part | Nội dung | Mục đích |
|------|----------|----------|
| **Part A** | Monte Carlo Simulation | So sánh power của **4 tests** phát hiện heteroscedasticity |
| **Part B** | QTE Analysis | Kiểm tra quan hệ nhân quả giữa 2 chuỗi thời gian qua các quantile |

### 4 Tests được so sánh trong Part A:

| Test | Loại | Ý tưởng chính |
|------|------|---------------|
| **Sup-GQ** | Bootstrap | Max F-statistic qua tất cả breakpoints |
| **Fisher's G** | Bootstrap | Kết hợp p-values từ tất cả breakpoints |
| **Breusch-Pagan** | Asymptotic | Variance phụ thuộc vào X? |
| **White** | Asymptotic | Variance phụ thuộc vào X và X²? |

### Heteroscedasticity là gì?

**Homoscedasticity** = Phương sai của error ĐỒNG NHẤT (cùng 1 giá trị)
**Heteroscedasticity** = Phương sai của error THAY ĐỔI theo thời gian/biến

```
Homoscedasticity:                    Heteroscedasticity:
    |   * * * * * * * * *                |                *   *
    | * * * * * * * * * *                |           * *    *
y   | * * * * * * * * * *            y   | * * *  *    *  *
    |   * * * * * * * * *                |  * *  *  *      *
    |________________________ x          |________________________ x
    Variance đều khắp nơi               Variance tăng dần
```

**Tại sao cần phát hiện?**
- OLS không còn efficient nếu có heteroscedasticity
- Standard errors bị sai → các test thống kê không đáng tin
- Cần điều chỉnh (robust SE, GLS, ...)

---

## 2. Cấu Hình (Configuration)

```r
# Dòng 21-28 trong file
RUN_PART_A <- TRUE    # Có chạy Part A không?
RUN_PART_B <- TRUE    # Có chạy Part B không?
R_SIMS <- 1000        # Số lần Monte Carlo
B_BOOT <- 499         # Số lần Bootstrap
```

### Ý nghĩa các tham số:

| Tham số | Giá trị | Ý nghĩa |
|---------|---------|---------|
| `R_SIMS` | 1000 | Lặp 1000 lần để ước lượng power chính xác |
| `B_BOOT` | 499 | 499 bootstrap samples để tính p-value |

### Thời gian chạy ước tính:
- **Part A:** 3-6 giờ (phụ thuộc số CPU cores)
- **Part B:** ~30 phút

---

## 3. Part A: Monte Carlo Power Comparison

### 3.1 Mục tiêu

So sánh **Power** (khả năng phát hiện) của 4 tests:

| Test | Tên đầy đủ | Package | Ý tưởng |
|------|------------|---------|---------|
| **Sup-GQ** | Supremum Goldfeld-Quandt | - | Tìm F-statistic LỚN NHẤT qua tất cả breakpoints |
| **Fisher's G** | Fisher's Combined Test | - | Kết hợp p-values: G = -2Σlog(p) |
| **Breusch-Pagan** | Breusch-Pagan Test | lmtest | Regress ε² ~ X |
| **White** | White Test | - | Regress ε² ~ X + X² |

### 3.2 Data Generating Process (DGP)

**Mô hình tạo dữ liệu giả lập:**

```r
# Dòng 111-117
generate_data <- function(n, delta) {
  x <- rnorm(n, mean = 0, sd = 1)
  var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))
  epsilon <- rnorm(n, mean = 0, sd = sqrt(var_structure))
  y <- 1 + x + epsilon
  return(data.frame(y = y, x = x))
}
```

**Giải thích bằng hình:**

```
Mô hình: y = 1 + x + ε

Với N = 100 observations:

Observation:  1  2  3  ...  50 | 51  52  ...  100
              ←── Nửa đầu ──→ | ←── Nửa sau ──→
Variance:          σ² = 1     |    σ² = 1 + delta

Delta = 0: σ² = 1 toàn bộ     → Homoscedastic (H₀ đúng)
Delta = 1: σ² = 1 | σ² = 2    → Variance gấp đôi ở nửa sau
Delta = 3: σ² = 1 | σ² = 4    → Variance gấp bốn ở nửa sau
```

### 3.3 Test 1 & 2: Sup-GQ và Fisher's G

```r
# Dòng 123-170
calc_statistics <- function(y, x, trim = 0.15) {
  ...
}
```

**Thuật toán chi tiết:**

```
BƯỚC 1: Xác định các điểm chia có thể (breakpoints)
        trim = 15% → τ chạy từ observation 15 đến 85 (với N=100)

BƯỚC 2: Với MỖI điểm chia τ:

        Dữ liệu: [obs 1, ..., obs τ] | [obs τ+1, ..., obs n]
                    Nhóm 1           |      Nhóm 2

        a) Chạy OLS riêng cho mỗi nhóm:
           Nhóm 1: y₁ = β₀ + β₁x₁ + ε₁  →  RSS₁ = Σ(residuals)²
           Nhóm 2: y₂ = β₀ + β₁x₂ + ε₂  →  RSS₂ = Σ(residuals)²

        b) Tính F-statistic:
           F = (RSS₂/df₂) / (RSS₁/df₁)

           Nếu variance nhóm 2 > nhóm 1 → F sẽ lớn

        c) Tính p-value của F

BƯỚC 3: Tổng hợp kết quả
        Sup-GQ = max(F)           ← Lấy F lớn nhất
        G = -2 × Σ log(p-values)  ← Kết hợp p-values
```

**Tại sao dùng `.lm.fit()` thay vì `lm()`?**

```r
# Dòng 145, 152
fit1 <- .lm.fit(X1, y1)  # Nhanh hơn và ổn định số học
```

- `.lm.fit()` là internal function của R, dùng QR decomposition
- Nhanh hơn `lm()` vì không parse formula
- Ổn định số học hơn `solve(t(X) %*% X) %*% t(X) %*% y`

### 3.4 Test 3: Breusch-Pagan Test

```r
# Dòng 172-178
breusch_pagan_test <- function(y, x) {
  model <- lm(y ~ x)
  bp_result <- bptest(model)  # từ package lmtest
  return(bp_result$p.value)
}
```

**Nguyên lý:**
1. Fit model: y = β₀ + β₁x + ε
2. Lấy residuals: ε̂
3. Test: Var(ε) có phụ thuộc vào X không?
4. Dùng LM (Lagrange Multiplier) test

**Ưu điểm:** Nhanh, dùng asymptotic distribution (không cần bootstrap)
**Nhược điểm:** Chỉ test linear dependence của variance

### 3.5 Test 4: White Test

```r
# Dòng 180-192
white_test <- function(y, x) {
  n <- length(y)
  model <- lm(y ~ x)
  resid_sq <- resid(model)^2      # ε̂²
  x_sq <- x^2
  aux_model <- lm(resid_sq ~ x + x_sq)  # Hồi quy phụ
  r_squared <- summary(aux_model)$r.squared
  white_stat <- n * r_squared     # nR² ~ χ²(2)
  p_value <- pchisq(white_stat, df = 2, lower.tail = FALSE)
  return(p_value)
}
```

**Nguyên lý:**
1. Fit model: y = β₀ + β₁x + ε
2. Lấy ε̂² (squared residuals)
3. Regress: ε̂² ~ x + x² (auxiliary regression)
4. Test statistic: nR² ~ χ²(2) under H₀

**Ưu điểm:** Phát hiện cả non-linear heteroscedasticity
**Nhược điểm:** Có thể kém power với small samples

### 3.6 Parametric Bootstrap cho Sup-GQ và G

**Vấn đề:** Sup-GQ và G không có phân phối lý thuyết chuẩn → không có critical value sẵn

**Giải pháp:** Dùng Bootstrap để xấp xỉ phân phối dưới H₀

```r
# Dòng 194-220
wild_bootstrap <- function(y_obs, x_obs, B, stats_obs) {
  # Fit model
  null_model <- lm(y_obs ~ x_obs)
  y_hat <- fitted(null_model)
  e_hat <- resid(null_model)

  # Ước lượng variance CHUNG dưới H₀
  sigma_hat <- sqrt(sum(e_hat^2) / (n - 2))

  for (b in 1:B) {
    # Tạo errors MỚI từ N(0, σ²) - ĐỒNG NHẤT!
    e_star <- rnorm(n, 0, sigma_hat)
    y_star <- y_hat + e_star

    # Tính test statistics trên bootstrap sample
    stats_star <- calc_statistics(y_star, x_obs)
  }

  # p-value = tỷ lệ bootstrap stat ≥ observed stat
  pval_sup <- mean(boot_sup_vals >= stats_obs$sup)
  pval_g <- mean(boot_g_vals >= stats_obs$g)
}
```

**QUAN TRỌNG - Tại sao dùng Parametric Bootstrap?**

```
H₀: Homoscedasticity (variance đồng nhất)

→ Dưới H₀, tất cả errors phải có CÙNG variance σ²
→ Bootstrap phải tạo errors từ N(0, σ²) với σ² ĐỒNG NHẤT
→ KHÔNG dùng wild bootstrap (e_hat * v) vì nó giữ nguyên pattern heteroscedasticity!
```

### 3.7 Monte Carlo Simulation với Parallel Processing

```r
# Dòng 229-290
n_cores <- detectCores() - 2  # Để lại 2 cores cho hệ thống
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(lmtest))  # Load lmtest trên mỗi worker
parLapply(cl, 1:R, function(r) { ... })  # Chạy song song
stopCluster(cl)
```

**Quy trình Monte Carlo:**

```
Với mỗi delta ∈ {0, 1, 3}:
    Lặp R = 1000 lần (song song trên nhiều cores):
        1. Generate data với delta
        2. Tính Sup-GQ và G observed
        3. Bootstrap 499 lần → p-values cho Sup-GQ, G
        4. Tính p-values cho BP và White (asymptotic)
        5. Reject nếu p-value < 0.05

    Power = (Số lần reject) / 1000
```

**Tại sao cần Parallel?**
- 1000 MC × 499 Bootstrap = ~500,000 iterations mỗi delta
- Không parallel: 10+ giờ
- 10 cores parallel: ~1-2 giờ

---

## 4. Part B: Quantile Transfer Entropy (QTE)

### 4.1 Quantile Regression là gì?

**OLS thông thường:** Ước lượng **TRUNG BÌNH** có điều kiện E[Y|X]

**Quantile Regression:** Ước lượng **PHÂN VỊ** có điều kiện Q_τ[Y|X]

```
τ = 0.10: Quantile 10% (đuôi trái - worst 10% outcomes)
τ = 0.50: Quantile 50% (median - trung vị)
τ = 0.90: Quantile 90% (đuôi phải - best 10% outcomes)
```

**Ví dụ thực tế:**
- τ = 0.10 của stock returns = worst 10% days → quan trọng cho risk management
- τ = 0.90 = best 10% days

### 4.2 Transfer Entropy là gì?

**Câu hỏi:** "X có chứa thông tin hữu ích để dự báo Y không?"

```
Restricted model (R):   Q_τ(Y_t | Y_{t-k})           ← Chỉ dùng Y quá khứ
Unrestricted model (U): Q_τ(Y_t | Y_{t-k}, X_{t-k})  ← Dùng cả X quá khứ

Nếu model U tốt hơn model R nhiều → X có causality lên Y
```

### 4.3 Cách tính QTE

```r
# Dòng 385-407
calculate_qte <- function(Y, X, tau, k) {
  # Chuẩn bị data với lag k
  Y_t <- Y[(k+1):n]      # Y hiện tại
  Y_lag <- Y[1:(n-k)]    # Y quá khứ
  X_lag <- X[1:(n-k)]    # X quá khứ

  # Restricted model
  model_r <- rq(Y_t ~ Y_lag, tau = tau)
  loss_r <- sum(check_function(residuals(model_r), tau))

  # Unrestricted model
  model_u <- rq(Y_t ~ Y_lag + X_lag, tau = tau)
  loss_u <- sum(check_function(residuals(model_u), tau))

  # QTE
  qte <- log(loss_r) - log(loss_u)
  return(qte)
}
```

**Công thức:**

```
QTE = log(Loss_R) - log(Loss_U)

Loss = Σ ρ_τ(residuals)  (quantile loss function)

ρ_τ(u) = u × (τ - I(u<0))
       = { τ × u      nếu u ≥ 0
       = { (τ-1) × u  nếu u < 0
```

**Interpretation:**
- QTE > 0: Loss_R > Loss_U → Model có X tốt hơn → X có causality
- QTE ≈ 0: Hai models tương đương → X không có causality

### 4.4 Stationary Bootstrap

**Vấn đề:** Time series có autocorrelation → không thể bootstrap i.i.d.

**Giải pháp:** Block bootstrap - lấy các blocks liên tiếp

```r
# Dòng 409-426
stationary_bootstrap_indices <- function(n, block_length) {
  # Block length ngẫu nhiên (geometric distribution)
  # Trung bình block length = n^(1/3) ≈ 12 với n=2000
}
```

```
Original: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ...]

Stationary Bootstrap:
  Chọn start = 3, block length = 4  → [3, 4, 5, 6]
  Chọn start = 8, block length = 2  → [8, 9]
  Chọn start = 1, block length = 3  → [1, 2, 3]
  ...

Result: [3, 4, 5, 6, 8, 9, 1, 2, 3, ...]
```

### 4.5 Bootstrap Test cho QTE

```r
# Dòng 428-452
bootstrap_qte_test <- function(Y, X, tau, k, B) {
  qte_obs <- calculate_qte(Y, X, tau, k)

  for (b in 1:B) {
    # Resample X và Y ĐỘC LẬP → phá vỡ dependence (impose H₀)
    idx_X <- stationary_bootstrap_indices(n, block_length)
    idx_Y <- stationary_bootstrap_indices(n, block_length)  # KHÁC idx_X!

    X_star <- X[idx_X]
    Y_star <- Y[idx_Y]

    qte_boot[b] <- calculate_qte(Y_star, X_star, tau, k)
  }

  p_value <- mean(qte_boot >= qte_obs)
}
```

**Key insight:** Resample X và Y ĐỘC LẬP để impose H₀ (X không cause Y)

---

## 5. Hướng Dẫn Chạy

### Trong RStudio:

1. Mở file `Coursework_Complete.R`
2. Kiểm tra configuration:
   ```r
   RUN_PART_A <- TRUE
   RUN_PART_B <- TRUE
   ```
3. Nhấn **Ctrl+Shift+Enter** (Source) hoặc nút **Source**

### Trong Terminal:

```bash
cd /path/to/monte_prj
Rscript Coursework_Complete.R
```

### Dừng chương trình:

- RStudio: Nhấn nút **Stop** (đỏ) hoặc **Esc**
- Terminal: **Ctrl+C**

### Tips:

1. **Chạy thử nhanh:** Tạm giảm R_SIMS = 50, B_BOOT = 99
2. **Theo dõi tiến độ:** Xem console - hiện delta đang chạy và thời gian
3. **Máy yếu:** Có thể chạy qua đêm

---

## 6. Đọc Hiểu Kết Quả

### Part A - Bảng Power (4 tests):

```
Table 1: Empirical Size (δ=0) and Power (δ=1,3) at 5% Significance Level
─────────────────────────────────────────────────────────────────────────────
   Delta       Sup-GQ   Fisher's G Breusch-Pagan        White
─────────────────────────────────────────────────────────────────────────────
       0        0.052        0.048        0.051        0.049   ← SIZE CHECK
       1        0.450        0.380        0.250        0.220   ← POWER (mild)
       3        0.850        0.920        0.650        0.580   ← POWER (strong)
─────────────────────────────────────────────────────────────────────────────
```

**Cách đọc:**

| Delta | Ý nghĩa | Kỳ vọng |
|-------|---------|---------|
| 0 | Size check (H₀ đúng) | ≈ 0.05 (nominal level) |
| 1 | Heteroscedasticity nhẹ | 0.2 - 0.5 |
| 3 | Heteroscedasticity mạnh | 0.6 - 0.95 |

**Interpretation:**
- **Delta = 0** phải ≈ 0.05. Nếu > 0.10 → test bị **over-sized**
- So sánh 4 tests: Test nào có power cao hơn ở delta = 1, 3 → test đó tốt hơn
- **Sup-GQ và G** thường mạnh hơn vì chúng được thiết kế cho structural break
- **BP và White** là general tests, có thể kém power với break cụ thể

### Part B - Bảng QTE:

```
Table 2: QTE Estimates (*** = significant at 5%)
─────────────────────────────────────────────────────────────────────────
Lag (k)  | tau = 0.10      | tau = 0.50      | tau = 0.90
─────────────────────────────────────────────────────────────────────────
1        | 0.0234***       | 0.0089          | 0.0312***
2        | 0.0198***       | 0.0045          | 0.0256**
3        | 0.0156          | 0.0032          | 0.0189
5        | 0.0098          | 0.0021          | 0.0134
─────────────────────────────────────────────────────────────────────────
```

**Cách đọc:**
- `***` = p-value < 0.05 → có causality có ý nghĩa thống kê
- QTE > 0 với `***` → X có ảnh hưởng đến Y tại quantile đó

**Patterns thường thấy:**
- Significant ở τ = 0.10, 0.90 nhưng không ở 0.50 → **Tail risk spillover**
- QTE giảm khi k tăng → Effect giảm dần theo thời gian

---

## Tổng Hợp Các Khái Niệm

| Thuật ngữ | Định nghĩa |
|-----------|------------|
| **Size** | P(Reject H₀ \| H₀ đúng) - tỷ lệ "báo động giả" |
| **Power** | P(Reject H₀ \| H₁ đúng) - khả năng phát hiện |
| **Sup-GQ** | max(F-statistics) qua tất cả breakpoints |
| **Fisher's G** | -2Σlog(p-values) - combined test |
| **Breusch-Pagan** | LM test cho Var(ε) ~ X |
| **White** | nR² test cho Var(ε) ~ X + X² |
| **Bootstrap** | Resampling để xấp xỉ phân phối |
| **QTE** | log(Loss_R) - log(Loss_U) - đo causality |
| **Quantile τ** | Điểm mà τ% observations nằm dưới |

---

## Packages Sử Dụng

| Package | Mục đích |
|---------|----------|
| `parallel` | Parallel processing (nhiều cores) |
| `lmtest` | Breusch-Pagan test - hàm `bptest()` |
| `quantreg` | Quantile regression - hàm `rq()` |
| `tseries` | Time series analysis |
| `zoo`, `xts` | Time series objects |
| `dplyr` | Data manipulation |
| `quantmod` | Financial data |

---

## Sơ Đồ Tóm Tắt

```
┌─────────────────────────────────────────────────────────────────┐
│                    COURSEWORK_COMPLETE.R                        │
├─────────────────────────────────────────────────────────────────┤
│  PART A: Monte Carlo Power Comparison (4 Tests)                 │
│  ├── generate_data(): Tạo data với heteroscedasticity           │
│  ├── calc_statistics(): Tính Sup-GQ và Fisher's G               │
│  ├── breusch_pagan_test(): BP test (lmtest package)             │
│  ├── white_test(): White test (nR² statistic)                   │
│  ├── wild_bootstrap(): Parametric bootstrap (impose H₀)         │
│  ├── parLapply(): Parallel processing                           │
│  └── Output: Power table cho Delta = 0, 1, 3                    │
├─────────────────────────────────────────────────────────────────┤
│  PART B: Quantile Transfer Entropy                              │
│  ├── calculate_qte(): Tính QTE cho mỗi (tau, k)                 │
│  ├── stationary_bootstrap(): Bootstrap cho time series          │
│  ├── bootstrap_qte_test(): Test significance                    │
│  └── Output: QTE estimates với p-values                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## So Sánh 4 Tests

| Tiêu chí | Sup-GQ | Fisher's G | Breusch-Pagan | White |
|----------|--------|------------|---------------|-------|
| **Loại** | Bootstrap | Bootstrap | Asymptotic | Asymptotic |
| **H₁ design** | Structural break | Structural break | General | General |
| **Speed** | Chậm | Chậm | Nhanh | Nhanh |
| **Power với break** | Cao | Cao | Trung bình | Trung bình |
| **Power general het.** | Trung bình | Trung bình | Tốt | Tốt |

---

*Cập nhật: 2026-01-05*
