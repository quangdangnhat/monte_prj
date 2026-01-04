# GIẢI THÍCH CHI TIẾT: Coursework_Complete.R

## Tổng Quan

File `Coursework_Complete.R` là bài tập môn **Financial Econometrics** năm học 2025-2026, bao gồm 2 phần chính:

| Phần | Nội dung | Mục đích |
|------|----------|----------|
| **Part A** | Monte Carlo Power Comparison | So sánh sức mạnh (power) của các bài kiểm tra heteroscedasticity |
| **Part B** | Quantile Transfer Entropy (QTE) | Phân tích quan hệ nhân quả giữa các chuỗi thời gian |

---

## [0] CONFIGURATION (Dòng 17-35)

### Các biến cấu hình:

```r
RUN_PART_A <- TRUE    # Chạy Part A? (TRUE/FALSE)
RUN_PART_B <- TRUE    # Chạy Part B? (TRUE/FALSE)
QUICK_TEST <- TRUE    # Chế độ test nhanh?
```

### Chế độ chạy:

| Chế độ | R (MC reps) | B (Bootstrap) | Thời gian ước tính |
|--------|-------------|---------------|-------------------|
| `QUICK_TEST = TRUE` | 50 | 99 | ~5 phút |
| `QUICK_TEST = FALSE` | 1000 | 499 | ~1 giờ |

**Giải thích:**
- `R_SIMS`: Số lần lặp Monte Carlo (càng lớn → kết quả càng ổn định)
- `B_BOOT`: Số lần bootstrap (càng lớn → p-value càng chính xác)

---

## [1] SETUP & PACKAGES (Dòng 37-76)

### Packages được sử dụng:

| Package | Mục đích |
|---------|----------|
| `lmtest` | Breusch-Pagan test cho heteroscedasticity |
| `quantreg` | Quantile regression (hồi quy phân vị) |
| `tseries` | Phân tích chuỗi thời gian |
| `zoo`, `xts` | Xử lý dữ liệu chuỗi thời gian |
| `dplyr` | Thao tác dữ liệu |
| `quantmod` | Tải dữ liệu tài chính |

---

# PART A: MONTE CARLO POWER COMPARISON

## Mục tiêu

So sánh **power** (sức mạnh) của 4 bài kiểm tra heteroscedasticity:

1. **Sup-GQ** (Supremum Goldfeld-Quandt)
2. **Fisher's G** (Combined G statistic)
3. **Breusch-Pagan**
4. **White Test**

## [A.1] PARAMETERS (Dòng 91-108)

```r
N <- 100              # Kích thước mẫu
R <- R_SIMS           # Số lần Monte Carlo
B <- B_BOOT           # Số lần Bootstrap
Delta <- c(0, 1, 3)   # Mức độ heteroscedasticity
trim <- 0.15          # Trimming 15% hai đầu
alpha <- 0.05         # Mức ý nghĩa 5%
```

### Ý nghĩa của Delta:

| Delta | Phương sai sau breakpoint | Ý nghĩa |
|-------|--------------------------|---------|
| 0 | σ² = 1 (không đổi) | **Size check**: H₀ đúng, tỷ lệ reject phải ≈ 0.05 |
| 1 | σ² = 2 (gấp đôi) | Heteroscedasticity nhẹ |
| 3 | σ² = 4 (gấp 4) | Heteroscedasticity mạnh |

---

## [A.2] FUNCTIONS (Dòng 110-211)

### 1. `generate_data(n, delta)` - Tạo dữ liệu mô phỏng

**Mô hình DGP (Data Generating Process):**

```
y_t = 1 + x_t + ε_t

Trong đó:
- x_t ~ N(0, 1)
- ε_t ~ N(0, σ_t²)
- σ_t² = 1           nếu t ≤ n/2 (nửa đầu)
- σ_t² = 1 + delta   nếu t > n/2 (nửa sau)
```

**Giải thích:**
- Đây là mô hình hồi quy đơn giản với **structural break** ở giữa mẫu
- Phương sai thay đổi đột ngột tại điểm τ = 0.5

```r
generate_data <- function(n, delta) {
  x <- rnorm(n, mean = 0, sd = 1)                    # Biến độc lập
  var_structure <- c(rep(1, n/2), rep(1 + delta, n/2))  # Cấu trúc phương sai
  epsilon <- rnorm(n, mean = 0, sd = sqrt(var_structure)) # Sai số
  y <- 1 + x + epsilon                               # Biến phụ thuộc
  return(data.frame(y = y, x = x))
}
```

---

### 2. `calc_statistics(y, x, trim)` - Tính Sup-GQ và Fisher's G

**Thuật toán:**

```
Bước 1: Xác định lưới breakpoint τ ∈ [15%, 85%]
Bước 2: Với mỗi τ:
    a) Chia mẫu thành 2 phần: [1, τ] và [τ+1, n]
    b) Chạy OLS riêng biệt cho mỗi phần
    c) Tính RSS₁, RSS₂ (Residual Sum of Squares)
    d) Tính F-statistic: F = (RSS₂/df₂) / (RSS₁/df₁)
    e) Tính p-value của F
Bước 3:
    - Sup-GQ = max(F) qua tất cả τ
    - G = -2 × Σlog(p_values) (Fisher's method)
```

**Công thức chi tiết:**

$$F_{GQ}(\tau) = \frac{RSS_2 / (n - \tau - 2)}{RSS_1 / (\tau - 2)}$$

$$G = -2 \sum_{i=1}^{M} \log(p_i)$$

Dưới H₀, G ~ χ²(2M) với M là số điểm trong lưới.

---

### 3. `breusch_pagan_test(y, x)` - Breusch-Pagan Test

**Nguyên lý:**
- Kiểm tra xem phương sai của sai số có phụ thuộc vào biến độc lập không
- Sử dụng package `lmtest`

```r
breusch_pagan_test <- function(y, x) {
  model <- lm(y ~ x)
  bp_test <- bptest(model)    # Từ package lmtest
  return(bp_test$p.value)
}
```

---

### 4. `white_test(y, x)` - White Test

**Nguyên lý:**
- Hồi quy ε̂² lên x và x²
- Kiểm tra R² của hồi quy phụ

```r
white_test <- function(y, x) {
  model <- lm(y ~ x)
  resid_sq <- resid(model)^2        # ε̂²
  x_sq <- x^2
  aux_model <- lm(resid_sq ~ x + x_sq)  # Hồi quy phụ
  r_squared <- summary(aux_model)$r.squared
  white_stat <- n * r_squared       # nR² ~ χ²(2)
  p_value <- pchisq(white_stat, df = 2, lower.tail = FALSE)
  return(p_value)
}
```

---

### 5. `wild_bootstrap(y_obs, x_obs, B, stats_obs)` - Wild Bootstrap

**Tại sao dùng Wild Bootstrap?**
- Phân phối của Sup-GQ và G không chuẩn → không có critical value lý thuyết
- Wild bootstrap bảo toàn cấu trúc heteroscedasticity

**Thuật toán:**

```
Bước 1: Fit model dưới H₀: y = β₀ + β₁x + ε
Bước 2: Lấy ŷ (fitted) và ê (residuals)
Bước 3: Lặp B lần:
    a) Sinh v ~ Rademacher: v ∈ {-1, +1} với xác suất bằng nhau
    b) Tạo y* = ŷ + ê × v (bootstrap sample)
    c) Tính T* (statistic trên bootstrap sample)
Bước 4: p-value = Tỷ lệ (T* ≥ T_obs)
```

```r
wild_bootstrap <- function(y_obs, x_obs, B, stats_obs) {
  null_model <- lm(y_obs ~ x_obs)
  y_hat <- fitted(null_model)      # ŷ
  e_hat <- resid(null_model)       # ê

  for (b in 1:B) {
    v <- sample(c(-1, 1), n, replace = TRUE)  # Rademacher
    y_star <- y_hat + e_hat * v               # Bootstrap sample
    stats_star <- calc_statistics(y_star, x_obs)
    # Lưu kết quả...
  }

  pval_sup <- mean(boot_sup_vals >= stats_obs$sup)
  pval_g <- mean(boot_g_vals >= stats_obs$g)
  return(list(pval_sup = pval_sup, pval_g = pval_g))
}
```

---

## [A.3] MONTE CARLO SIMULATION (Dòng 214-273)

**Quy trình:**

```
Với mỗi delta ∈ {0, 1, 3}:
    Lặp R lần (Monte Carlo):
        1. Sinh dữ liệu với delta
        2. Tính Sup-GQ và G observed
        3. Wild bootstrap → p-values cho Sup-GQ, G
        4. Tính p-values cho BP và White
        5. Đếm số lần reject H₀ (p < 0.05)
    Power = (số reject) / R
```

---

## [A.4] OUTPUT - Bảng kết quả Part A

```
Table 1: Empirical Size (δ=0) and Power (δ=1,3) at 5% Significance Level
─────────────────────────────────────────────────────────────────────────
Delta      | Sup-GQ     | G (Fisher) | Breusch-Pagan | White
─────────────────────────────────────────────────────────────────────────
0 (Size)   | 0.050      | 0.048      | 0.052         | 0.048
1          | 0.320      | 0.380      | 0.250         | 0.220
3          | 0.850      | 0.920      | 0.650         | 0.580
─────────────────────────────────────────────────────────────────────────
```

**Cách đọc kết quả:**
- **Delta = 0 (Size):** Phải ≈ 0.05 (nominal level). Nếu > 0.05 → test bị **over-sized**
- **Delta = 1, 3 (Power):** Càng cao càng tốt. Test nào có power cao hơn → phát hiện heteroscedasticity tốt hơn

---

# PART B: QUANTILE TRANSFER ENTROPY (QTE)

## Mục tiêu

Kiểm tra **quan hệ nhân quả Granger theo phân vị** từ X (uncertainty) đến Y (returns):

> "X có chứa thông tin hữu ích để dự báo phân vị τ của Y không?"

## [B.1] PARAMETERS (Dòng 325-338)

```r
tau_values <- c(0.10, 0.50, 0.90)  # Các phân vị
k_values <- c(1, 2, 3, 5)          # Các độ trễ
```

**Ý nghĩa các phân vị:**

| τ | Ý nghĩa | Ứng dụng |
|---|---------|----------|
| 0.10 | Đuôi trái (10% thấp nhất) | Extreme losses, risk management |
| 0.50 | Trung vị | Xu hướng trung bình |
| 0.90 | Đuôi phải (10% cao nhất) | Extreme gains |

---

## [B.2] DATA PREPARATION (Dòng 340-374)

### Hàm sinh dữ liệu mô phỏng:

```r
generate_qte_data <- function(n = 2000, causality = 0.15) {
  # X: AR(1) process (mô phỏng uncertainty)
  x[t] = 0.3 × x[t-1] + ε_t

  # Y: Phụ thuộc vào X với causality
  y[t] = 0.2 × y[t-1] + causality × x[t-1] + η_t

  # η_t có phương sai phụ thuộc vào |x[t-1]| (GARCH-like)
}
```

**Lưu ý:** Trong thực tế, thay thế bằng dữ liệu thật (EPU và S&P500 returns).

---

## [B.3] QTE FUNCTIONS (Dòng 376-456)

### 1. `check_function(u, tau)` - Hàm kiểm tra phân vị

$$\rho_\tau(u) = u \cdot (\tau - \mathbf{1}_{u < 0})$$

Đây là hàm mất mát của quantile regression:
- Nếu u ≥ 0: ρ(u) = τ × u
- Nếu u < 0: ρ(u) = (τ - 1) × u

---

### 2. `calculate_qte(Y, X, tau, k)` - Tính QTE

**Công thức QTE:**

$$QTE_{\tau,k} = \log\left(\sum \rho_\tau(\hat{u}_R)\right) - \log\left(\sum \rho_\tau(\hat{u}_U)\right)$$

Trong đó:
- **Restricted model (R):** $Q_\tau(Y_t | Y_{t-k})$ - chỉ dùng lag của Y
- **Unrestricted model (U):** $Q_\tau(Y_t | Y_{t-k}, X_{t-k})$ - dùng cả lag của X

**Diễn giải:**
- QTE > 0: X chứa thông tin hữu ích cho dự báo Y
- QTE ≈ 0: X không có quan hệ nhân quả với Y

```r
calculate_qte <- function(Y, X, tau, k) {
  # Chuẩn bị dữ liệu với lag k
  Y_t <- Y[(k+1):n]
  Y_lag <- Y[1:(n-k)]
  X_lag <- X[1:(n-k)]

  # Restricted: Q(Y_t | Y_lag)
  model_r <- rq(Y_t ~ Y_lag, tau = tau)
  loss_r <- sum(check_function(residuals(model_r), tau))

  # Unrestricted: Q(Y_t | Y_lag, X_lag)
  model_u <- rq(Y_t ~ Y_lag + X_lag, tau = tau)
  loss_u <- sum(check_function(residuals(model_u), tau))

  qte <- log(loss_r) - log(loss_u)
  return(qte)
}
```

---

### 3. `stationary_bootstrap_indices(n, block_length)` - Stationary Bootstrap

**Tại sao dùng Stationary Bootstrap?**
- Dữ liệu chuỗi thời gian có autocorrelation
- Block bootstrap bảo toàn cấu trúc phụ thuộc

**Thuật toán:**
- Block length ngẫu nhiên theo phân phối geometric
- Trung bình block length = n^(1/3)

---

### 4. `bootstrap_qte_test(Y, X, tau, k, B)` - Bootstrap Test

**Quy trình:**

```
1. Tính QTE_obs từ dữ liệu gốc
2. Lặp B lần:
   a) Resample X và Y độc lập (impose H₀: X không gây ra Y)
   b) Tính QTE* từ bootstrap sample
3. P-value = Tỷ lệ (QTE* ≥ QTE_obs)
```

---

## [B.4-B.5] OUTPUT - Bảng kết quả Part B

```
Table 2: QTE Estimates (*** = significant at 5%)
─────────────────────────────────────────────────────────────────────────
Lag (k)  | tau = 0.10      | tau = 0.50      | tau = 0.90
─────────────────────────────────────────────────────────────────────────
1        | 0.0234***       | 0.0089          | 0.0312***
2        | 0.0198**        | 0.0045          | 0.0256**
3        | 0.0156          | 0.0032          | 0.0189
5        | 0.0098          | 0.0021          | 0.0134
─────────────────────────────────────────────────────────────────────────
```

**Cách đọc kết quả:**
- `***` = p-value < 0.05 → có quan hệ nhân quả có ý nghĩa thống kê
- QTE cao ở τ = 0.10 và 0.90 → X ảnh hưởng đến đuôi phân phối của Y
- QTE giảm khi k tăng → hiệu ứng giảm dần theo thời gian

---

# SƠ ĐỒ LUỒNG CHƯƠNG TRÌNH

```
┌─────────────────────────────────────────────────────────────┐
│                    START                                    │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  [0] CONFIGURATION                                          │
│  - Chọn QUICK_TEST = TRUE/FALSE                             │
│  - Chọn RUN_PART_A, RUN_PART_B                              │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  [1] SETUP & PACKAGES                                       │
│  - Load: lmtest, quantreg, tseries, zoo, xts, dplyr         │
└─────────────────────────────────────────────────────────────┘
                           │
           ┌───────────────┴───────────────┐
           ▼                               ▼
┌──────────────────────┐       ┌──────────────────────┐
│      PART A          │       │      PART B          │
│   (if RUN_PART_A)    │       │   (if RUN_PART_B)    │
└──────────────────────┘       └──────────────────────┘
           │                               │
           ▼                               ▼
┌──────────────────────┐       ┌──────────────────────┐
│ For delta in {0,1,3} │       │ For tau in {.1,.5,.9}│
│   For r in 1:R       │       │   For k in {1,2,3,5} │
│     - Generate data  │       │     - Calculate QTE  │
│     - Wild bootstrap │       │     - Bootstrap test │
│     - Count rejects  │       │     - Store results  │
└──────────────────────┘       └──────────────────────┘
           │                               │
           ▼                               ▼
┌──────────────────────┐       ┌──────────────────────┐
│    Output Table 1    │       │  Output Tables 2, 3  │
│  Size & Power Table  │       │  QTE & P-values      │
└──────────────────────┘       └──────────────────────┘
           │                               │
           └───────────────┬───────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                    FINAL SUMMARY                            │
│  - Print kết quả Part A                                     │
│  - Print significant QTE results Part B                     │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                         END                                 │
└─────────────────────────────────────────────────────────────┘
```

---

# BẢNG TỔNG HỢP CÁC KHÁI NIỆM

| Thuật ngữ | Định nghĩa | Công thức/Ý nghĩa |
|-----------|------------|-------------------|
| **Size** | Tỷ lệ reject khi H₀ đúng | Phải ≈ α (nominal level) |
| **Power** | Tỷ lệ reject khi H₁ đúng | Càng cao càng tốt |
| **Sup-GQ** | Maximum F-test across breakpoints | max(F(τ)) |
| **Fisher's G** | Combined p-values | -2Σlog(p) |
| **Wild Bootstrap** | Resampling với Rademacher | y* = ŷ + ê×v |
| **QTE** | Quantile Transfer Entropy | log(Loss_R) - log(Loss_U) |
| **Quantile τ** | Phân vị thứ τ | P(Y ≤ Q_τ) = τ |

---

# HƯỚNG DẪN SỬ DỤNG

## Chạy test nhanh (5 phút):

```r
QUICK_TEST <- TRUE
RUN_PART_A <- TRUE
RUN_PART_B <- TRUE
source("Coursework_Complete.R")
```

## Chạy đầy đủ cho submission (~1 giờ):

```r
QUICK_TEST <- FALSE
RUN_PART_A <- TRUE
RUN_PART_B <- TRUE
source("Coursework_Complete.R")
```

## Chỉ chạy Part A:

```r
RUN_PART_A <- TRUE
RUN_PART_B <- FALSE
source("Coursework_Complete.R")
```

---

# LƯU Ý QUAN TRỌNG

1. **Seed:** `set.seed(123)` đảm bảo kết quả reproducible
2. **Part B Data:** Đang dùng dữ liệu mô phỏng → thay bằng dữ liệu thật (EPU, S&P500)
3. **Thời gian chạy:** Full simulation rất lâu → chạy trên máy mạnh hoặc cluster
4. **Interpretation:** Size check (δ=0) quan trọng để đánh giá test có reliable không
