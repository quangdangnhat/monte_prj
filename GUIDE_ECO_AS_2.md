# HƯỚNG DẪN CHI TIẾT: Phân Tích Hồi Quy Lượng Tử (Quantile Regression)

## Giới Thiệu

File `ECO_AS_2.R` thực hiện phân tích **hồi quy lượng tử** để nghiên cứu mối quan hệ giữa:
- **EPU (Economic Policy Uncertainty)**: Chỉ số bất ổn chính sách kinh tế
- **S&P500 Returns**: Lợi suất của chỉ số chứng khoán S&P500

**Mục tiêu**: Xem EPU có ảnh hưởng đến lợi suất chứng khoán hay không, đặc biệt ở các trạng thái thị trường khác nhau (tăng mạnh, giảm mạnh, trung bình).

---

## PHẦN 1: CÀI ĐẶT VÀ IMPORT THƯ VIỆN

### 1.1 Cài Đặt Packages (Dòng 1-28)

```r
install.packages("quantmod")      # Lấy dữ liệu tài chính từ Yahoo/FRED
install.packages("quantreg")      # Hồi quy lượng tử
install.packages("tseries")       # Kiểm định tính dừng (ADF test)
install.packages("boot")          # Bootstrap
install.packages("dplyr")         # Xử lý dữ liệu
```

**Giải thích cho newbie:**
- `quantmod`: Giống như "tải dữ liệu chứng khoán tự động"
- `quantreg`: Thư viện chính để chạy hồi quy lượng tử
- `tseries`: Kiểm tra xem dữ liệu có "ổn định" không
- `boot`: Thực hiện Bootstrap (kiểm tra độ tin cậy)
- `dplyr`: Công cụ làm sạch và biến đổi dữ liệu

### 1.2 Load Thư Viện (Dòng 30-44)

```r
library(quantreg)      # Hồi quy lượng tử
library(tseries)       # Kiểm định tính dừng
library(dplyr)         # Xử lý dữ liệu
library(zoo)           # Chuỗi thời gian
```

---

## PHẦN 2: LẤY VÀ XỬ LÝ DỮ LIỆU (BƯỚC 1)

### 2.1 Tải Dữ Liệu S&P500 (Dòng 56-62)

```r
getSymbols("^GSPC", from = "2015-01-01", to = "2025-12-31")
```

**Giải thích:**
- `^GSPC`: Mã chứng khoán của S&P500
- Lấy dữ liệu từ 2015 đến 2025
- Dữ liệu bao gồm: Open, High, Low, Close, Volume

### 2.2 Tính Log Returns (Dòng 70-74)

```r
returns <- diff(log(data$Close))
returns <- na.omit(returns)
```

**Tại sao dùng Log Returns?**

```
                    Giá ngày 1    Giá ngày 2
Công thức đơn giản: (P2 - P1)/P1 = (110 - 100)/100 = 10%
Log returns:        ln(P2/P1) = ln(110/100) = 9.53%
```

**Ưu điểm Log Returns:**
1. **Tính cộng được**: Log returns có thể cộng qua thời gian
2. **Phân phối chuẩn hơn**: Gần với phân phối chuẩn (Gaussian)
3. **Đối xứng**: Tăng 10% và giảm 10% có giá trị tuyệt đối tương đương

---

## PHẦN 3: KIỂM TRA TÍNH DỪNG (BƯỚC 2)

### 3.1 ADF Test (Dòng 80-88)

```r
adf.test(returns)      # Test cho Returns
adf.test(data$EPU)     # Test cho EPU
```

**Tính dừng (Stationarity) là gì?**

```
Chuỗi DỪNG:                    Chuỗi KHÔNG DỪNG:
    ___    ___                       /
   /   \  /   \                     /
  /     \/     \                   /
 /               \               /
(dao động quanh trung bình)    (có xu hướng tăng/giảm)
```

**Đọc kết quả ADF Test:**
- **H0 (Giả thuyết không)**: Chuỗi KHÔNG dừng (có unit root)
- **p-value < 0.05**: Bác bỏ H0 → Chuỗi DỪNG ✓
- **p-value ≥ 0.05**: Không bác bỏ H0 → Chuỗi KHÔNG dừng ✗

**Kết quả trong code:**
```r
# Returns: p-value = 0.01 → DỪNG ✓
# EPU:     p-value = 0.2755 → KHÔNG DỪNG ✗
```

### 3.2 Lấy Sai Phân Để Làm Dừng (Dòng 90-96)

```r
diff_EPU <- diff(data$EPU)
data$EPU_diff <- c(NA, diff(data$EPU))
```

**Sai phân là gì?**
```
Dữ liệu gốc:  100, 105, 108, 115
Sai phân:     NA,   5,   3,   7  (chênh lệch giữa các ngày)
```

---

## PHẦN 4: HỒI QUY LƯỢNG TỬ (BƯỚC 3)

### 4.1 Tại Sao Dùng Hồi Quy Lượng Tử?

**So sánh OLS và Quantile Regression:**

```
OLS (Hồi quy thường):
                    •
               •    •   •
          •  •  •  -------- đường hồi quy (trung bình)
        •     •       •
      •         •

Chỉ cho 1 đường thẳng ở TRUNG BÌNH

Quantile Regression:
                    •
               •    •   • -------- τ = 0.90 (phân vị 90%)
          •  •  •  -------- τ = 0.50 (trung vị)
        •     •       • -------- τ = 0.25 (phân vị 25%)
      •         • -------- τ = 0.10 (phân vị 10%)

Cho NHIỀU đường thẳng ở các phân vị khác nhau!
```

### 4.2 Ý Nghĩa Các Phân Vị (τ - tau)

| τ (tau) | Ý nghĩa thị trường | Mô tả |
|---------|-------------------|-------|
| 0.10 | Giảm mạnh | Thị trường sụt giảm nghiêm trọng (10% ngày tệ nhất) |
| 0.25 | Giảm nhẹ | Thị trường đang yếu |
| 0.50 | Bình thường | Trung vị - trạng thái trung bình |
| 0.75 | Tăng nhẹ | Thị trường đang khá |
| 0.90 | Tăng mạnh | Thị trường bùng nổ (10% ngày tốt nhất) |

### 4.3 Tạo Biến Độ Trễ (Dòng 102-104)

```r
data_clean$Return_lag <- dplyr::lag(data_clean$Return, 1)
data_clean$dEPU_lag <- dplyr::lag(data_clean$EPU, 1)
```

**Tại sao cần độ trễ (lag)?**

```
Ngày        EPU    Return
Jan-2020    120    ← EPU tháng 1 có thể ảnh hưởng
Feb-2020    130    -2%   ← Return tháng 2
Mar-2020    150    -5%

Mô hình: Return_t = f(EPU_{t-1}, Return_{t-1})
```

Chúng ta dùng **EPU của tháng trước** để dự đoán **Return của tháng này** vì:
- EPU công bố với độ trễ
- Nhà đầu tư cần thời gian phản ứng

### 4.4 Mô Hình Giới Hạn (Restricted) - Dòng 117-121

```r
model_restricted <- rq(Return ~ Return_lag,
                       data = data_clean,
                       tau = 0.25)
```

**Giải thích:**
```
Return_t = α + β₁ × Return_{t-1} + ε

Trong đó:
- α (Intercept): Hệ số chặn
- β₁: Hệ số của Return tháng trước
- KHÔNG có EPU trong mô hình này
```

### 4.5 Mô Hình Không Giới Hạn (Unrestricted) - Dòng 129-133

```r
model_unrestricted <- rq(Return ~ Return_lag + dEPU_lag,
                         data = data_clean,
                         tau = 0.25)
```

**Giải thích:**
```
Return_t = α + β₁ × Return_{t-1} + β₂ × EPU_{t-1} + ε

Trong đó:
- α (Intercept): Hệ số chặn
- β₁: Hệ số của Return tháng trước
- β₂: Hệ số của EPU → ĐÂY LÀ HỆ SỐ QUAN TRỌNG!
```

---

## PHẦN 5: TÍNH QTE - QUANTILE TREATMENT EFFECT

### 5.1 QTE Là Gì?

**QTE đo lường**: EPU có "truyền thông tin" đến Returns hay không?

```
QTE = log(Mất mát mô hình giới hạn) - log(Mất mát mô hình không giới hạn)
    = log(ρ₁) - log(ρ₂)
```

**Ý nghĩa:**
- **QTE > 0**: Thêm EPU vào mô hình làm GIẢM mất mát → EPU CÓ ích
- **QTE ≈ 0**: EPU không thêm thông tin gì
- **QTE < 0**: (Hiếm khi xảy ra) Thêm EPU làm mô hình tệ hơn

### 5.2 Hàm Mất Mát Phân Vị (Check Function)

```r
# Phần dư
u1 <- residuals(model_restricted)
u2 <- residuals(model_unrestricted)

# Tách phần dư dương và âm
neg_u1 <- u1[u1 < 0]
pos_u1 <- u1[u1 >= 0]

# Hàm mất mát phân vị
sum_pos1 <- sum(pos_u1 * tau)           # u ≥ 0: τ × u
sum_n1 <- sum(neg_u1 * (tau - 1))       # u < 0: (τ-1) × u
```

**Minh họa hàm mất mát:**

```
Mất mát ρ_τ(u)
    ↑
    |     /
    |    /  ← slope = τ (với u > 0)
    |   /
    |  /
----+------→ u (phần dư)
   /|
  / |  ← slope = τ-1 (với u < 0)
 /  |
```

Với τ = 0.25:
- Dự đoán CAO hơn thực tế (u < 0): Phạt 75%
- Dự đoán THẤP hơn thực tế (u > 0): Phạt 25%

### 5.3 Công Thức QTE Cuối Cùng (Dòng 169)

```r
QTE_exact <- log(sum_p1) - log(sum_p2)
```

---

## PHẦN 6: VÒNG LẶP CHO TẤT CẢ PHÂN VỊ (Dòng 183-284)

### 6.1 Cấu Trúc Vòng Lặp

```r
tau_vec <- c(0.10, 0.25, 0.50, 0.75, 0.90)

for (i in 1:length(tau_vec)) {
  tau <- tau_vec[i]

  # Bước 1: Chạy mô hình giới hạn
  model_restricted <- rq(Return ~ Return_lag, data = data_model, tau = tau)

  # Bước 2: Chạy mô hình không giới hạn
  model_unrestricted <- rq(Return ~ Return_lag + dEPU_lag, data = data_model, tau = tau)

  # Bước 3: Tính QTE
  # ... (như đã giải thích ở trên)

  # Bước 4: Lưu kết quả
  results <- rbind(results, data.frame(...))
}
```

### 6.2 Khoảng Tin Cậy (Dòng 239-248)

```r
# Khoảng tin cậy 90%
ci_lower_90 <- theta - 1.645 * std_error
ci_upper_90 <- theta + 1.645 * std_error

# Khoảng tin cậy 95%
ci_lower_95 <- theta - 1.96 * std_error
ci_upper_95 <- theta + 1.96 * std_error
```

**Cách đọc khoảng tin cậy:**
```
Nếu CI 95% là [0.001, 0.003]:
→ Cả hai giá trị đều > 0
→ Hệ số có ý nghĩa thống kê (khác 0)

Nếu CI 95% là [-0.002, 0.003]:
→ Bao gồm cả giá trị âm và dương
→ Có thể hệ số thực sự = 0
→ KHÔNG có ý nghĩa thống kê
```

---

## PHẦN 7: BOOTSTRAP (BƯỚC 7 - Dòng 330-500)

### 7.1 Bootstrap Là Gì?

**Vấn đề**: Chúng ta chỉ có 1 bộ dữ liệu, làm sao biết kết quả có đáng tin?

**Giải pháp Bootstrap**: Tạo ra nhiều bộ dữ liệu "giả" bằng cách lấy mẫu có hoàn lại

```
Dữ liệu gốc:     [A, B, C, D, E]

Bootstrap 1:     [A, A, C, D, E]  ← lấy ngẫu nhiên có hoàn lại
Bootstrap 2:     [B, C, C, D, D]
Bootstrap 3:     [A, B, D, D, E]
...
Bootstrap 499:   [A, C, C, E, E]
```

### 7.2 Tham Số Bootstrap (Dòng 336-342)

```r
B <- 499              # Số lần bootstrap
tau_vec <- c(0.10, 0.50, 0.90)  # 3 phân vị cần kiểm định
n <- nrow(data_model)
l <- round(n^(1/3))   # Kích thước khối
```

**Tại sao B = 499?**
- Số lẻ để tính p-value chính xác hơn
- Đủ lớn để có kết quả ổn định
- Không quá lớn để tiết kiệm thời gian

### 7.3 Hàm Tính QTE Nhanh (Dòng 358-393)

```r
tinh_QTE_nhanh <- function(data, tau) {
  model1 <- rq(Return ~ Return_lag, data = data, tau = tau)
  model2 <- rq(Return ~ Return_lag + dEPU_lag, data = data, tau = tau)

  # Tính mất mát và QTE như trên...

  return(QTE_exact)
}
```

### 7.4 Hàm Bootstrap (Dòng 396-409)

```r
bootstrap_don_gian <- function() {
  # Tạo HAI chỉ mục ĐỘC LẬP
  boot_y <- sample(1:n, n, replace = TRUE)  # Cho Y (Return)
  boot_x <- sample(1:n, n, replace = TRUE)  # Cho X (EPU) - ĐỘC LẬP!

  data_boot <- data.frame(
    Return = data_model$Return[boot_y],
    Return_lag = data_model$Return_lag[boot_y],
    dEPU_lag = data_model$dEPU_lag[boot_x]  # Độc lập!
  )

  return(data_boot)
}
```

**Tại sao lấy mẫu ĐỘC LẬP cho X và Y?**

Đây là kỹ thuật kiểm định **Granger causality** bằng bootstrap:
- H0: EPU KHÔNG gây ra Returns (không có quan hệ)
- Nếu H0 đúng → X và Y độc lập
- Bằng cách xáo trộn độc lập, ta mô phỏng H0

### 7.5 Vòng Lặp Bootstrap Chính (Dòng 423-440)

```r
set.seed(123)  # Để kết quả có thể tái lập

for (b in 1:B) {
  # Tạo dữ liệu bootstrap
  data_boot <- bootstrap_don_gian()

  # Tính QTE cho 3 phân vị
  for (j in 1:3) {
    QTE_boot[b, j] <- tinh_QTE_nhanh(data_boot, tau_vec[j])
  }

  # Hiển thị tiến độ
  if (b %% 100 == 0) cat("Đã xong", b, "lần\n")
}
```

### 7.6 Tính P-value (Dòng 455-476)

```r
for (i in 1:3) {
  # Đếm số lần QTE_boot >= QTE_thuc
  count <- sum(QTE_boot[, i] >= QTE_thuc[i], na.rm = TRUE)

  # Tính p-value một phía
  p_val <- (1 + count) / (B + 1)

  ket_qua$p_value[i] <- p_val
  ket_qua$Y_nghia[i] <- ifelse(p_val < 0.05, "CÓ", "KHÔNG")
}
```

**Giải thích p-value:**

```
Nếu H0 đúng (không có quan hệ):
→ QTE thực tế nên nằm trong phân phối bootstrap
→ Có nhiều QTE_boot >= QTE_thực

Nếu H0 sai (CÓ quan hệ):
→ QTE thực tế lớn hơn hầu hết QTE_boot
→ Ít QTE_boot >= QTE_thực
→ p-value nhỏ
```

**Công thức:**
```
p-value = (1 + số lần QTE_boot >= QTE_thực) / (B + 1)
        = (1 + count) / 500

Nếu count = 2:  p-value = 3/500 = 0.006 → Bác bỏ H0
Nếu count = 290: p-value = 291/500 = 0.582 → Không bác bỏ H0
```

---

## PHẦN 8: ĐỌC KẾT QUẢ

### 8.1 Bảng Kết Quả Cuối Cùng (Dòng 495-499)

```
Phan_vi  Trang_thai      QTE_thuc   p_value  Y_nghia
τ=0.1    Lỗ thấp         0.0039     0.582    KHÔNG
τ=0.5    Trung bình      0.0436     0.004    CÓ
τ=0.9    Lợi nhuận cao   0.1319     0.002    CÓ
```

### 8.2 Cách Diễn Giải

**τ = 0.10 (Thị trường giảm mạnh):**
- p-value = 0.582 > 0.05
- **Kết luận**: KHÔNG có bằng chứng EPU ảnh hưởng đến returns
- **Ý nghĩa**: Khi thị trường đang lỗ nặng, EPU không phải yếu tố quan trọng

**τ = 0.50 (Thị trường bình thường):**
- p-value = 0.004 < 0.05
- **Kết luận**: CÓ bằng chứng EPU ảnh hưởng đến returns
- **Ý nghĩa**: Trong điều kiện bình thường, EPU truyền thông tin đến thị trường

**τ = 0.90 (Thị trường tăng mạnh):**
- p-value = 0.002 < 0.05
- **Kết luận**: CÓ bằng chứng EPU ảnh hưởng đến returns
- **Ý nghĩa**: Khi thị trường đang tăng mạnh, EPU có ảnh hưởng

---

## PHẦN 9: TÓM TẮT QUY TRÌNH

```
┌─────────────────────────────────────────────────────────────┐
│  BƯỚC 1: LẤY DỮ LIỆU                                        │
│  • S&P500 từ Yahoo Finance                                  │
│  • EPU từ website/file CSV                                  │
│  • Tính log returns                                         │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  BƯỚC 2: KIỂM TRA TÍNH DỪNG                                 │
│  • ADF Test                                                 │
│  • Lấy sai phân nếu cần                                     │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  BƯỚC 3: HỒI QUY LƯỢNG TỬ                                   │
│  • Mô hình giới hạn: Return ~ Return_lag                    │
│  • Mô hình không giới hạn: Return ~ Return_lag + EPU_lag    │
│  • Chạy cho các τ = 0.10, 0.25, 0.50, 0.75, 0.90            │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  BƯỚC 4: TÍNH QTE                                           │
│  • QTE = log(ρ₁) - log(ρ₂)                                  │
│  • Đo lường ảnh hưởng của EPU                               │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  BƯỚC 5: BOOTSTRAP                                          │
│  • Tạo 499 bộ dữ liệu giả (H0: không có quan hệ)            │
│  • Tính QTE cho mỗi bộ                                      │
│  • Tính p-value                                             │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  BƯỚC 6: KẾT LUẬN                                           │
│  • p-value < 0.05: BÁC BỎ H0 → CÓ ảnh hưởng                 │
│  • p-value ≥ 0.05: Không bác bỏ H0 → KHÔNG có bằng chứng    │
└─────────────────────────────────────────────────────────────┘
```

---

## PHẦN 10: CÁC LỖI THƯỜNG GẶP VÀ CÁCH KHẮC PHỤC

### 10.1 Lỗi "NA values"

```r
# Vấn đề: Có NA trong dữ liệu
data_clean <- na.omit(data)  # Xóa các hàng có NA
```

### 10.2 Lỗi Hội Tụ Trong Bootstrap

```r
# Dùng tryCatch để bắt lỗi
val <- tryCatch({
  tinh_QTE_nhanh(data_boot, tau_vec[j])
}, error = function(e) return(NA))
```

### 10.3 Kết Quả Không Ổn Định

```r
# Đặt seed để kết quả có thể tái lập
set.seed(123)
```

---

## THUẬT NGỮ QUAN TRỌNG

| Thuật ngữ | Tiếng Anh | Ý nghĩa |
|-----------|-----------|---------|
| Phân vị | Quantile (τ) | Điểm chia dữ liệu thành các phần |
| Tính dừng | Stationarity | Dữ liệu dao động quanh giá trị trung bình |
| Hồi quy lượng tử | Quantile Regression | Hồi quy ở các phân vị khác nhau |
| QTE | Quantile Treatment Effect | Đo lường ảnh hưởng của biến X lên Y |
| Bootstrap | Bootstrap | Phương pháp lấy mẫu lại |
| p-value | p-value | Xác suất để kết quả xảy ra ngẫu nhiên |
| EPU | Economic Policy Uncertainty | Chỉ số bất ổn chính sách kinh tế |
| Log returns | Logarithmic Returns | Lợi suất logarit |

---

## KẾT LUẬN NGHIÊN CỨU

Từ kết quả phân tích:

1. **Khi thị trường giảm mạnh (τ=0.10)**: EPU không ảnh hưởng đáng kể
   - Nhà đầu tư đã hoảng loạn, không quan tâm đến EPU

2. **Khi thị trường bình thường (τ=0.50)**: EPU CÓ ảnh hưởng
   - Nhà đầu tư lý trí hơn, phản ứng với thông tin

3. **Khi thị trường tăng mạnh (τ=0.90)**: EPU CÓ ảnh hưởng mạnh
   - Có thể EPU cao → nhà đầu tư tìm tài sản an toàn → đẩy giá lên

---

*Hướng dẫn này được tạo để giúp người mới học hiểu rõ về phân tích hồi quy lượng tử và bootstrap trong tài chính.*
