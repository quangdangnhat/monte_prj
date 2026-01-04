install.packages("quantmod")
install.packages("dplyr")
install.packages("quantreg")        # Ước lượng mô hình hồi quy lượng tử
install.packages("quantregGrowth")  # Hỗ trợ thêm cho mô hình lượng tử

# 2. Chuỗi thời gian và kiểm định
install.packages("tseries")         # Kiểm định ADF, KPSS
install.packages("urca")           # Kiểm định tính dừng nâng cao (ur.df)
install.packages("forecast")       # Công cụ chuỗi thời gian, kiểm tra ACF/PACF
install.packages("vars")           # Phân tích VAR (nếu cần phân tích đa biến)
# 3. Bootstrap và tính toán
install.packages("boot")           # Hàm bootstrap cơ bản
install.packages("np")             # Ước lượng phi tham số (kernel density)
install.packages("sandwich")       # Ước lượng phương sai mạnh mẽ

# 4. Xử lý dữ liệu
install.packages("dplyr")          # Biến đổi và làm sạch dữ liệu
install.packages("tidyr")          # Định dạng lại dữ liệu
install.packages("zoo")           # Xử lý dữ liệu chuỗi thời gian
install.packages("xts")           # Đối tượng chuỗi thời gian nâng cao



# Cài đặt các package cần thiết
install.packages(c("rmarkdown", "knitr", "tinytex", "kableExtra", "ggplot2"))

# Cài đặt LaTeX (cần cho xuất file PDF)
tinytex::install_tinytex()
 # Cần cho ARCH-LM Test
library(tseries)
library(quantmod)
library(dplyr)
library(quantreg)      # Hồi quy lượng tử
library(quantregGrowth)# Mở rộng hồi quy lượng tử
library(tseries)       # Kiểm định tính dừng
library(urca)          # Kiểm định tính dừng nâng cao
library(forecast)      # Công cụ chuỗi thời gian
library(boot)          # Bootstrap
library(np)            # Nonparametric estimation
library(sandwich)      # Robust variance
library(dplyr)         # Data manipulation
library(tidyr)         # Data tidying
library(zoo)           # Time series objects
library(xts)           # Extended time series

# GIẢ ĐỊNH: returns là vector/chuỗi thời gian lợi suất của bạn

#-------------------------------------------------
# BƯỚC 1: LẤY DATA & XƯ LÝ
#-------------------------------------------------

# Data EPU : web site
# Lấy dữ liệu S&P500:
## FRED:
## YAHOO:
getSymbols("^GSPC", from = "2015-01-01", to = "2025-12-31")

sp500_df <- as.data.frame(GSPC) %>%
  mutate(Date = index(GSPC), Ticker = "S&P500") %>%
  rename(Open = GSPC.Open, High = GSPC.High, Low = GSPC.Low, 
         Close = GSPC.Close, Volume = GSPC.Volume, Adjusted = GSPC.Adjusted) %>%
  select(Date, Ticker, Open, High, Low, Close, Volume, Adjusted)
# XUẤT ĐẺ XƯ LÝ

write.csv(sp500_df, "D:/eco_ass.csv", row.names = TRUE)

# NHẬP LẠI EPU-SP500
data <- read.csv("C:/Users/ADMIN/Downloads/EPU_SP_COMBINE.csv")

# Tính log returns
returns <- diff(log(Cl(GSPC)))
returns <- diff(log(data$Close))

returns <- na.omit(returns)

#-------------------------------------------------
# BƯỚC 2:  Kiểm tra tính dừng
#-------------------------------------------------

adf.test(returns) 
#Dickey-Fuller = -5.105, Lag order = 4, p-value = 0.01
adf.test(data$EPU)
#Dickey-Fuller = -2.7249, Lag order = 4, p-value = 0.2755


# Augmented Dickey-Fuller test
# H0: Chuỗi không dừng (có unit root)
# p-value < 0.05 → Bác bỏ H0 → Chuỗi DỪNG
# Lấy sai phân bậc 1
diff_returns <- diff(returns)
diff_EPU <- diff(data$EPU)
data$Return <- c(NA, diff(log(data$Close)))

data$EPU_diff <- c(NA, diff(data$EPU))
#bỏ NA
data_clean <- na.omit(data)

#-------------------------------------------------
# BƯỚC 3: CHẠY QR 
#-------------------------------------------------
# Tạo biến độ trễ cho Y (Return) và X (đã xử lý)
data_clean$Return_lag <- dplyr::lag(data_clean$Return, 1)
data_clean$dEPU_lag <- dplyr::lag(data_clean$EPU, 1)


# Kiểm tra độ lệch chuẩn
sd(data_clean$Return)      # Nên ~0.04 (4% hàng tháng)
sd(data_clean$dEPU_lag)    # ???
# Nên làm trước khi chạy mô hình
data_clean$Return_std <- scale(data_clean$Return)
data_clean$dEPU_std <- scale(data_clean$dEPU_lag)

# Loại bỏ NA do tạo độ trễ
data_clean <- na.omit(data_clean)
## CASE 1: restrichted

# 1. Mô hình GIỚI HẠN (chỉ có Return_lag)
model_restricted <- rq(Return ~ Return_lag, 
                       data = data_clean, 
                       tau = 0.25)
summary(model_restricted)
#           coefficients lower bd upper bd
#(Intercept) -0.00803     -0.01188  0.00012
#Return_lag   0.25220      0.04909  0.36130

## CASE 2: un-restrichted


# 2. Mô hình KHÔNG GIỚI HẠN (có cả Return_lag và dEPU_lag)
model_unrestricted <- rq(Return ~ Return_lag + dEPU_lag, 
                         data = data_clean, 
                         tau = 0.25)
summary(model_unrestricted)
#           coefficients lower bd upper bd
#(Intercept) -0.02215     -0.03915  0.00688
#Return_lag   0.18094     -0.05283  0.34710
#dEPU_lag     0.00012     -0.00015  0.00024

# Tính thủ công cho 1 phân vị
tau <- 0.25; h <- 0.05

# Phần dư
u1 <- residuals(model_restricted)
u2 <- residuals(model_unrestricted)
print (u1)
print (u2)

neg_u1<- u1[u1<0]
neg_u2<- u2[u2<0]

pos_u1 <- u1[u1 >= 0]     # Phần không âm (bao gồm cả 0)
pos_u2 <- u2[u2 >= 0]     # Phần không âm (bao gồm cả 0)

print(neg_u1)
print(neg_u2)

# Hàm mất mát phân vị
sum_pos1 <- sum(pos_u1 * tau, na.rm = TRUE)   
sum_pos2 <- sum(pos_u2 * tau, na.rm = TRUE)        # u≥0: τ*u
# u≥0: τ*u
sum_n1 <- sum(neg_u1 * (tau-1), na.rm = TRUE)
sum_n2 <- sum(neg_u2 * (tau-1), na.rm = TRUE)

sum_p1 <- sum_n1+sum_pos1
sum_p2 <- sum_n2+sum_pos2


# QTE
QTE_exact <- log(sum_p1) - log(sum_p2)
print(QTE_exact)
# Các phân vị cần phân tích
tau_vec <- c(0.10, 0.25, 0.50, 0.75, 0.90)
tau_labels <- c("Giảm mạnh (τ=0.10)", "Giảm nhẹ (τ=0.25)", 
                "Bình thường (τ=0.50)", "Tăng nhẹ (τ=0.75)", 
                "Tăng mạnh (τ=0.90)")
data_model<-data_clean
## CASE1+2:
# Lưu trữ kết quả
results <- data.frame()

cat("\n=== KẾT QUẢ HỒI QUY LƯỢNG TỬ ===\n")

for (i in 1:length(tau_vec)) {
  tau <- tau_vec[i]
  label <- tau_labels[i]
  
  # 1. Mô hình GIỚI HẠN (chỉ có Return_lag)
  model_restricted <- rq(Return ~ Return_lag, 
                         data = data_model, 
                         tau = tau)
  
  # 2. Mô hình KHÔNG GIỚI HẠN (có cả Return_lag và dEPU_lag)
  model_unrestricted <- rq(Return ~ Return_lag + dEPU_lag, 
                           data = data_model, 
                           tau = tau)
  
  # Lấy hệ số
  coef_restricted <- coef(model_restricted)
  coef_unrestricted <- coef(model_unrestricted)
  
  # Tính R² cho mỗi mô hình (pseudo R²)
  r2_restricted <- 1 - model_restricted$rho / model_restricted$rho[1]
  r2_unrestricted <- 1 - model_unrestricted$rho / model_unrestricted$rho[1]
  # 1. PHIÊN BẢN ĐƠN GIẢN: Giả định f_τ(·) ≈ 1
  theta<-coef_unrestricted[2]
  QTE_simple <- theta
  
  # 2. PHIÊN BẢN CHÍNH XÁC HƠN: Ước lượng f_τ(·) nhanh
  u1 <- residuals(model_restricted)
  u2 <- residuals(model_unrestricted)
  
  neg_u1<- u1[u1<0]
  neg_u2<- u2[u2<0]
  pos_u1 <- u1[u1 >= 0]     # Phần không âm (bao gồm cả 0)
  pos_u2 <- u2[u2 >= 0]     # Phần không âm (bao gồm cả 0)
  
  print(neg_u1)
  print(neg_u2)
  
  # Hàm mất mát phân vị
  sum_pos1 <- sum(pos_u1 * tau, na.rm = TRUE)   
  sum_pos2 <- sum(pos_u2 * tau, na.rm = TRUE)        # u≥0: τ*u
  # u≥0: τ*u
  sum_n1 <- sum(neg_u1 * (tau-1), na.rm = TRUE)
  sum_n2 <- sum(neg_u2 * (tau-1), na.rm = TRUE)
  
  sum_p1 <- sum_n1+sum_pos1
  sum_p2 <- sum_n2+sum_pos2
  
  

  # QTE
  QTE_exact <- log(sum_p1) - log(sum_p2)
  print(QTE_exact)
  
  # Trong thực tế, nên lấy từ summary(model) trong Bước 4
  std_error <- abs(theta) / 2  # Ước lượng đơn giản
  
  ci_lower_90 <- theta - 1.645 * std_error
  ci_upper_90 <- theta + 1.645 * std_error
  
  ci_lower_95 <- theta - 1.96 * std_error
  ci_upper_95 <- theta + 1.96 * std_error
  
  # Kiểm tra ý nghĩa thống kê
  significant_90 <- ifelse(ci_lower_90 > 0 || ci_upper_90 < 0, "CÓ (90%)", "KHÔNG")
  significant_95 <- ifelse(ci_lower_95 > 0 || ci_upper_95 < 0, "CÓ (95%)", "KHÔNG")
  
  # Lưu kết quả
  results <- rbind(results, 
                   data.frame(
                     Phan_vi = tau,
                     Trang_thai = label,
                     # Mô hình giới hạn
                     Restricted_Intercept = coef_restricted[1],
                     Restricted_ReturnLag = coef_restricted[2],
                     Restricted_R2 = r2_restricted,
                     # Mô hình không giới hạn
                     Unrestricted_Intercept = coef_unrestricted[1],
                     Unrestricted_ReturnLag = coef_unrestricted[2],
                     Unrestricted_EPULag = coef_unrestricted[3],
                     Unrestricted_R2 = r2_unrestricted,
                     QTE = QTE_exact,
                     CI_90 = sprintf("[%.4f, %.4f]", ci_lower_90, ci_upper_90),
                     Y_nghia_90 = significant_90,
                     CI_95 = sprintf("[%.4f, %.4f]", ci_lower_95, ci_upper_95),
                     Y_nghia_95 = significant_95
                   )
  )
  
  # Hiển thị kết quả
  cat("\n", label, ":\n")
  cat("1. Mô hình GIỚI HẠN (chỉ có Return_lag):\n")
  print(summary(model_restricted))
  
  cat("\n2. Mô hình KHÔNG GIỚI HẠN (có cả Return_lag và dEPU_lag):\n")
  print(summary(model_unrestricted))
  
  # So sánh R²
  r2_diff <- r2_unrestricted - r2_restricted
  cat(sprintf("\n→ Cải thiện R² khi thêm EPU: %.4f\n", r2_diff))
  
  cat("----------------------------------------\n")
}

# Hiển thị bảng tổng hợp
cat("\n=== BẢNG TỔNG HỢP HỆ SỐ ===\n")
print(results)

# ====================
#-------------------------------------------------
# BƯỚC 5: TÍNH QTE TỪ KẾT QUẢ CÓ SẴ LẤY DATA & XƯ LÝ
#-------------------------------------------------
# Tạo dataframe lưu kết quả QTE cuối cùng
QTE_final <- data.frame()

for (i in 1:nrow(results)) {
  # Lấy thông tin từ kết quả đã tính ở Bước 4
  tau <- results$Phan_vi[i]
  label <- results$Trang_thai[i]
  theta <- results$Unrestricted_EPULag[i]  # θ_τ đã có sẵn!
  QTE_loop<-results$QTE[i]
  print(QTE_loop)
  # 4. Lưu kết quả
  QTE_final <- rbind(QTE_final,
                     data.frame(
                       Phan_vi = tau,
                       Trang_thai = label,
                       Theta_τ = round(theta, 6),
                       QTE_chínhxác = QTE_loop,
                       
                       R2_cảithiện = round(results$Unrestricted_R2[i] - results$Restricted_R2[i], 4)
                     )
  )
  
  # Hiển thị từng phân vị
  cat(sprintf("\n%s:\n", label))
  cat(sprintf("  θ_τ = %.6f\n", theta))
  cat(sprintf("  QTE ≈ %.6f\n", ifelse(is.na(QTE_exact), QTE_simple, QTE_loop)))
  cat(sprintf("  Khoảng tin cậy 90%%: [%.4f, %.4f] → %s\n", 
              ci_lower_90, ci_upper_90, significant_90))
  cat(sprintf("  R² cải thiện khi thêm EPU: %.4f\n", 
              results$Unrestricted_R2[i] - results$Restricted_R2[i]))
}

# Hiển thị bảng tổng hợp
cat("=== BẢNG TỔNG HỢP KẾT QUẢ QTE CHO TẤT CẢ PHÂN VỊ ===\n")
print(QTE_final)
#-------------------------------------------------
# BƯỚC 7: LẤY BOOSTRAPING
#-------------------------------------------------
# ==================== CÀI ĐẶT GÓI ====================
library(quantreg)

# ==================== THAM SỐ BOOTSTRAP ====================
B <- 499                     # Số lần bootstrap
tau_vec <- c(0.10, 0.50, 0.90)  # 3 phân vị cần kiểm định
k <- 1                       # Độ trễ 1 tháng

# Kích thước khối (dựa trên số quan sát)
n <- nrow(data_model)
l <- round(n^(1/3))          # Quy tắc n^(1/3)

cat("=== THAM SỐ BOOTSTRAP ===\n")
cat("Số quan sát:", n, "\n")
cat("Kích thước khối:", l, "month\n")
cat("Số lần bootstrap:", B, "\n\n")

# ==================== QTE THỰC TẾ (ĐÃ TÍNH) ====================
# Sử dụng kết quả đã tính trước đó-QTE chinh xác
QTE_thuc <- c(0.003926615   , 0.043611075 , 0.131923898       )  # τ=0.10, 0.50, 0.90
names(QTE_thuc) <- paste0("τ=", tau_vec)

cat("=== QTE THỰC TẾ ===\n")


# ==================== HÀM TÍNH QTE NHANH ====================
tinh_QTE_nhanh <- function(data, tau) {
  # Ước lượng mô hình không giới hạn
  model1 <- rq(Return ~ Return_lag , data = data, tau = tau, method = "fn")
  model2 <- rq(Return ~ Return_lag + dEPU_lag, data = data, tau = tau,method = "fn")

  
  
  # 2. PHIÊN BẢN CHÍNH XÁC HƠN: Ước lượng f_τ(·) nhanh
  u1 <- residuals(model1)
  u2 <- residuals(model2)
  neg_u1<- u1[u1<0]
  neg_u2<- u2[u2<0]
  pos_u1 <- u1[u1 >= 0]     # Phần không âm (bao gồm cả 0)
  pos_u2 <- u2[u2 >= 0]     # Phần không âm (bao gồm cả 0)
  
  print(neg_u1)
  print(neg_u2)
  
  # Hàm mất mát phân vị
  sum_pos1 <- sum(pos_u1 * tau, na.rm = TRUE)   
  sum_pos2 <- sum(pos_u2 * tau, na.rm = TRUE)        # u≥0: τ*u
  # u≥0: τ*u
  sum_n1 <- sum(neg_u1 * (tau-1), na.rm = TRUE)
  sum_n2 <- sum(neg_u2 * (tau-1), na.rm = TRUE)
  
  sum_p1 <- sum_n1+sum_pos1
  sum_p2 <- sum_n2+sum_pos2
  
  
  
  # QTE
  QTE_exact <- log(sum_p1) - log(sum_p2)
  
  
  return(QTE_exact)
}

# ==================== HÀM BOOTSTRAP KHỐI ĐƠN GIẢN ====================
bootstrap_don_gian <- function() {
  # Tạo HAI chỉ mục độc lập
  boot_y <- sample(1:n, n, replace = TRUE)  # Đơn giản: bootstrap thường
  boot_x <- sample(1:n, n, replace = TRUE)  # Độc lập với boot_y
  
  # Tạo dữ liệu bootstrap
  data_boot <- data.frame(
    Return = data_model$Return[boot_y],
    Return_lag = data_model$Return_lag[boot_y],
    dEPU_lag = data_model$dEPU_lag[boot_x]  # Độc lập!
  )
  
  return(data_boot)
}

# ==================== THỰC HIỆN BOOTSTRAP ====================
cat("=== BẮT ĐẦU BOOTSTRAP ===\n")

# Ma trận lưu kết quả bootstrap
QTE_boot <- matrix(NA, nrow = B, ncol = 3)
colnames(QTE_boot) <- paste0("τ=", tau_vec)

# Bootstrap loop

#View(data_boot)
summary(data_boot)
set.seed(123)  # Để tái lập kết quả
for (b in 1:B) {
  # Tạo dữ liệu bootstrap
  data_boot <- bootstrap_don_gian()
  
  # Tính QTE cho 3 phân vị
  for (j in 1:3) {
    val <- tryCatch({tinh_QTE_nhanh(data_boot, tau_vec[j])}, error = function(e) return(NA))
    
    QTE_boot[b, j] <- val
    
    # Đếm số lần QTE_boot >= QTE_thuc
    
    }

  
  # Hiển thị tiến độ
  if (b %% 100 == 0) cat("Đã xong", b, "lần\n")
}

# ==================== TÍNH GIÁ TRỊ P ====================
cat("\n=== KẾT QUẢ KIỂM ĐỊNH ===\n")

# Bảng kết quả
QTE_boot
ket_qua <- data.frame(
  Phan_vi = tau_vec,
  Trang_thai = c("Lỗ thấp", "Trung bình", "Lợi nhuận cao"),
  QTE_thuc = QTE_thuc,
  p_value = NA,
  Y_nghia = NA
)
# Tính p-value cho từng phân vị
for (i in 1:3) {
  # Đếm số lần QTE_boot >= QTE_thuc
  count <- sum(QTE_boot[, i] >= QTE_thuc[i], na.rm = TRUE)
  
  # Tính p-value một phía
  p_val <- (1 + count) / (B + 1)
  
  ket_qua$p_value[i] <- round(p_val, 4)
  ket_qua$Y_nghia[i] <- ifelse(p_val < 0.05, "CÓ", "KHÔNG")
  # Hiển thị chi tiết
  cat(sprintf("\nτ=%.2f (%s):\n", tau_vec[i], ket_qua$Trang_thai[i]))
  cat(sprintf("  QTE thực: %.6f\n", QTE_thuc[i]))
  cat(sprintf("  số lần QTE boot > QTE thực: %.6f\n", count))
  
  cat(sprintf("  p-value: %.4f\n", p_val))
  cat(sprintf("  Ý nghĩa (α=0.05): %s\n", ket_qua$Y_nghia[i]))
  
  # Phân phối bootstrap
  cat(sprintf("  Phân vị 95%% bootstrap: [%.6f, %.6f]\n",
              quantile(QTE_boot[, i], 0.025, na.rm = TRUE),
              quantile(QTE_boot[, i], 0.975, na.rm = TRUE)))
}

# ==================== HIỂN THỊ BẢNG KẾT QUẢ ====================
cat("BẢNG KẾT QUẢ BOOTSTRAP\n")
print(ket_qua)


# ==================== KẾT LUẬN ====================
cat("KẾT LUẬN:\n")

for (i in 1:3) {
  if (ket_qua$p_value[i] < 0.05) {
    cat(sprintf("• τ=%.2f: BÁC BỎ H₀ (p=%.4f < 0.05). CÓ bằng chứng truyền thông tin.\n",
                tau_vec[i], ket_qua$p_value[i]))
  } else {
    cat(sprintf("• τ=%.2f: KHÔNG bác bỏ H₀ (p=%.4f ≥ 0.05). KHÔNG có bằng chứng.\n",
                tau_vec[i], ket_qua$p_value[i]))
  }
}
#       Phan_vi    Trang_thai    QTE_thuc p_value Y_nghia
#τ=0.1     0.1       Lỗ thấp 0.003926615   0.582   KHÔNG Trong thực tế không có chuyện truyền dẫn tích cực và từ biến động tới lợi nhuận 
#τ=0.5     0.5    Trung bình 0.043611075   0.004   CÓ trong thực tế có khả năng truyền dẫn từ biến động tới lợi nhuận 
#τ=0.9     0.9 Lợi nhuận cao 0.131923898   0.002      CÓ Trong thực tế có khả năng bị bất động sẽ đẩy lợi nhuận tăng.
# Lỗ thấp: Boostap Khẳng định rằng các kết quả của mô hình trên chỉ có tính chính xác tới dữ liệu đó chứ không có tính bao quát cho tất cả.
