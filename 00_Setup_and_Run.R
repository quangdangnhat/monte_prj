# =============================================================================
# HƯỚNG DẪN CHẠY PROJECT TRÊN RSTUDIO
# Financial Econometrics Coursework 2025-2026
# =============================================================================

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BƯỚC 1: CÀI ĐẶT CÁC PACKAGES CẦN THIẾT                                   ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Chạy đoạn code này TRƯỚC TIÊN để cài đặt tất cả packages cần thiết
# Bôi đen đoạn code và nhấn Ctrl+Enter (Windows) hoặc Cmd+Enter (Mac)

packages_needed <- c(
  "lmtest",      # Breusch-Pagan test
  "quantreg",    # Quantile regression cho QTE
  "tseries",     # ADF test cho stationarity
  "zoo",         # Time series manipulation
  "xts",         # Extended time series
  "dplyr",       # Data manipulation
  "ggplot2",     # Plotting
  "quantmod",    # Download financial data
  "readr",       # Read CSV files
  "rmarkdown",   # Generate PDF report
  "knitr"        # Knitr for reports
)

# Cài đặt packages chưa có
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "http://cran.r-project.org")
  } else {
    cat(sprintf("%s: OK\n", pkg))
  }
}

cat("Checking and installing packages...\n")
cat("=====================================\n")
for (pkg in packages_needed) {
  install_if_missing(pkg)
}
cat("=====================================\n")
cat("All packages installed!\n\n")

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BƯỚC 2 (TÙY CHỌN): CÀI TINYTEX ĐỂ TẠO PDF                               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Nếu bạn muốn tạo báo cáo PDF, cần cài TinyTeX
# Bỏ comment dòng dưới và chạy (mất khoảng 5-10 phút)

# install.packages("tinytex")
# tinytex::install_tinytex()

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BƯỚC 3: KIỂM TRA WORKING DIRECTORY                                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Kiểm tra thư mục hiện tại
cat("Current working directory:\n")
cat(getwd(), "\n\n")

# Liệt kê các file R trong thư mục
cat("R files in directory:\n")
print(list.files(pattern = "\\.R$"))

# Nếu không thấy các file Part_A_Monte_Carlo.R, Part_B_QTE.R...
# Hãy set working directory bằng một trong các cách:
#
# Cách 1: Dùng menu RStudio
#   Session -> Set Working Directory -> Choose Directory
#   Chọn thư mục chứa các file R
#
# Cách 2: Dùng code (thay đường dẫn phù hợp)
#   setwd("C:/Users/YourName/Documents/monte_prj")  # Windows
#   setwd("~/Documents/monte_prj")                   # Mac/Linux

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BƯỚC 4: CHẠY PART A - MONTE CARLO SIMULATION                            ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# OPTION A: Chạy Quick Test (~3-5 phút) - Để kiểm tra code hoạt động
cat("\n=== ĐỂ CHẠY QUICK TEST, BỎ COMMENT DÒNG DƯỚI ===\n")
# source("Part_A_Quick_Test.R")

# OPTION B: Chạy Full Simulation (~45-60 phút) - Cho kết quả chính xác
cat("=== ĐỂ CHẠY FULL SIMULATION, BỎ COMMENT DÒNG DƯỚI ===\n")
# source("Part_A_Monte_Carlo.R")

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BƯỚC 5: CHẠY PART B - QTE ANALYSIS                                      ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Bước 5a: Chuẩn bị dữ liệu (download hoặc dùng simulated data)
cat("\n=== ĐỂ CHUẨN BỊ DỮ LIỆU, BỎ COMMENT DÒNG DƯỚI ===\n")
# source("Part_B_Data_Download.R")

# Bước 5b: Chạy QTE Analysis (~10-15 phút)
cat("=== ĐỂ CHẠY QTE ANALYSIS, BỎ COMMENT DÒNG DƯỚI ===\n")
# source("Part_B_QTE.R")

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BƯỚC 6: TẠO BÁO CÁO PDF                                                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Sau khi chạy Part A và Part B, tạo báo cáo PDF
cat("\n=== ĐỂ TẠO PDF REPORT, BỎ COMMENT DÒNG DƯỚI ===\n")
# rmarkdown::render("Coursework_Report.Rmd", output_format = "pdf_document")

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  HƯỚNG DẪN NHANH                                                          ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                    HƯỚNG DẪN NHANH                                        ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                           ║\n")
cat("║  1. Chạy file này trước để cài packages (Ctrl+Shift+Enter)               ║\n")
cat("║                                                                           ║\n")
cat("║  2. Mở Part_A_Quick_Test.R và chạy toàn bộ (Ctrl+Shift+Enter)            ║\n")
cat("║     -> Kết quả hiện trong Console                                        ║\n")
cat("║                                                                           ║\n")
cat("║  3. Mở Part_B_Data_Download.R và chạy toàn bộ                            ║\n")
cat("║     -> Tạo file dữ liệu                                                  ║\n")
cat("║                                                                           ║\n")
cat("║  4. Mở Part_B_QTE.R và chạy toàn bộ                                      ║\n")
cat("║     -> Kết quả QTE hiện trong Console                                    ║\n")
cat("║                                                                           ║\n")
cat("║  5. Mở Coursework_Report.Rmd và nhấn 'Knit' để tạo PDF                   ║\n")
cat("║                                                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════════╝\n")
