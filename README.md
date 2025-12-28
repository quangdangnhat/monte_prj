# Financial Econometrics Coursework 2025-2026

## Project Structure

```
monte_prj/
├── Part_A_Monte_Carlo.R      # Full Monte Carlo simulation (R=1000, B=499)
├── Part_A_Quick_Test.R       # Quick test version (R=50, B=99)
├── Part_B_Data_Download.R    # Data preparation for QTE analysis
├── Part_B_QTE.R              # Quantile Transfer Entropy analysis
├── Coursework_Report.Rmd     # Complete report with theory
├── run_all.R                 # Master script to run everything
├── monte.R                   # Original code (reference)
└── README.md                 # This file
```

## Requirements

```r
install.packages(c("lmtest", "quantreg", "tseries", "zoo",
                   "xts", "readr", "dplyr", "ggplot2",
                   "quantmod", "rmarkdown"))
```

For PDF generation:
```r
install.packages("tinytex")
tinytex::install_tinytex()
```

## How to Run

### Option 1: Run All (Interactive)
```r
source("run_all.R")
```

### Option 2: Run Parts Separately

**Part A - Monte Carlo:**
```r
# Quick test (~5 min)
source("Part_A_Quick_Test.R")

# Full simulation (~1 hour)
source("Part_A_Monte_Carlo.R")
```

**Part B - QTE Analysis:**
```r
# Step 1: Prepare data
source("Part_B_Data_Download.R")

# Step 2: Run QTE analysis
source("Part_B_QTE.R")
```

**Generate Report:**
```r
rmarkdown::render("Coursework_Report.Rmd")
```

## Output Files

- `Part_A_Results.csv` - Monte Carlo power comparison results
- `Part_B_QTE_Results.csv` - QTE analysis results
- `Coursework_Report.pdf` - Complete report

## Part A: Tests Implemented

1. **Sup-GQ Test** - Maximum F-statistic across split points
2. **Fisher's G Test** - Combined p-value test: G = -2Σlog(p)
3. **Breusch-Pagan Test** - Variance depends on regressors
4. **White Test** - General heteroscedasticity test

## Part B: QTE Analysis

- **Quantiles tested**: τ = 0.10, 0.50, 0.90
- **Lags tested**: k = 1, 2, 3, 5
- **Bootstrap method**: Stationary Block Bootstrap (B=499)
- **Block length**: n^(1/3) rule of thumb

## Notes

- For real data analysis, download EPU from https://www.policyuncertainty.com/
- If data files not found, simulated data will be used for demonstration
- Full simulation (Part A) takes approximately 1 hour
