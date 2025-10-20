---
editor_options: 
  markdown: 
    wrap: sentence
---

# 🪙 Gold Time Series Forecasting

## Overview

This project focuses on modeling the volatility process in gold returns and mapping volatility regimes over time.
I use various statistical & econometric methods (outlined below)—to (i) detect when volatility shifts, (ii) characterize how it behaves within each regime (clustering, persistence, leverage effects), and (iii) generate regime-specific forecasts.
Analyses are run across three windows (1990–present, 2010–present, 2018–Oct 2025).

I then link volatility regimes and pricing dynamics to macroeconomic and geopolitical uncertainty indices (EPU, GPR).

## Practical Significance

Gold serves as a **safe-haven asset** and **inflation hedge**, especially during periods of macroeconomic instability and geopolitical tension—such as right now (October 2025)\
This project seeks to:

1.  Identify episodes of instability or regime change in the gold market.\
2.  Test whether gold prices follow a **random walk** (they do).\
3.  Quantify how **economic policy uncertainty (EPU)** and **geopolitical risk (GPR)** relate to changes in gold volatility.\
4.  Build a foundation for **ARIMA, EGARCH, and structural-break–adjusted models** to improve short-term forecasting accuracy.

------------------------------------------------------------------------

## Statistical Methods

1.  **Descriptive Analysis:** Price, log-price, and log-return visualization.

2.  **Stationarity Tests:** ADF; variance-ratio for random-walk behavior.

3.  **ACF/PACF Diagnostics:** Identify order / type of process.

4.  **White-Noise/Independence:** Ljung–Box on levels and (standardized) residuals.

5.  **Structural Breaks:** Bai–Perron (means/parameters) and regime segmentation.

6.  **Variance Changepoints:** PELT on returns/return-squares for variance shifts.

7.  **Volatility Models:** ARCH/GARCH (t-errors), **EGARCH** for asymmetry; rolling/segmented re-fits for regime sensitivity.

8.  **Forecasting:** Multi-step variance forecasts; half-life of shocks; comparison by window.

## Repository Structure

``` text
gold-ts-forecast/
├── data/
│   ├── raw/
│   │   ├── 01_raw_economic_data.RData
│   │   ├── 02_raw_uncertainty.Rdata
│   │   ├── _data_note.Rmd
│   │   └── _data_note.html
│   └── clean/
│       ├── 01_clean_economic_data.RData
│       ├── 02_clean_uncertainty_data.Rdata
│       └── 03_lb_and_varbreaks.RData
│
├── scripts/
│   ├── 01_get_economic_data.R
│   ├── 02_clean_transform_economic.R
│   ├── 03_get_uncertainty_data.R
│   ├── 04_uncertainty_clean_data.R
│   ├── 05_explore_gold_acf_pacf.R
│   ├── 06_stationarity_tests.R
│   ├── 07_struc_break.R
│   ├── 07b_struc_break.R
│   ├── 08_white_noise_test.R
│   └── 09_changepoint_var_test.R
│
├── figures/
│   ├── 01_overview/
│   ├── 02_diagnostics/
│   │   ├── 01_acf_pacf/
│   │   ├── 02_adf/
│   │   └── 03_white_noise/
│   ├── 03_breaks/
│   ├── 04_models/
│   ├── 05_forecasts/
│   └── 06_macro/
│
├── results/
│   └── diagnostics/
│       └── 02_diagnostics/
│           ├── 01_acf_pacf/
│           ├── 02_adf/
│           └── 03_white_noise/
│
├── paper/
├
```

## Folder Guide

| Path | Purpose |
|----|----|
| `data/raw/` | Unprocessed datasets and the data note (`_data_note.*`). |
| `data/clean/` | Cleaned datasets and derived objects used downstream. |
| `scripts/` | Analysis scripts in numerical run order (01 → 09). |
| `figures/` | Plots grouped by stage (overview, diagnostics, breaks, etc.). |
| `results/diagnostics/` | Tabular outputs from tests (ADF, Ljung–Box, white noise). |
| `paper/` | Manuscript/report materials. |
| `renv/`, `renv.lock` | Reproducible R environment and locked package versions. |

## Author

**Cy Coldiron**\
B.A. Economics; Statistics & Data Science — **UC Santa Barbara**\
Visiting Student — **University of Cape Town**\
*Time-Series Econometrics Project · 2025*
