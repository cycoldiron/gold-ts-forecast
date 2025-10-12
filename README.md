# 🪙 Gold Time Series Forecasting

## Overview
This project analyzes the **dynamic behavior and predictability of gold prices** through a comprehensive time-series econometric framework.  
The primary goal is to understand the **stability, volatility, and structural evolution** of gold returns over three times frames (1990 / 2010 / 2018 - Oct 2025) — providing insights relevant to **macroeconomic uncertainty, market efficiency, and financial risk management**.

The analysis combines exploratory data visualization, unit-root testing, structural-break detection, and diagnostic evaluation to identify the statistical properties of gold prices and returns prior to model-based forecasting.

---

## Practical Significance
Gold serves as a **safe-haven asset** and **inflation hedge**, especially during periods of macroeconomic instability and geopolitical tension.  
This project seeks to:

1. Identify episodes of instability or regime change in the gold market.  
2. Test whether gold prices follow a **random walk** or exhibit **mean reversion**.  
3. Quantify how **economic policy uncertainty (EPU)** and **geopolitical risk (GPR)** relate to changes in gold volatility.  
4. Build a foundation for **ARIMA, GARCH, and structural-break–adjusted models** to improve short-term forecasting accuracy.

---

## Statistical Methods
1. **Descriptive Time-Series Analysis** — visualization of price levels, log-prices, and log-returns.  
2. **Stationarity & Unit Root Tests** — Augmented Dickey-Fuller (ADF), variance-ratio tests.  
3. **ACF/PACF Diagnostics** — identifying autoregressive (AR) and moving-average (MA) components.  
4. **White-Noise & Independence Tests** — Ljung-Box and residual autocorrelation checks.  
5. **Structural Break Tests** — Bai–Perron multiple breakpoint detection and regime segmentation.  
6. **Variance Change Detection** — rolling variance and changepoint variance analysis.  
7. *(Later stages)* — ARIMA and GARCH modeling for forecasting.

## Repository Structure

```text
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

| Path                         | Purpose                                                        |
|-----------------------------|----------------------------------------------------------------|
| `data/raw/`                 | Unprocessed datasets and the data note (`_data_note.*`).      |
| `data/clean/`               | Cleaned datasets and derived objects used downstream.         |
| `scripts/`                  | Analysis scripts in numerical run order (01 → 09).            |
| `figures/`                  | Plots grouped by stage (overview, diagnostics, breaks, etc.). |
| `results/diagnostics/`      | Tabular outputs from tests (ADF, Ljung–Box, white noise).     |
| `paper/`                    | Manuscript/report materials.                                   |
| `renv/`, `renv.lock`        | Reproducible R environment and locked package versions.        |


## Author
**Cy Coldiron**  
B.A. Economics; Statistics & Data Science — **UC Santa Barbara**  
Visiting Student — **University of Cape Town**  
*Time-Series Econometrics Project · 2025*

