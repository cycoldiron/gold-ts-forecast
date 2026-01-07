# Gold Time Series Forecasting

## Overview

This project focuses on modeling the volatility process in gold returns and mapping volatility regimes over time. I use various statistical & econometric methods to (i) identify unique volatility regimes (ii) characterize how volatility behaves within each regime (clustering, persistence, leverage effects), and (iii) generate regime-specific forecasts. Analyses are run across three windows (1990–present, 2010–present, 2018–Oct 2025).

``` text
## Repository structure

```text
gold-ts-forecast/
├── .gitignore
├── .RData
├── .Rhistory
├── .Rprofile
├── data/
├── drafts/
│   ├── 01_draft_files/
│   ├── 01_draft.html
│   ├── gold_volatility_paper_2026-01-07.qmd
│   └── images/
├── gold_volatility_paper_2026-01-07.pdf
├── gold-ts-forecast.Rproj
├── LICENSE
├── README.html
├── README.md
├── renv/
├── renv.lock
└── scripts/

└── scripts/
    ├── 01_get_economic_data.R.R
    ├── 02_clean_transform_economic.R
    ├── 03_get_uncertainty_data.R
    ├── 04_uncertainty_clean_data.R
    ├── 05_explore_gold_acf_pacf.R
    ├── 06_stationarity_tests.R
    ├── 07_struc_break.R
    ├── 07b_struc_break.R
    ├── 08_white_noise_test.R
    ├── 08b_white_noise_test.R
    ├── 09_changepoint_var_test.R
    ├── 10_garch_workflow.R
    ├── 10b_pacf_variance.R
    ├── 10c_revised_model_testing.R
    ├── 10d_summary_table_revision.R
    ├── 10e_egarch_diagnostics_revised.R
    ├── 11_variance_break_attribution.R
    ├── 12_break_aware_egarch_regime_dumm.R
    ├── 12_variance_model_testing.R
    ├── 13_break_aware_egarch_regime_dumm.R
    ├── 14_variance_diagnostics.R
    ├── 14b_residual_updates.R
    ├── 15_build_egarch_volatility.R
    ├── 16_uncertainty_links_figures.R
    ├── 16b_uncertainty_links_figures.R
    ├── 16c_overlay_plots.R
    └── 17_break_data.R
```

|     |
|-----|
|     |
