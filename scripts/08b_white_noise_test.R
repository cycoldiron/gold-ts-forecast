# —— packages ——————————————————————————————————————————————
library(dplyr)
library(gt)
library(purrr)
library(tibble)
library(fs)

# —— paths ————————————————————————————————————————————————
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/03_white_noise"
dir_create(out_dir, recurse = TRUE)

# —— data: daily log returns ——————————————————————————————
# assumes `gold_full_cleaned` is already in memory as you showed
r <- na.omit(gold_full_cleaned$close_log_return)

# keep/dismiss mean in the mean model
t_mu    <- t.test(r)
keep_mu <- (t_mu$p.value < 0.05)

# —— fit mean candidates ————————————————————————————————
m00 <- arima(r, order = c(0,0,0), include.mean = keep_mu) # ARMA(0,0)
m10 <- arima(r, order = c(1,0,0), include.mean = keep_mu) # AR(1)
m01 <- arima(r, order = c(0,0,1), include.mean = keep_mu) # MA(1)

fits <- list(
  "ARMA(0,0)" = m00,
  "AR(1)"     = m10,
  "MA(1)"     = m01
)

# —— table 1: mean model comparison (AIC/BIC) ————————————
tab_ic <- tibble(
  Model = names(fits),
  AIC   = sapply(fits, AIC),
  BIC   = sapply(fits, BIC)
)

gt_ic <- tab_ic |>
  mutate(
    AIC = round(AIC, 2),
    BIC = round(BIC, 2)
  ) |>
  gt(rowname_col = "Model") |>
  tab_header(
    title    = md("**Mean Model Comparison — Daily Gold Log Returns**"),
    subtitle = paste0(
      "N = ", length(r),
      " | Mean included: ", ifelse(keep_mu, "Yes (t-test p < 0.05)", "No (demeaned)")
    )
  ) |>
  fmt_number(columns = c(AIC, BIC), decimals = 2) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_body(columns = AIC, rows = AIC == min(AIC)),
      cells_body(columns = BIC, rows = BIC == min(BIC))
    )
  ) |>
  opt_row_striping() |>
  tab_source_note(md("Source: Stooq (XAUUSD), author’s calculations."))

gt::gtsave(gt_ic, filename = file.path(out_dir, "mean_model_comparison.png"))

# —— table 2: Ljung–Box on mean-model residuals ————————————
lb_df <- imap_dfr(fits, ~{
  fitdf <- ifelse(.y == "ARMA(0,0)", 0, 1)  # number of mean parameters for LB
  lb <- Box.test(residuals(.x), lag = 20, type = "Ljung-Box", fitdf = fitdf)
  tibble(
    Model   = .y,
    Lags    = 20,
    `X^2 (Q)` = unname(lb$statistic),
    df      = unname(lb$parameter),
    `p-value` = lb$p.value
  )
})

gt_lb <- lb_df |>
  mutate(
    `X^2 (Q)` = round(`X^2 (Q)`, 3),
    `p-value` = signif(`p-value`, 3)
  ) |>
  gt() |>
  tab_header(
    title = md("**Ljung–Box Test on Mean-Model Residuals**"),
    subtitle = "H₀: no autocorrelation up to lag L"
  ) |>
  fmt_number(columns = c(`X^2 (Q)`), decimals = 3) |>
  opt_row_striping() |>
  tab_source_note(md("Note: Q ≈ χ²(df) under H₀.  Source: author’s calculations."))

gt::gtsave(gt_lb, filename = file.path(out_dir, "ljung_box_residuals.png"))

# (Optional) also save HTML versions
# gtsave(gt_ic, file.path(out_dir, "mean_model_comparison.html"))
# gtsave(gt_lb, file.path(out_dir, "ljung_box_residuals.html"))
