# =============================================================================
# 06_garch_workflow_FINAL.R
# Volatility diagnostics + GARCH family fit + diagnostics (polished)
# Outputs: PNG figures/tables + TXT summary (no HTML)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(gt); library(FinTS); library(rugarch); library(fs)
})

# --------- Paths --------------------------------------------------------------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/04_arch_garch"
dir_create(out_dir, recurse = TRUE)

# --------- Data ---------------------------------------------------------------
# Assumes `gold_full_cleaned` has: date, close_log_return
r_all  <- gold_full_cleaned$close_log_return
dates  <- gold_full_cleaned$date[!is.na(r_all)]
r      <- na.omit(r_all)
n_obs  <- length(r)

# --------- Mean decision ------------------------------------------------------
t_mu    <- t.test(r)
keep_mu <- (t_mu$p.value < 0.05)   # TRUE => include mean; FALSE => demean

# --------- Helpers ------------------------------------------------------------
save_gt_png <- function(gt_obj, basename, width = 1200, height = 650) {
  gtsave(gt_obj, file.path(out_dir, paste0(basename, ".png")),
         vwidth = width, vheight = height)
}
fmt_p <- function(x) ifelse(x < 1e-3, "< 0.001", sprintf("%.3f", x))
acf_df <- function(x, lag.max = 40) {
  a <- acf(x, lag.max = lag.max, plot = FALSE)
  tibble(lag = 1:lag.max, acf = as.numeric(a$acf)[-1])
}
get_ic <- function(fit, n) {
  ll <- tryCatch(as.numeric(fit@fit$LLH), error = function(e) NA_real_)
  k  <- tryCatch(length(coef(fit)),       error = function(e) NA_integer_)
  if (is.na(ll) || is.na(k)) return(c(AIC = NA_real_, BIC = NA_real_))
  c(AIC = -2*ll + 2*k, BIC = -2*ll + k*log(n))
}
checkmark <- function(x) ifelse(isTRUE(x), "\u2713", "\u2717")  # ✓ / ✗
diag_tables <- function(std_resid, lags_arch = 12, lag_lb = 20) {
  lb_r   <- Box.test(std_resid,   lag = lag_lb, type = "Ljung-Box")
  lb_r2  <- Box.test(std_resid^2, lag = lag_lb, type = "Ljung-Box")
  arch_r <- FinTS::ArchTest(std_resid, lags = lags_arch)
  tibble(
    Test      = c("Ljung–Box on standardized residuals",
                  "Ljung–Box on squared standardized residuals",
                  paste0("ARCH–LM on standardized residuals (lags = ", lags_arch, ")")),
    Statistic = c(unname(lb_r$statistic), unname(lb_r2$statistic), unname(arch_r$statistic)),
    `df / Lags` = c(unname(lb_r$parameter), unname(lb_r2$parameter), lags_arch),
    `p-value` = c(lb_r$p.value, lb_r2$p.value, arch_r$p.value)
  )
}

# =============================================================================
# 1) ARCH DIAGNOSTICS
# =============================================================================

# McLeod–Li via Ljung–Box on r^2 for lags 1..L
max_lag <- 20
ml_tbl <- bind_rows(lapply(1:max_lag, function(L) {
  bt <- Box.test(r^2, lag = L, type = "Ljung-Box")
  tibble(Lag = L, Q = unname(bt$statistic), p = bt$p.value)
}))
gt_ml <- ml_tbl |>
  mutate(Q = round(Q, 3), `p-value` = fmt_p(p)) |>
  rename(`Q (LB on r²)` = Q) |>
  gt() |>
  tab_header(
    title    = md("**McLeod–Li Portmanteau Test (Squared Returns)**"),
    subtitle = "H₀: no ARCH effects (no autocorrelation in r² up to lag L)"
  ) |>
  opt_row_striping() |>
  cols_width(Lag ~ px(90), `Q (LB on r²)` ~ px(180), `p-value` ~ px(140))
save_gt_png(gt_ml, "01_mcleod_li_arch_test")

# Engle’s ARCH–LM at lags 6/12/24
arch_lags <- c(6, 12, 24)
arch_rows <- map_dfr(arch_lags, function(L) {
  a <- FinTS::ArchTest(r, lags = L)
  tibble(Lags = L, `F-Statistic` = unname(a$statistic), `p-value` = fmt_p(a$p.value))
})
gt_arch <- arch_rows |>
  mutate(`F-Statistic` = round(`F-Statistic`, 3)) |>
  gt() |>
  tab_header(
    title    = md("**Engle ARCH–LM Test (Returns)**"),
    subtitle = "H₀: no ARCH effects"
  ) |>
  opt_row_striping() |>
  cols_width(Lags ~ px(90), `F-Statistic` ~ px(160), `p-value` ~ px(140))
save_gt_png(gt_arch, "02_arch_lm_test")

# ACF(|r|) and ACF(r^2) — polished ggplots
ci <- 1.96 / sqrt(n_obs)

p_acf_abs <- acf_df(abs(r), 40) |>
  ggplot(aes(lag, acf)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = -ci, ymax = ci), alpha = 0.08) +
  geom_hline(yintercept = c(-ci, ci), linetype = "dashed") +
  geom_col(width = 0.8) +
  labs(title = "ACF of Absolute Returns (|r|) — Daily Gold Log Returns",
       x = "Lag (trading days)", y = "ACF") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(out_dir, "03_acf_absr.png"), p_acf_abs, width = 10, height = 7, dpi = 300)

p_acf_r2 <- acf_df(r^2, 40) |>
  ggplot(aes(lag, acf)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = -ci, ymax = ci), alpha = 0.08) +
  geom_hline(yintercept = c(-ci, ci), linetype = "dashed") +
  geom_col(width = 0.8) +
  labs(title = expression(paste("ACF of Absolute Returns Squared (", r^2, ") — Daily Gold Log Returns")),
       x = "Lag (trading days)", y = "ACF") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(out_dir, "04_acf_rsq.png"), p_acf_r2, width = 10, height = 7, dpi = 300)

# =============================================================================
# 2) FIT CANDIDATE VARIANCE MODELS + IC
# =============================================================================

spec_s11_norm <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = keep_mu),
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  distribution.model = "norm"
)
spec_s11_t <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = keep_mu),
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)
spec_gjr_t <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = keep_mu),
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)
spec_egarch_t <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = keep_mu),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)

fit_safe <- function(spec, x) tryCatch(ugarchfit(spec, x), error = function(e) NULL)
fit_list <- list(
  "sGARCH(1,1) — Normal"    = fit_safe(spec_s11_norm, r),
  "sGARCH(1,1) — Student-t" = fit_safe(spec_s11_t,    r),
  "GJR(1,1) — Student-t"    = fit_safe(spec_gjr_t,    r),
  "EGARCH(1,1) — Student-t" = fit_safe(spec_egarch_t, r)
)

ic_rows <- imap_dfr(fit_list, function(fit, lab) {
  if (is.null(fit)) tibble(Model = lab, AIC = NA_real_, BIC = NA_real_, Converged = FALSE)
  else {
    ic <- get_ic(fit, n_obs)
    tibble(Model = lab, AIC = ic["AIC"], BIC = ic["BIC"], Converged = (fit@fit$convergence == 0))
  }
})

gt_ic <- ic_rows |>
  mutate(AIC = round(AIC, 2), BIC = round(BIC, 2), Converged = checkmark(Converged)) |>
  gt(rowname_col = "Model") |>
  tab_header(
    title    = md("**Variance Model Comparison — GARCH Family**"),
    subtitle = paste0("N = ", n_obs, " | Mean included: ",
                      ifelse(keep_mu, "Yes (t-test p < 0.05)", "No (demeaned)"))
  ) |>
  opt_row_striping() |>
  cols_width(AIC ~ px(120), BIC ~ px(120), Converged ~ px(120)) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_body(columns = AIC, rows = AIC == min(AIC, na.rm = TRUE)),
      cells_body(columns = BIC, rows = BIC == min(BIC, na.rm = TRUE))
    )
  )
save_gt_png(gt_ic, "05_variance_model_ic")

valid_idx <- which(ic_rows$Converged == TRUE & !is.na(ic_rows$AIC) & is.finite(ic_rows$AIC))
best_idx  <- if (length(valid_idx)) valid_idx[which.min(ic_rows$AIC[valid_idx])] else which(names(fit_list)=="sGARCH(1,1) — Student-t")
best_lab  <- ic_rows$Model[best_idx]
best_fit  <- fit_list[[best_lab]]

# Persistence & half-life (when alpha/beta exist)
cf <- tryCatch(coef(best_fit), error = function(e) NULL)
alpha <- tryCatch(unname(cf["alpha1"]), error = function(e) NA_real_)
beta  <- tryCatch(unname(cf["beta1"]),  error = function(e) NA_real_)
persist <- if (is.na(alpha) || is.na(beta)) NA_real_ else alpha + beta
half_life <- if (!is.na(persist) && persist > 0 && persist < 1) log(0.5)/log(persist) else NA_real_

subtitle_persist <- if (is.na(persist)) {
  "EGARCH/GJR do not always yield α, β; persistence interpreted via their own parameters."
} else if (persist >= 1) {
  "α+β ≥ 1 ⇒ near-unit persistence; EGARCH can remain stationary in log-variance even then."
} else {
  "Persistence near 1 ⇒ long memory; half-life shows decay speed."
}

gt_persist <- tibble(
  Model = best_lab, alpha = alpha, beta = beta,
  `alpha + beta` = persist, `Half-life (days)` = half_life
) |>
  mutate(across(where(is.numeric), ~round(.x, 4)),
         `Half-life (days)` = ifelse(is.na(`Half-life (days)`), NA, round(`Half-life (days)`, 2))) |>
  gt() |>
  tab_header(title = md("**Volatility Persistence — Selected Model**"),
             subtitle = subtitle_persist) |>
  opt_row_striping() |>
  cols_width(everything() ~ px(180))
save_gt_png(gt_persist, "06_persistence_table")

# =============================================================================
# 3) RESIDUAL DIAGNOSTICS + SIGN-BIAS + SIGMA PLOT (with dates)
# =============================================================================

res_std <- residuals(best_fit, standardize = TRUE)
gt_diag <- diag_tables(res_std, lags_arch = 12, lag_lb = 20) |>
  mutate(Statistic = round(Statistic, 3), `p-value` = fmt_p(`p-value`)) |>
  gt() |>
  tab_header(
    title    = md("**Residual Diagnostics — Selected Model**"),
    subtitle = "Large p-values are desirable (no remaining structure)."
  ) |>
  opt_row_striping() |>
  cols_width(Test ~ px(520), `df / Lags` ~ px(120), Statistic ~ px(140), `p-value` ~ px(140))
save_gt_png(gt_diag, "07_residual_diagnostics")

sb_raw <- tryCatch(rugarch::signbias(best_fit), error = function(e) NULL)
if (!is.null(sb_raw)) {
  sb_df <- as.data.frame(sb_raw) |> tibble::rownames_to_column("Component")
  nm <- tolower(gsub("[^a-z]", "", names(sb_df)))
  idx_t <- if (any(nm %in% c("tstat","tstatistics","tvalue","tvalues"))) which(nm %in% c("tstat","tstatistics","tvalue","tvalues")) else grep("^t", nm)
  idx_p <- if (any(nm %in% c("prob","pvalue","pvalues","pval","p"))) which(nm %in% c("prob","pvalue","pvalues","pval","p")) else grep("^p", nm)
  if (length(idx_t) && length(idx_p)) {
    names(sb_df)[idx_t[1]] <- "t_stat"; names(sb_df)[idx_p[1]] <- "p_value"
    gt_sb <- sb_df |>
      select(Component, t_stat, p_value) |>
      mutate(t_stat = round(as.numeric(t_stat), 3),
             p_value = fmt_p(as.numeric(p_value))) |>
      gt() |>
      tab_header(
        title    = md("**Sign-Bias (Asymmetry) Test — Selected Model**"),
        subtitle = "Significant negative-bias ⇒ negative shocks raise volatility more."
      ) |>
      opt_row_striping() |>
      cols_width(Component ~ px(260), t_stat ~ px(140), p_value ~ px(140)) |>
      cols_label(t_stat = "t-stat", p_value = "p-value")
    save_gt_png(gt_sb, "08_sign_bias_test")
  }
}

# Sigma vs DATE
sig <- sigma(best_fit)
df_sig <- tibble(date = dates, sigma = as.numeric(sig))
p_sig <- ggplot(df_sig, aes(date, sigma)) +
  geom_line(linewidth = 0.4) +
  labs(title = expression(paste("Estimated Conditional Volatility (", sigma[t], ")")),
       subtitle = best_lab, x = "Date", y = expression(sigma[t])) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(out_dir, "11_sigma_series.png"), p_sig, width = 11, height = 5.8, dpi = 300)

# Save model summary
capture.output(show(best_fit), file = file.path(out_dir, "10_fit_summary.txt"))

message("Done. PNGs and TXT saved to: ", out_dir)
print(dir_ls(out_dir))
