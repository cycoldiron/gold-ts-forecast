# ================================================================
# EGARCH(1,1)–t: Baseline vs Break-Aware Variants (OOS variance forecasts)
# Variants (all are EGARCH(1,1) with Student-t errors):
#   1) Baseline (single recursive roll)
#   2) Break-Aware: Regime-refit, short (PELT segments; many small refits)
#   3) Break-Aware: Regime-refit, merged ≥500 obs (stabilized refits)
#   4) Variance-Dummy: single roll with break dummy in variance equation
#
# Scoring: RMSE/MAE on e_t = r_t^2 - sigma_hat_t^2, plus QLIKE (lower is better)
# DM test: QLIKE (H1: Alternative < Baseline) using HAC variance of loss diff
# Output: PNG + HTML GT tables, tracker CSVs
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(xts); library(zoo)
  library(changepoint)
  library(rugarch)
  library(gt); library(fs); library(scales)
  library(tidyr); library(purrr); library(readr); library(stringr)
})

# -------------------------------
# 0) User options
# -------------------------------
use_parallel       <- FALSE   # placeholder; not used here
out_dir            <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/05_forecasts"

# Baseline rolling setup
refit_every_base   <- 500L
refit_window_base  <- "recursive"

# Regime rolling (Option A: short-regime)
refit_every_reg    <- 250L
refit_window_reg   <- "moving"
warmup_frac        <- 0.35
warmup_min         <- 50L

# Merged-regime minimum length (in observations)
min_reg_len        <- 500L

# Solver / fit controls (more robust)
solver_ctrl <- list(trace = 0, rel.tol = 1e-8, x.tol = 1e-8, eval.max = 5000, iter.max = 5000)
fit_ctrl    <- list(stationarity = 1, fixed.se = 0)

# -------------------------------
# 1) Data (expects gold_full_cleaned with date, close_log_return)
# -------------------------------
stopifnot(exists("gold_full_cleaned"))
df <- gold_full_cleaned %>%
  arrange(date) %>%
  filter(!is.na(close_log_return)) %>%
  transmute(Date = as.Date(date), r = close_log_return)

stopifnot(nrow(df) > 500)
r_xts <- xts(df$r, order.by = df$Date)

# Safe slicer: subset xts by observed dates (no "start/end" string)
slice_xts_range <- function(x, start_date, end_date) {
  idx <- zoo::index(x)
  x[idx >= as.Date(start_date) & idx <= as.Date(end_date)]
}

# -------------------------------
# 2) Variance changepoints (PELT) -> regimes
# -------------------------------
cp_idx <- tryCatch({
  fit_cp <- cpt.var(df$r, method = "PELT")
  cpts(fit_cp)
}, error = function(e) integer(0))

break_dates <- if (length(cp_idx)) df$Date[cp_idx] else as.Date(character(0))
message("Variance break count: ", length(break_dates))
if (length(break_dates)) message("Breaks at: ", paste(format(break_dates), collapse = ", "))

make_regime_id <- function(dates, breaks) {
  if (!length(breaks)) return(rep(1L, length(dates)))
  cut(dates,
      breaks = c(min(dates) - 1, sort(unique(breaks)), max(dates)),
      labels = FALSE, include.lowest = TRUE)
}

reg_idx <- tibble(
  date      = as.Date(index(r_xts)),
  regime_id = make_regime_id(as.Date(index(r_xts)), break_dates)
)

regimes <- reg_idx %>%
  group_by(regime_id) %>%
  summarise(start = min(date), end = max(date), .groups = "drop") %>%
  arrange(regime_id)

# Sanity
stopifnot(all(regimes$start <= regimes$end),
          all(regimes$start >= min(df$Date)),
          all(regimes$end   <= max(df$Date)))

# Regime lengths (obs) — useful to print/inspect
regime_lengths <- regimes %>%
  rowwise() %>%
  mutate(n = NROW(slice_xts_range(r_xts, start, end))) %>%
  ungroup()

# -------------------------------
# 3) EGARCH(1,1)–t spec
# -------------------------------
spec_eg <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)

# -------------------------------
# 4) Robust OOS date mapping helpers (n.ahead = 1)
# -------------------------------
oos_dates_simple <- function(roll_obj, r_idx) {
  dens <- roll_obj@forecast$density
  if (is.null(dens) || NROW(dens) == 0L) return(as.Date(character()))
  as.Date(tail(r_idx, NROW(dens)))
}

extract_oos_full <- function(roll_obj, r_xts) {
  dens  <- roll_obj@forecast$density
  stopifnot(!is.null(dens), NROW(dens) > 0L)
  dates <- oos_dates_simple(roll_obj, index(r_xts))
  tibble(
    date = dates,
    r    = as.numeric(dens[, "Realized"]),
    r2   = r^2,
    sig2 = as.numeric(dens[, "Sigma"])^2
  )
}

# -------------------------------
# 5) Baseline: single recursive roll
# -------------------------------
n_start_base <- round(0.35 * NROW(r_xts))  # warmup = 35% of full sample
message("\n[Baseline] Rolling…")
roll_baseline <- ugarchroll(
  spec_eg, data = r_xts,
  n.ahead = 1, n.start = n_start_base,
  refit.every = refit_every_base, refit.window = refit_window_base,
  solver = "nlminb", solver.control = solver_ctrl,
  fit.control = fit_ctrl,
  calculate.VaR = FALSE, keep.coef = TRUE
)
OOS_base <- extract_oos_full(roll_baseline, r_xts)
message("Baseline OOS window: ", min(OOS_base$date, na.rm=TRUE), " – ", max(OOS_base$date, na.rm=TRUE))
message("Baseline OOS n = ", nrow(OOS_base))

# -------------------------------
# 6) Break-Aware (SHORT-REGIME): per-regime rolls and splice OOS
# -------------------------------
safe_refit_every <- function(n_seg, n_start_seg, target_every) {
  rem <- n_seg - n_start_seg
  if (rem <= 1) return(NA_integer_)
  max(1L, min(target_every, rem - 1L))
}

roll_one_regime <- function(start_date, end_date, regime_id, refit_window = "moving") {
  seg_xts <- slice_xts_range(r_xts, start_date, end_date)
  n_seg   <- NROW(seg_xts)
  if (n_seg <= 2) return(NULL)
  
  n_start_seg <- max(warmup_min, round(warmup_frac * n_seg))
  if (n_seg <= (n_start_seg + 1)) return(NULL)
  
  refit_every_local <- safe_refit_every(n_seg, n_start_seg, refit_every_reg)
  if (is.na(refit_every_local)) return(NULL)
  
  message(sprintf("Regime %d | %s..%s | n=%d | warmup=%d | refit.every=%d",
                  regime_id, start_date, end_date, n_seg, n_start_seg, refit_every_local))
  
  out <- tryCatch({
    roll_seg <- ugarchroll(
      spec_eg, data = seg_xts,
      n.ahead = 1, n.start = n_start_seg,
      refit.every = refit_every_local, refit.window = refit_window,
      solver = "nlminb", solver.control = solver_ctrl,
      fit.control = fit_ctrl,
      calculate.VaR = FALSE, keep.coef = TRUE
    )
    dens <- roll_seg@forecast$density
    if (is.null(dens) || NROW(dens) == 0L) return(NULL)
    
    tibble(
      date        = oos_dates_simple(roll_seg, index(seg_xts)),
      r           = as.numeric(dens[, "Realized"]),
      r2          = r^2,
      sig2        = as.numeric(dens[, "Sigma"])^2,
      regime_id   = regime_id,
      break_start = start_date,
      break_dummy = as.integer(regime_id > 1L)
    )
  }, error = function(e) {
    message(sprintf("Regime %d failed: %s", regime_id, conditionMessage(e)))
    NULL
  })
  out
}

message("\n[Break-Aware: short-regime] Rolling across regimes…")
OOS_regime_short <- pmap_dfr(
  list(regimes$start, regimes$end, regimes$regime_id),
  roll_one_regime
)
if (!nrow(OOS_regime_short)) stop("No OOS data produced in short-regime rolls.")
message("Break-aware (short) OOS window: ", min(OOS_regime_short$date), " – ", max(OOS_regime_short$date))
message("Break-aware (short) OOS n = ", nrow(OOS_regime_short))

# -------------------------------
# 7) Break-Aware (MERGED-REGIME): enforce minimum regime length
# -------------------------------
merge_regimes_minlen <- function(reg_tbl, min_len, r_xts) {
  reg2 <- reg_tbl %>%
    rowwise() %>%
    mutate(n = NROW(slice_xts_range(r_xts, start, end))) %>%
    ungroup()
  
  out <- list(); i <- 1L
  while (i <= nrow(reg2)) {
    j <- i; nsum <- reg2$n[i]
    start_i <- reg2$start[i]; end_j <- reg2$end[i]
    while (nsum < min_len && j < nrow(reg2)) {
      j <- j + 1L; nsum <- nsum + reg2$n[j]; end_j <- reg2$end[j]
    }
    out[[length(out) + 1L]] <- tibble(start = start_i, end = end_j)
    i <- j + 1L
  }
  bind_rows(out) %>% mutate(regime_id = row_number()) %>% select(regime_id, start, end)
}

regimes_merged <- merge_regimes_minlen(regimes, min_reg_len, r_xts)

message("\n[Merged-Regime] Rolling across merged regimes…")
OOS_regime_merged <- pmap_dfr(
  list(regimes_merged$start, regimes_merged$end, regimes_merged$regime_id),
  ~ roll_one_regime(..1, ..2, ..3, refit_window = "moving")
)
if (!nrow(OOS_regime_merged)) stop("No OOS data produced in merged-regime rolls.")
message("Break-aware (merged) OOS window: ", min(OOS_regime_merged$date), " – ", max(OOS_regime_merged$date))
message("Break-aware (merged) OOS n = ", nrow(OOS_regime_merged))

# -------------------------------
# 8) Variance-Dummy (single roll; external regressor in variance)
#     D_t = 1 post-first-break, else 0
# -------------------------------
break_dummy_series <- if (length(break_dates)) {
  as.integer(index(r_xts) >= min(break_dates))
} else {
  rep(0L, NROW(r_xts))
}
Xvar <- matrix(break_dummy_series, ncol = 1)

spec_eg_dummy <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1),
                        external.regressors = Xvar),
  distribution.model = "std"
)

message("\n[Variance-Dummy] Rolling single model with variance regressor…")
roll_dummy <- ugarchroll(
  spec_eg_dummy, data = r_xts,
  n.ahead = 1, n.start = n_start_base,
  refit.every = refit_every_base, refit.window = refit_window_base,
  solver = "nlminb", solver.control = solver_ctrl,
  fit.control = fit_ctrl,
  calculate.VaR = FALSE, keep.coef = TRUE
)
OOS_dummy <- extract_oos_full(roll_dummy, r_xts)
message("Variance-Dummy OOS window: ", min(OOS_dummy$date), " – ", max(OOS_dummy$date))
message("Variance-Dummy OOS n = ", nrow(OOS_dummy))

# -------------------------------
# 9) Align overlap & score
# -------------------------------
mse   <- function(x) mean(x^2, na.rm = TRUE)
mae   <- function(x) mean(abs(x),  na.rm = TRUE)
qlike <- function(s2_hat, r2) {
  s2 <- pmax(s2_hat, .Machine$double.eps)
  log(s2) + r2/s2
}

Z <- OOS_base %>%
  select(date, r2_base = r2, sig2_base = sig2) %>%
  inner_join(select(OOS_regime_short,  date, r2_short = r2, sig2_short = sig2), by = "date") %>%
  inner_join(select(OOS_regime_merged, date, r2_merge = r2, sig2_merge = sig2), by = "date") %>%
  inner_join(select(OOS_dummy,         date, r2_dummy = r2, sig2_dummy = sig2), by = "date")

stopifnot(nrow(Z) > 25)

Z <- Z %>%
  mutate(
    e_base   = r2_base  - sig2_base,
    e_short  = r2_short - sig2_short,
    e_merge  = r2_merge - sig2_merge,
    e_dummy  = r2_dummy - sig2_dummy,
    L_base   = qlike(sig2_base,  r2_base),
    L_short  = qlike(sig2_short, r2_short),
    L_merge  = qlike(sig2_merge, r2_merge),
    L_dummy  = qlike(sig2_dummy, r2_dummy)
  )

# Accuracy table (Option B names)
acc_tbl <- tibble(
  Model = c("EGARCH(1,1)–t — Baseline (single recursive roll)",
            "EGARCH(1,1)–t — Break-Aware (Regime-refit, short)",
            "EGARCH(1,1)–t — Break-Aware (Regime-refit, merged ≥500)",
            "EGARCH(1,1)–t — Variance-Dummy (single roll)"),
  RMSE  = c(sqrt(mse(Z$e_base)),  sqrt(mse(Z$e_short)), sqrt(mse(Z$e_merge)), sqrt(mse(Z$e_dummy))),
  MAE   = c(mae(Z$e_base),        mae(Z$e_short),       mae(Z$e_merge),       mae(Z$e_dummy)),
  QLIKE = c(mean(Z$L_base),        mean(Z$L_short),      mean(Z$L_merge),      mean(Z$L_dummy))
) %>%
  mutate(`ΔQLIKE vs Baseline` = QLIKE - QLIKE[Model == "EGARCH(1,1)–t — Baseline (single recursive roll)"]) %>%
  relocate(`ΔQLIKE vs Baseline`, .after = QLIKE)

# -------------------------------
# 10) Proper DM test on QLIKE loss (H1: Alt < Base)
# -------------------------------
dm_on_losses <- function(L_alt, L_base, alternative = c("less","two.sided","greater")) {
  alternative <- match.arg(alternative)
  d <- L_alt - L_base
  d <- d[is.finite(d)]
  Tn <- length(d)
  if (Tn < 25 || !is.finite(var(d)) || var(d) <= .Machine$double.eps) {
    return(list(stat = 0, p = 1, mean_diff = mean(d)))
  }
  m <- mean(d)
  dc <- d - m
  # Newey–West bandwidth
  h <- floor(Tn^(1/3))
  # autocovariance at lag k
  acvf <- function(x, k) mean(x[1:(length(x)-k)] * x[(1+k):length(x)])
  gamma0 <- acvf(dc, 0)
  S <- gamma0
  if (h >= 1) {
    for (k in 1:h) {
      w <- 1 - k/(h+1)
      S <- S + 2 * w * acvf(dc, k)
    }
  }
  var_mean <- S / Tn
  tstat <- m / sqrt(var_mean)
  pval <- switch(alternative,
                 less      = pnorm(tstat),               # H1: mean(d) < 0
                 greater   = 1 - pnorm(tstat),
                 two.sided = 2 * min(pnorm(tstat), 1 - pnorm(tstat)))
  list(stat = as.numeric(tstat), p = as.numeric(pval), mean_diff = m)
}

dm_df <- bind_rows(
  {
    res <- dm_on_losses(Z$L_short, Z$L_base, "less")
    tibble(`Alt vs Baseline` = "Short-Regime vs Baseline (QLIKE)",
           `DM stat` = res$stat, `p-value` = res$p,
           `Mean ΔQLIKE (Alt−Base)` = res$mean_diff,
           `Decision (5%)` = ifelse(res$p < 0.05 & res$mean_diff < 0, "Alt better",
                                    ifelse(res$p < 0.05 & res$mean_diff > 0, "Baseline better", "No diff")))
  },
  {
    res <- dm_on_losses(Z$L_merge, Z$L_base, "less")
    tibble(`Alt vs Baseline` = "Merged-Regime vs Baseline (QLIKE)",
           `DM stat` = res$stat, `p-value` = res$p,
           `Mean ΔQLIKE (Alt−Base)` = res$mean_diff,
           `Decision (5%)` = ifelse(res$p < 0.05 & res$mean_diff < 0, "Alt better",
                                    ifelse(res$p < 0.05 & res$mean_diff > 0, "Baseline better", "No diff")))
  },
  {
    res <- dm_on_losses(Z$L_dummy, Z$L_base, "less")
    tibble(`Alt vs Baseline` = "Variance-Dummy vs Baseline (QLIKE)",
           `DM stat` = res$stat, `p-value` = res$p,
           `Mean ΔQLIKE (Alt−Base)` = res$mean_diff,
           `Decision (5%)` = ifelse(res$p < 0.05 & res$mean_diff < 0, "Alt better",
                                    ifelse(res$p < 0.05 & res$mean_diff > 0, "Baseline better", "No diff")))
  }
)

# -------------------------------
# 11) Save trackers & pretty tables
# -------------------------------
dir_create(out_dir, recurse = TRUE)

# Tracker CSVs for audit (per-regime variants)
tracker_short  <- file.path(out_dir, "break_aware_oos_tracker_short.csv")
tracker_merged <- file.path(out_dir, "break_aware_oos_tracker_merged.csv")
write_csv(select(OOS_regime_short,  date, regime_id, break_start, r2, sig2), tracker_short)
write_csv(select(OOS_regime_merged, date, regime_id, break_start, r2, sig2), tracker_merged)

date_range <- paste0(min(Z$date), " – ", max(Z$date))
n_eval     <- nrow(Z)

# Accuracy table UI (Option B names already baked into acc_tbl)
acc_gt <- acc_tbl |>
  gt(rowname_col = "Model") |>
  fmt_number(columns = c(RMSE, MAE, QLIKE, `ΔQLIKE vs Baseline`),
             decimals = 6, drop_trailing_zeros = TRUE) |>
  tab_header(
    title = md("**OOS Variance Forecast Accuracy — EGARCH(1,1)–t (all variants)**"),
    subtitle = md(paste0("Evaluation window: ", date_range, " &nbsp;&nbsp;|&nbsp;&nbsp; n = ", n_eval))
  ) |>
  tab_source_note(md("Lower is better for all metrics. QLIKE is a log-loss and can be negative.")) |>
  tab_source_note(md("Break-aware (Regime-refit) models are per-regime **ugarchroll** re-estimates using PELT breaks.")) |>
  tab_source_note(md("Variance-Dummy uses a single **ugarchroll** with a break dummy in the **variance** equation.")) |>
  tab_options(table.font.size = px(13), data_row.padding = px(5))

# DM table UI
dm_gt <- dm_df |>
  gt() |>
  fmt_number(columns = c(`DM stat`, `p-value`, `Mean ΔQLIKE (Alt−Base)`), decimals = 4) |>
  tab_header(
    title = md("**Diebold–Mariano (QLIKE): Alternatives vs Baseline**"),
    subtitle = md("H1: Alternative has **lower** expected QLIKE loss than Baseline")
  ) |>
  tab_options(table.font.size = px(13), data_row.padding = px(5))

# Safer PNG saves (falls back to HTML if PNG hiccups)
save_gt_pair <- function(gt_obj, stem, out_dir) {
  png_path  <- file.path(out_dir, paste0(stem, ".png"))
  html_path <- file.path(out_dir, paste0(stem, ".html"))
  ok <- TRUE
  tryCatch({ gtsave(gt_obj, png_path) }, error = function(e) { ok <<- FALSE })
  if (!ok) message("[WARN] PNG save failed for ", stem, " — saving HTML only.")
  gtsave(gt_obj, html_path)
  list(png = png_path, html = html_path)
}

acc_paths <- save_gt_pair(acc_gt, "oos_variance_accuracy_all", out_dir)
dm_paths  <- save_gt_pair(dm_gt,  "dm_qlike_all", out_dir)

message("Saved tables:\n  ", acc_paths$html, "\n  ", dm_paths$html,
        "\nTracker CSVs:\n  ", tracker_short, "\n  ", tracker_merged)

# -------------------------------
# 12) Console summary
# -------------------------------
message("\n=== Overlap window: ", min(Z$date), " – ", max(Z$date), " | n = ", n_eval, " ===")
print(acc_tbl)
print(dm_df)

# -------------------------------
# (Optional) Cache snapshot to skip modeling next time
# -------------------------------
# dir_create(file.path(out_dir, "cache"), recurse = TRUE)
# saveRDS(list(
#   OOS_base = OOS_base,
#   OOS_regime_short = OOS_regime_short,
#   OOS_regime_merged = OOS_regime_merged,
#   OOS_dummy = OOS_dummy,
#   Z = Z,
#   acc_tbl = acc_tbl
# ), file = file.path(out_dir, "cache/oos_snap.rds"))
