# ================================================================
# EGARCH(1,1)–t (variance dummy) — Diagnostics (Greek labels)
# - Keeps your estimation exactly as-is
# - Aesthetic/labeling upgrades:
#     * Parameter table uses "Parameter" (capitalized) and Greek symbols
#     * Ljung–Box table: spec in title, columns aligned, green highlight for p > 0.05
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(xts); library(zoo)
  library(changepoint)
  library(rugarch)
  library(gt); library(fs); library(scales)
  library(readr)
})

# -------------------------------
# Paths
# -------------------------------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/06_variance_model_diagnostics"
dir_create(out_dir, recurse = TRUE)

# -------------------------------
# Data
# -------------------------------
stopifnot(exists("gold_full_cleaned"))
df <- gold_full_cleaned %>%
  arrange(date) %>%
  filter(!is.na(close_log_return)) %>%
  transmute(Date = as.Date(date), r = close_log_return)
stopifnot(nrow(df) > 500)
r_xts <- xts(df$r, order.by = df$Date)

# -------------------------------
# Variance dummy (after first variance break)
# -------------------------------
first_break <- tryCatch({
  cp <- changepoint::cpt.var(df$r, method = "PELT")
  b  <- changepoint::cpts(cp)
  if (length(b)) df$Date[min(b)] else NA
}, error = function(e) NA)

break_dummy <- if (is.na(first_break)) rep(0L, NROW(r_xts)) else as.integer(index(r_xts) >= first_break)
Xvar <- matrix(break_dummy, ncol = 1)

# -------------------------------
# EGARCH(1,1)–t spec and fit
# -------------------------------
spec_vd <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1),
                        external.regressors = Xvar),
  distribution.model = "std"
)

solver_ctrl <- list(trace = 0, rel.tol = 1e-8, x.tol = 1e-8, eval.max = 5000, iter.max = 5000)
fit_ctrl    <- list(stationarity = 1, fixed.se = 0)

fit_full <- ugarchfit(
  spec = spec_vd, data = r_xts,
  solver = "nlminb",
  solver.control = solver_ctrl,
  fit.control = fit_ctrl
)

# -------------------------------
# Helpers
# -------------------------------
round2 <- function(x) if (is.numeric(x)) round(x, 2) else x

pretty_param <- function(x) {
  dplyr::recode(x,
                "omega"  = "\u03C9",            # ω
                "alpha1" = "\u03B1",            # α
                "beta1"  = "\u03B2",            # β
                "gamma1" = "\u03B3",            # γ
                "vxreg1" = "\u03B4 (dummy)",    # δ (variance dummy)
                "shape"  = "\u03BD"             # ν (Student-t df)
  )
}

# -------------------------------
# Styled Parameter Table (EGARCH(1,1)–t, Variance Dummy)
# -------------------------------
get_param_table <- function(fit) {
  mc <- as.data.frame(fit@fit$matcoef)
  colnames(mc) <- c("Estimate","Std_Error","t_value","p_value")
  
  pretty_param <- function(x) dplyr::recode(
    x,
    "omega"  = "\u03C9",            # ω
    "alpha1" = "\u03B1",            # α
    "beta1"  = "\u03B2",            # β
    "gamma1" = "\u03B3",            # γ
    "vxreg1" = "\u03B4 (dummy)",    # δ
    "shape"  = "\u03BD"             # ν
  )
  
  tibble(Param_raw = rownames(mc)) |>
    dplyr::bind_cols(tibble::as_tibble(mc)) |>
    dplyr::mutate(
      Parameter = pretty_param(Param_raw),
      sig = p_value < 0.05
    ) |>
    dplyr::select(Parameter, Estimate, Std_Error, t_value, p_value, sig)
}

param_tbl <- get_param_table(fit_full)

param_gt <- param_tbl |>
  gt::gt(rowname_col = "Parameter") |>
  # nice human labels + 3dp for clearer tiny effects
  gt::fmt_number(columns = c(Estimate, Std_Error, t_value, p_value), decimals = 3) |>
  gt::cols_label(
    Estimate  = "Estimate",
    Std_Error = "Std. Error",
    t_value   = "t-value",
    p_value   = "p-value"
  ) |>
  # table header like your baseline example
  gt::tab_header(
    title    = gt::md("**EGARCH(1,1)–t (Variance dummy) — Estimated Parameters**"),
    subtitle = gt::md("Standard errors and significance")
  ) |>
  # center numeric, left-align parameter names
  gt::cols_align(columns = everything(), align = "center") |>
  gt::cols_align(columns = 1, align = "left") |>
  # p-value "badge": red for p<0.05, light gray otherwise
  gt::tab_style(
    style = list(gt::cell_fill(color = "#EA4335"), gt::cell_text(color = "white", weight = "bold")),
    locations = gt::cells_body(columns = "p_value", rows = sig == TRUE)
  ) |>
  gt::tab_style(
    style = list(gt::cell_fill(color = "#F1F3F4"), gt::cell_text(color = "black")),
    locations = gt::cells_body(columns = "p_value", rows = sig == FALSE)
  ) |>
  # row striping for readability
  gt::opt_row_striping() |>
  # slightly larger font + padding, bold column headers
  gt::tab_options(
    table.font.size = gt::px(13),
    data_row.padding = gt::px(6),
    column_labels.font.weight = "bold"
  ) |>
  # hide helper column
  gt::cols_hide("sig")

gtsave(param_gt, file.path(out_dir, "params_vardummy.html"))
try(gtsave(param_gt, file.path(out_dir, "params_vardummy.png")), silent = TRUE)

# -------------------------------
# Numeric checks (2 dp)
# -------------------------------
cf <- coef(fit_full)
val_beta1 <- unname(cf["beta1"])
val_shape <- unname(cf["shape"])

sanity_tbl <- tibble(
  Check     = c("EGARCH persistence |β| < 1", "Student-t ν > 2 (finite variance)"),
  Value     = c(round(abs(val_beta1), 2), round(val_shape, 2)),
  Threshold = c(1.00, 2.00),
  Status    = c(ifelse(abs(val_beta1) < 1, "OK", "Flag"),
                ifelse(val_shape > 2, "OK", "Flag"))
)

sanity_gt <- sanity_tbl %>%
  gt() %>%
  fmt_number(columns = c(Value, Threshold), decimals = 2) %>%
  tab_header(title = md("**Parameter diagnostics**")) %>%
  tab_options(table.font.size = px(13), data_row.padding = px(5))

gtsave(sanity_gt, file.path(out_dir, "sanity_numeric.html"))
try(gtsave(sanity_gt, file.path(out_dir, "sanity_numeric.png")), silent = TRUE)

# -------------------------------
# Sub-sample stability (Greek; 2 dp)
# -------------------------------
fit_subset <- function(frac) {
  n <- NROW(r_xts); n_cut <- max(100L, floor(frac * n))
  Xsub <- matrix(break_dummy[1:n_cut], ncol = 1)
  spec_sub <- ugarchspec(
    mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
    variance.model = list(model = "eGARCH", garchOrder = c(1,1),
                          external.regressors = Xsub),
    distribution.model = "std"
  )
  suppressWarnings(
    ugarchfit(spec = spec_sub, data = r_xts[1:n_cut],
              solver = "nlminb", solver.control = solver_ctrl, fit.control = fit_ctrl)
  )
}
fit_60 <- fit_subset(0.60)
fit_80 <- fit_subset(0.80)

extract_keys <- function(fit) {
  cf <- tryCatch(coef(fit), error = function(e) NULL)
  pick <- c("omega","alpha1","beta1","gamma1","vxreg1","shape")
  out <- setNames(rep(NA_real_, length(pick)), pick)
  if (!is.null(cf)) out[names(cf)[names(cf) %in% pick]] <- cf[names(cf) %in% pick]
  as.list(out)
}

stab_tbl_raw <- bind_rows(
  tibble(Sample = "Full",      !!!extract_keys(fit_full)),
  tibble(Sample = "First 60%", !!!extract_keys(fit_60)),
  tibble(Sample = "First 80%", !!!extract_keys(fit_80))
) %>% mutate(across(-Sample, round2))

stab_tbl <- stab_tbl_raw %>%
  rename(
    "\u03C9" = omega,
    "\u03B1" = alpha1,
    "\u03B2" = beta1,
    "\u03B3" = gamma1,
    "\u03B4 (dummy)" = vxreg1,
    "\u03BD" = shape
  )

stab_gt <- stab_tbl %>%
  gt() %>%
  fmt_number(columns = -Sample, decimals = 2) %>%
  tab_header(title = md("**Sub-sample stability (key parameters)**")) %>%
  tab_options(table.font.size = px(13), data_row.padding = px(5))

gtsave(stab_gt, file.path(out_dir, "subsample_stability.html"))
try(gtsave(stab_gt, file.path(out_dir, "subsample_stability.png")), silent = TRUE)

# -------------------------------
# Residual diagnostics prep
# -------------------------------
z_vec   <- as.numeric(residuals(fit_full, standardize = TRUE))
sig_vec <- as.numeric(sigma(fit_full))
idx     <- index(r_xts)
z_xts   <- xts(z_vec, order.by = idx)
sig_xts <- xts(sig_vec, order.by = idx)

# -------------------------------
# Ljung–Box table (spec in title; green if p > 0.05)
# -------------------------------
.lb_p <- function(x, lag) Box.test(x, lag = lag, type = "Ljung-Box")$p.value

lb_tbl <- tibble(
  Series   = c("Standardized residuals", "Squared residuals"),
  `Lag 10` = c(.lb_p(z_vec, 10), .lb_p(z_vec^2, 10)),
  `Lag 20` = c(.lb_p(z_vec, 20), .lb_p(z_vec^2, 20))
)

lb_gt <- lb_tbl %>%
  gt() %>%
  fmt_number(columns = c(`Lag 10`, `Lag 20`), decimals = 2) %>%
  tab_header(
    title    = md("**Ljung–Box tests — EGARCH(1,1)–t (Variance dummy)**"),
    subtitle = md("p-values (H\u2080: no serial correlation)")
  ) %>%
  cols_align(columns = everything(), align = "center") %>%
  cols_align(columns = Series, align = "left") %>%
  # highlight cells that FAIL to reject H0 at 5% (p > 0.05)
  tab_style(
    style = list(cell_fill(color = "#E6F4EA"), cell_text(color = "#137333", weight = "bold")),
    locations = cells_body(columns = `Lag 10`, rows = `Lag 10` > 0.05)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#E6F4EA"), cell_text(color = "#137333", weight = "bold")),
    locations = cells_body(columns = `Lag 20`, rows = `Lag 20` > 0.05)
  ) %>%
  tab_source_note(md("Green cells fail to reject H\u2080 at 5% (no evidence of autocorrelation).")) %>%
  tab_options(table.font.size = px(13), data_row.padding = px(6))

gtsave(lb_gt, file.path(out_dir, "ljung_box.html"))
try(gtsave(lb_gt, file.path(out_dir, "ljung_box.png")), silent = TRUE)

# -------------------------------
# Optional: diagnostics grid image (unchanged, just saving again)
# -------------------------------
png(file.path(out_dir, "diagnostics_grid.png"), width = 2000, height = 1600, res = 220)

# layout: first row is title (spans 3 cols), then 3x3 panels
mat <- matrix(c(1,1,1,
                2,3,4,
                5,6,7,
                8,9,10), nrow = 4, byrow = TRUE)
layout(mat, heights = c(0.12, 0.96, 0.96, 0.96))

# (1) Title panel
par(mar = c(0,0,0,0))
plot.new()
overall_title <- paste0(
  "Diagnostics — EGARCH(1,1)–t with variance dummy   |   ",
  format(min(index(r_xts)), "%Y-%m-%d"), " to ", format(max(index(r_xts)), "%Y-%m-%d")
)
text(x = 0.5, y = 0.5, labels = overall_title, cex = 1.1, font = 2)

# Common panel params
panel_par <- function() par(mar = c(3.4, 3.6, 3.2, 1.2), mgp = c(2.0, 0.6, 0), cex.main = 1.0, cex.lab = 0.95, cex.axis = 0.85)
plot_xts_line <- function(x, main) plot(zoo::as.zoo(x), type = "l", xlab = "", ylab = "", main = main)

# (2) Standardized residuals
panel_par(); plot_xts_line(z_xts, main = "Standardized residuals")

# (3) ACF: z
panel_par(); acf(z_vec, lag.max = 40, main = "ACF: standardized residuals")

# (4) ACF: z^2
panel_par(); acf(z_vec^2, lag.max = 40, main = "ACF: squared residuals")

# (5) Histogram
panel_par(); hist(z_vec, breaks = 60, main = "Histogram: standardized residuals", xlab = "")

# (6) QQ
panel_par(); qqnorm(z_vec, main = "Normal Q–Q plot"); qqline(z_vec)

# (7) σ_t
panel_par(); plot_xts_line(sig_xts, main = "Estimated conditional volatility (σ)")

# (8) |z|
panel_par(); plot_xts_line(abs(z_xts), main = "Absolute residuals")

# (9) z^2
panel_par(); plot_xts_line(z_xts^2, main = "Squared residuals")

# (10) cumsum z
panel_par(); plot_xts_line(xts(cumsum(z_vec), order.by = idx), main = "Cumulative standardized residuals")

dev.off()

message("Saved to: ", out_dir)
