# 07_egarch11_diagnostics_table_COMPACT.R
# Compact PNG table for EGARCH(1,1) diagnostics, with regime highlighting
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(gt); library(FinTS); library(rugarch); library(fs); library(scales)
})

# -------- Paths --------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/04_models/01_arch_garch"
dir_create(out_dir, recurse = TRUE)

# -------- Data --------
df <- gold_full_cleaned |> select(date, r = close_log_return) |> tidyr::drop_na()
r  <- df$r
keep_mu <- (t.test(r)$p.value < 0.05)

# -------- EGARCH(1,1) fit (Student-t) --------
spec_e11_t <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = keep_mu),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)
fit_safe <- function(spec, x) tryCatch(ugarchfit(spec, x), error = function(e) NULL)

fit_full <- fit_safe(spec_e11_t, r)
stopifnot(!is.null(fit_full), fit_full@fit$convergence == 0)
res_full <- residuals(fit_full, standardize = TRUE)

post_2010 <- df |> filter(date >= as.Date("2010-01-01"))
post_2018 <- df |> filter(date >= as.Date("2018-01-01"))

fit_2010 <- fit_safe(spec_e11_t, post_2010$r)
fit_2018 <- fit_safe(spec_e11_t, post_2018$r)

res_2010 <- if (!is.null(fit_2010) && fit_2010@fit$convergence==0) residuals(fit_2010, standardize=TRUE) else NULL
res_2018 <- if (!is.null(fit_2018) && fit_2018@fit$convergence==0) residuals(fit_2018, standardize=TRUE) else NULL

# -------- Helpers --------
acf_max5 <- function(x) {
  a <- acf(x^2, lag.max = 5, plot = FALSE)
  round(max(abs(as.numeric(a$acf)[-1])), 3)
}

# Keep ONLY: LB on squared residuals + ARCH–LM
diag_block_compact <- function(std_resid, label, L_lb = 10, L_arch = 12) {
  lb_r2  <- Box.test(std_resid^2, lag = L_lb, type = "Ljung-Box")
  arch   <- FinTS::ArchTest(std_resid, lags = L_arch)
  mx_acf <- acf_max5(std_resid)
  
  tibble(
    Sample = label,
    Test   = c("LB on squared standardized residuals",
               "ARCH–LM on standardized residuals"),
    `df / Lags` = c(L_lb, L_arch),
    Statistic   = round(c(unname(lb_r2$statistic),
                          unname(arch$statistic)), 3),
    `p-value`   = c(lb_r2$p.value, arch$p.value),  # numeric for coloring
    `max |ACF(r²)| (1–5)` = c(mx_acf, NA_real_)    # once per sample
  )
}

diag_all <- bind_rows(
  diag_block_compact(res_full, "Full sample"),
  if (!is.null(res_2010)) diag_block_compact(res_2010, "Post-2010"),
  if (!is.null(res_2018)) diag_block_compact(res_2018, "Post-2018")
)

# -------- Compact GT table --------
gt_tbl <- diag_all |>
  gt(groupname_col = "Sample", row_group_as_column = FALSE) |>
  tab_header(
    title = md("**Residual Diagnostics — EGARCH(1,1), Student-t**")
  ) |>
  row_group_order(c("Full sample", "Post-2010", "Post-2018")) |>
  cols_width(
    Test ~ px(360),
    `df / Lags` ~ px(95),
    Statistic ~ px(110),
    `p-value` ~ px(95),
    `max |ACF(r²)| (1–5)` ~ px(160)
  ) |>
  cols_align(
    align = "center",
    columns = c(`df / Lags`, Statistic, `p-value`, `max |ACF(r²)| (1–5)`)
  ) |>
  # Color p-values: red (reject) for small, green (fail to reject) for large
  data_color(
    columns = `p-value`,
    colors = col_bin(
      palette = c("#e74c3c", "#2ecc71"),
      bins = c(-Inf, 0.05, Inf),
      right = FALSE
    ),
    apply_to = "fill"
  ) |>
  # Format p-values (no bold)
  fmt(
    columns = `p-value`,
    fns = function(x) ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))
  ) |>
  # Hide missing values as blanks
  sub_missing(
    columns = `max |ACF(r²)| (1–5)`,
    missing_text = ""
  ) |>
  opt_row_striping() |>
  tab_options(
    table.font.size   = px(12),
    data_row.padding  = px(3),
    heading.padding   = px(6),
    column_labels.padding = px(4)
  ) |>
  # Regime highlighting (row-group headers)
  tab_style(
    style = list(cell_fill(color = "#f2f2f2"), cell_text(weight = "bold")),
    locations = cells_row_groups(groups = "Full sample")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#f2f2f2"), cell_text(weight = "bold")),
    locations = cells_row_groups(groups = "Post-2010")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#f2f2f2"), cell_text(weight = "bold")),
    locations = cells_row_groups(groups = "Post-2018")
  )

gtsave(
  gt_tbl,
  filename = file.path(out_dir, "egarch11_residual_diagnostics_COMPACT.png"),
  vwidth = 980, vheight = 620
)

message("Saved: ", file.path(out_dir, "egarch11_residual_diagnostics_COMPACT.png"))
