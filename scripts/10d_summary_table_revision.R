# --------------------------------------------
# 07_garch_fit_summary_gt.R  (updated)
# --------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(gt); library(fs)
})

# ---------- Paths ----------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/04_models/01_arch_garch"
dir_create(out_dir, recurse = TRUE)

# ---------- Coefficients (from your summary) ----------
params <- tibble::tribble(
  ~Parameter, ~Estimate, ~`Robust SE`, ~`t value`, ~`p-value`,
  "μ",          0.000198, 0.000066,   3.0295, 0.00245,
  "ω",         -0.065860, 0.002070, -31.8230, 0.00000,
  "α₁",         0.036760, 0.006525,   5.6335, 0.00000,
  "β₁",         0.992998, 0.000182, 5470.7789, 0.00000,
  "γ₁",         0.115676, 0.003996,  28.9471, 0.00000,
  "ν (shape)",  4.394630, 0.232122,  18.9324, 0.00000
)

# ---------- gt table ----------
tbl_params <- params |>
  mutate(
    Estimate  = round(Estimate, 4),
    `Robust SE` = round(`Robust SE`, 4),
    `t value` = round(`t value`, 2)        # <-- round t-values to 2 dp
  ) |>
  gt() |>
  tab_header(
    title = md("**EGARCH(1,1)-t Model — Estimated Parameters**"),
    subtitle = md("*Robust standard errors and significance*")
  ) |>
  cols_label(
    Parameter = md("**Parameter**"),
    Estimate  = md("**Estimate**"),
    `Robust SE` = md("**Robust SE**"),
    `t value` = md("**t-value**"),
    `p-value`  = md("**p-value**")
  ) |>
  fmt_number(columns = c(Estimate, `Robust SE`), decimals = 4) |>
  fmt_number(columns = `t value`, decimals = 2) |>
  fmt_number(columns = `p-value`, decimals = 5) |>
  data_color(
    columns = `p-value`,
    colors = scales::col_bin(
      palette = c("#D32F2F", "#FBC02D", "#388E3C"),
      domain  = c(0, 0.1),
      bins    = c(0, 0.01, 0.05, 0.1)
    )
  ) |>
  opt_row_striping()

# ---------- Save ----------
gtsave(tbl_params,
       file.path(out_dir, "06_fit_summary_params.png"),
       vwidth = 1400, vheight = 600)
