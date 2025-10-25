# 06_garch_model_comparison_PNG.R — one PNG with model AIC/BIC (BIC-best bold)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr)
  library(gt); library(rugarch); library(fs)
})

# ---------- Paths ----------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/04_models/01_arch_garch"
dir_create(out_dir, recurse = TRUE)

# ---------- Data ----------
r <- na.omit(gold_full_cleaned$close_log_return)
n_obs  <- length(r)
keep_mu <- (t.test(r)$p.value < 0.05)

save_gt_png <- function(gt_obj, basename, width = 1400, height = 520) {
  gtsave(gt_obj, file.path(out_dir, paste0(basename, ".png")),
         vwidth = width, vheight = height)
}

# ---------- Robust fitting (minimal) ----------
.fit <- function(spec, x, solver = "hybrid") {
  tryCatch(
    ugarchfit(spec, x, solver = solver,
              fit.control = list(scale = 1),
              solver.control = list(trace = 0)),
    error = function(e) NULL
  )
}
fit_safe <- function(spec, x) {
  f <- .fit(spec, x, "hybrid")
  if (is.null(f) || f@fit$convergence != 0 || length(f@fit$LLH) != 1) {
    xz <- as.numeric(scale(x))
    for (s in c("hybrid","gosolnp","nlminb")) {
      f <- .fit(spec, xz, s)
      if (!is.null(f) && f@fit$convergence == 0 && length(f@fit$LLH) == 1) break
    }
  }
  f
}

spec_build <- function(model, p, q, dist) {
  ugarchspec(
    mean.model     = list(armaOrder = c(0,0), include.mean = keep_mu),
    variance.model = list(model = model, garchOrder = c(p,q)),
    distribution.model = dist
  )
}

label_for <- function(model, p, q, dist) {
  base <- switch(
    model,
    "sGARCH"  = if (q == 0) paste0("ARCH(", p, ")") else paste0("GARCH(", p, ",", q, ")"),
    "gjrGARCH"= paste0("GJR(",  p, ",", q, ")"),
    "eGARCH"  = paste0("EGARCH(",p, ",", q, ")"),
    "iGARCH"  = paste0("IGARCH(",p, ",", q, ")"),
    model
  )
  paste0(
    base,
    if (dist == "norm") " — normal"
    else if (dist == "std") " — t"
    else if (dist == "sstd") " — skew-t" else ""
  )
}


# ---------- Candidate set (baseline-first) ----------
grid <- tribble(
  ~model,     ~p, ~q, ~dist,
  "sGARCH",    1,  0, "std",    # ARCH(1)
  "sGARCH",    2,  0, "std",    # ARCH(2)
  "sGARCH",    3,  0, "std",    # ARCH(3)
  "sGARCH",    1,  1, "std",    # GARCH(1,1) — t
  "sGARCH",    1,  1, "sstd",   # GARCH(1,1) — skew-t
  "sGARCH",    1,  2, "std",    # GARCH(1,2)
  "sGARCH",    2,  1, "std",    # GARCH(2,1)
  "iGARCH",    1,  1, "std",    # IGARCH(1,1)
  "gjrGARCH",  1,  1, "std",    # GJR(1,1) — t
  "gjrGARCH",  1,  1, "sstd",   # GJR(1,1) — skew-t
  "eGARCH",    1,  1, "std",    # EGARCH(1,1) — t
  "eGARCH",    1,  1, "sstd"    # EGARCH(1,1) — skew-t
)

# ---------- Fit -> one tibble row per spec ----------
fit_one <- function(model, p, q, dist) {
  lab  <- label_for(model, p, q, dist)
  fit  <- fit_safe(spec_build(model, p, q, dist), r)
  
  row <- tibble(Model = lab, AIC = NA_real_, BIC = NA_real_, Converged = FALSE)
  
  if (is.null(fit) || fit@fit$convergence != 0 || length(fit@fit$LLH) != 1) return(row)
  ll <- suppressWarnings(as.numeric(fit@fit$LLH))
  k  <- tryCatch(length(coef(fit)), error = function(e) NA_integer_)
  if (!is.finite(ll) || !is.finite(k)) return(row)
  
  row$Converged <- TRUE
  row$AIC <- -2*ll + 2*k
  row$BIC <- -2*ll + k*log(n_obs)
  row
}

rows_list <- pmap(grid, fit_one)
fit_rows  <- bind_rows(rows_list)

# ---------- Table (bold best BIC, highlight best row) ----------
best_idx <- which(is.finite(fit_rows$BIC))
best_idx <- if (length(best_idx)) best_idx[which.min(fit_rows$BIC[best_idx])] else NA_integer_

tbl <- fit_rows |>
  mutate(
    AIC = round(AIC, 2),
    BIC = round(BIC, 2),
    Converged = ifelse(Converged, "\u2713", "\u2717")
  )

# ---------- Table (bold best BIC, highlight best row) ----------
gt_tbl <- gt(tbl, rowname_col = "Model") |>
  tab_header(
    title = md("**Variance Model Comparison — Baselines & ARCH/GARCH Family**")
    # (optional title change done here; no subtitle)
  ) |>
  opt_row_striping() |>
  cols_width(AIC ~ px(140), BIC ~ px(140), Converged ~ px(110))


if (!is.na(best_idx)) {
  gt_tbl <- gt_tbl |>
    # bold the best BIC cell
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = BIC, rows = best_idx)
    ) |>
    # highlight the entire best row in yellow
    tab_style(
      style = cell_fill(color = "#FFF59D"),
      locations = cells_body(columns = everything(), rows = best_idx)
    )
}

save_gt_png(gt_tbl, "05_variance_model_ic")
message("Saved PNG: ", file.path(out_dir, "05_variance_model_ic.png"))
