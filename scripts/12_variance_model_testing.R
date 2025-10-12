# =============================================================================
# 06c_egarch_grid_and_compare.R
# Richer EGARCH grid, clean diagnostics (Student-t only), LR test, and sigma plots
# Outputs: PNG gt() tables + PNG plots to out_dir
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(gt); library(FinTS); library(rugarch); library(fs)
})

# ---------- Paths -------------------------------------------------------------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/05_egarch_grid"
dir_create(out_dir, recurse = TRUE)

# ---------- Data --------------------------------------------------------------
# expects: gold_full_cleaned$date, gold_full_cleaned$close_log_return
r_all <- gold_full_cleaned$close_log_return
dates <- gold_full_cleaned$date[!is.na(r_all)]
r     <- na.omit(r_all)
n_obs <- length(r)

# ---------- Helpers -----------------------------------------------------------
save_gt_png <- function(gt_obj, basename, width = 1200, height = 650) {
  gtsave(gt_obj, file.path(out_dir, paste0(basename, ".png")),
         vwidth = width, vheight = height)
}
fmt_p <- function(x) ifelse(is.na(x), NA_character_,
                            ifelse(x < 1e-3, "< 0.001", sprintf("%.3f", x)))
get_ic <- function(fit, n) {
  if (is.null(fit)) return(c(AIC = NA_real_, BIC = NA_real_))
  ll <- tryCatch(as.numeric(fit@fit$LLH), error = function(e) NA_real_)
  k  <- tryCatch(length(coef(fit)),       error = function(e) NA_integer_)
  if (is.na(ll) || is.na(k)) return(c(AIC = NA_real_, BIC = NA_real_))
  c(AIC = -2*ll + 2*k, BIC = -2*ll + k*log(n))
}
checkmark <- function(x) ifelse(isTRUE(x), "\u2713", "\u2717")  # ✓ / ✗

# mean decision (keep intercept if significant)
t_mu    <- t.test(r)
keep_mu <- (t_mu$p.value < 0.05)

# ---------- EGARCH spec + safe fit -------------------------------------------
spec_egarch <- function(p, q, dist = "std", include_mean = keep_mu) {
  ugarchspec(
    mean.model     = list(armaOrder = c(0,0), include.mean = include_mean),
    variance.model = list(model = "eGARCH", garchOrder = c(p,q)),
    distribution.model = dist
  )
}
fit_safe <- function(p, q, dist) {
  sp <- spec_egarch(p, q, dist)
  tryCatch(
    ugarchfit(sp, r, solver = "hybrid", solver.control = list(trace = 0)),
    error = function(e) NULL
  )
}

# ---------- EGARCH(p,q) GRID (IC table uses both t and skew-t) ----------------
pq <- expand_grid(p = 1:3, q = 1:3, dist = c("std","sstd"))
eg_fits <- pq %>%
  mutate(fit = pmap(list(p, q, dist), ~ fit_safe(..1, ..2, ..3))) %>%
  mutate(Converged = map_lgl(fit, ~ !is.null(.x) && .x@fit$convergence == 0),
         AIC = map_dbl(fit, ~ get_ic(.x, n_obs)["AIC"]),
         BIC = map_dbl(fit, ~ get_ic(.x, n_obs)["BIC"])) %>%
  mutate(Model = paste0("EGARCH(", p, ",", q, ")"),
         Innovations = ifelse(dist=="std","Student-t","Skew-t")) %>%
  select(Model, p, q, Innovations, fit, Converged, AIC, BIC) %>%
  filter(!is.na(AIC) & !is.na(BIC))  # drop any NA fits (non-converged)

# --- IC table: keep ONLY Student-t rows --------------------------------------
ic_tbl <- eg_fits %>%
  filter(Innovations == "Student-t") %>%            # <<< add this line
  transmute(Model,
            Innovations,
            Converged = checkmark(Converged),
            AIC_raw = AIC,
            BIC_raw = BIC,
            AIC = round(AIC, 2),
            BIC = round(BIC, 2))

minAIC <- min(ic_tbl$AIC_raw, na.rm = TRUE)
minBIC <- min(ic_tbl$BIC_raw, na.rm = TRUE)

gt_ic <- ic_tbl %>%
  select(Model, Innovations, Converged, AIC, BIC, AIC_raw, BIC_raw) %>%
  gt() %>%
  tab_header(
    title = md("**Variance Model Comparison — EGARCH(p,q) Grid (Student-t)**"),
    subtitle = paste0("N = ", n_obs, " | Mean included: ",
                      ifelse(keep_mu, "Yes (t-test p < 0.05)", "No (demeaned)"))
  ) %>%
  opt_row_striping() %>%
  cols_width(Model ~ px(200), Innovations ~ px(140),
             AIC ~ px(120), BIC ~ px(120), Converged ~ px(120)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = AIC, rows = AIC_raw == minAIC)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = BIC, rows = BIC_raw == minBIC)) %>%
  cols_hide(columns = c(AIC_raw, BIC_raw))

save_gt_png(gt_ic, "egarch_grid_ic")

# ---------- Diagnostics: (Student-t models ONLY) ------------------------------
diag_one <- function(f) {
  if (is.null(f)) return(tibble(Test=character(), p_num=double()))
  rs <- residuals(f, standardize = TRUE)
  tibble(
    Test  = c("LB(ε̂)", "LB(ε̂²)", "ARCH-LM (12)"),
    p_num = c(Box.test(rs,   lag = 20, type = "Ljung-Box")$p.value,
              Box.test(rs^2, lag = 20, type = "Ljung-Box")$p.value,
              FinTS::ArchTest(rs, lags = 12)$p.value)
  )
}

diag_w_all <- eg_fits %>%
  mutate(diag = map(fit, diag_one)) %>%
  select(Model, Innovations, diag) %>%
  unnest(diag) %>%
  pivot_wider(names_from = Test, values_from = p_num) %>%
  mutate(
    `LB(ε̂) p (num)`       = `LB(ε̂)`,
    `LB(ε̂²) p (num)`      = `LB(ε̂²)`,
    `ARCH-LM (12) p (num)` = `ARCH-LM (12)`,
    `LB(ε̂) p-value`       = fmt_p(`LB(ε̂)`),
    `LB(ε̂²) p-value`      = fmt_p(`LB(ε̂²)`),
    `ARCH-LM (12) p-value` = fmt_p(`ARCH-LM (12)`)
  )

# Keep ONLY Student-t rows for the diagnostics display
diag_w <- diag_w_all %>%
  filter(Innovations == "Student-t") %>%
  arrange(Model)

# Build display with a spacer column between LB(ε̂) and LB(ε̂²)
diag_disp <- diag_w %>%
  transmute(
    Model, Innovations,
    `LB(ε̂) p-value`,
    ` ` = "",                        # spacer
    `LB(ε̂²) p-value`,
    `ARCH-LM (12) p-value`
  )

gt_diag <- diag_disp %>%
  gt() %>%
  tab_header(
    title = md("**Residual Diagnostics — EGARCH(p,q) (Student-t only)**"),
    subtitle = "LB(ε̂): mean autocorr (want large p). LB(ε̂²)/ARCH-LM: leftover ARCH (want large p)."
  ) %>%
  opt_row_striping() %>%
  cols_width(
    Model ~ px(200),
    Innovations ~ px(160),
    `LB(ε̂) p-value` ~ px(170),
    ` ` ~ px(28),                    # visual gap
    `LB(ε̂²) p-value` ~ px(230),
    `ARCH-LM (12) p-value` ~ px(180)
  ) %>%
  cols_align(columns = ` `, align = "center") %>%
  # faint vertical rule after LB(ε̂)
  tab_style(
    style = cell_borders(sides = "right", color = "#D9D9D9", weight = px(2)),
    locations = cells_body(columns = `LB(ε̂) p-value`)
  ) %>%
  # bold where underlying numeric p >= 0.05
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `LB(ε̂) p-value`,
      rows = diag_w$`LB(ε̂) p (num)` >= 0.05
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `LB(ε̂²) p-value`,
      rows = diag_w$`LB(ε̂²) p (num)` >= 0.05
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `ARCH-LM (12) p-value`,
      rows = diag_w$`ARCH-LM (12) p (num)` >= 0.05
    )
  )
save_gt_png(gt_diag, "egarch_grid_diagnostics")

# ---------- Likelihood-Ratio Test: EGARCH(1,1) vs EGARCH(2,1) — Student-t -----
fit_e11_std <- eg_fits %>% filter(Model=="EGARCH(1,1)", Innovations=="Student-t") %>% pull(fit) %>% .[[1]]
fit_e21_std <- eg_fits %>% filter(Model=="EGARCH(2,1)", Innovations=="Student-t") %>% pull(fit) %>% .[[1]]

lr_row <- tryCatch({
  ll1 <- as.numeric(fit_e11_std@fit$LLH)
  ll2 <- as.numeric(fit_e21_std@fit$LLH)
  k1  <- length(coef(fit_e11_std))
  k2  <- length(coef(fit_e21_std))
  LR  <- 2 * (ll2 - ll1)
  df  <- k2 - k1
  pLR <- 1 - pchisq(LR, df)
  tibble(Comparison = "EGARCH(2,1) vs EGARCH(1,1) — Student-t",
         LR = round(LR, 3), df = df, `p-value` = fmt_p(pLR))
}, error = function(e) tibble(Comparison = "EGARCH(2,1) vs EGARCH(1,1) — Student-t",
                              LR = NA_real_, df = NA_integer_, `p-value` = NA_character_))

gt_lr <- lr_row %>%
  gt() %>%
  tab_header(
    title = md("**Likelihood-Ratio Test (Nested Variance Orders)**"),
    subtitle = "H₀: extra EGARCH term(s) do not improve fit; reject H₀ if p < 0.05."
  ) %>%
  cols_width(Comparison ~ px(520), LR ~ px(140), df ~ px(100), `p-value` ~ px(140)) %>%
  opt_row_striping()
save_gt_png(gt_lr, "egarch_lr_test")

# ---------- σ_t plots: EGARCH(1,1) — Student-t + best AIC model ---------------
# EGARCH(1,1) — Student-t
sig11 <- sigma(fit_e11_std)
df_sig11 <- tibble(date = dates, sigma = as.numeric(sig11))
p_sig11 <- ggplot(df_sig11, aes(date, sigma)) +
  geom_line(linewidth = 0.4) +
  labs(title = expression(paste("Estimated Conditional Volatility (", sigma[t], ")")),
       subtitle = "EGARCH(1,1) — Student-t",
       x = "Date", y = expression(sigma[t])) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(out_dir, "egarch11_sigma.png"), p_sig11, width = 11, height = 5.8, dpi = 300)

# Best AIC model (for reference)
best_idx <- which.min(eg_fits$AIC)
best_fit <- eg_fits$fit[[best_idx]]
best_lab <- paste0(eg_fits$Model[[best_idx]], " — ", eg_fits$Innovations[[best_idx]])
sig_best <- sigma(best_fit)
df_sig_best <- tibble(date = dates, sigma = as.numeric(sig_best))
p_sig_best <- ggplot(df_sig_best, aes(date, sigma)) +
  geom_line(linewidth = 0.4) +
  labs(title = expression(paste("Estimated Conditional Volatility (", sigma[t], ")")),
       subtitle = best_lab, x = "Date", y = expression(sigma[t])) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(out_dir, "egarch_best_sigma.png"), p_sig_best, width = 11, height = 5.8, dpi = 300)

message("Done. Saved IC, Student-t diagnostics, LR test, and sigma plots to: ", out_dir)
