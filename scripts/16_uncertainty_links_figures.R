# ================================================================
# 15_uncertainty_links_figures_REDUCED.R
# High-signal visuals linking gold volatility (EGARCH σ_t) to VIX & EPU
# Saves PNGs to: /Users/cycoldiron/Desktop/gold-ts-forecast/figures/07_uncertainty_figures
#
# Inputs:
#   - data/clean/04_egarch_data.RData  (sigma_series, rv_series, break_dates, fit_*)
#   - data/clean/02_clean_uncertainty_data.Rdata (if VIX/EPU not in memory)
#     Objects expected:
#       * vix_full_clean: tibble(date, vix_close)  [daily]
#       * epu_full_clean: tibble(date_m, epu)      [monthly]
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(lubridate); library(ggplot2); library(scales)
  library(slider); library(fs)
})

# -------------------------------
# Paths / output
# -------------------------------
root   <- "/Users/cycoldiron/Desktop/gold-ts-forecast"
in_vol <- file.path(root, "data/clean/04_egarch_data.RData")
in_unc <- file.path(root, "data/clean/02_clean_uncertainty_data.Rdata")
outdir <- file.path(root, "figures/07_uncertainty_figures")
fs::dir_create(outdir, recurse = TRUE)

# -------------------------------
# Load volatility data
# -------------------------------
stopifnot(file.exists(in_vol))
load(in_vol)  # sigma_series, rv_series, break_dates, fit_base, fit_dummy
stopifnot(all(c("date","r","sigma_dummy") %in% names(sigma_series)))

# -------------------------------
# Ensure uncertainty data available
# -------------------------------
if (!exists("vix_full_clean") || !exists("epu_full_clean")) {
  if (file.exists(in_unc)) load(in_unc) else
    stop("Missing VIX/EPU data. Load `vix_full_clean` and `epu_full_clean`, or place 02_clean_uncertainty_data.Rdata.")
}
stopifnot(all(c("date","vix_close") %in% names(vix_full_clean)))
stopifnot(all(c("date_m","epu")     %in% names(epu_full_clean)))

# -------------------------------
# Build daily & monthly merges
# -------------------------------
vix_day <- vix_full_clean %>%
  transmute(date = as.Date(date), vix_close = as.numeric(vix_close))

daily <- sigma_series %>%
  select(date, r, sigma_dummy) %>%
  left_join(vix_day, by = "date") %>%
  arrange(date)

monthly <- daily %>%
  mutate(date_m = floor_date(date, "month")) %>%
  group_by(date_m) %>%
  summarise(
    m_sigma = mean(sigma_dummy, na.rm = TRUE),
    m_vix   = mean(vix_close,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(epu_full_clean %>% transmute(date_m = as.Date(date_m), epu = as.numeric(epu)),
            by = "date_m")

# -------------------------------
# Helpers
# -------------------------------
save_png <- function(plot, filename, width = 12, height = 6, dpi = 300) {
  ggsave(file.path(outdir, filename), plot, width = width, height = height, dpi = dpi)
}

theme_clean <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(face = "bold"))

# ================================================================
# FIGURE 1 (KEPT, UPGRADED)
# Scatter-density: σ_t vs VIX (+ LOESS, binned means, reference lines)
# ================================================================
p1_df <- daily %>% filter(is.finite(vix_close), is.finite(sigma_dummy))

# binned means (equal-width bins on VIX)
n_bins <- 30
vix_brks <- pretty(range(p1_df$vix_close, na.rm = TRUE), n = n_bins)
bin_means <- p1_df %>%
  mutate(vix_bin = cut(vix_close, breaks = vix_brks, include.lowest = TRUE)) %>%
  group_by(vix_bin) %>%
  summarise(
    vix_mid = mean(range(vix_close)[1] + (as.numeric(vix_bin) - 0.5) *
                     diff(range(vix_close)) / (length(unique(vix_bin)))),
    sigma_mean = mean(sigma_dummy, na.rm = TRUE),
    n = dplyr::n(), .groups = "drop"
  ) %>% filter(is.finite(vix_mid), is.finite(sigma_mean))

p1 <- ggplot(p1_df, aes(vix_close, sigma_dummy)) +
  geom_bin2d(bins = 40) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  geom_point(data = bin_means, aes(vix_mid, sigma_mean),
             inherit.aes = FALSE, size = 1.5, alpha = 0.9) +
  scale_fill_viridis_c(name = "count") +
  geom_vline(xintercept = c(20, 30), linetype = "dashed", color = "grey40") +
  labs(title = "Gold conditional volatility vs VIX",
       subtitle = "2D counts, LOESS trend, and binned means; dashed lines at VIX≈20/30",
       x = "VIX (close, daily)", y = expression(hat(sigma)[t])) +
  theme_clean
save_png(p1, "01_scatter_sigma_vs_vix_hi.png")

# ================================================================
# FIGURE 2 (KEPT, UPGRADED)
# Lead–lag correlations: EPU leads/lags vs monthly σ_t with bootstrap CIs
#   Convention: +k => EPU leads σ_t by k months
# ================================================================
shift_vec <- function(x, k) {
  if (k > 0) dplyr::lag(x,  n = k)      # EPU leads
  else if (k < 0) dplyr::lead(x, n = -k) # EPU lags
  else x
}

# z-scores for comparability
p2_base <- monthly %>%
  transmute(date_m,
            z_sigma = as.numeric(scale(m_sigma)),
            z_epu   = as.numeric(scale(epu)))

lags <- -6:6

# bootstrap CI for correlation at each lag
boot_corr <- function(x, y, B = 500L, seed = 123) {
  set.seed(seed)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 10) return(c(NA, NA, NA))
  r_hat <- suppressWarnings(cor(x, y))
  n <- length(x)
  rb <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx]))
  })
  q <- quantile(rb, c(0.025, 0.975), na.rm = TRUE)
  c(r = r_hat, lo = q[[1]], hi = q[[2]])
}

p2_df <- map_dfr(lags, function(L) {
  y <- p2_base$z_sigma
  x <- shift_vec(p2_base$z_epu, L)
  out <- boot_corr(y, x, B = 600)
  tibble(lag = L, corr = out[1], lo = out[2], hi = out[3])
})

# identify max |corr|
imax <- which.max(abs(p2_df$corr %>% replace_na(0)))
lab_point <- p2_df[imax, ]

p2 <- ggplot(p2_df, aes(lag, corr)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
  geom_col(width = 0.8, fill = "#1f77b4") +
  geom_point(data = lab_point, aes(lag, corr), color = "#d62728", size = 2) +
  geom_text(data = lab_point,
            aes(lag, corr, label = sprintf("max |r|=%.2f @ %s", corr,
                                           ifelse(lag>=0, paste0("+",lag), lag))),
            vjust = -0.8, color = "#d62728", fontface = "bold", size = 3.3) +
  labs(title = "Lead–lag correlations: EPU vs monthly σ_t",
       subtitle = "+k: EPU leads by k months; bars with bootstrap 95% CIs (ribbons)",
       x = "Lag (months)", y = "Correlation") +
  theme_clean
save_png(p2, "02_leadlag_epu_sigma_boot.png")

# ================================================================
# FIGURE 3 (KEPT, UPGRADED)
# Joint-stress heatmap: VIX decile × EPU decile → median monthly σ_t
#   - masks sparse cells (n < n_min)
#   - overlays counts
# ================================================================
n_min <- 8

heat_df <- monthly %>%
  filter(is.finite(m_vix), is.finite(epu), is.finite(m_sigma)) %>%
  mutate(vix_d = ntile(m_vix, 10),
         epu_d = ntile(epu, 10)) %>%
  group_by(vix_d, epu_d) %>%
  summarise(med_sigma = median(m_sigma, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop") %>%
  mutate(show = n >= n_min)

p3 <- ggplot(heat_df, aes(vix_d, epu_d, fill = ifelse(show, med_sigma, NA))) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(show, n, "")), size = 3, color = "white") +
  scale_fill_viridis_c(name = "Median σ_t", na.value = "grey85") +
  scale_x_continuous(breaks = seq(1,10,2)) +
  scale_y_continuous(breaks = seq(1,10,2)) +
  labs(title = "Median monthly σ_t by VIX and EPU deciles",
       subtitle = paste0("Cells annotated with counts (masked where n<", n_min, ")"),
       x = "VIX decile (monthly mean)", y = "EPU decile") +
  theme_clean
save_png(p3, "03_heat_sigma_vix_epu_deciles_masked.png", width = 7, height = 6)


# ---- Put near the top (after theme_clean) ----
pal <- c("σ_t (EGARCH)" = "#1f77b4",   # blue
         "VIX"           = "#d62728")  # red
band_fill <- "grey92"

# ---- Colors & helpers (put once, above the plot) ----
pal <- c("σ_t (EGARCH)" = "#1f77b4",  # blue
         "VIX"           = "#d62728") # red
band_fill <- "grey92"

# Build long-form data for colored overlay
p4_df <- daily %>%
  mutate(z_sig = as.numeric(scale(sigma_dummy)),
         z_vix = as.numeric(scale(vix_close))) %>%
  tidyr::pivot_longer(c(z_sig, z_vix),
                      names_to = "series", values_to = "z") %>%
  mutate(series = dplyr::recode(series,
                                z_sig = "σ_t (EGARCH)",
                                z_vix = "VIX"))

# Plot
p4 <- ggplot(p4_df, aes(date, z, color = series)) +
  { if (nrow(bands) > 0)
    geom_rect(data = bands,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, alpha = 0.55) } +
  geom_line(linewidth = 0.45) +
  scale_color_manual(values = pal, guide = guide_legend(title = NULL)) +
  labs(title = "Gold volatility vs VIX (z-scores, daily)",
       subtitle = if (nrow(bands) > 0) "Shaded bands: variance regimes (PELT)" else NULL,
       x = "Date", y = "z-score") +
  theme_clean +
  theme(legend.position = "top")

save_png(p4, "04_overlay_sigma_vix_daily_colored.png", width = 14)




# -------------------------------
# Done
# -------------------------------
message("Saved high-signal PNGs to: ", outdir)
