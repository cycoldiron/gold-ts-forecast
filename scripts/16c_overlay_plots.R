# ================================================================
# 15_uncertainty_links_OVERLAYS_SMOOTHED_WITH_BANDS.R
# Smoothed overlays with correct variance-regime bands (PELT)
# Outputs -> .../figures/07_uncertainty_figures/overlay_plots
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(lubridate); library(ggplot2); library(scales)
  library(slider); library(fs)
})

# -------------------------------
# Paths / output
# -------------------------------
root    <- "/Users/cycoldiron/Desktop/gold-ts-forecast"
in_vol  <- file.path(root, "data/clean/04_egarch_data.RData")
in_unc  <- file.path(root, "data/clean/02_clean_uncertainty_data.Rdata")
outmain <- file.path(root, "figures/07_uncertainty_figures")
outdir  <- file.path(outmain, "overlay_plots")
fs::dir_create(outdir, recurse = TRUE)

# -------------------------------
# Load data
# -------------------------------
stopifnot(file.exists(in_vol))
load(in_vol)  # sigma_series, rv_series, break_dates, fit_base, fit_dummy
stopifnot(all(c("date","r","sigma_dummy") %in% names(sigma_series)))

if (!exists("vix_full_clean") || !exists("epu_full_clean")) {
  if (file.exists(in_unc)) load(in_unc) else
    stop("Missing VIX/EPU data. Ensure `vix_full_clean` and `epu_full_clean` exist.")
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
  left_join(epu_full_clean %>%
              transmute(date_m = as.Date(date_m), epu = as.numeric(epu)),
            by = "date_m")

# -------------------------------
# Helpers & theme (centered titles)
# -------------------------------
save_png <- function(plot, filename, width = 12, height = 6, dpi = 300) {
  ggsave(file.path(outdir, filename), plot, width = width, height = height, dpi = dpi)
}

theme_clean <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(face = "bold"))

pal_two    <- c("σ_t (EGARCH)" = "#1f77b4", "VIX" = "#d62728")
pal_epu    <- c("σ_t (EGARCH)" = "#1f77b4", "EPU" = "#2ca02c")
pal_three  <- c("σ_t (EGARCH)" = "#1f77b4", "VIX" = "#d62728", "EPU" = "#2ca02c")
band_fill  <- "grey92"

roll_mean <- function(x, k) slider::slide_dbl(x, mean, .before = k-1, .complete = TRUE, na.rm = TRUE)
zscore    <- function(x) as.numeric(scale(x))

# -------------------------------
# Variance-regime bands from break_dates (PELT on gold)
# Supports paired ranges OR single changepoints
# -------------------------------
make_bands <- function(break_dates, x_min, x_max) {
  if (is.null(break_dates) || length(break_dates) == 0) {
    return(tibble(start = as.Date(character()), end = as.Date(character())))
  }
  bd <- sort(as.Date(break_dates))
  use_pairs <- (length(bd) %% 2 == 0)
  if (use_pairs) {
    starts <- bd[seq(1, length(bd), by = 2)]
    ends   <- bd[seq(2, length(bd), by = 2)]
    bands  <- tibble(start = starts, end = ends)
  } else {
    starts <- bd
    ends   <- c(bd[-1] - 1L, x_max)
    bands  <- tibble(start = starts, end = ends)
  }
  bands %>%
    mutate(start = pmax(start, x_min),
           end   = pmin(end,   x_max)) %>%
    filter(is.finite(start), is.finite(end), end >= start)
}

# Daily & monthly bands (monthly uses the same Date limits)
bands_daily <- make_bands(break_dates, min(daily$date, na.rm = TRUE), max(daily$date, na.rm = TRUE))
bands_month <- bands_daily

# ================================================================
# DAILY OVERLAYS (σ_t & VIX) — with regime bands
# ================================================================
daily_z <- daily %>%
  transmute(date,
            z_sig = zscore(sigma_dummy),
            z_vix = zscore(vix_close)) %>%
  mutate(z_sig_21 = roll_mean(z_sig, 21),
         z_vix_21 = roll_mean(z_vix, 21))

ov_daily_long <- daily_z %>%
  select(date,
         `σ_t (EGARCH)` = z_sig_21,
         VIX            = z_vix_21) %>%
  pivot_longer(-date, names_to = "series", values_to = "z_21")

ov_daily_raw <- daily_z %>%
  select(date,
         `σ_t (EGARCH)` = z_sig,
         VIX            = z_vix) %>%
  pivot_longer(-date, names_to = "series", values_to = "z_raw")

plt_overlay_daily_bands <- ggplot() +
  { if (nrow(bands_daily) > 0)
    geom_rect(data = bands_daily,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, alpha = 0.50) } +
  geom_line(data = ov_daily_raw, aes(date, z_raw, color = series),
            alpha = 0.25, linewidth = 0.25) +
  geom_line(data = ov_daily_long, aes(date, z_21, color = series),
            linewidth = 0.9) +
  scale_color_manual(values = pal_two, guide = guide_legend(title = NULL)) +
  labs(title = "Gold volatility and VIX (z-scores, daily; 21-day smoothing)",
       subtitle = "Dark = 21-day rolling mean; light = raw daily values. Shaded bands = variance regimes in gold (PELT).",
       x = "Date", y = "z-score") +
  theme_clean + theme(legend.position = "top")
save_png(plt_overlay_daily_bands, "overlay_daily_sigma_vix_21d_with_bands.png", width = 14)

plt_overlay_daily_smooth_only <- ggplot() +
  { if (nrow(bands_daily) > 0)
    geom_rect(data = bands_daily,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, alpha = 0.50) } +
  geom_line(data = ov_daily_long, aes(date, z_21, color = series), linewidth = 0.9) +
  scale_color_manual(values = pal_two, guide = guide_legend(title = NULL)) +
  labs(title = "Gold volatility and VIX (z-scores, daily; 21-day smoothing only)",
       subtitle = "Shaded bands = variance regimes in gold (PELT).",
       x = "Date", y = "z-score") +
  theme_clean + theme(legend.position = "top")
save_png(plt_overlay_daily_smooth_only, "overlay_daily_sigma_vix_21d_only_with_bands.png", width = 14)

# Windowed
plot_overlay_window_two <- function(start_date, end_date, fname, title_suffix) {
  bands_win <- bands_daily %>%
    mutate(start2 = pmax(start, start_date),
           end2   = pmin(end,   end_date)) %>%
    filter(end2 >= start_date & start2 <= end_date) %>%
    transmute(start = start2, end = end2)
  
  g <- ov_daily_long %>%
    filter(date >= start_date, date <= end_date) %>%
    ggplot(aes(date, z_21, color = series)) +
    { if (nrow(bands_win) > 0)
      geom_rect(data = bands_win,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = band_fill, alpha = 0.50) } +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal_two, guide = guide_legend(title = NULL)) +
    labs(title = paste0("Gold volatility and VIX (z-scores; 21-day smoothing) — ", title_suffix),
         subtitle = "Shaded bands = variance regimes in gold (PELT).",
         x = "Date", y = "z-score") +
    theme_clean + theme(legend.position = "top")
  save_png(g, fname, width = 14)
}

plot_overlay_window_two(as.Date("1990-01-01"), as.Date("2025-12-31"),
                        "overlay_daily_sigma_vix_1990_2025_21d.png", "1990–2025")
plot_overlay_window_two(as.Date("2005-01-01"), as.Date("2015-12-31"),
                        "overlay_daily_sigma_vix_2005_2015_21d.png", "2005–2015")
plot_overlay_window_two(as.Date("2015-01-01"), as.Date("2025-12-31"),
                        "overlay_daily_sigma_vix_2015_2025_21d.png", "2015–2025")

# ================================================================
# MONTHLY OVERLAYS (σ_t, VIX, EPU) — with regime bands
# ================================================================
monthly_z <- monthly %>%
  transmute(date_m,
            z_sigma = zscore(m_sigma),
            z_vix   = zscore(m_vix),
            z_epu   = zscore(epu)) %>%
  mutate(z_sigma_3m = roll_mean(z_sigma, 3),
         z_vix_3m   = roll_mean(z_vix, 3),
         z_epu_3m   = roll_mean(z_epu, 3))

# 3-series (σ_t, VIX, EPU)
ov_monthly_three <- monthly_z %>%
  select(date_m,
         `σ_t (EGARCH)` = z_sigma_3m,
         VIX            = z_vix_3m,
         EPU            = z_epu_3m) %>%
  pivot_longer(-date_m, names_to = "series", values_to = "z_3m")

plt_overlay_monthly_three <- ggplot() +
  { if (nrow(bands_month) > 0)
    geom_rect(data = bands_month,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, alpha = 0.50) } +
  geom_line(data = ov_monthly_three, aes(date_m, z_3m, color = series), linewidth = 0.9) +
  scale_color_manual(values = pal_three, guide = guide_legend(title = NULL)) +
  labs(title = "Gold volatility, VIX, and EPU (z-scores, monthly; 3-month smoothing)",
       subtitle = "Shaded bands = variance regimes in gold (PELT). Series are monthly means, then 3-month rolling averages.",
       x = "Month", y = "z-score") +
  theme_clean + theme(legend.position = "top")
save_png(plt_overlay_monthly_three, "overlay_monthly_sigma_vix_epu_3m_with_bands.png", width = 14)

# σ_t & EPU only
ov_monthly_epu <- monthly_z %>%
  select(date_m,
         `σ_t (EGARCH)` = z_sigma_3m,
         EPU            = z_epu_3m) %>%
  pivot_longer(-date_m, names_to = "series", values_to = "z_3m")

plt_overlay_monthly_epu_only <- ggplot() +
  { if (nrow(bands_month) > 0)
    geom_rect(data = bands_month,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, alpha = 0.50) } +
  geom_line(data = ov_monthly_epu, aes(date_m, z_3m, color = series), linewidth = 0.9) +
  scale_color_manual(values = pal_epu, guide = guide_legend(title = NULL)) +
  labs(title = "Gold volatility and EPU (z-scores, monthly; 3-month smoothing)",
       subtitle = "Shaded bands = variance regimes in gold (PELT). Series are monthly means, then 3-month rolling averages.",
       x = "Month", y = "z-score") +
  theme_clean + theme(legend.position = "top")
save_png(plt_overlay_monthly_epu_only, "overlay_monthly_sigma_epu_3m_with_bands.png", width = 14)

# -------------------------------
message("Saved overlays (with proper PELT bands) to: ", outdir)
# -------------------------------
