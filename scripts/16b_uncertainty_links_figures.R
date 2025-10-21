# ================================================================
# 15_uncertainty_links_CORE_AND_EVENTS.R
# Core figures linking gold volatility (EGARCH σ_t), gold returns, VIX & EPU
# Outputs -> .../figures/07_uncertainty_figures/figures_b
#
# Requires:
#   - data/clean/04_egarch_data.RData  (sigma_series, rv_series, break_dates, fit_*)
#   - data/clean/02_clean_uncertainty_data.Rdata  (if vix/epu not already in env)
#     vix_full_clean: tibble(date, vix_close)  [daily]
#     epu_full_clean: tibble(date_m, epu)      [monthly]
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
outdir  <- file.path(outmain, "summary_figures")
fs::dir_create(outdir, recurse = TRUE)

# -------------------------------
# Load data
# -------------------------------
stopifnot(file.exists(in_vol))
load(in_vol)  # sigma_series, rv_series, break_dates, fit_base, fit_dummy
stopifnot(all(c("date","r","sigma_dummy") %in% names(sigma_series)))

if (!exists("vix_full_clean") || !exists("epu_full_clean")) {
  if (file.exists(in_unc)) load(in_unc) else
    stop("Missing VIX/EPU data. Ensure `vix_full_clean` and `epu_full_clean` are available or save 02_clean_uncertainty_data.Rdata.")
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

# Regime bands from PELT breaks (if any)
if (length(break_dates)) {
  bd <- sort(unique(as.Date(break_dates)))
  bands <- tibble(start = bd, end = c(bd[-1] - 1L, max(daily$date))) %>% filter(start <= end)
} else {
  bands <- tibble(start = as.Date(character()), end = as.Date(character()))
}

# Palettes
pal_overlay_two <- c("σ_t (EGARCH)" = "#1f77b4", "VIX" = "#d62728")
pal_overlay_three <- c("VIX (z)" = "#d62728", "σ_t (EGARCH, z)" = "#1f77b4", "Gold returns (z)" = "#2ca02c")
band_fill <- "grey92"

# ================================================================
# CORE FIGURE 1 — Daily overlay: z(σ_t) vs z(VIX) with regime bands
# ================================================================
overlay_two_df <- daily %>%
  mutate(z_sig = as.numeric(scale(sigma_dummy)),
         z_vix = as.numeric(scale(vix_close))) %>%
  pivot_longer(c(z_sig, z_vix), names_to = "series", values_to = "z") %>%
  mutate(series = recode(series, z_sig = "σ_t (EGARCH)", z_vix = "VIX"))

plt_core_overlay_sigma_vix <- ggplot(overlay_two_df, aes(date, z, color = series)) +
  { if (nrow(bands) > 0)
    geom_rect(data = bands,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, alpha = 0.55) } +
  geom_line(linewidth = 0.45) +
  scale_color_manual(values = pal_overlay_two, guide = guide_legend(title = NULL)) +
  labs(title = "Gold volatility and VIX (z-scores, daily)",
       subtitle = if (nrow(bands) > 0) "Shaded bands: variance regimes (PELT)" else NULL,
       x = "Date", y = "z-score") +
  theme_clean + theme(legend.position = "top")
save_png(plt_core_overlay_sigma_vix, "core_01_overlay_sigma_vix_daily.png", width = 14)

# ================================================================
# CORE FIGURE 2 — σ_t vs VIX (daily): 2D counts + LOESS + binned means
# ================================================================
scatter_df <- daily %>% filter(is.finite(vix_close), is.finite(sigma_dummy))
n_bins <- 30
rng     <- range(scatter_df$vix_close, na.rm = TRUE)
vix_seq <- seq(rng[1], rng[2], length.out = n_bins + 1)
binned_means <- scatter_df %>%
  mutate(vix_bin = cut(vix_close, breaks = vix_seq, include.lowest = TRUE)) %>%
  group_by(vix_bin) %>%
  summarise(vix_mid = mean(c(min(vix_close), max(vix_close)), na.rm = TRUE),
            sigma_mean = mean(sigma_dummy, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(is.finite(vix_mid), is.finite(sigma_mean))

plt_core_scatter_sigma_vix <- ggplot(scatter_df, aes(vix_close, sigma_dummy)) +
  geom_bin2d(bins = 40) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  geom_point(data = binned_means, aes(vix_mid, sigma_mean), inherit.aes = FALSE,
             size = 1.5, alpha = 0.9) +
  scale_fill_viridis_c(name = "count") +
  geom_vline(xintercept = c(20, 30), linetype = "dashed", color = "grey40") +
  labs(title = "Gold volatility vs VIX (daily)",
       subtitle = "2D counts with LOESS trend and binned means; dashed lines at VIX≈20/30",
       x = "VIX (close, daily)", y = expression(hat(sigma)[t])) +
  theme_clean
save_png(plt_core_scatter_sigma_vix, "core_02_scatter_sigma_vs_vix.png")

# ================================================================
# CORE FIGURE 3 — VIX×EPU deciles → median monthly σ_t  (masked counts)
# ================================================================
n_min <- 8
heatmap_df <- monthly %>%
  filter(is.finite(m_vix), is.finite(epu), is.finite(m_sigma)) %>%
  mutate(vix_d = ntile(m_vix, 10), epu_d = ntile(epu, 10)) %>%
  group_by(vix_d, epu_d) %>%
  summarise(med_sigma = median(m_sigma, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(show = n >= n_min)

plt_core_heat_vix_epu <- ggplot(heatmap_df, aes(vix_d, epu_d, fill = ifelse(show, med_sigma, NA))) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(show, n, "")), size = 3, color = "white") +
  scale_fill_viridis_c(name = "Median σ_t", na.value = "grey85") +
  scale_x_continuous(breaks = seq(1,10,2)) +
  scale_y_continuous(breaks = seq(1,10,2)) +
  labs(title = "Median monthly gold volatility by VIX and EPU deciles",
       subtitle = paste0("Counts annotated; cells masked when n<", n_min),
       x = "VIX decile (monthly mean)", y = "EPU decile") +
  theme_clean
save_png(plt_core_heat_vix_epu, "core_03_heat_sigma_vix_epu_deciles.png", width = 7, height = 6)

# ================================================================
# NEW SET A — 3-series overlays: VIX, σ_t, and gold log returns (z)
# Four windows: 1990–2025, 1997–2005, 2005–2015, 2015–2025
# ================================================================
overlay_three_df <- daily %>%
  transmute(date,
            `σ_t (EGARCH, z)` = as.numeric(scale(sigma_dummy)),
            `VIX (z)`         = as.numeric(scale(vix_close)),
            `Gold returns (z)`= as.numeric(scale(r))) %>%
  pivot_longer(-date, names_to = "series", values_to = "z")

plot_overlay_window <- function(start_date, end_date, fname, title_suffix) {
  # clip regime bands to window
  bands_win <- if (nrow(bands)) {
    bands %>%
      mutate(start2 = pmax(start, start_date),
             end2   = pmin(end,   end_date)) %>%
      filter(end2 >= start_date & start2 <= end_date) %>%
      transmute(start = start2, end = end2)
  } else tibble(start = as.Date(character()), end = as.Date(character()))
  
  g <- overlay_three_df %>%
    filter(date >= start_date, date <= end_date) %>%
    ggplot(aes(date, z, color = series)) +
    { if (nrow(bands_win) > 0)
      geom_rect(data = bands_win,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = band_fill, alpha = 0.55) } +
    geom_line(linewidth = 0.45) +
    scale_color_manual(values = pal_overlay_three, guide = guide_legend(title = NULL)) +
    labs(title = paste0("VIX, gold volatility, and gold returns (z-scores) — ", title_suffix),
         subtitle = if (nrow(bands_win) > 0) "Shaded bands: variance regimes (PELT)" else NULL,
         x = "Date", y = "z-score") +
    theme_clean + theme(legend.position = "top")
  save_png(g, fname, width = 14)
}

plot_overlay_window(as.Date("1990-01-01"), as.Date("2025-12-31"),
                    "panelA_overlay3_1990_2025.png",  "1990–2025")
plot_overlay_window(as.Date("1997-01-01"), as.Date("2005-12-31"),
                    "panelA_overlay3_1997_2005.png",  "1997–2005")
plot_overlay_window(as.Date("2005-01-01"), as.Date("2015-12-31"),
                    "panelA_overlay3_2005_2015.png",  "2005–2015")
plot_overlay_window(as.Date("2015-01-01"), as.Date("2025-12-31"),
                    "panelA_overlay3_2015_2025.png",  "2015–2025")

# ================================================================
# NEW SET B — Event window around VIX spikes (±30d)
# Two panels: mean returns & mean σ_t across all top-1% spikes,
# with 2008 and 2020 spike episodes overlaid.
# ================================================================
window <- 30

# all top 1% VIX days
spike_dates <- daily %>%
  filter(is.finite(vix_close)) %>%
  filter(vix_close >= quantile(vix_close, 0.99, na.rm = TRUE)) %>%
  pull(date)

# representative crisis peaks
crisis_2008 <- daily %>%
  filter(date >= as.Date("2008-08-01"), date <= as.Date("2009-06-30")) %>%
  slice_max(vix_close, n = 1, with_ties = FALSE) %>% pull(date)
crisis_2020 <- daily %>%
  filter(date >= as.Date("2020-02-01"), date <= as.Date("2020-05-31")) %>%
  slice_max(vix_close, n = 1, with_ties = FALSE) %>% pull(date)

build_event_df <- function(center_date, label) {
  i <- which(daily$date == center_date)
  if (!length(i)) return(tibble())
  lo <- max(1, i - window); hi <- min(nrow(daily), i + window)
  tibble(rel_day = (lo:hi) - i,
         r = daily$r[lo:hi],
         sigma = daily$sigma_dummy[lo:hi],
         which = label)
}

avg_df <- map_dfr(spike_dates, ~ build_event_df(.x, "Average top 1%")) %>%
  group_by(rel_day, which) %>%
  summarise(mean_r = mean(r, na.rm = TRUE),
            mean_sigma = mean(sigma, na.rm = TRUE), .groups = "drop")

crisis_df <- bind_rows(
  build_event_df(crisis_2008, "2008 spike"),
  build_event_df(crisis_2020, "2020 spike")
)

palette_event <- c("Average top 1%" = "#1f77b4", "2008 spike" = "#ff7f0e", "2020 spike" = "#2ca02c")

# returns panel
plt_event_returns <- ggplot() +
  geom_line(data = avg_df, aes(rel_day, mean_r, color = which), linewidth = 0.9) +
  geom_line(data = crisis_df, aes(rel_day, r, color = which), alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = palette_event) +
  labs(title = "Gold returns around VIX spikes (±30 trading days)",
       subtitle = "Average across top-1% VIX days vs 2008 and 2020 episodes",
       x = NULL, y = "Return") +
  theme_clean + theme(legend.position = "top")

# volatility panel
plt_event_sigma <- ggplot() +
  geom_line(data = avg_df, aes(rel_day, mean_sigma, color = which), linewidth = 0.9) +
  geom_line(data = crisis_df, aes(rel_day, sigma, color = which), alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = palette_event, guide = "none") +
  labs(title = "Gold volatility (σ_t) around VIX spikes (±30 trading days)",
       subtitle = "Average across top-1% VIX days vs 2008 and 2020 episodes",
       x = "Days from spike", y = expression(hat(sigma)[t])) +
  theme_clean

# combine & save (uses patchwork if available; otherwise save separate panels)
if (requireNamespace("patchwork", quietly = TRUE)) {
  ev_panel <- plt_event_returns / plt_event_sigma + patchwork::plot_layout(heights = c(1,1.05))
  save_png(ev_panel, "panelB_event_window_returns_and_sigma.png", width = 14, height = 8)
} else {
  message("Package 'patchwork' not found — saving event panels separately.")
  save_png(plt_event_returns, "panelB_event_window_returns.png", width = 14, height = 4.5)
  save_png(plt_event_sigma,   "panelB_event_window_sigma.png",   width = 14, height = 4.5)
}

# -------------------------------
message("Saved core & event figures to: ", outdir)
# -------------------------------
