# 03_explore_gold_acf_pacf.R  (v5)
# Full-sample ACF/PACF with improved aesthetics + unified y-scales
# Saves PNGs to .../figures/02_diagnostics/01_acf_pacf

library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
library(glue)

# ---------- Data ----------
# gold <- readr::read_csv("data/clean/gold_daily.csv", show_col_types = FALSE)
gold <- gold %>%
  arrange(date) %>%
  mutate(
    close_log        = as.numeric(close_log),
    close_log_return = as.numeric(close_log_return)
  )

end_date <- max(gold$date, na.rm = TRUE)

# ---------- Helpers ----------
tidy_acf <- function(x, max_lag = 60, type = c("acf","pacf")) {
  type <- match.arg(type); x <- as.numeric(x); n_eff <- sum(!is.na(x))
  if (type == "acf") {
    obj  <- stats::acf(x, lag.max = max_lag, plot = FALSE, na.action = na.pass)
    vals <- as.numeric(obj$acf)[-1]; lags <- seq_along(vals)
  } else {
    obj  <- stats::pacf(x, lag.max = max_lag, plot = FALSE, na.action = na.pass)
    vals <- as.numeric(obj$acf);     lags <- seq_along(vals)
  }
  tibble(lag = as.integer(lags), value = vals, ci = 1.96/sqrt(n_eff), type = type)
}

lb_table <- function(x, max_lag = 20, fitdf = 0) {
  x <- as.numeric(na.omit(x))
  tibble(lag = 1:max_lag) |>
    mutate(
      Q       = map_dbl(lag, ~ Box.test(x, lag = .x, type = "Ljung-Box", fitdf = fitdf)$statistic %>% as.numeric()),
      p_value = map_dbl(lag, ~ Box.test(x, lag = .x, type = "Ljung-Box", fitdf = fitdf)$p.value   %>% as.numeric())
    )
}

# Aesthetic tweaks: colored CI lines, clearer zero line, unified y-limits, CENTERED TITLES
plot_corr <- function(df, title, ylab, y_lim = c(-1, 1),
                      ci_color = "#2B6CB0", zero_color = "#9AA0A6", bar_fill = "#3C4043") {
  ggplot(df, aes(lag, value)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = zero_color) +
    geom_hline(aes(yintercept =  ci), linetype = "22", linewidth = 0.4, color = ci_color) +
    geom_hline(aes(yintercept = -ci), linetype = "22", linewidth = 0.4, color = ci_color) +
    geom_col(width = 0.88, fill = bar_fill) +
    coord_cartesian(ylim = y_lim, expand = FALSE) +
    labs(title = title, x = "Lag (trading days)", y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),  # centered
      panel.grid.minor  = element_blank(),
      panel.grid.major.x= element_blank(),
      axis.title.x      = element_text(margin = margin(t = 6)),
      axis.title.y      = element_text(margin = margin(r = 6))
    )
}

make_diagnostics <- function(df, label, max_lag = 60, save = TRUE,
                             out_dir = "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/01_acf_pacf") {
  
  x_logp <- df$close_log
  x_ret  <- df$close_log_return %>% stats::na.omit()
  
  # Unified y-scale across all panels
  y_lim <- c(-1, 1)
  
  # ----- Log price -----
  acf_logp  <- tidy_acf(x_logp, max_lag, "acf")
  pacf_logp <- tidy_acf(x_logp, max_lag, "pacf")
  g1 <- plot_corr(acf_logp,  "ACF — log price",  "ACF",  y_lim = y_lim)
  g2 <- plot_corr(pacf_logp, "PACF — log price", "PACF", y_lim = y_lim)
  fig_lp <- (g1 | g2)
  
  # ----- Returns -----
  acf_ret  <- tidy_acf(x_ret, max_lag, "acf")
  pacf_ret <- tidy_acf(x_ret, max_lag, "pacf")
  g3 <- plot_corr(acf_ret,  "ACF — daily log return",  "ACF",  y_lim = y_lim)
  g4 <- plot_corr(pacf_ret, "PACF — daily log return", "PACF", y_lim = y_lim)
  fig_ret <- (g3 | g4)
  
  # Optional: Q-stats (full-sample only)
  dir.create("results/diagnostics", showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(lb_table(x_logp, 20, 0), glue("results/diagnostics/lbq_logprice_{label}.csv"))
  readr::write_csv(lb_table(x_ret,  20, 0), glue("results/diagnostics/lbq_logreturns_{label}.csv"))
  
  # ----- Save -----
  if (save) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, glue("acf_pacf_logprice_{label}.png")),  fig_lp,  width = 11, height = 4.2, dpi = 300)
    ggsave(file.path(out_dir, glue("acf_pacf_logreturns_{label}.png")), fig_ret, width = 11, height = 4.2, dpi = 300)
  }
  
  list(fig_logprice = fig_lp, fig_returns = fig_ret)
}

# ---------- Run ONLY the full window & Save ----------
full_range <- c(as.Date("1990-01-01"), end_date)
sub <- gold %>% filter(date >= full_range[1], date <= full_range[2])
make_diagnostics(sub, label = "full", max_lag = 60, save = TRUE)

message("\nSaved ACF/PACF PNGs (full sample only) to:\n",
        "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/01_acf_pacf\n",
        "Lag units = trading days; y-axis fixed to [-1, 1] across all panels.\n")
