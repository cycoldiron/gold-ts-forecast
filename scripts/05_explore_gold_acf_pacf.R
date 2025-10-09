# 03_explore_gold_acf_pacf.R  (v2)
# Plots overview + ACF/PACF with Q-stats; saves to figures/01_overview and 02_diagnostics

library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
library(glue)

# ---------- Source caption (update to your exact provenance) ----------
SOURCE_CAPTION <- "Source: Stooq (XAUUSD daily close, USD/oz, 1990–present)."

# ---------- Data ----------
# gold <- readr::read_csv("data/clean/gold_daily.csv", show_col_types = FALSE)
gold <- gold %>% arrange(date) %>%
  mutate(close_log = as.numeric(close_log),
         close_log_return = as.numeric(close_log_return))

end_date <- max(gold$date, na.rm = TRUE)

# ---------- Helpers ----------
PANEL_TITLE <- function(suffix, .start, .end) {
  glue("Gold (daily, {year(.start)}–{year(.end)}): {suffix}")
}

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

plot_corr <- function(df, title, ylab) {
  ggplot(df, aes(lag, value)) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_hline(aes(yintercept =  ci), linetype = 2, linewidth = 0.2) +
    geom_hline(aes(yintercept = -ci), linetype = 2, linewidth = 0.2) +
    geom_col(width = 0.9) +
    labs(title = title, x = "Lag (trading days)", y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    )
}

make_ts_overview <- function(df, label, save = TRUE) {
  start_ <- min(df$date, na.rm = TRUE); end_ <- max(df$date, na.rm = TRUE)
  
  p_logprice <- ggplot(df, aes(date, close_log)) +
    geom_line(linewidth = 0.4) +
    labs(title = PANEL_TITLE("log price", start_, end_), x = NULL, y = "log price (USD/oz)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14), panel.grid.minor = element_blank())
  
  p_returns <- ggplot(df, aes(date, close_log_return)) +
    geom_line(linewidth = 0.35) +
    labs(title = PANEL_TITLE("daily log return", start_, end_), x = NULL, y = "Δ log price (daily)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14), panel.grid.minor = element_blank())
  
  fig <- (p_logprice / p_returns) +
    plot_annotation(caption = SOURCE_CAPTION) &
    theme(plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 6)))
  
  if (save) {
    dir.create("figures/01_overview", showWarnings = FALSE, recursive = TRUE)
    ggsave(glue("figures/01_overview/ts_overview_{label}.png"), fig, width = 11, height = 6.2, dpi = 300)
  }
  fig
}

make_diagnostics <- function(df, label, max_lag = 60, save = TRUE) {
  start_ <- min(df$date, na.rm = TRUE); end_ <- max(df$date, na.rm = TRUE)
  x_logp <- df$close_log
  x_ret  <- df$close_log_return %>% stats::na.omit()
  
  # Log price
  acf_logp  <- tidy_acf(x_logp, max_lag, "acf")
  pacf_logp <- tidy_acf(x_logp, max_lag, "pacf")
  g1 <- plot_corr(acf_logp,  PANEL_TITLE("ACF — log price", start_, end_),  "ACF")
  g2 <- plot_corr(pacf_logp, PANEL_TITLE("PACF — log price", start_, end_), "PACF")
  fig_lp <- (g1 | g2) + plot_annotation(caption = SOURCE_CAPTION) &
    theme(plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 6)))
  
  # Returns
  acf_ret  <- tidy_acf(x_ret, max_lag, "acf")
  pacf_ret <- tidy_acf(x_ret, max_lag, "pacf")
  g3 <- plot_corr(acf_ret,  PANEL_TITLE("ACF — daily log return", start_, end_),  "ACF")
  g4 <- plot_corr(pacf_ret, PANEL_TITLE("PACF — daily log return", start_, end_), "PACF")
  fig_ret <- (g3 | g4) + plot_annotation(caption = SOURCE_CAPTION) &
    theme(plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 6)))
  
  # Q-stats (saved to results/diagnostics)
  dir.create("results/diagnostics", showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(lb_table(x_logp, 20, 0), glue("results/diagnostics/lbq_logprice_{label}.csv"))
  readr::write_csv(lb_table(x_ret,  20, 0), glue("results/diagnostics/lbq_logreturns_{label}.csv"))
  
  if (save) {
    dir.create("figures/02_diagnostics", showWarnings = FALSE, recursive = TRUE)
    ggsave(glue("figures/02_diagnostics/acf_pacf_logprice_{label}.png"),  fig_lp,  width = 11, height = 4.2, dpi = 300)
    ggsave(glue("figures/02_diagnostics/acf_pacf_logreturns_{label}.png"), fig_ret, width = 11, height = 4.2, dpi = 300)
  }
  list(fig_logprice = fig_lp, fig_returns = fig_ret)
}

# ---------- Windows to analyze ----------
windows <- list(
  full     = c(as.Date("1990-01-01"), end_date),
  post_gfc = c(as.Date("2010-01-01"), end_date),
  recent   = c(as.Date("2018-01-01"), end_date),
  last1y   = c(end_date - 365,        end_date),
  last30d  = c(end_date - 30,         end_date)
)

# ---------- Run & Save ----------
walk(names(windows), function(nm) {
  rng <- windows[[nm]]
  sub <- gold %>% filter(date >= rng[1], date <= rng[2])
  make_ts_overview(sub, label = nm)
  make_diagnostics(sub, label = nm, max_lag = 60)
})

# Console hint
message("\nSaved overview PNGs to figures/01_overview/ and diagnostics PNGs to figures/02_diagnostics/.\n",
        "Lag units = trading days for these daily series.\n")
