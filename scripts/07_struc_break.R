# scripts/03_breaks_bai_perron.R
# Runs Bai–Perron multiple-break tests on gold returns across windows
# Saves figures to figures/03_breaks

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(stringr)
  library(strucchange)  # Bai–Perron tests
  library(purrr)
  library(readr)
})

# ── Config ─────────────────────────────────────────────────────────────────────
fig_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/03_breaks"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Toggle for a quick pass (skips CI calc and uses bigger minimum segment)
FAST_PASS <- FALSE

# Windows (consistent with your outline)
windows <- list(
  FULL      = c(as.Date("1990-01-01"), NA),  # NA = max date in data
  POST_GFC  = c(as.Date("2010-01-01"), NA),
  RECENT    = c(as.Date("2018-01-01"), NA)
)

# Minimum regime size as a fraction of window length
min_frac_by_window <- c(FULL = 0.10, POST_GFC = 0.10, RECENT = 0.10)
if (FAST_PASS) {
  min_frac_by_window <- c(FULL = 0.20, POST_GFC = 0.15, RECENT = 0.15)
}

# ── Input ──────────────────────────────────────────────────────────────────────
# Uses the tibble `gold_full_cleaned` already in memory
req_cols <- c("date","close_log","close_log_return")
stopifnot(all(req_cols %in% names(gold_full_cleaned)))

gold0 <- gold_full_cleaned %>%
  transmute(
    date = as.Date(date),
    close_log = close_log,
    r = close_log_return
  ) %>%
  arrange(date) %>%
  filter(!is.na(r))

# ── Helpers ────────────────────────────────────────────────────────────────────
run_bai_perron <- function(df, min_frac = 0.10, with_ci = TRUE) {
  n <- nrow(df)
  if (n < 30) {
    return(list(
      fsup = NA, bp_all = NULL, bp_fit = NULL, m_hat = 0,
      breaks_tbl = tibble()
    ))
  }
  h <- max(ceiling(n * min_frac), 5L)
  
  # supF test
  t <- Sys.time()
  fsup <- tryCatch(
    sctest(r ~ 1, type = "supF", data = df, from = h, to = n - h),
    error = function(e) NA
  )
  message("  supF: ", round(difftime(Sys.time(), t, "secs"), 2), "s")
  
  # dynamic-programming path and selected breaks by BIC
  t <- Sys.time()
  bp_all <- breakpoints(r ~ 1, data = df, h = h)
  m_hat  <- which.min(BIC(bp_all)) - 1L
  bp_fit <- breakpoints(r ~ 1, data = df, h = h, breaks = m_hat)
  message("  breakpoints (path+fit): ", round(difftime(Sys.time(), t, "secs"), 2), "s; m_hat=", m_hat)
  
  # confidence intervals (optional / may be slow)
  breaks_tbl <- tibble()
  if (with_ci && m_hat > 0) {
    t <- Sys.time()
    ci <- tryCatch(confint(bp_fit), error = function(e) NULL)
    if (!is.null(ci) && !is.null(ci$confint)) {
      idx <- ci$confint
      if (!is.null(dim(idx)) && nrow(idx) > 0) {
        breaks_tbl <- tibble(
          break_no   = seq_len(nrow(idx)),
          i_lower    = idx[, "lower"],
          i_break    = idx[, "breakpoints"],
          i_upper    = idx[, "upper"],
          date_lower = df$date[pmax(1, pmin(n, idx[, "lower"]))],
          date_break = df$date[pmax(1, pmin(n, idx[, "breakpoints"]))],
          date_upper = df$date[pmax(1, pmin(n, idx[, "upper"]))]
        )
      }
    }
    message("  confint: ", round(difftime(Sys.time(), t, "secs"), 2), "s")
  } else if (m_hat > 0) {
    bk <- bp_fit$breakpoints
    breaks_tbl <- tibble(
      break_no   = seq_along(bk),
      i_break    = bk,
      date_break = df$date[bk]
    )
  }
  
  list(fsup = fsup, bp_all = bp_all, bp_fit = bp_fit, m_hat = m_hat,
       breaks_tbl = breaks_tbl)
}

theme_ts <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6))
    )
}

plot_returns_with_breaks <- function(df, brk, window_name, m_hat) {
  ggplot(df, aes(date, r)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
    { if (nrow(brk) > 0 && all(c("date_lower","date_upper") %in% names(brk)))
      geom_rect(data = brk,
                aes(xmin = date_lower, xmax = date_upper, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, alpha = 0.08)
    } +
    geom_line(linewidth = 0.4) +
    { if (nrow(brk) > 0)
      geom_vline(aes(xintercept = date_break), data = brk,
                 linetype = "dashed", linewidth = 0.4)
    } +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = sprintf("Gold log returns — Bai–Perron mean breaks (%s)", window_name),
      subtitle = sprintf("Selected breaks m = %d", m_hat),
      x = "Date",
      y = "Log return",
      caption = "Source: Stooq (XAUUSD daily close, USD/oz). Shaded bands: 95% CI (if shown)."
    ) +
    theme_ts()
}

plot_price_with_breaks <- function(df, brk, window_name) {
  ggplot(df, aes(date, close_log)) +
    { if (nrow(brk) > 0 && all(c("date_lower","date_upper") %in% names(brk)))
      geom_rect(data = brk,
                aes(xmin = date_lower, xmax = date_upper, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, alpha = 0.08)
    } +
    geom_line(linewidth = 0.5) +
    { if (nrow(brk) > 0)
      geom_vline(aes(xintercept = date_break), data = brk,
                 linetype = "dashed", linewidth = 0.4)
    } +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = sprintf("Gold log price — break lines from returns (%s)", window_name),
      x = "Date",
      y = "log(Price, USD/oz)",
      caption = "Source: Stooq (XAUUSD daily close, USD/oz). Shaded bands: 95% CI (if shown)."
    ) +
    theme_ts()
}

save_figs <- function(p_ret, p_price, nm) {
  ggsave(file.path(fig_dir, paste0("returns_breaks_", tolower(nm), ".png")),
         p_ret, width = 11, height = 4.6, dpi = 300)
  ggsave(file.path(fig_dir, paste0("logprice_with_returns_breaks_", tolower(nm), ".png")),
         p_price, width = 11, height = 4.6, dpi = 300)
}

# ── Run and save ───────────────────────────────────────────────────────────────
all_results <- imap(windows, function(lims, nm) {
  message("== ", nm, " ==")
  t_win <- Sys.time()
  
  lim_start <- lims[1]
  lim_end   <- if (is.na(lims[2])) max(gold0$date) else lims[2]
  
  dfw <- gold0 %>%
    filter(date >= lim_start, date <= lim_end) %>%
    arrange(date)
  
  message("  window: ", lim_start, " to ", lim_end, " (n=", nrow(dfw), ")")
  
  # choose min_frac and whether to compute CIs
  min_frac <- unname(min_frac_by_window[[nm]])
  with_ci  <- !FAST_PASS
  
  t0 <- Sys.time()
  res <- run_bai_perron(dfw, min_frac = min_frac, with_ci = with_ci)
  message("  total window time: ", round(difftime(Sys.time(), t0, "secs"), 2), "s")
  
  # plots
  p_ret   <- plot_returns_with_breaks(dfw, res$breaks_tbl, nm, res$m_hat)
  p_price <- plot_price_with_breaks(
    df = dfw %>% select(date, close_log),
    brk = res$breaks_tbl,
    window_name = nm
  )
  save_figs(p_ret, p_price, nm)
  
  # compact tibble for tables
  out <- tibble(window = nm, m_hat = res$m_hat) %>%
    bind_cols(res$breaks_tbl)
  
  message("== ", nm, " done in ", round(difftime(Sys.time(), t_win, "secs"), 2), "s")
  out
})

breaks_summary <- bind_rows(all_results)

# CSV output (for your Table 2)
readr::write_csv(breaks_summary,
                 file.path(fig_dir, if (FAST_PASS)
                   "bai_perron_break_dates_fast.csv"
                   else
                     "bai_perron_break_dates_with_ci.csv")
)

message("Done. Figures + CSV saved in: ", fig_dir,
        if (FAST_PASS) " [FAST PASS]" else "")
