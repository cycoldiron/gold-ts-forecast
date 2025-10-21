# ================================================================
# ADF Summary (from gold_full_cleaned) -> PNG
# - Window labels show only dates (e.g., "Jan 1990 – Oct 2025")
# - ADF stat shows τ (urca) or p-value (tseries fallback)
# - All numeric displays rounded to 2 decimals
# - Output: /Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/02_adf/adf_summary.png
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(purrr)
  library(tidyr); library(rlang)
  library(urca);  library(tseries)
  library(gt);    library(scales); library(fs)
})

# ---------- Helpers ----------
schwert_max <- function(T) floor(12 * (T/100)^(1/4))
friendly    <- c(close_log = "Gold: Log Price", close_log_ret = "Gold: Daily Log Return")
tau_key     <- function(type) c(trend="tau3", drift="tau2", none="tau1")[type]

# label_style: "bY" -> "Jan 1990", "mY" -> "01/1990"
window_label <- function(start, end, style = c("bY","mY")){
  style <- match.arg(style)
  fmt <- if (style == "bY") "%b %Y" else "%m/%Y"
  paste0(format(as.Date(start), fmt), " \u2013 ", format(as.Date(end), fmt))
}

# ADF with robust fallbacks:
# 1) urca::ur.df with selectlags="AIC" (cap = Schwert)
# 2) urca fixed-lag step-down if needed
# 3) tseries::adf.test fallback (use p-value)
adf_row <- function(x, type){
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 25) return(NULL)
  
  maxlag <- schwert_max(length(x)); tk <- tau_key(type)
  
  try_urdf <- function(lags, mode = c("AIC","Fixed")){
    mode <- match.arg(mode)
    o <- try(ur.df(x, type = type, lags = lags, selectlags = mode), silent = TRUE)
    if (inherits(o, "try-error")) return(NULL)
    tau <- suppressWarnings(as.numeric(unname(o@teststat[tk])))
    cv5 <- suppressWarnings(as.numeric(unname(o@cval[tk, "5pct"])))
    if (!is.finite(tau) || !is.finite(cv5)) return(NULL)
    tibble(source = "urca", k = o@lags, tau = tau, cv5 = cv5, pval = NA_real_)
  }
  
  # (1) urca, AIC-selected
  res <- try_urdf(maxlag, "AIC")
  if (!is.null(res)) return(res)
  
  # (2) urca, fixed lags step-down
  for (kk in seq.int(min(maxlag, 12L), 0L, by = -1L)) {
    res2 <- try_urdf(kk, "Fixed")
    if (!is.null(res2)) return(res2)
  }
  
  # (3) tseries fallback — report p-value, and k ≈ T^(1/3)
  k_fallback <- floor((length(x) - 1)^(1/3))
  at <- try(adf.test(x, k = k_fallback), silent = TRUE)
  if (inherits(at, "try-error")) return(NULL)
  tibble(source = "tseries", k = k_fallback, tau = NA_real_, cv5 = NA_real_, pval = at$p.value)
}

normalize_windows <- function(df, windows = NULL) {
  if (is.null(windows)) {
    end_date <- max(df$date, na.rm = TRUE)
    windows <- list(
      full     = c(as.Date("1990-01-01"), end_date),
      post_gfc = c(as.Date("2010-01-01"), end_date),
      recent   = c(as.Date("2018-01-01"), end_date),
      last1y   = c(end_date - 365,        end_date),
      last30d  = c(end_date - 30,         end_date)
    )
  }
  windows <- lapply(windows, as.Date)
  windows <- windows[vapply(windows, function(x) length(x)==2 && all(!is.na(x)), logical(1))]
  if (!length(windows)) abort("No valid windows supplied/built.")
  windows
}

# Clip each window to data range; require >=25 finite obs
summarise_window <- function(df, winlab, series_key, col, det, windows, label_style = "bY"){
  data_min <- min(df$date[is.finite(df[[col]])], na.rm = TRUE)
  data_max <- max(df$date[is.finite(df[[col]])], na.rm = TRUE)
  rng_in   <- windows[[winlab]]
  rng      <- c(max(rng_in[1], data_min), min(rng_in[2], data_max))
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] > rng[2]) return(NULL)
  
  sub <- df %>% dplyr::filter(date >= rng[1], date <= rng[2])
  if (sum(is.finite(sub[[col]])) < 25) return(NULL)
  
  res <- adf_row(sub[[col]], det); if (is.null(res)) return(NULL)
  
  # Decision + display stat (τ for urca; p-value for tseries), both to 2 decimals
  decision <- ifelse(res$source == "urca",
                     ifelse(res$tau < res$cv5, "Reject unit root", "Fail to reject"),
                     ifelse(res$pval < 0.05,   "Reject unit root", "Fail to reject"))
  stat_disp <- ifelse(res$source == "urca",
                      sprintf("%.2f", res$tau),
                      sprintf("p=%.2f", res$pval))
  
  tibble(
    Series   = friendly[[series_key]] %||% series_key,
    Window   = window_label(rng[1], rng[2], style = label_style),  # dates only
    Det      = str_to_title(det),
    Stat     = stat_disp,
    k        = sprintf("%.2f", as.numeric(res$k)),                 # round k
    Decision = decision
  )
}

make_adf_compact <- function(df, windows = NULL, label_style = "bY"){
  stopifnot(all(c("date","close_log","close_log_return") %in% names(df)))
  df <- df |> dplyr::mutate(date = as.Date(date)) |> dplyr::arrange(date)
  windows <- normalize_windows(df, windows)
  
  order_pref <- c("full","post_gfc","recent","last1y","last30d")
  win_names  <- intersect(order_pref, names(windows))
  if (!length(win_names)) win_names <- names(windows)
  
  out <- dplyr::bind_rows(
    purrr::map(win_names, ~ summarise_window(df, .x, "close_log_ret", "close_log_return", "drift", windows, label_style)),
    purrr::map(win_names, ~ summarise_window(df, .x, "close_log",      "close_log",         "trend", windows, label_style))
  ) |> dplyr::bind_rows()
  
  if (!nrow(out)) rlang::abort("All windows too short or ADF failed across all windows/series.")
  
  out |>
    gt::gt(groupname_col = "Series") |>
    gt::cols_label(
      Window   = "",
      Det      = "Deterministic",
      Stat     = html("<i>ADF stat</i>"),
      k        = "Lag k",
      Decision = "Decision @ 5%"
    ) |>
    gt::tab_header(title = gt::md("Stationarity summary — **Gold prices vs returns**")) |>
    gt::opt_row_striping() |>
    gt::data_color(columns = "Decision",
                   colors  = scales::col_factor(c("#0a7f3f","#c9252d"),
                                                levels = c("Reject unit root","Fail to reject"))) |>
    gt::fmt_markdown(columns = "Stat") |>
    gt::tab_style(style = list(gt::cell_fill(color = "#f5f7ff"), gt::cell_text(weight = "bold")),
                  locations = gt::cells_row_groups()) |>
    gt::tab_footnote(
      footnote = gt::md("Stat shows **τ** (urca) or **p-value** (tseries) when urca fails. Lags: Schwert cap with AIC for urca; k≈T^(1/3) for tseries."),
      locations = gt::cells_title("title")
    ) |>
    gt::tab_options(table.width = gt::pct(66))
}

# ---------- Build df from gold_full_cleaned (robust to column names) ----------
stopifnot(exists("gold_full_cleaned"))
g <- gold_full_cleaned %>% mutate(date = as.Date(date)) %>% arrange(date)

price_cols_log   <- intersect(c("close_log","log_price","p_log"), names(g))
price_cols_level <- intersect(c("close","price","close_usd","p"), names(g))
if (length(price_cols_log)) {
  close_log_vec <- suppressWarnings(as.numeric(g[[price_cols_log[1]]]))
} else if (length(price_cols_level)) {
  close_log_vec <- log(suppressWarnings(as.numeric(g[[price_cols_level[1]]])))
} else {
  close_log_vec <- rep(NA_real_, nrow(g))
}
ret_cols <- intersect(c("close_log_return","log_return","log_ret","r","ret"), names(g))
if (length(ret_cols)) {
  ret_vec <- suppressWarnings(as.numeric(g[[ret_cols[1]]]))
} else if (any(is.finite(close_log_vec))) {
  ret_vec <- c(NA_real_, diff(close_log_vec))
} else {
  ret_vec <- rep(NA_real_, nrow(g))
}

df <- tibble(date = g$date, close_log = close_log_vec, close_log_return = ret_vec)

# ---------- Windows (Option 2) ----------
end_date <- max(df$date[is.finite(df$close_log) | is.finite(df$close_log_return)], na.rm = TRUE)
windows  <- list(
  full     = c(as.Date("1990-01-01"), end_date),
  post_gfc = c(as.Date("2010-01-01"), end_date),
  recent   = c(as.Date("2018-01-01"), end_date),
  last1y   = c(end_date - 365,        end_date),
  last30d  = c(end_date - 30,         end_date)
)

# ---------- Build & save ----------
# label_style: "bY" => "Jan 1990 – Oct 2025"; use "mY" for "01/1990 – 10/2025"
gt_adf <- make_adf_compact(df, windows = windows, label_style = "bY")

outdir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/02_adf"
fs::dir_create(outdir, recurse = TRUE)
out_png <- file.path(outdir, "adf_summary.png")
gt::gtsave(gt_adf, filename = out_png)

message("Saved ADF table to: ", out_png)
# If PNG rendering errors: install.packages("webshot2"); webshot2::install_phantomjs()
