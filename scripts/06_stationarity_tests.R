# --- deps
library(dplyr); library(stringr); library(purrr)
library(urca);  library(tseries); library(gt)

# pretty window label
label_window <- function(name, start, end) {
  paste0(str_to_title(gsub("_"," ", name)),
         " (", format(start), " to ", format(end), ")")
}

# friendly series names for titles
friendly_name <- function(series_name) {
  c(close_log     = "Gold: Log Price",
    close_log_ret = "Gold: Daily Log Return")[series_name] %||% series_name
}

# map deterministic part -> intended τ row
tau_row_name <- function(type) c(trend="tau3", drift="tau2", none="tau1")[type]

# robust extractor for τ stat + CVs (+ k and a convenience p-value)
adf_urca_row <- function(x, type){
  x <- stats::na.omit(as.numeric(x))
  if (length(x) < 25) return(NULL)
  
  fit     <- ur.df(x, type = type, selectlags = "AIC")
  tau_key <- tau_row_name(type)
  
  # ----- critical values (always in @cval)
  cv_1 <- unname(fit@cval[tau_key, "1pct"])
  cv_5 <- unname(fit@cval[tau_key, "5pct"])
  cv_10<- unname(fit@cval[tau_key, "10pct"])
  
  # ----- τ statistic (robust to structure/naming quirks)
  ts_vec <- fit@teststat
  # if names missing, try to backfill from cval rownames
  if (is.null(names(ts_vec)) || !tau_key %in% names(ts_vec)) {
    nms <- names(ts_vec)
    if (is.null(nms)) {
      # fallback: assume first τ is the one present; align by cval order
      nms <- rownames(fit@cval)
      if (!is.null(nms) && length(nms) >= length(ts_vec)) names(ts_vec) <- nms[seq_along(ts_vec)]
    }
  }
  # final pick; if still missing, choose the first name that starts with "tau"
  tau_stat <- suppressWarnings(unname(ts_vec[tau_key]))
  if (is.na(tau_stat)) {
    cand <- names(ts_vec)[grepl("^tau", names(ts_vec))]
    if (length(cand)) tau_stat <- unname(ts_vec[cand[1]])
  }
  
  # convenience p-value from tseries (decision should rely on CVs)
  pval <- suppressWarnings(tryCatch(as.numeric(adf.test(x, k = fit@lags)$p.value),
                                    error = function(e) NA_real_))
  
  tibble(
    k        = fit@lags,
    tau_stat = as.numeric(tau_stat),
    cv_1pct  = cv_1, cv_5pct = cv_5, cv_10pct = cv_10,
    p_value  = pval
  )
}

# build a vertical 11-row block for ONE window
one_window_block <- function(wname, series_name, col, det, data, windows){
  rng <- windows[[wname]]
  sub <- data %>% filter(date >= rng[1], date <= rng[2])
  res <- adf_urca_row(sub[[col]], det)
  if (is.null(res)) return(NULL)
  
  decision <- ifelse(res$tau_stat < res$cv_5pct, "Reject unit root", "Fail to reject")
  
  tibble(
    window_label = label_window(wname, rng[1], rng[2]),
    Metric = c("Series","Deterministic","Sample","Observations",
               "k (lagged Δ terms)","ADF τ statistic",
               "Critical value (1%)","Critical value (5%)","Critical value (10%)",
               "p-value (tseries, approx.)","Decision @ 5%"),
    Value  = c(friendly_name(series_name),
               str_to_title(det),
               paste(format(rng[1]), "to", format(rng[2])),
               format(nrow(sub), big.mark = ","),
               as.character(res$k),
               sprintf("%.3f", res$tau_stat),
               sprintf("%.2f",  res$cv_1pct),
               sprintf("%.2f",  res$cv_5pct),
               sprintf("%.2f",  res$cv_10pct),
               ifelse(is.na(res$p_value), "—", sprintf("%.3f", res$p_value)),
               decision)
  )
}

# main wrapper: vertical gt for a series across 3 windows (highlight group headers)
make_adf_gt_for_series <- function(series_name, col, det,
                                   window_names = c("full","post_gfc","recent"),
                                   data = df, windows = windows){
  
  blocks <- map_dfr(window_names,
                    ~ one_window_block(.x, series_name, col, det, data, windows))
  
  gt(blocks, groupname_col = "window_label") |>
    cols_label(Metric = "", Value = "") |>
    tab_header(title = md(paste0("ADF Summary — **", friendly_name(series_name), "**"))) |>
    opt_row_striping() |>
    # highlight each window header
    tab_style(
      style = list(cell_fill(color = "#f5f7ff"),
                   cell_text(weight = "bold")),
      locations = cells_row_groups()
    ) |>
    tab_footnote(
      footnote = md("Decision uses **τ vs 5% critical value** from **urca**; the `tseries` p-value is shown only for reference."),
      locations = cells_title("title")
    ) |>
    tab_options(table.width = pct(66))
}

# ---- build the two tables (kept windows: Full, Post-GFC, Recent)
gt_price  <- make_adf_gt_for_series("close_log",     "close_log",        det = "trend",
                                    window_names = c("full","post_gfc","recent"),
                                    data = df, windows = windows)
gt_return <- make_adf_gt_for_series("close_log_ret", "close_log_return", det = "drift",
                                    window_names = c("full","post_gfc","recent"),
                                    data = df, windows = windows)

gt_price
gt_return


# Ensure a PNG-capable backend is available
if (!requireNamespace("webshot2", quietly = TRUE) &&
    !requireNamespace("chromote", quietly = TRUE)) {
  stop("Please install either {webshot2} or {chromote} to save gt tables as PNGs.")
}

fig_dir <- "figures/02_diagnostics"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Tweak width/height for crisp PNGs; adjust if you want larger/smaller images
gt::gtsave(gt_price,
           filename = file.path(fig_dir, "adf_log_price.png"),
           vwidth   = 1100, vheight = 900)

gt::gtsave(gt_return,
           filename = file.path(fig_dir, "adf_log_returns.png"),
           vwidth   = 1100, vheight = 900)






