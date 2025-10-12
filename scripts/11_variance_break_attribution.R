# =============================================================================
# 03_breaks_attribution.R
# Classify variance changepoints as Positive-led or Negative-led
# (keeping your PELT-on-returns approach used in the figures)
# Outputs: one GT table + one bar chart per window (PNG only)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(changepoint); library(slider); library(gt); library(fs)
})

# ---------- Paths -------------------------------------------------------------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/03_breaks"
dir_create(out_dir, recurse = TRUE)

# ---------- Colors (match your figures) --------------------------------------
col_bar_pos <- "#2ca02c"  # green
col_bar_neg <- "#d62728"  # red
col_bar_unc <- "#7f7f7f"  # grey

# ---------- Your helper (unchanged) ------------------------------------------
get_window_df <- function(start_date) {
  gold_full_cleaned %>%
    filter(date >= as.Date(start_date)) %>%
    arrange(date) %>%
    filter(!is.na(close_log_return)) %>%
    transmute(
      Date = as.Date(date),
      r    = close_log_return,
      r2   = close_log_return^2
    ) |>
    mutate(r2_roll30 = slide_dbl(r2, mean, .before = 29, .complete = TRUE))
}

# Two-year window start (matches your plotting script)
ret_idx  <- gold_full_cleaned %>% arrange(date) %>% filter(!is.na(close_log_return)) %>% pull(date)
start_2y <- if (length(ret_idx) >= 504) ret_idx[length(ret_idx) - 504 + 1] else max(gold_full_cleaned$date, na.rm = TRUE) - 730

# ---------- Detect breaks (same as your plots) -------------------------------
detect_breaks <- function(r) {
  # You used: cpt.var(df$r, method = "PELT")
  # Keep defaults for penalty to mirror the figure.
  obj <- tryCatch(cpt.var(r, method = "PELT"), error = function(e) NULL)
  if (is.null(obj)) integer(0) else cpts(obj)
}

# ---------- Classify each break ----------------------------------------------
classify_breaks <- function(df, cp_idx, h = 5) {
  if (length(cp_idx) == 0L) {
    return(tibble())
  }
  r <- df$r; dates <- df$Date
  pmap_dfr(list(k = seq_along(cp_idx), i = cp_idx), function(k, i) {
    i0 <- max(1, i - h); i1 <- min(length(r), i + h)
    rw <- r[i0:i1]; dw <- dates[i0:i1]
    
    # largest absolute return in window
    j_rel   <- which.max(abs(rw))
    r_star  <- rw[j_rel]; date_star <- dw[j_rel]
    sign_max <- ifelse(r_star > 0, "Positive", ifelse(r_star < 0, "Negative", "Zero"))
    
    # share of local variance by sign
    pos_var <- sum((rw[rw > 0])^2); neg_var <- sum((rw[rw < 0])^2)
    tot_var <- pos_var + neg_var
    pct_pos <- ifelse(tot_var > 0, pos_var / tot_var, NA_real_)
    
    # break-day return and sign
    r_break   <- r[i]; sign_break <- ifelse(r_break > 0, "Positive", ifelse(r_break < 0, "Negative", "Zero"))
    
    label <- dplyr::case_when(
      sign_max == "Positive" ~ "Positive-led",
      sign_max == "Negative" ~ "Negative-led",
      TRUE                   ~ "Unclear"
    )
    
    tibble(
      Break = k,
      Break_Date = dates[i],
      Break_Return = r_break,
      Break_Sign = sign_break,
      MaxAbs_Return = r_star,
      MaxAbs_Date = date_star,
      MaxAbs_Sign = sign_max,
      Pct_Variance_Positive = pct_pos,
      Classification = label
    
    )
  })
}

# ---------- Summarize & plot helpers -----------------------------------------
summarize_class <- function(tbl) {
  if (!nrow(tbl)) return(list(summary = tibble(n_breaks = 0, n_pos_led = 0, share_pos_led = NA_real_, p_binom = NA_real_), bar = NULL))
  sumrow <- tbl %>%
    mutate(is_pos = Classification == "Positive-led") %>%
    summarise(n_breaks = n(),
              n_pos_led = sum(is_pos),
              share_pos_led = n_pos_led / n_breaks)
  test <- binom.test(sumrow$n_pos_led, sumrow$n_breaks, p = 0.5)
  sumrow$p_binom <- unname(test$p.value)
  list(summary = sumrow, bar = NULL)
}

save_gt_tbl <- function(tbl, title, subtitle, filebase, add_summary = NULL) {
  g <- tbl %>%
    mutate(
      Break_Return = round(Break_Return, 4),
      MaxAbs_Return = round(MaxAbs_Return, 4),
      Pct_Variance_Positive = ifelse(is.na(Pct_Variance_Positive), NA, round(100 * Pct_Variance_Positive, 1))
    ) %>%
    gt() %>%
    tab_header(title = md(paste0("**", title, "**")),
               subtitle = subtitle) %>%
    # label for display only
    cols_label(
      Break                  = "Break #",
      Break_Date             = "Break date",
      Break_Return           = "Return on break day",
      Break_Sign             = "Break sign",
      MaxAbs_Return          = "Max |return| in window",
      MaxAbs_Date            = "Date of max |return|",
      MaxAbs_Sign            = "Sign of max |return|",
      Pct_Variance_Positive  = "Positive variance share (%)",
      Classification         = "Classification"
    ) %>%
    opt_row_striping() %>%
    # IMPORTANT: use data column names here, not the labels
    cols_width(
      Break ~ px(90),
      Break_Date ~ px(140),
      Break_Return ~ px(160),
      Break_Sign ~ px(130),
      MaxAbs_Return ~ px(180),
      MaxAbs_Date ~ px(160),
      MaxAbs_Sign ~ px(160),
      Pct_Variance_Positive ~ px(220),
      Classification ~ px(160)
    )
  
  if (!is.null(add_summary)) {
    note_line <- sprintf(
      "Positive-led: %d / %d (%.1f%%), Binomial p = %.3f",
      add_summary$n_pos_led, add_summary$n_breaks,
      100 * add_summary$share_pos_led, add_summary$p_binom
    )
    g <- g %>%
      tab_footnote(footnote = note_line, locations = cells_title(groups = "title"))
  }
  
  gtsave(g, file.path(out_dir, paste0(filebase, ".png")), vwidth = 1400, vheight = 800)
}


save_bar <- function(tbl, title, filebase) {
  if (!nrow(tbl)) return(invisible(NULL))
  p <- tbl %>%
    count(Classification) %>%
    mutate(Classification = factor(Classification, levels = c("Positive-led","Negative-led","Unclear"))) %>%
    ggplot(aes(Classification, n, fill = Classification)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Positive-led" = col_bar_pos,
                                 "Negative-led" = col_bar_neg,
                                 "Unclear" = col_bar_unc), guide = "none") +
    labs(title = title, x = "Classification", y = "Count") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(out_dir, paste0(filebase, "_bar.png")), p, width = 8.5, height = 5.2, dpi = 300, bg = "white")
}

# ---------- Windows to analyze (align with your figs) ------------------------
windows <- list(
  `1990-present` = min(gold_full_cleaned$date, na.rm = TRUE),
  `2010-present` = as.Date("2010-01-01"),
  `2018-present` = as.Date("2018-01-01"),
  `Past-2y`      = start_2y
)

# ---------- Run attribution per window ---------------------------------------
for (nm in names(windows)) {
  df <- get_window_df(windows[[nm]])
  if (nrow(df) < 30) next
  
  cp_idx <- detect_breaks(df$r)
  attr_tbl <- classify_breaks(df, cp_idx, h = 5)  # change h for sensitivity
  
  # summary + binomial
  sumlist <- summarize_class(attr_tbl)
  summ    <- sumlist$summary
  
  # save table + bar chart
  ttl <- paste0("Attribution of Variance Breaks — ", nm)
  sub <- paste0("Classification uses sign of the largest |return| in ±5 trading days around each break")
  base <- paste0("break_attribution_", gsub("[^A-Za-z0-9]+", "_", nm))
  
  save_gt_tbl(attr_tbl, ttl, sub, base, add_summary = summ)
  save_bar(attr_tbl, ttl, base)
}

message("Saved attribution tables and bars in: ", out_dir)
