suppressPackageStartupMessages({
  library(dplyr)
  library(gt)
  library(purrr)
  library(scales)
})

# Reuse your summary; build it if missing (from `results`)
if (!exists("headers_summary")) {
  if (!exists("results")) stop("Need `results` or `headers_summary` in memory.")
  headers_summary <- bind_rows(lapply(results, `[[`, "header"))
}

# ---- Robust date coercion for window labels ----
.to_date <- function(x) {
  if (length(x) == 0 || is.na(x)) return(NA_Date_)
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
  # try character
  out <- suppressWarnings(as.Date(x))
  if (is.na(out)) NA_Date_ else out
}

pretty_window_label <- function(lims) {
  start <- .to_date(lims[1])
  end   <- .to_date(lims[2])
  # Use max(gold0$date) if end is NA (i.e., "present")
  if (is.na(end)) sprintf("%s–Pres", format(start, "%Y"))
  else sprintf("%s–%s", format(start, "%Y"), format(end, "%Y"))
}

# Build nice labels in the order of your windows list
win_labels <- map_chr(windows, pretty_window_label)

# ---- Assemble display table with your requested edits ----
tbl_disp <- headers_summary %>%
  mutate(
    Window = factor(names(windows)[match(window, names(windows))],
                    levels = names(windows),
                    labels = win_labels),
    `n (obs)`          = n,
    `h (min seg, obs)` = h,        # h = minimal segment length used in search
    `supF p`           = supF_p,
    `Num. breaks`      = m_hat,
    Decision           = ifelse(m_hat > 0, "break(s)", "no break")
  ) %>%
  select(Window, `n (obs)`, `h (min seg, obs)`, `supF p`, `Num. breaks`, Decision) %>%
  arrange(Window)

gt_tbl <- tbl_disp %>%
  gt() %>%
  tab_header(
    title = md("**Bai–Perron mean-break test on gold log returns**"),
    subtitle = md("Windows as defined in outline. **h** = minimal segment length (observations) used in the search.")
  ) %>%
  fmt_number(columns = c(`n (obs)`, `h (min seg, obs)`), decimals = 0, use_seps = TRUE) %>%
  fmt_number(columns = `supF p`, decimals = 3) %>%
  cols_align(columns = everything(), align = "center") %>%
  # Wider columns
  cols_width(
    Window ~ px(180),
    `n (obs)` ~ px(130),
    `h (min seg, obs)` ~ px(190),
    `supF p` ~ px(110),
    `Num. breaks` ~ px(140),
    Decision ~ px(140)
  ) %>%
  data_color(
    columns = Decision,
    colors = col_factor(
      palette = c("no break" = "#FCA5A5", "break(s)" = "#86EFAC"),
      domain  = c("no break", "break(s)")
    )
  ) %>%
  tab_footnote(footnote = "Num. breaks selected by BIC.") %>%
  tab_source_note(source_note = "Source: Stooq (XAUUSD daily close, USD/oz). Bai–Perron via strucchange.") %>%
  tab_options(table.width = px(980))

# Save to your figures directory
if (!exists("fig_dir")) {
  fig_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/03_breaks"
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
}
gt_out_path <- file.path(fig_dir, "table_bai_perron_returns_gt.png")
gtsave(gt_tbl, gt_out_path, vwidth = 1300, vheight = 520, expand = 10)

message("Saved updated gt table to: ", gt_out_path)

