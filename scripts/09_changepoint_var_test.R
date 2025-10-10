# ---------------- Two-year figure fix: robust y-cap and explicit limits ---------------
library(dplyr)
library(ggplot2)
library(changepoint)
library(slider)
library(fs)

out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/03_breaks"
dir_create(out_dir, recurse = TRUE)

col_spikes <- "#B8B8B8"
col_roll   <- "#1F77B4"
col_break  <- "#D62728"

get_window_df <- function(start_date) {
  gold_full_cleaned %>%
    filter(date >= as.Date(start_date)) %>%
    arrange(date) %>%
    filter(!is.na(close_log_return)) %>%
    transmute(Date = as.Date(date),
              r = close_log_return,
              r2 = close_log_return^2) |>
    mutate(r2_roll30 = slide_dbl(r2, mean, .before = 29, .complete = TRUE))
}

# --- Start for past two years (already good in your last run)
ret_idx  <- gold_full_cleaned %>% arrange(date) %>% filter(!is.na(close_log_return)) %>% pull(date)
start_2y <- if (length(ret_idx) >= 504) ret_idx[length(ret_idx) - 504 + 1] else max(gold_full_cleaned$date, na.rm = TRUE) - 730

# --- Build data
df <- get_window_df(start_2y)

# ---- DEBUG (optional): print a few quantiles so you can see what's happening
qs <- quantile(df$r2, probs = c(0, .9, .95, .99, .999), na.rm = TRUE)
print(qs)

# ---- Robust cap (never 0/NA) + explicit limits
# Prefer the 99th percentile of strictly-positive r2; fall back to 99.9th; ensure epsilon floor
q_pos <- tryCatch(quantile(df$r2[df$r2 > 0], 0.99, na.rm = TRUE), error = function(e) NA_real_)
y_cap <- max(q_pos, quantile(df$r2, 0.999, na.rm = TRUE), 1e-8, na.rm = TRUE)

df <- df %>% mutate(y_plot = pmin(r2, y_cap))

y_lim_upper <- y_cap * 1.10  # give a little headroom so grid shows

# ---- Change-points (safe)
brk_dates <- as.Date(character())
if (nrow(df) >= 20) {
  cp <- tryCatch(cpt.var(df$r, method = "PELT"), error = function(e) NULL)
  if (!is.null(cp) && length(cpts(cp))) brk_dates <- df$Date[cpts(cp)]
}
bl <- data.frame(Date = brk_dates)

# ---- Plot
p_2y <- ggplot(df) +
  geom_segment(aes(x = Date, xend = Date, y = 0, yend = y_plot),
               linewidth = 0.22, alpha = 0.7, colour = col_spikes) +
  geom_line(aes(x = Date, y = r2_roll30),
            linewidth = 1.0, colour = col_roll, na.rm = TRUE) +
  { if (nrow(bl) > 0)
    geom_vline(data = bl, aes(xintercept = Date),
               linetype = "22", linewidth = 1.0, colour = col_break) } +
  scale_y_continuous(labels = scales::label_scientific(digits = 2),
                     limits = c(0, y_lim_upper), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "Variance breaks in gold log returns - Past Two Years",
    subtitle = "Vertical dashed lines mark estimated variance change-points (PELT) from daily returns",
    x        = "Date",
    y        = expression(r[t]^2~"(Squared daily log return)"),
    caption  = "Source: Stooq (XAUUSD)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5),
    plot.caption     = element_text(hjust = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

pngf <- file.path(out_dir, "variance-breaks-r2_past-2y.png")
ggsave(pngf, p_2y, width = 11, height = 4.2, dpi = 300, bg = "white")
message("Saved: ", pngf)
