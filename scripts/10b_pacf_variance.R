# =============================================================================
# 06a_pacf_rsq_only.R
# PACF of squared returns to guide ARCH order
# Output: 03_pacf_rsq.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(fs)
})

# --------- Paths --------------------------------------------------------------
out_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/04_models/01_arch_garch"
dir_create(out_dir, recurse = TRUE)

# --------- Data ---------------------------------------------------------------
# Assumes `gold_full_cleaned` has: date, close_log_return
r <- na.omit(gold_full_cleaned$close_log_return)
n <- length(r)

# --------- PACF helper --------------------------------------------------------
pacf_df <- function(x, lag.max = 40) {
  p <- pacf(x, lag.max = lag.max, plot = FALSE)
  tibble(lag = 1:lag.max, pacf = as.numeric(p$acf))
}

# 95% CI (large-sample)
ci <- 1.96 / sqrt(n)

# --------- Plot: PACF of squared returns -------------------------------------
df_pacf_r2 <- pacf_df(r^2, 40)

g <- ggplot(df_pacf_r2, aes(lag, pacf)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = -ci, ymax = ci), alpha = 0.08) +
  geom_hline(yintercept = c(-ci, ci), linetype = "dashed") +
  geom_col(width = 0.8) +
  labs(
    title = "PACF — Squared Returns (ARCH Order Heuristic)",
    x = "Lag (trading days)",
    y = expression(PACF(r^2))
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "03_pacf_rsq.png"), g, width = 10, height = 7, dpi = 300)
message("Saved: ", file.path(out_dir, "03_pacf_rsq.png"))
