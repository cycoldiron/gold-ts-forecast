# ================= Regime table from variance changepoints ====================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(lubridate); library(changepoint); library(fs)
})

# ---- Paths -------------------------------------------------------------------
save_path <- "/Users/cycoldiron/Desktop/gold-ts-forecast/data/clean/05_breaks_data.Rdata"
dir_create(path_dir(save_path), recurse = TRUE)

# ---- Expect `gold_full_cleaned` to exist with: date, close_log_return --------
stopifnot(all(c("date","close_log_return") %in% names(gold_full_cleaned)))

# ---- Prep full sample (drop NAs, order) --------------------------------------
returns_df <- gold_full_cleaned %>%
  arrange(date) %>%
  filter(!is.na(close_log_return)) %>%
  transmute(
    Date = as.Date(date),
    r    = close_log_return,
    r2   = close_log_return^2
  )

n_obs <- nrow(returns_df)
if (n_obs < 5) stop("Not enough non-missing returns to estimate changepoints.")

# ---- Changepoint detection on variance (PELT) --------------------------------
cp <- tryCatch(cpt.var(returns_df$r, method = "PELT", class = TRUE), error = function(e) NULL)
if (is.null(cp)) stop("cpt.var() failed; cannot construct regimes.")

cp_idx <- cpts(cp)                    # integer indices of break *endpoints*
# Add boundaries to convert to contiguous regimes
idx_starts <- c(1L, cp_idx + 1L)
idx_ends   <- c(cp_idx, n_obs)

# ---- Build regime tibble -----------------------------------------------------
regimes_tbl <- map2_dfr(idx_starts, idx_ends, function(i1, i2) {
  seg <- returns_df[i1:i2, , drop = FALSE]
  tibble(
    regime_id   = length(which(idx_ends <= i2)),  # sequential 1..K
    start_date  = seg$Date[1],
    end_date    = seg$Date[nrow(seg)],
    n_days      = nrow(seg),
    mean_variance = mean(seg$r2, na.rm = TRUE),   # E[r_t^2] within regime
    avg_volatility = sd(seg$r, na.rm = TRUE)      # sd(r_t) within regime
  )
})

# Optional: list of exact break dates (endpoints of preceding regime)
break_dates <- tibble(break_date = returns_df$Date[cp_idx])

# ---- Save --------------------------------------------------------------------
# Save both the regimes table and the raw break dates for convenience
save(regimes_tbl, break_dates, file = save_path)
message("Saved regimes table to: ", save_path)

# ---- (Optional) quick peek ---------------------------------------------------
# print(regimes_tbl)
