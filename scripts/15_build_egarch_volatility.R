# ================================================================
# 04_egarch_data: Build & save volatility series (no uncertainty data)
# Output: /Users/cycoldiron/Desktop/gold-ts-forecast/data/clean/04_egarch_data.RData
# Requires: gold_full_cleaned with columns date, close_log_return
# ================================================================
suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(xts); library(zoo)
  library(changepoint); library(rugarch); library(slider); library(fs); library(tidyr)
})

stopifnot(exists("gold_full_cleaned"))
df <- gold_full_cleaned %>%
  arrange(date) %>%
  transmute(date = as.Date(date),
            r    = close_log_return) %>%
  filter(is.finite(r))

stopifnot(nrow(df) > 500)
r_xts <- xts::xts(df$r, order.by = df$date)

# --- 1) Find variance changepoints & build 0/1 variance dummy (post first break)
cp_idx      <- tryCatch(cpts(changepoint::cpt.var(df$r, method = "PELT")), error = function(e) integer(0))
break_dates <- if (length(cp_idx)) df$date[cp_idx] else as.Date(character(0))
break_dummy <- if (length(break_dates)) as.integer(df$date >= min(break_dates)) else rep(0L, nrow(df))
Xvar        <- matrix(break_dummy, ncol = 1)

# --- 2) EGARCH(1,1)–t specs and fits (full sample for complete σ_t)
solver_ctrl <- list(trace = 0, rel.tol = 1e-8, x.tol = 1e-8, eval.max = 5000, iter.max = 5000)
fit_ctrl    <- list(stationarity = 1, fixed.se = 0)

spec_base <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)

spec_dummy <- ugarchspec(
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "eGARCH", garchOrder = c(1,1),
                        external.regressors = Xvar),
  distribution.model = "std"
)

fit_base  <- ugarchfit(spec_base,  data = r_xts,
                       solver = "nlminb", solver.control = solver_ctrl, fit.control = fit_ctrl)
fit_dummy <- ugarchfit(spec_dummy, data = r_xts,
                       solver = "nlminb", solver.control = solver_ctrl, fit.control = fit_ctrl)

sigma_series <- tibble(
  date        = df$date,
  r           = df$r,
  sigma_base  = as.numeric(sigma(fit_base)),
  sigma_dummy = as.numeric(sigma(fit_dummy)),
  break_dummy = break_dummy
)

# --- 3) Model-free realized volatility (21-day)
rv_series <- sigma_series %>%
  mutate(rv21 = sqrt(slider::slide_dbl(r^2, sum, .before = 20, .complete = TRUE)))

# --- 4) Save only volatility objects (no VIX/EPU)
out_dir  <- "/Users/cycoldiron/Desktop/gold-ts-forecast/data/clean"
fs::dir_create(out_dir, recurse = TRUE)
out_path <- file.path(out_dir, "04_egarch_data.RData")

save(sigma_series, rv_series, break_dates, fit_base, fit_dummy, file = out_path)
message("Saved volatility data to: ", out_path)
