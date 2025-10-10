library(dplyr)
library(lubridate)
library(zoo)
library(janitor)
library(rlang)

# --- EPU (monthly) ------------------------------------------------------------
# epu_raw columns seen: Year, Month, News_Based_Policy_Uncert_Index
epu_full_clean <- epu_raw %>%
  clean_names() %>%
  transmute(
    year   = as.integer(year),
    month  = as.integer(month),
    date_m = as.Date(as.yearmon(sprintf("%04d-%02d", year, month))),
    epu    = as.numeric(news_based_policy_uncert_index)
  ) %>%
  filter(date_m >= as.Date("1990-01-01")) %>%
  arrange(date_m)

# --- GPR (monthly) ------------------------------------------------------------
# gpr_raw columns vary by file; we coalesce "new" (GPR/GPRT/GPRA) and "historical" (GPRH/GPRHT/GPRHA)
# 'month' in your str() is POSIXct monthly; we convert to yearmon -> Date (month end)
gpr_full_clean <- gpr_raw %>%
  clean_names() %>%
  mutate(date_m = as.Date(as.yearmon(month))) %>%
  transmute(
    date_m,
    year        = year(date_m),
    month       = month(date_m),
    gpr         = coalesce(as.numeric(.data[["gpr"]]),  as.numeric(.data[["gprh"]])),
    gpr_threats = coalesce(as.numeric(.data[["gprt"]]), as.numeric(.data[["gprht"]])),
    gpr_acts    = coalesce(as.numeric(.data[["gpra"]]), as.numeric(.data[["gprha"]]))
  ) %>%
  filter(date_m >= as.Date("1990-01-01")) %>%
  arrange(date_m)

# --- OPTIONAL: VIX (daily) ----------------------------------------------------
# If you have vix_raw loaded, this standardizes to snake case with a clean 'date' and 'vix' column
# (Rename the source columns below if yours differ.)
# vix_clean <- vix_raw %>%
#   clean_names() %>%
#   mutate(date = as.Date(date)) %>%
#   transmute(
#     date,
#     vix = as.numeric(vix)   # adjust the source column name if needed
#   ) %>%
#   filter(date >= as.Date("1990-01-01")) %>%
#   arrange(date)

# --- (Optional) Save cleaned monthly proxies ----------------------------------
# dir.create("data/clean", recursive = TRUE, showWarnings = FALSE)
# save(epu_clean, gpr_clean, file = file.path("data", "clean", "02_clean_uncertainty.Rdata"), compress = "bzip2")
