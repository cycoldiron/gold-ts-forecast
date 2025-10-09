# Packages
library(readr)
library(dplyr)
library(janitor)
library(stringr)
library(lubridate)
library(fs)

# Paths
proj_dir  <- "/Users/cycoldiron/Desktop/gold-ts-forecast"
raw_dir   <- file.path(proj_dir, "data", "raw")
clean_dir <- file.path(proj_dir, "data", "cleaned")
dir_create(clean_dir, recurse = TRUE)

# Helper: basic cleaning
basic_clean <- function(df) {
  # standardize names, drop empty rows, trim blanks to NA
  df <- df |>
    janitor::clean_names() |>
    janitor::remove_empty("rows") |>
    mutate(across(where(is.character), ~ na_if(trimws(.x), "")))
  
  # find a date-like column if present
  nms <- names(df)
  date_guess <- nms[str_detect(nms, "(^|_)(date|observation_date|time|timestamp)($|_)")]
  
  if (length(date_guess) > 0) {
    date_col <- date_guess[1]
    # parse common formats; keep it simple/robust
    parsed <- parse_date_time(df[[date_col]],
                              orders = c("Ymd","ymd","Y-m-d","mdy","dmy","m/d/Y","d/m/Y"),
                              tz = "UTC")
    df$date <- as.Date(parsed)
    if (!identical(date_col, "date")) df[[date_col]] <- NULL
    df <- df |>
      filter(!is.na(date)) |>
      arrange(date) |>
      relocate(date, .before = 1)
  }
  
  # drop rows that are entirely NA
  df <- df |> filter(if_any(everything(), ~ !is.na(.x)))
  # de-duplicate exact duplicate rows
  df <- df |> distinct()
  df
}

# Read raw (as-is), then clean
gold_full_cleaned <- read_csv(file.path(raw_dir, "gold_spot_1990present.csv"), show_col_types = FALSE) |> basic_clean()
tips_full_cleaned <- read_csv(file.path(raw_dir, "tips_10y_real_yield_fred.csv"), show_col_types = FALSE) |> basic_clean()
usd_full_cleaned  <- read_csv(file.path(raw_dir, "usd_broad_index_fred.csv"), show_col_types = FALSE) |> basic_clean()
vix_full_cleaned  <- read_csv(file.path(raw_dir, "vix_index_fred.csv"), show_col_types = FALSE) |> basic_clean()

gold_full_cleaned <- gold_full_cleaned %>%
  arrange(date) %>%
  # guard against nonpositive prices
  mutate(across(c(open, high, low, close), ~ ifelse(.x > 0, .x, NA_real_))) %>%
  # log prices for all OHLC
  mutate(across(c(open, high, low, close), log, .names = "{.col}_log")) %>%
  # log return using close (Δ log price); first obs = NA
  mutate(close_log_return = close_log - dplyr::lag(close_log),
         close_return_pct = (close / dplyr::lag(close)) - 1)

# 2) Overwrite the RData file with updated data frames
save_path <- "/Users/cycoldiron/Desktop/gold-ts-forecast/data/clean/01_clean_economic_data.RData"

save(
  gold_full_cleaned,
  tips_full_cleaned,
  usd_full_cleaned,
  file = save_path
)