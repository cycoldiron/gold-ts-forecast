# scripts/01_get_data.R
# Purpose: Create folders, pull daily series from FRED, and save raw CSVs + a short data note.

# ---- Setup ----
needed <- c("quantmod", "tidyverse", "lubridate", "glue")
to_install <- needed[!needed %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
library(quantmod); library(tidyverse); library(lubridate); library(glue)

# Repo folders
dirs <- c("data/raw", "data/clean", "scripts", "figures", "results", "paper")
dir.create("data", showWarnings = FALSE)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Series map (FRED IDs) ----
fred_series <- c(
  gold_pm = "GOLDPMGBD228NLBM",   # LBMA Gold Price PM (USD/oz)
  gold_am = "GOLDAMGBD228NLBM",   # LBMA Gold Price AM (USD/oz)  [aux/robustness]
  vix     = "VIXCLS",             # CBOE VIX close
  usd_brd = "DTWEXBGS",           # Fed Nominal Broad Dollar Index (daily)
  tips10  = "DFII10"              # 10y TIPS real yield (percent)
)

# ---- Download helper ----
get_fred <- function(sym) {
  getSymbols(sym, src = "FRED", auto.assign = FALSE)
}

# Download all, quietly
raw_xts <- lapply(fred_series, get_fred)

# Convert to tibble with 'date' + column name
to_tibble <- function(x, name) {
  tibble(date = as_date(index(x)),
         !!name := as.numeric(x[,1]))
}
raw_tbls <- Map(to_tibble, raw_xts, names(fred_series))

# Save each series as its own CSV in data/raw
iwalk(raw_tbls, function(tb, nm) {
  out <- file.path("data/raw", paste0(nm, "_fred_", fred_series[[nm]], ".csv"))
  write_csv(tb, out)
})

# Also save a single joined “wide” raw file (by inner join to keep common dates only)
raw_wide <- reduce(raw_tbls, full_join, by = "date") |> arrange(date)
write_csv(raw_wide, "data/raw/raw_all_wide.csv")

# ---- Data note (markdown) ----
note <- glue(
  "# Data Note (raw)\n\n",
  "Pulled on: {Sys.Date()}\n\n",
  "**Gold (PM)**: FRED mirror of LBMA/ICE, ID: {fred_series['gold_pm']}, units: USD/oz, freq: daily (London ~15:00).  \n",
  "**Gold (AM)**: FRED mirror of LBMA/ICE, ID: {fred_series['gold_am']}, units: USD/oz, freq: daily (London ~10:30).  \n",
  "**VIX**: CBOE VIX close, ID: {fred_series['vix']}, units: index, freq: daily.  \n",
  "**USD Broad**: Fed Nominal Broad U.S. Dollar Index, ID: {fred_series['usd_brd']}, units: index (2006=100), freq: daily.  \n",
  "**10y TIPS**: Real yield, ID: {fred_series['tips10']}, units: percent, freq: daily.  \n\n",
  "Notes: PM is the primary gold series; AM used only for robustness. Use USD Broad (DTWEXBGS) for macro co-movement; DXY optional.\n"
)
writeLines(note, "data/raw/_data_note.md")

message("Done: raw series saved to data/raw/ and note written.")
