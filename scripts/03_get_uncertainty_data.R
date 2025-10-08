library(readxl)
library(dplyr)
library(lubridate)
library(zoo)



epu_raw <- read_csv("https://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.csv",
                show_col_types = FALSE) 
# 2) GPR (Global)
gpr_raw <- read_xls("/Users/cycoldiron/Desktop/gold-ts-forecast/data/raw/gpr_index_raw.xls")

# --- Save uncertainty tibbles to a single .Rdata file -------------------------
# Ensure target folder exists
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# Path for the save file (match your requested name & case)
save_path <- file.path("data", "raw", "02_raw_uncertainty.Rdata")

# Save all three objects in one R data file
# (assumes epu_raw, gpr_raw, and vix_raw exist in the environment)
save(epu_raw, gpr_raw, vix_raw, file = save_path, compress = "bzip2")

# Optional: message to confirm
message("Saved: ", normalizePath(save_path))
