# ---- White-noise diagnostics table + figure (fixed) ------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(gt)
library(fs)

fig_dir <- "/Users/cycoldiron/Desktop/gold-ts-forecast/figures/02_diagnostics/03_white_noise"
dir_create(fig_dir, recurse = TRUE)

# lb_tbl already exists from your previous run:
# columns: Window, N, Lag, Stat, Pvalue

# ---------- (A) Course-style table (Q-stat & p-values at L = 10, 20, 40) ----
lb_table_wide <- lb_tbl %>%
  mutate(LagLab = paste0("L=", Lag)) %>%
  select(Window, N, LagLab, Stat, Pvalue) %>%
  pivot_wider(
    id_cols = c(Window, N),
    names_from = LagLab,
    values_from = c(Stat, Pvalue),
    names_glue = "{.value} ({LagLab})"
  ) %>%
  arrange(factor(Window, levels = c("1990-Pres","2010-Pres","2018-Pres")))

gt_lb <- lb_table_wide %>%
  gt(rowname_col = "Window") %>%
  fmt_number(columns = starts_with("Stat"), decimals = 1) %>%
  fmt_number(columns = starts_with("Pvalue"), decimals = 3) %>%
  tab_spanner(label = "Ljung–Box Q statistic", columns = starts_with("Stat")) %>%
  tab_spanner(label = "p-value", columns = starts_with("Pvalue")) %>%
  cols_label(N = "N obs.") %>%
  tab_header(title = md("**White-noise tests (daily log returns)**")) %>%
  tab_source_note(md("H0: no autocorrelation up to lag L. Cells with p < 0.05 reject H0."))

# Shade any p < 0.05 in light red
for (nm in names(lb_table_wide)) {
  if (grepl("^Pvalue", nm)) {
    gt_lb <- gt_lb %>%
      data_color(
        columns = all_of(nm),
        colors  = scales::col_bin(palette = c("#fde0dd", "white"),
                                  domain  = c(0,1),
                                  bins    = c(0, 0.05, 1.01))
      )
  }
}

tbl_png <- file.path(fig_dir, "lb_white_noise_table.png")
tbl_csv <- file.path(fig_dir, "lb_white_noise_table.csv")
gtsave(gt_lb, filename = tbl_png, vwidth = 1200, vheight = 360)
write.csv(lb_table_wide, tbl_csv, row.names = FALSE)

# ---------- (B) Visualization: p-values by lag (per window) ------------------
p_lb <- lb_tbl %>%
  mutate(Window = factor(Window, levels = c("1990-Pres","2010-Pres","2018-Pres"))) %>%
  ggplot(aes(x = Lag, y = Pvalue)) +
  geom_col(width = 6, alpha = 0.9) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  facet_wrap(~ Window, ncol = 1, scales = "free_x") +
  labs(
    title = "Ljung–Box p-values by lag (daily log returns)",
    subtitle = "Dashed line = 5% level; bars above the line ⇒ fail to reject white noise up to that lag.",
    x = "Lag (L)",
    y = "p-value"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig_png <- file.path(fig_dir, "lb_white_noise_pvalues.png")
ggsave(fig_png, plot = p_lb, width = 9, height = 8, dpi = 300)

message("Saved:\n  - ", tbl_png, "\n  - ", tbl_csv, "\n  - ", fig_png)
