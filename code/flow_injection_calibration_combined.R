# Flow Injection Calibration Combined.R
#   Generates merged Figure 4: two-panel calibration (rate + count).
#   User functions called: none
rm(list = ls())

library(ggplot2)
library(patchwork)

# data from "Flow Injection Plateau Rates.docx" and "Flow Injection Plateau Counts.docx"
conc <- c(50, 250, 500, 1000, 1500)

rate96 <- c(6503, 36351, 60374, 116243, 177676)
rate21 <- c(5343, 34423, 58888, 114807, 171630)
rate06 <- c(4461, 32907, 57233, 111963, 167024)

count96 <- c(624.3, 3490, 5796, 11159, 17057)
count21 <- c(448.8, 2892, 4947, 9644, 14417)
count06 <- c(267.6, 1974, 3434, 6718, 10021)

dwell_labels <- c("1x96 ms", "4x21 ms", "10x6 ms")

rate_df <- data.frame(
  conc  = rep(conc, 3),
  y     = c(rate96, rate21, rate06),
  dwell = factor(rep(dwell_labels, each = length(conc)), levels = dwell_labels)
)

count_df <- data.frame(
  conc  = rep(conc, 3),
  y     = c(count96, count21, count06),
  dwell = factor(rep(dwell_labels, each = length(conc)), levels = dwell_labels)
)

shared_theme <- theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.tag          = element_text(face = "bold", size = 10)
  )

p_rate <- ggplot(rate_df, aes(conc, y, shape = dwell)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "black") +
  scale_shape_manual(values = c(1, 16, 2)) +
  scale_x_continuous(limits = c(0, 1500)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Conc (ng/mL)", y = "Dwell Time Rate Average") +
  shared_theme

p_count <- ggplot(count_df, aes(conc, y, shape = dwell)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "black") +
  scale_shape_manual(values = c(1, 16, 2)) +
  scale_x_continuous(limits = c(0, 1500)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Conc (ng/mL)", y = "Dwell Time Count Sum") +
  shared_theme

fig4 <- p_rate + p_count +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title    = element_blank())

print(fig4)

# publication PDF
if (FALSE) {
  fig_dir <- file.path("manuscript", "figures")
  pdf(file.path(fig_dir, "figure_04.pdf"), width = 7, height = 3.5)
  print(fig4)
  dev.off()
}
