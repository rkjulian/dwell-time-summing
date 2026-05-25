# Timing Figure.R
#   Generates the temporal ion distribution in Figure 1.
#   User functions called: none
# Author: Fred Lytle, Randall Julian
# Last edit: 5/19/2026

rm(list=ls())

set.seed(1)
dt <- rexp(50,0.5)
t <- cumsum(dt)

# -- base R version (partial figure, no annotations) ---------------------

# set plot window to w = 5", h = 2.5"
par(yaxt='n')
# par(xaxt='n')
plot(t,rep(2.4,length(t)),type='h',xlim=c(0,100),ylim=c(0,3),
     ylab='',xlab='',xaxp=c(0,100,4),lwd=1.5)
lines(x=c(0,25,50,75,100),y=rep(3.5,5),type='h',lty=2,lwd=1)
abline(h=0,lwd=1)

# -- ggplot version (complete figure with all annotations) ---------------

library(ggplot2)

td <- 21
tp <- 4
tc <- 100
n_seg <- as.integer(tc / (td + tp))
seg_starts <- (0:(n_seg - 1)) * (td + tp)

t_events <- t[t <= tc]
ion_df <- data.frame(x = t_events)

count_labels <- sapply(seg_starts, function(s) {
  in_seg   <- sum(t_events >= s & t_events < s + td + tp)
  in_dwell <- sum(t_events >= s & t_events < s + td)
  paste0(in_seg, "/", in_dwell)
})

pause_df <- data.frame(
  xmin = seg_starts + td,
  xmax = seg_starts + td + tp,
  ymin = 0, ymax = 1
)

bounds <- seq(0, tc, by = td + tp)

box_pad <- 3
y_arrow <- 1.07
y_label <- 1.19

tp_left  <- seg_starts[2] + td
tp_right <- seg_starts[3]
tp_mid   <- (tp_left + tp_right) / 2

fig1 <- ggplot() +
  annotate("rect", xmin = -box_pad, xmax = tc + box_pad,
           ymin = 0, ymax = 1,
           fill = NA, color = "black", linewidth = 0.5) +
  geom_rect(data = pause_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey80") +
  geom_segment(data = ion_df,
               aes(x = x, xend = x, y = 0, yend = 0.88),
               linewidth = 0.5) +
  annotate("segment", x = bounds, xend = bounds,
           y = 0, yend = 1,
           linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = seg_starts + td / 2, y = 0.95,
           label = count_labels, size = 3.8) +
  # t_d: label above, double-headed arrow below
  annotate("text", x = td / 2, y = y_label,
           label = expression(italic(t)[italic(d)] ~ "=" ~ 21),
           size = 4) +
  annotate("segment", x = 0, xend = td, y = y_arrow, yend = y_arrow,
           arrow = arrow(ends = "both", length = unit(0.07, "inches"),
                         type = "closed"),
           linewidth = 0.4) +
  # t_p: label above, two inward-pointing arrows below
  annotate("text", x = tp_mid, y = y_label,
           label = expression(italic(t)[italic(p)] ~ "=" ~ 4),
           size = 4) +
  annotate("segment", x = tp_left - 8, xend = tp_left,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(ends = "last", length = unit(0.07, "inches"),
                         type = "closed"),
           linewidth = 0.4) +
  annotate("segment", x = tp_right + 8, xend = tp_right,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(ends = "last", length = unit(0.07, "inches"),
                         type = "closed"),
           linewidth = 0.4) +
  # x-axis
  annotate("segment", x = seq(25, 75, 25), xend = seq(25, 75, 25),
           y = 0, yend = -0.04, linewidth = 0.4) +
  annotate("text", x = seq(0, tc, 25), y = -0.10,
           label = seq(0, tc, 25), size = 3.5) +
  # t_c: arrow above, label below
  annotate("segment", x = 0, xend = tc, y = -0.22, yend = -0.22,
           arrow = arrow(ends = "both", length = unit(0.07, "inches"),
                         type = "closed"),
           linewidth = 0.4) +
  annotate("text", x = tc / 2, y = -0.32,
           label = expression(italic(t)[italic(c)] ~ "=" ~ 100),
           size = 4) +
  coord_cartesian(xlim = c(-7, 107), ylim = c(-0.38, 1.26), clip = "off") +
  theme_void()

print(fig1)

# publication PDF
if (FALSE) {
  fig_dir <- file.path("manuscript", "figures")
  pdf(file.path(fig_dir, "figure_01.pdf"), width = 5, height = 2.5)
  print(fig1)
  dev.off()
}
