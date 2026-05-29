# allan_deviation_analysis.R
#   Overlapping Allan deviation of flow injection plateau rates.
#   Produces Figure 6: log-log Allan deviation vs averaging time for
#   all five concentrations, with chi-squared confidence intervals.
#   Produces Figure S1: raw vs detrended sensitivity comparison.
#   Reports log-log slopes with CIs for Supplementary Table S4.
#
#   Sources: code/allan_deviation.R, code/functions/getIndexSelectedFilenameFlow.R
#   Data:    data/flow_injection_data/ (96ms, 2-channel, naloxone 40eV)
#   Output:  manuscript/figures/figure_06.pdf, manuscript/figures/figure_S1.pdf
rm(list=ls())

library(dplyr)
library(tibble)

dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

for(file in list.files(funcPath, full.names = TRUE)) {
  source(file = file, echo = FALSE)
}
source(file.path("code", "allan_deviation.R"))

# --- Parameters ---

concentrations <- c('1500ng', '1000ng', '500ng', '250ng', '50ng')
conc_numeric   <- c(1500, 1000, 500, 250, 50)
Chan      <- '2'
Dwell     <- '96ms'
Pause     <- '4ms'
Energy    <- '40eV'
Trans     <- 'any'
tau_0     <- 0.2
start_time <- 3.0
edge_trim  <- 5
n_reps     <- 3L
edf_min    <- 10L   # truncation: minimum edf per replicate

dataFiles <- list.files(dataPath)

# --- Load all 15 files (5 conc x 3 rep), store raw and detrended ---

data_detrended <- list()
data_raw       <- list()

cat("Loading 96 ms plateau data:\n")

for(ci in seq_along(concentrations)) {
  Conc <- concentrations[ci]
  conc_val <- conc_numeric[ci]

  for(rep in 1:n_reps) {
    Replicate <- as.character(rep)
    params <- list(Conc=Conc, Chan=Chan, Dwell=Dwell, Pause=Pause,
                   Replicate=Replicate, Trans=Trans, Energy=Energy)
    indSelected <- getIndexSelectedFilenameFlow(dataFiles, params)
    filename <- dataFiles[indSelected[1]]

    raw <- read.delim(paste0(dataPath, filename), header = FALSE)
    plateau_idx <- which(raw$V1 >= start_time)
    time_pts <- raw$V1[plateau_idx]
    rate <- raw$V2[plateau_idx]

    fit <- lm(rate ~ time_pts)
    resid <- fit$residuals

    n_pts <- length(resid)
    keep <- (edge_trim + 1):(n_pts - edge_trim)

    data_detrended[[length(data_detrended) + 1]] <- tibble(
      concentration = conc_val,
      replicate = rep,
      rate = resid[keep]
    )

    data_raw[[length(data_raw) + 1]] <- tibble(
      concentration = conc_val,
      replicate = rep,
      rate = rate[keep]
    )

    cat(sprintf("  %7s rep %d: %d plateau pts, %d after trim, mean = %.0f counts/s\n",
                Conc, rep, length(plateau_idx), length(keep), mean(rate)))
  }
}

plateau_detrended <- bind_rows(data_detrended)
plateau_raw       <- bind_rows(data_raw)

# --- Compute Allan deviation per concentration x replicate ---

cat("\nComputing overlapping Allan deviation...\n")

allan_detrended <- compute_allan_deviation_grouped(
  data = plateau_detrended,
  value_col = "rate",
  tau_0 = tau_0,
  group_cols = c("concentration", "replicate")
)

allan_raw <- compute_allan_deviation_grouped(
  data = plateau_raw,
  value_col = "rate",
  tau_0 = tau_0,
  group_cols = c("concentration", "replicate")
)

# --- Average Allan variances across replicates, then take sqrt ---

average_allan <- function(allan_results) {
  allan_results |>
    group_by(concentration, tau, m) |>
    summarise(
      n_samples = first(n_samples),
      allan_variance = mean(allan_variance),
      allan_deviation = sqrt(mean(allan_variance)),
      n_pairs = first(n_pairs),
      .groups = "drop"
    )
}

allan_avg <- average_allan(allan_detrended)
allan_avg_raw <- average_allan(allan_raw)

# --- Compute CIs and edf for detrended data ---

ci_result <- compute_allan_ci(
  allan_variance = allan_avg$allan_variance,
  n_samples = allan_avg$n_samples,
  m = allan_avg$m,
  n_reps = n_reps
)

allan_avg <- bind_cols(allan_avg, ci_result)

# --- Truncation based on edf threshold ---

cat(sprintf("\nTruncation criterion: edf >= %d per replicate (pooled >= %d)\n",
            edf_min, n_reps * edf_min))

tau_max_by_conc <- allan_avg |>
  filter(edf >= n_reps * edf_min) |>
  group_by(concentration) |>
  summarise(tau_max = max(tau), .groups = "drop")

cat("Maximum tau meeting edf threshold:\n")
for(i in seq_len(nrow(tau_max_by_conc))) {
  cat(sprintf("  %4d ng/mL: tau_max = %.1f s\n",
              tau_max_by_conc$concentration[i],
              tau_max_by_conc$tau_max[i]))
}

tau_cutoff <- 3.0
allan_avg_plot <- filter(allan_avg, tau <= tau_cutoff)
allan_avg_raw_plot <- filter(allan_avg_raw, tau <= tau_cutoff)

# --- Log-log slopes with CIs ---

cat("\nAllan deviation log-log slopes (white noise = -0.500):\n")
cat(sprintf("  %10s  %6s  %6s  %5s  %7s  %19s  %s\n",
            "Conc", "tau_lo", "tau_hi", "n_pts", "slope",
            "95% CI", "includes -0.5"))

slope_results <- list()
for(cv in sort(unique(allan_avg_plot$concentration), decreasing = TRUE)) {
  sub <- filter(allan_avg_plot, concentration == cv)
  log_fit <- lm(log10(allan_deviation) ~ log10(tau), data = sub)
  beta <- coef(log_fit)[2]
  ci <- confint(log_fit, "log10(tau)", level = 0.95)
  includes_half <- ci[1] <= -0.5 && ci[2] >= -0.5

  slope_results[[length(slope_results) + 1]] <- tibble(
    concentration = cv,
    tau_min = min(sub$tau),
    tau_max = max(sub$tau),
    n_points = nrow(sub),
    beta = beta,
    ci_lower = ci[1],
    ci_upper = ci[2],
    includes_minus_half = includes_half
  )

  cat(sprintf("  %4d ng/mL  %5.2f  %5.1f  %5d  %+.4f  [%+.4f, %+.4f]  %s\n",
              cv, min(sub$tau), max(sub$tau), nrow(sub), beta,
              ci[1], ci[2], ifelse(includes_half, "yes", "NO")))
}

slope_table <- bind_rows(slope_results)

# --- Figure 6: detrended Allan deviation with error bars ---

conc_levels <- sort(unique(allan_avg_plot$concentration), decreasing = TRUE)
pch_vec <- c(1, 16, 2, 0, 6)
conc_labels <- paste0(conc_levels, " ng/mL")

anchor <- filter(allan_avg_plot, concentration == conc_levels[1])
anchor_tau <- anchor$tau[1]
anchor_sigma <- anchor$allan_deviation[1]
ref_tau <- 10^seq(log10(min(allan_avg_plot$tau) * 0.8),
                  log10(max(allan_avg_plot$tau) * 1.2), length.out = 100)
ref_sigma <- anchor_sigma * (ref_tau / anchor_tau)^(-0.5)

plot_figure_6 <- function() {
  par(mar = c(5, 4, 2, 2) + 0.1)
  sub1 <- filter(allan_avg_plot, concentration == conc_levels[1])
  plot(sub1$tau, sub1$allan_deviation, log = "xy",
       pch = pch_vec[1], cex = 0.8,
       xlab = expression(paste("Averaging time ", tau, " (s)")),
       ylab = expression(paste("Allan deviation ", sigma[A](tau))),
       xlim = c(min(allan_avg_plot$tau), tau_cutoff),
       ylim = range(c(allan_avg_plot$ci_lower, allan_avg_plot$ci_upper)),
       xaxt = "n")
  axis(1, at = c(0.2, 0.5, 1.0, 2.0, 3.0),
       labels = c("0.2", "0.5", "1.0", "2.0", "3.0"))

  arrows(sub1$tau, sub1$ci_lower, sub1$tau, sub1$ci_upper,
         code = 3, angle = 90, length = 0.03, lwd = 0.7)
  lines(sub1$tau, sub1$allan_deviation, lwd = 1)

  for(i in 2:length(conc_levels)) {
    sub <- filter(allan_avg_plot, concentration == conc_levels[i])
    points(sub$tau, sub$allan_deviation, pch = pch_vec[i], cex = 0.8)
    arrows(sub$tau, sub$ci_lower, sub$tau, sub$ci_upper,
           code = 3, angle = 90, length = 0.03, lwd = 0.7)
    lines(sub$tau, sub$allan_deviation, lwd = 1)
  }

  lines(ref_tau, ref_sigma, lty = 2, lwd = 1.5)
  legend("bottomleft", legend = conc_labels, pch = pch_vec, bty = 'n',
         cex = 0.8, pt.cex = 1.2)
  par(mar = c(5, 4, 4, 2) + 0.1)
}

plot_figure_6()

pdf(file.path("manuscript", "figures", "figure_06.pdf"), width = 6.5, height = 5)
plot_figure_6()
dev.off()
cat("\nFigure 6 saved to manuscript/figures/figure_06.pdf\n")

# --- Figure S1: raw vs detrended sensitivity analysis ---

plot_figure_s1 <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1) + 0.1)

  # Panel A: raw (un-detrended)
  sub1_raw <- filter(allan_avg_raw_plot, concentration == conc_levels[1])
  plot(sub1_raw$tau, sub1_raw$allan_deviation, log = "xy",
       pch = pch_vec[1], cex = 1.2,
       xlab = expression(paste("Averaging time ", tau, " (s)")),
       ylab = expression(paste("Allan deviation ", sigma[A](tau))),
       xlim = c(min(allan_avg_raw_plot$tau), tau_cutoff),
       ylim = range(allan_avg_raw_plot$allan_deviation),
       xaxt = "n", main = "(A) Raw plateau rates")
  axis(1, at = c(0.2, 0.5, 1.0, 2.0, 3.0),
       labels = c("0.2", "0.5", "1.0", "2.0", "3.0"))
  lines(sub1_raw$tau, sub1_raw$allan_deviation, lwd = 1)
  for(i in 2:length(conc_levels)) {
    sub <- filter(allan_avg_raw_plot, concentration == conc_levels[i])
    points(sub$tau, sub$allan_deviation, pch = pch_vec[i], cex = 1.2)
    lines(sub$tau, sub$allan_deviation, lwd = 1)
  }
  lines(ref_tau, ref_sigma, lty = 2, lwd = 1.5)

  # Panel B: detrended
  sub1_det <- filter(allan_avg_plot, concentration == conc_levels[1])
  plot(sub1_det$tau, sub1_det$allan_deviation, log = "xy",
       pch = pch_vec[1], cex = 1.2,
       xlab = expression(paste("Averaging time ", tau, " (s)")),
       ylab = expression(paste("Allan deviation ", sigma[A](tau))),
       xlim = c(min(allan_avg_plot$tau), tau_cutoff),
       ylim = range(allan_avg_plot$allan_deviation),
       xaxt = "n", main = "(B) Detrended residuals")
  axis(1, at = c(0.2, 0.5, 1.0, 2.0, 3.0),
       labels = c("0.2", "0.5", "1.0", "2.0", "3.0"))
  lines(sub1_det$tau, sub1_det$allan_deviation, lwd = 1)
  for(i in 2:length(conc_levels)) {
    sub <- filter(allan_avg_plot, concentration == conc_levels[i])
    points(sub$tau, sub$allan_deviation, pch = pch_vec[i], cex = 1.2)
    lines(sub$tau, sub$allan_deviation, lwd = 1)
  }
  lines(ref_tau, ref_sigma, lty = 2, lwd = 1.5)
  legend("bottomleft", legend = conc_labels, pch = pch_vec, bty = 'n',
         cex = 0.6, pt.cex = 0.8)

  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

plot_figure_s1()

pdf(file.path("manuscript", "figures", "figure_S1.pdf"), width = 10, height = 5)
plot_figure_s1()
dev.off()
cat("Figure S1 saved to manuscript/figures/figure_S1.pdf\n")
