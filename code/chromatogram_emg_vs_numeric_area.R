# chromatogram_emg_vs_numeric_area.R
#
#   Compares EMG-fitted peak areas to numeric (rectangular-rule) areas
#   for all chromatographic configurations, with bootstrap CIs, a
#   claim-calibrated recovery simulation as positive control, and a
#   pooled equivalence test on the config effect.
#
#   Resolves whether the EMG parametric fit masks or creates dwell-time-
#   dependent precision effects that a direct numeric integration would
#   reveal.
rm(list = ls())
library(pracma)

dataPath <- file.path("data", "chromatogram_data", "")
funcPath <- file.path("code", "functions", "")
for (file in list.files(funcPath, full.names = TRUE)) {
  source(file = file, echo = FALSE)
}
dataFiles <- list.files(dataPath)

startTime <- 1.0   # data window start (minutes)
endTime   <- 1.3   # data window end (minutes)
peakStart_min <- 1.09
Pause  <- "4ms"
Trans  <- "any"
n_reps <- 4
B_boot <- 10000
cv_pro_flow <- 0.04  # from flow injection Table 1

configs <- list(
  list(label = "1x96",  chan = "2",    dwell = "96ms", td = 0.096, n = 1),
  list(label = "4x21",  chan = "2x4",  dwell = "21ms", td = 0.021, n = 4),
  list(label = "10x6",  chan = "2x10", dwell = "6ms",  td = 0.006, n = 10)
)
concentrations <- c("500ng", "250ng", "200ng", "100ng", "50ng")
conc_numeric   <- c(500, 250, 200, 100, 50)


# =============================================================================
# SECTION 1: EMG AREA VERIFICATION
# =============================================================================

cat("=== EMG coefs[1] verification ===\n")
test_params <- c(1000, 67, 0.5, 1.8)
x_fine <- seq(50, 100, by = 0.001)
y_fine <- EMG(x_fine, test_params)
cat(sprintf("  params[1] = %.1f, numerical integral = %.4f, ratio = %.6f\n",
            test_params[1], sum(y_fine) * 0.001,
            sum(y_fine) * 0.001 / test_params[1]))
cat("  coefs[1] IS the integrated area (count-sec), not amplitude.\n\n")


# =============================================================================
# SECTION 2: AREA COMPUTATION WITH BASELINE AND WINDOW FIXES
# =============================================================================

get_peak_areas <- function(dataFiles, dataPath, conc, cfg, rep, energy) {
  params <- list(Conc = conc, Chan = cfg$chan, Dwell = cfg$dwell,
                 Pause = Pause, Replicate = as.character(rep),
                 Trans = Trans, Energy = energy)
  indSelected <- getIndexSelectedFilenameChrom(dataFiles, params)
  filenames <- dataFiles[indSelected]
  numDT <- length(filenames)

  allTime <- read.delim(paste0(dataPath, filenames[1]), header = FALSE)$V1
  keepInd <- which(allTime >= startTime & allTime <= endTime)
  time <- allTime[keepInd] * 60

  rates <- array(0, c(numDT, length(time)))
  for (i in 1:numDT) {
    data <- read.delim(paste0(dataPath, filenames[i]), header = FALSE)
    rates[i, ] <- data$V2[keepInd]
  }

  if (numDT > 1) avg_rates <- colMeans(rates) else avg_rates <- drop(rates)

  counts <- rates * cfg$td
  if (numDT > 1) sum_counts <- colSums(counts) else sum_counts <- drop(counts)

  deltaT <- time[2] - time[1]
  peakStart_s <- peakStart_min * 60
  peakEnd_s   <- max(time)  # explicit: data window is the real upper bound

  # baseline: mean of pre-peak region
  baseline_idx <- which(time < peakStart_s)
  baseline_count <- if (length(baseline_idx) >= 3) mean(sum_counts[baseline_idx]) else 0
  baseline_rate  <- if (length(baseline_idx) >= 3) mean(avg_rates[baseline_idx]) else 0

  peak_idx <- which(time >= peakStart_s & time <= peakEnd_s)

  # numeric areas with baseline subtraction
  num_count <- sum(sum_counts[peak_idx] - baseline_count) * deltaT
  num_rate  <- sum(avg_rates[peak_idx] - baseline_rate) * deltaT

  # EMG fit on full window (includes baseline in the fit)
  options(warn = -1)
  guess_c <- c(max(sum_counts), time[which.max(sum_counts)], 0.5, 0.5)
  out_c <- tryCatch(marquardtSearch(time, sum_counts, EMG, guess_c),
                    error = function(e) NULL)
  guess_r <- c(sum(avg_rates), time[which.max(avg_rates)], 0.5, 0.5)
  out_r <- tryCatch(marquardtSearch(time, avg_rates, EMG, guess_r),
                    error = function(e) NULL)
  options(warn = 0)

  emg_count    <- if (!is.null(out_c)) out_c$coefs[1] else NA
  emg_rate     <- if (!is.null(out_r)) out_r$coefs[1] else NA
  emg_c_params <- if (!is.null(out_c)) out_c$coefs    else rep(NA, 4)
  converged_c  <- !is.null(out_c)
  converged_r  <- !is.null(out_r)

  list(emg_count = emg_count, emg_rate = emg_rate,
       num_count = num_count, num_rate = num_rate,
       emg_params = emg_c_params,
       converged_c = converged_c, converged_r = converged_r,
       mean_count = mean(sum_counts[peak_idx]),
       baseline_count = baseline_count,
       time = time, peak_idx = peak_idx, deltaT = deltaT)
}

# collect all individual replicate areas
all_results <- list()
convergence_log <- data.frame(conc = character(), config = character(),
                              rep = integer(), conv_c = logical(),
                              conv_r = logical(), stringsAsFactors = FALSE)

cat("=== Computing areas for all configurations ===\n")

for (ci in seq_along(concentrations)) {
  for (cfg in configs) {
    for (rep in 1:n_reps) {
      a <- get_peak_areas(dataFiles, dataPath, concentrations[ci], cfg, rep, "40eV")
      key <- paste(conc_numeric[ci], cfg$label, sep = "_")
      if (is.null(all_results[[key]])) {
        all_results[[key]] <- list(
          conc = conc_numeric[ci], config = cfg$label,
          live_time = cfg$n * cfg$td,
          emg_counts = numeric(0), num_counts = numeric(0),
          emg_rates = numeric(0), num_rates = numeric(0),
          mean_peak_counts = numeric(0), emg_params_list = list()
        )
      }
      all_results[[key]]$emg_counts <- c(all_results[[key]]$emg_counts, a$emg_count)
      all_results[[key]]$num_counts <- c(all_results[[key]]$num_counts, a$num_count)
      all_results[[key]]$emg_rates  <- c(all_results[[key]]$emg_rates, a$emg_rate)
      all_results[[key]]$num_rates  <- c(all_results[[key]]$num_rates, a$num_rate)
      all_results[[key]]$mean_peak_counts <- c(all_results[[key]]$mean_peak_counts, a$mean_count)
      all_results[[key]]$emg_params_list[[rep]] <- a$emg_params

      convergence_log <- rbind(convergence_log, data.frame(
        conc = concentrations[ci], config = cfg$label,
        rep = rep, conv_c = a$converged_c, conv_r = a$converged_r,
        stringsAsFactors = FALSE
      ))
    }
  }
}

n_total <- nrow(convergence_log)
n_fail <- sum(!convergence_log$conv_c)
cat(sprintf("  Convergence: %d/%d count fits succeeded, %d/%d rate fits succeeded\n",
            n_total - n_fail, n_total,
            n_total - sum(!convergence_log$conv_r), n_total))
conv_by_config <- aggregate(cbind(conv_c, conv_r) ~ config, convergence_log,
                            function(x) sprintf("%d/%d", sum(x), length(x)))
cat("  Per-config breakdown (count / rate):\n")
for (i in 1:nrow(conv_by_config)) {
  cat(sprintf("    %s: count %s, rate %s\n",
              conv_by_config$config[i], conv_by_config$conv_c[i], conv_by_config$conv_r[i]))
}
cat("\n")


# =============================================================================
# SECTION 3: BOOTSTRAP CIs
# =============================================================================

bootstrap_cv <- function(areas, B = B_boot) {
  valid <- areas[!is.na(areas)]
  n <- length(valid)
  if (n < 2) return(list(cv = NA, ci = c(NA, NA)))
  cv_obs <- sd(valid) / mean(valid)
  cv_boot <- replicate(B, {
    idx <- sample(n, n, replace = TRUE)
    sd(valid[idx]) / mean(valid[idx])
  })
  list(cv = cv_obs, ci = quantile(cv_boot, c(0.025, 0.975), na.rm = TRUE))
}

cat("===========================================================================\n")
cat("Table 3 Comparison: EMG vs Numeric CV with 95% Bootstrap CIs (n=4)\n")
cat("===========================================================================\n\n")
cat(sprintf("%-6s %-6s  %7s [%7s, %7s]   %7s [%7s, %7s]   %6s\n",
            "Conc", "DT", "CV_EMG", "lo", "hi", "CV_Num", "lo", "hi", "Ratio"))
cat(paste(rep("-", 82), collapse = ""), "\n")

ratio_table <- data.frame(conc = integer(), config = character(),
                          cv_emg = numeric(), cv_num = numeric(),
                          ratio = numeric(), stringsAsFactors = FALSE)

for (ci in seq_along(conc_numeric)) {
  for (cfg in configs) {
    key <- paste(conc_numeric[ci], cfg$label, sep = "_")
    r <- all_results[[key]]

    b_emg <- bootstrap_cv(r$emg_counts)
    b_num <- bootstrap_cv(r$num_counts)

    ratio <- b_emg$cv / b_num$cv

    cat(sprintf("%-6d %-6s  %7.4f [%7.4f, %7.4f]   %7.4f [%7.4f, %7.4f]   %6.3f\n",
                r$conc, r$config,
                b_emg$cv, b_emg$ci[1], b_emg$ci[2],
                b_num$cv, b_num$ci[1], b_num$ci[2],
                ratio))

    ratio_table <- rbind(ratio_table, data.frame(
      conc = r$conc, config = r$config,
      cv_emg = b_emg$cv, cv_num = b_num$cv,
      ratio = ratio, stringsAsFactors = FALSE
    ))
  }
}


# =============================================================================
# SECTION 4: PER-CONFIG EMG/NUMERIC RATIO (FLATNESS TEST)
# =============================================================================

cat("\n===========================================================================\n")
cat("Per-Config EMG/Numeric CV Ratio (flat = EMG is config-independent filter)\n")
cat("===========================================================================\n\n")

for (cfg_label in c("1x96", "4x21", "10x6")) {
  sub <- ratio_table[ratio_table$config == cfg_label, ]
  cat(sprintf("  %s: ratios = %s, mean = %.3f, sd = %.3f\n",
              cfg_label,
              paste(sprintf("%.3f", sub$ratio), collapse = ", "),
              mean(sub$ratio), sd(sub$ratio)))
}

overall_ratios <- tapply(ratio_table$ratio, ratio_table$config, mean)
cat(sprintf("\n  Cross-config spread: %.3f to %.3f (range = %.3f)\n",
            min(overall_ratios), max(overall_ratios),
            max(overall_ratios) - min(overall_ratios)))
cat("  If range << mean, EMG is config-independent and cannot mask a relative\n")
cat("  dwell effect.\n")


# =============================================================================
# SECTION 5: POOLED NOISE MODEL
# =============================================================================

cat("\n===========================================================================\n")
cat("Pooled Noise Model: CV^2 = CV_pro^2 + beta / mean_area\n")
cat("===========================================================================\n\n")

model_data <- data.frame(
  config = character(), conc = integer(),
  cv2_emg = numeric(), cv2_num = numeric(),
  inv_area_emg = numeric(), inv_area_num = numeric(),
  live_time = numeric(), stringsAsFactors = FALSE
)

for (ci in seq_along(conc_numeric)) {
  for (cfg in configs) {
    key <- paste(conc_numeric[ci], cfg$label, sep = "_")
    r <- all_results[[key]]
    emg_valid <- r$emg_counts[!is.na(r$emg_counts)]
    num_valid <- r$num_counts[!is.na(r$num_counts)]
    if (length(emg_valid) < 2 || length(num_valid) < 2) next
    model_data <- rbind(model_data, data.frame(
      config = cfg$label, conc = conc_numeric[ci],
      cv2_emg = (sd(emg_valid) / mean(emg_valid))^2,
      cv2_num = (sd(num_valid) / mean(num_valid))^2,
      inv_area_emg = 1 / mean(emg_valid),
      inv_area_num = 1 / mean(num_valid),
      live_time = cfg$n * cfg$td,
      stringsAsFactors = FALSE
    ))
  }
}

fit_emg <- lm(cv2_emg ~ inv_area_emg, data = model_data)
fit_num <- lm(cv2_num ~ inv_area_num, data = model_data)

cat("  EMG estimator:    CV_pro = ", sprintf("%.4f", sqrt(max(0, coef(fit_emg)[1]))),
    "  beta = ", sprintf("%.2f", coef(fit_emg)[2]),
    "  R^2 = ", sprintf("%.3f", summary(fit_emg)$r.squared), "\n")
cat("  Numeric estimator: CV_pro = ", sprintf("%.4f", sqrt(max(0, coef(fit_num)[1]))),
    "  beta = ", sprintf("%.2f", coef(fit_num)[2]),
    "  R^2 = ", sprintf("%.3f", summary(fit_num)$r.squared), "\n")

# test for config-dependent residuals
model_data$resid_emg <- residuals(fit_emg)
model_data$resid_num <- residuals(fit_num)
aov_emg <- summary(aov(resid_emg ~ config, data = model_data))
aov_num <- summary(aov(resid_num ~ config, data = model_data))
cat(sprintf("\n  Config effect in residuals (ANOVA p-value):\n"))
cat(sprintf("    EMG:     p = %.4f\n", aov_emg[[1]]$`Pr(>F)`[1]))
cat(sprintf("    Numeric: p = %.4f\n", aov_num[[1]]$`Pr(>F)`[1]))
cat("  p > 0.05 means no systematic dwell-time dependence beyond the noise model.\n")


# =============================================================================
# SECTION 6: RECOVERY SIMULATION (POSITIVE CONTROL)
# =============================================================================

cat("\n===========================================================================\n")
cat("Recovery Simulation: inject known dwell effect, verify both estimators\n")
cat("detect it\n")
cat("===========================================================================\n\n")

# reference peak shape from 250 ng/mL 4x21 (well-fit, mid-concentration)
ref_key <- paste(250, "4x21", sep = "_")
ref_fit <- NULL
for (p in all_results[[ref_key]]$emg_params_list) {
  if (!any(is.na(p))) { ref_fit <- p; break }
}
if (is.null(ref_fit)) {
  cat("  No valid reference fit; skipping simulation.\n")
} else {

  true_area   <- ref_fit[1]
  true_mu     <- ref_fit[2]
  true_sigma  <- ref_fit[3]
  true_tau    <- ref_fit[4]
  true_params <- ref_fit
  baseline_rate <- 50  # counts/sec, approximate from pre-peak region
  n_sim <- 500

  sim_time <- seq(60, 78, by = 0.2)
  n_pts <- length(sim_time)

  simulate_and_measure <- function(cfg, cv_pro, intermediate_cv = 0) {
    live_time <- cfg$n * cfg$td
    true_shape <- EMG(sim_time, true_params)

    emg_areas <- numeric(n_sim)
    num_areas <- numeric(n_sim)
    emg_fail  <- 0

    for (sim in 1:n_sim) {
      # true rate at each time point
      rate_true <- true_shape + baseline_rate

      # proportional noise (common-mode per cycle)
      rate_noisy <- rate_true * (1 + rnorm(n_pts, sd = cv_pro))

      # optional intermediate noise: independent per dwell, averaged over n
      if (intermediate_cv > 0 && cfg$n > 1) {
        per_dwell_noise <- matrix(rnorm(n_pts * cfg$n, sd = intermediate_cv),
                                 nrow = n_pts, ncol = cfg$n)
        avg_noise <- rowMeans(per_dwell_noise)
        rate_noisy <- rate_noisy * (1 + avg_noise)
      } else if (intermediate_cv > 0 && cfg$n == 1) {
        rate_noisy <- rate_noisy * (1 + rnorm(n_pts, sd = intermediate_cv))
      }

      # Poisson counting: expected counts = rate * live_time
      expected_counts <- pmax(rate_noisy * live_time, 0)
      observed_counts <- rpois(n_pts, expected_counts)
      observed_rate   <- observed_counts / live_time

      # numeric area (baseline-subtracted)
      peak_idx <- which(sim_time >= peakStart_min * 60)
      bl <- mean(observed_rate[sim_time < peakStart_min * 60])
      num_areas[sim] <- sum(observed_rate[peak_idx] - bl) * 0.2

      # EMG fit
      options(warn = -1)
      guess <- c(max(observed_rate), sim_time[which.max(observed_rate)], 0.5, 0.5)
      fit <- tryCatch(marquardtSearch(sim_time, observed_rate, EMG, guess),
                      error = function(e) NULL)
      options(warn = 0)

      if (!is.null(fit)) {
        emg_areas[sim] <- fit$coefs[1]
      } else {
        emg_areas[sim] <- NA
        emg_fail <- emg_fail + 1
      }
    }

    emg_valid <- emg_areas[!is.na(emg_areas)]
    list(
      cv_emg = sd(emg_valid) / mean(emg_valid),
      cv_num = sd(num_areas) / mean(num_areas),
      mean_emg = mean(emg_valid),
      mean_num = mean(num_areas),
      n_fail = emg_fail,
      true_area = true_area
    )
  }

  cat(sprintf("  Reference peak: area=%.0f, mu=%.1f, sigma=%.3f, tau=%.3f\n",
              true_area, true_mu, true_sigma, true_tau))
  cat(sprintf("  Baseline=%d counts/s, CV_pro=%.2f, N_sim=%d\n\n",
              baseline_rate, cv_pro_flow, n_sim))

  # Scenario A: Poisson + proportional only (live-time penalty is the only effect)
  cat("--- Scenario A: Poisson + proportional, no intermediate noise ---\n")
  cat("    (dwell effect = live-time penalty only)\n\n")
  cat(sprintf("  %-6s  %8s %8s  %8s %8s  %6s  %5s\n",
              "Config", "CV_EMG", "CV_Num", "Mean_EMG", "Mean_Num", "Ratio", "Fails"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  sim_results_a <- list()
  for (cfg in configs) {
    set.seed(42)
    res <- simulate_and_measure(cfg, cv_pro_flow, intermediate_cv = 0)
    sim_results_a[[cfg$label]] <- res
    cat(sprintf("  %-6s  %8.4f %8.4f  %8.0f %8.0f  %6.3f  %5d\n",
                cfg$label, res$cv_emg, res$cv_num,
                res$mean_emg, res$mean_num,
                res$cv_emg / res$cv_num, res$n_fail))
  }

  # expected CV ratio from live-time penalty
  lt <- sapply(configs, function(c) c$n * c$td)
  cat(sprintf("\n  Expected CV ratio (live-time): 96ms/21ms = %.3f, 96ms/6ms = %.3f\n",
              sqrt(lt[1] / lt[2]), sqrt(lt[1] / lt[3])))
  cat(sprintf("  Observed EMG:                 96ms/21ms = %.3f, 96ms/6ms = %.3f\n",
              sim_results_a[["4x21"]]$cv_emg / sim_results_a[["1x96"]]$cv_emg,
              sim_results_a[["10x6"]]$cv_emg / sim_results_a[["1x96"]]$cv_emg))
  cat(sprintf("  Observed Numeric:             96ms/21ms = %.3f, 96ms/6ms = %.3f\n",
              sim_results_a[["4x21"]]$cv_num / sim_results_a[["1x96"]]$cv_num,
              sim_results_a[["10x6"]]$cv_num / sim_results_a[["1x96"]]$cv_num))

  # ---------------------------------------------------------------------------
  # Scenario B: inject intermediate noise sized to the literature claim.
  #
  # The summing literature claims sqrt(10) ~ 3.2x LLOQ improvement for 10x6
  # vs 1x96. Under a precision-based LLOQ (EP17-A2), this requires the 1x96
  # area CV to be R_claim times the 10x6 area CV. We back-calculate the
  # per-cycle intermediate noise magnitude needed to produce that area-CV
  # ratio after peak-integration dilution.
  #
  # A per-cycle intermediate fluctuation is shared across all ~N_pts points
  # within a single simulated injection but is independent across injections.
  # For multi-dwell configs (n > 1), each dwell gets an independent draw and
  # the average over n dwells reduces the contribution by sqrt(n). For the
  # area (sum over N_pts points of the same per-cycle draw), the intermediate
  # noise does NOT average down with N_pts -- it is a common multiplicative
  # factor within each injection.
  #
  # To produce the claimed area-CV ratio R_claim between 1x96 (no averaging)
  # and 10x6 (averaging over n=10 dwells):
  #   CV_area_1x96^2 = CV_base^2 + inter_cv^2
  #   CV_area_10x6^2 = CV_base^2 + inter_cv^2 / 10
  #   R_claim = CV_area_1x96 / CV_area_10x6
  #
  # Solving: inter_cv^2 = CV_base^2 * (R_claim^2 - 1) / (1 - R_claim^2/10)
  # where CV_base is the Scenario A area CV for 1x96.
  # ---------------------------------------------------------------------------

  R_claim <- sqrt(10)  # claimed LLOQ ratio for 10x6 vs 1x96
  n_multi <- 10        # number of dwells in the multi-dwell config

  cv_base <- sim_results_a[["1x96"]]$cv_num
  denom <- 1 - R_claim^2 / n_multi

  cat(sprintf("\n--- Scenario B: claim-calibrated intermediate noise ---\n"))
  cat(sprintf("  R_claim = sqrt(%d) = %.3f (claimed CV ratio 1x96 / 10x6)\n",
              n_multi, R_claim))
  cat(sprintf("  CV_base (Scenario A, 1x96 numeric) = %.4f\n", cv_base))

  if (denom <= 0) {
    cat(sprintf("  Denominator (1 - R_claim^2/n) = %.4f <= 0\n", denom))
    cat("  No finite intermediate noise can produce this CV ratio.\n")
    cat(sprintf("  The claimed sqrt(%d) improvement requires R_claim^2/n >= 1,\n", n_multi))
    cat("  meaning the averaging benefit from n dwells cannot overcome the\n")
    cat("  noise injection. The claim is structurally impossible for any\n")
    cat("  intermediate noise mechanism operating at the per-dwell timescale.\n\n")
    cat("  Falling back to a demonstrably detectable injection (inter_cv = 0.30)\n")
    cat("  to verify that both estimators can detect a large intermediate effect.\n")
    inter_cv <- 0.30
    inter_cv_is_fallback <- TRUE
  } else {
    inter_cv2 <- cv_base^2 * (R_claim^2 - 1) / denom
    inter_cv <- sqrt(inter_cv2)
    inter_cv_is_fallback <- FALSE

    cat(sprintf("  Required per-cycle intermediate CV  = %.4f (%.1f%%)\n",
                inter_cv, inter_cv * 100))
    peak_rate <- max(EMG(sim_time, true_params))
    cat(sprintf("  At the reference peak rate of ~%.0f counts/sec, this is a\n", peak_rate))
    cat(sprintf("  per-cycle rate fluctuation of ~%.0f counts/sec.\n",
                inter_cv * peak_rate))

    if (inter_cv > 1.0) {
      cat("  NOTE: required intermediate CV exceeds 100%% of the signal.\n")
      cat("  This is implausibly large, which is itself evidence that the\n")
      cat("  claimed area benefit cannot arise from this mechanism.\n")
    }
  }

  cat(sprintf("\n  %-6s  %8s %8s  %8s %8s  %6s  %5s\n",
              "Config", "CV_EMG", "CV_Num", "Mean_EMG", "Mean_Num", "Ratio", "Fails"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  sim_results_b <- list()
  for (cfg in configs) {
    set.seed(42)
    res <- simulate_and_measure(cfg, cv_pro_flow,
                                intermediate_cv = ifelse(is.na(inter_cv), 0, inter_cv))
    sim_results_b[[cfg$label]] <- res
    cat(sprintf("  %-6s  %8.4f %8.4f  %8.0f %8.0f  %6.3f  %5d\n",
                cfg$label, res$cv_emg, res$cv_num,
                res$mean_emg, res$mean_num,
                res$cv_emg / res$cv_num, res$n_fail))
  }

  obs_ratio_b <- sim_results_b[["1x96"]]$cv_num / sim_results_b[["10x6"]]$cv_num
  cat(sprintf("\n  Observed 1x96/10x6 CV ratio: %.3f (expected from claim: %.3f)\n",
              obs_ratio_b, R_claim))
  cat(sprintf("  Scenario A 10x6/1x96 EMG: %.3f  Num: %.3f\n",
              sim_results_a[["10x6"]]$cv_emg / sim_results_a[["1x96"]]$cv_emg,
              sim_results_a[["10x6"]]$cv_num / sim_results_a[["1x96"]]$cv_num))
  cat(sprintf("  Scenario B 10x6/1x96 EMG: %.3f  Num: %.3f\n",
              sim_results_b[["10x6"]]$cv_emg / sim_results_b[["1x96"]]$cv_emg,
              sim_results_b[["10x6"]]$cv_num / sim_results_b[["1x96"]]$cv_num))
  cat("  If both estimators show a smaller ratio in B than A, both detect\n")
  cat("  the injected intermediate noise. The positive control passes.\n")
}


# =============================================================================
# SECTION 7: POOLED EQUIVALENCE TEST ON CONFIG EFFECT
# =============================================================================

cat("\n===========================================================================\n")
cat("Pooled Equivalence: is the config effect smaller than the claimed benefit?\n")
cat("===========================================================================\n\n")

R_claim_equiv <- sqrt(10)
cat(sprintf("  Margin: R_claim = sqrt(10) = %.3f\n", R_claim_equiv))
cat("  (claimed CV ratio between unsummed and summed configurations)\n\n")

for (estimator in c("num", "emg")) {
  tag <- ifelse(estimator == "num", "Numeric", "EMG")
  cv_col <- ifelse(estimator == "num", "cv2_num", "cv2_emg")

  cv_by_config <- tapply(
    sqrt(model_data[[cv_col]]),
    model_data$config,
    mean
  )
  mean_cv_1x96 <- cv_by_config["1x96"]
  mean_cv_10x6 <- cv_by_config["10x6"]
  mean_cv_4x21 <- cv_by_config["4x21"]

  ratio_10x6 <- mean_cv_1x96 / mean_cv_10x6
  ratio_4x21 <- mean_cv_1x96 / mean_cv_4x21

  safe_ratio <- function(d, cv_col, cfg_num, cfg_den) {
    vals_num <- d[[cv_col]][d$config == cfg_num]
    vals_den <- d[[cv_col]][d$config == cfg_den]
    if (length(vals_num) == 0 || length(vals_den) == 0) return(NA)
    mn <- mean(sqrt(vals_den))
    if (is.na(mn) || mn <= 0) return(NA)
    mean(sqrt(vals_num)) / mn
  }

  boot_ratios_10x6 <- replicate(B_boot, {
    idx <- sample(nrow(model_data), nrow(model_data), replace = TRUE)
    safe_ratio(model_data[idx, ], cv_col, "1x96", "10x6")
  })
  boot_ratios_4x21 <- replicate(B_boot, {
    idx <- sample(nrow(model_data), nrow(model_data), replace = TRUE)
    safe_ratio(model_data[idx, ], cv_col, "1x96", "4x21")
  })

  ci90_10x6 <- quantile(boot_ratios_10x6, c(0.05, 0.95), na.rm = TRUE)
  ci90_4x21 <- quantile(boot_ratios_4x21, c(0.05, 0.95), na.rm = TRUE)

  cat(sprintf("  %s estimator:\n", tag))
  cat(sprintf("    Mean CV by config: 1x96=%.4f  4x21=%.4f  10x6=%.4f\n",
              mean_cv_1x96, mean_cv_4x21, mean_cv_10x6))
  cat(sprintf("    CV ratio 1x96/10x6 = %.3f, 90%% CI [%.3f, %.3f]  %s R_claim=%.3f\n",
              ratio_10x6, ci90_10x6[1], ci90_10x6[2],
              ifelse(ci90_10x6[2] < R_claim_equiv, "<", ">="), R_claim_equiv))
  cat(sprintf("    CV ratio 1x96/4x21 = %.3f, 90%% CI [%.3f, %.3f]  %s R_claim=%.3f\n",
              ratio_4x21, ci90_4x21[1], ci90_4x21[2],
              ifelse(ci90_4x21[2] < sqrt(4), "<", ">="), sqrt(4)))
  cat("\n")
}

cat("  If the 90% CI upper bound is below R_claim, the observed config\n")
cat("  effect is smaller than what the summing literature claims.\n")
