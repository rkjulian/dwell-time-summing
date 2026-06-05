# chromatogram_numeric_areas.R
#
#   Primary chromatographic area extraction by rectangular-rule numeric
#   integration over the peak window, with pre-peak baseline subtraction.
#   Matches the Simpson's/trapezoidal integration used by vendor software,
#   replacing the EMG parametric fit as the source of Table 3 areas and the
#   four chromatographic calibration scripts.
#
#   Iterates 5 concentrations x 3 configs x 4 replicates x 2 energies
#   (40eV analyte, 36eV naloxone-d5). For each combination it computes:
#     - count-sum area  : sum(sum_counts[peak] - baseline) * deltaT
#     - rate-average area: sum(avg_rates[peak] - baseline) * deltaT
#
#   Output is structured for copy-paste into the four calibration scripts
#   (per-replicate areas in the variable layout they already use) plus a
#   Table 3 summary (mean, sd, CV per conc x config).
#
#   User functions called: getIndexSelectedFilenameChrom()
rm(list = ls())

# directory structure
dataPath <- file.path("data", "chromatogram_data", "")
funcPath <- file.path("code", "functions", "")
for (file in list.files(funcPath, full.names = TRUE)) {
  source(file = file, echo = FALSE)
}
dataFiles <- list.files(dataPath)

# fixed parameters (data window and integration window)
startTime     <- 1.0    # data window start (minutes)
endTime       <- 1.3    # data window end   (minutes)
peakStart_min <- 1.09   # integration window start (minutes) = 65.4 s
Pause <- "4ms"
Trans <- "any"
n_reps <- 4

# configurations: label, channels token, dwell token, dwell time (s), n dwells
configs <- list(
  list(label = "96", chan = "2",    dwell = "96ms", td = 0.096, n = 1),
  list(label = "21", chan = "2x4",  dwell = "21ms", td = 0.021, n = 4),
  list(label = "06", chan = "2x10", dwell = "6ms",  td = 0.006, n = 10)
)

# concentrations ordered to match the calibration scripts: 50,100,200,250,500
conc_strings <- c("50ng", "100ng", "200ng", "250ng", "500ng")
conc_numeric <- c(50, 100, 200, 250, 500)

# energies: analyte (naloxone) at 40 eV, internal standard (d5) at 36 eV
energy_analyte <- "40eV"
energy_d5      <- "36eV"


# =============================================================================
# Numeric area for one conc x config x replicate x energy
# =============================================================================

get_numeric_areas <- function(dataFiles, dataPath, conc, cfg, rep, energy) {
  params <- list(Conc = conc, Chan = cfg$chan, Dwell = cfg$dwell,
                 Pause = Pause, Replicate = as.character(rep),
                 Trans = Trans, Energy = energy)
  indSelected <- getIndexSelectedFilenameChrom(dataFiles, params)
  filenames <- dataFiles[indSelected]
  numDT <- length(filenames)

  allTime <- read.delim(paste0(dataPath, filenames[1]), header = FALSE)$V1
  keepInd <- which(allTime >= startTime & allTime <= endTime)
  time <- allTime[keepInd] * 60   # minutes -> seconds

  rates <- array(0, c(numDT, length(time)))
  for (i in 1:numDT) {
    data <- read.delim(paste0(dataPath, filenames[i]), header = FALSE)
    rates[i, ] <- data$V2[keepInd]   # instrument data is a rate, counts/sec
  }

  # averaged rate over the n dwells, and summed counts over the n dwells
  if (numDT > 1) avg_rates <- colMeans(rates) else avg_rates <- drop(rates)
  counts <- rates * cfg$td
  if (numDT > 1) sum_counts <- colSums(counts) else sum_counts <- drop(counts)

  deltaT      <- time[2] - time[1]
  peakStart_s <- peakStart_min * 60
  peakEnd_s   <- max(time)   # data window upper bound (78.0 s)

  # baseline: mean of the pre-peak region
  baseline_idx   <- which(time < peakStart_s)
  baseline_count <- mean(sum_counts[baseline_idx])
  baseline_rate  <- mean(avg_rates[baseline_idx])

  peak_idx <- which(time >= peakStart_s & time <= peakEnd_s)

  # rectangular-rule numeric areas with baseline subtraction
  num_count <- sum(sum_counts[peak_idx] - baseline_count) * deltaT
  num_rate  <- sum(avg_rates[peak_idx]  - baseline_rate)  * deltaT

  list(num_count = num_count, num_rate = num_rate,
       live_time = cfg$n * cfg$td, deltaT = deltaT,
       n_peak = length(peak_idx), n_base = length(baseline_idx))
}


# =============================================================================
# Collect per-replicate areas for every combination
# =============================================================================

# areas[[config_label]]$count[[energy_tag]] is a 5 x 4 matrix (conc x replicate)
areas <- list()
for (cfg in configs) {
  areas[[cfg$label]] <- list(
    live_time = cfg$n * cfg$td,
    count = list(analyte = matrix(NA, 5, 4), d5 = matrix(NA, 5, 4)),
    rate  = list(analyte = matrix(NA, 5, 4), d5 = matrix(NA, 5, 4))
  )
}

ratio_check <- data.frame(config = character(), conc = integer(), rep = integer(),
                          ratio = numeric(), stringsAsFactors = FALSE)

for (ci in seq_along(conc_numeric)) {
  for (cfg in configs) {
    for (rep in 1:n_reps) {
      a40 <- get_numeric_areas(dataFiles, dataPath, conc_strings[ci], cfg, rep, energy_analyte)
      a36 <- get_numeric_areas(dataFiles, dataPath, conc_strings[ci], cfg, rep, energy_d5)

      areas[[cfg$label]]$count$analyte[ci, rep] <- a40$num_count
      areas[[cfg$label]]$rate$analyte[ci, rep]  <- a40$num_rate
      areas[[cfg$label]]$count$d5[ci, rep]      <- a36$num_count
      areas[[cfg$label]]$rate$d5[ci, rep]       <- a36$num_rate

      ratio_check <- rbind(ratio_check, data.frame(
        config = cfg$label, conc = conc_numeric[ci], rep = rep,
        ratio = a40$num_count / a40$num_rate, stringsAsFactors = FALSE
      ))
    }
  }
}


# =============================================================================
# VERIFICATION: count/rate area ratio = live_time (n * td) for each config
# =============================================================================

cat("===========================================================================\n")
cat("Verification: count-area / rate-area ratio = live_time (n * td)\n")
cat("===========================================================================\n")
for (cfg in configs) {
  sub <- ratio_check[ratio_check$config == cfg$label, ]
  cat(sprintf("  config %s: live_time = %.3f s, observed ratio range = [%.6f, %.6f]\n",
              cfg$label, cfg$n * cfg$td, min(sub$ratio), max(sub$ratio)))
}
cat("  (Each observed ratio should equal the config live_time to machine precision,\n")
cat("   which guarantees the count and rate CVs are identical per row.)\n\n")


# =============================================================================
# TABLE 3 SUMMARY: mean, sd, CV per conc x config (analyte, count and rate)
# =============================================================================

summarize <- function(M) {
  data.frame(
    mean = apply(M, 1, mean),
    sd   = apply(M, 1, sd),
    cv   = apply(M, 1, function(x) sd(x) / mean(x))
  )
}

cat("===========================================================================\n")
cat("Table 3 summary (analyte, 40eV) -- mean, sd, CV across 4 replicate areas\n")
cat("===========================================================================\n")
for (cfg in configs) {
  cat(sprintf("\n--- Config 1x%s (live_time %.3f s) ---\n", cfg$label, cfg$n * cfg$td))
  cnt <- summarize(areas[[cfg$label]]$count$analyte)
  rat <- summarize(areas[[cfg$label]]$rate$analyte)
  cat(sprintf("  %-6s  %12s %10s %8s   %12s %10s %8s\n",
              "Conc", "Count_mean", "Count_sd", "Count_CV",
              "Rate_mean", "Rate_sd", "Rate_CV"))
  for (ci in seq_along(conc_numeric)) {
    cat(sprintf("  %-6d  %12.2f %10.3f %8.4f   %12.2f %10.3f %8.4f\n",
                conc_numeric[ci],
                cnt$mean[ci], cnt$sd[ci], cnt$cv[ci],
                rat$mean[ci], rat$sd[ci], rat$cv[ci]))
  }
}
cat("\n")


# =============================================================================
# COPY-PASTE BLOCKS FOR THE FOUR CALIBRATION SCRIPTS
# =============================================================================

vec_c <- function(x) paste(sprintf("%.4f", x), collapse = ",")

# ---- non-normalized count calibration (avg/sd of count areas) --------------
cat("===========================================================================\n")
cat("PASTE -> chromatogram_calibration_peak_counts.R\n")
cat("===========================================================================\n")
cat("conc <- c(50,100,200,250,500)\n")
for (lbl in c("96", "21", "06")) {
  m <- areas[[lbl]]$count$analyte
  cat(sprintf("avg%s <- c(%s)\n", lbl, vec_c(apply(m, 1, mean))))
}
for (lbl in c("96", "21", "06")) {
  m <- areas[[lbl]]$count$analyte
  cat(sprintf("sd%s <- c(%s)\n", lbl, vec_c(apply(m, 1, sd))))
}
cat("\n")

# ---- non-normalized rate calibration (avg/sd of rate areas) ----------------
cat("===========================================================================\n")
cat("PASTE -> chromatogram_calibration_peak_rates.R\n")
cat("===========================================================================\n")
cat("conc <- c(50,100,200,250,500)\n")
for (lbl in c("96", "21", "06")) {
  m <- areas[[lbl]]$rate$analyte
  cat(sprintf("avg%s <- c(%s)\n", lbl, vec_c(apply(m, 1, mean))))
}
for (lbl in c("96", "21", "06")) {
  m <- areas[[lbl]]$rate$analyte
  cat(sprintf("sd%s <- c(%s)\n", lbl, vec_c(apply(m, 1, sd))))
}
cat("\n")

# ---- normalized count calibration (per-replicate analyte & d5 count areas) -
cat("===========================================================================\n")
cat("PASTE -> chromatogram_calibration_peak_counts_normalized.R\n")
cat("===========================================================================\n")
print_matrix_block <- function(varname, M) {
  cat(sprintf("%s <- array(0,dim=c(5,4))\n", varname))
  for (ci in 1:5) {
    cat(sprintf("%s[%d,] <- c(%s)\n", varname, ci, vec_c(M[ci, ])))
  }
}
print_matrix_block("cnt96", areas[["96"]]$count$analyte)
print_matrix_block("cnt21", areas[["21"]]$count$analyte)
print_matrix_block("cnt06", areas[["06"]]$count$analyte)
print_matrix_block("d5cnt96", areas[["96"]]$count$d5)
print_matrix_block("d5cnt21", areas[["21"]]$count$d5)
print_matrix_block("d5cnt06", areas[["06"]]$count$d5)
cat("\n")

# ---- normalized rate calibration (per-replicate analyte & d5 rate areas) ---
cat("===========================================================================\n")
cat("PASTE -> chromatogram_calibration_peak_rates_normalized.R\n")
cat("===========================================================================\n")
print_matrix_block("rt96", areas[["96"]]$rate$analyte)
print_matrix_block("rt21", areas[["21"]]$rate$analyte)
print_matrix_block("rt06", areas[["06"]]$rate$analyte)
print_matrix_block("d5rt96", areas[["96"]]$rate$d5)
print_matrix_block("d5rt21", areas[["21"]]$rate$d5)
print_matrix_block("d5rt06", areas[["06"]]$rate$d5)
cat("\n")
