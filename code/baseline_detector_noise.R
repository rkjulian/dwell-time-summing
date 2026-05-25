# Baseline_Detector_Noise.R
#   Measures pre-injection baseline noise (t <= 2.2 min) across all
#   concentration / dwell-time configurations.  Uses NADA::ros()
#   (Regression on Order Statistics) to handle left-censored zeros:
#   the instrument cannot report fractional ions, so zero-count
#   observations are treated as left-censored at a detection limit
#   of 0.5 counts.
#
#   Individual dwell-time counts (not averages across dwells) are
#   used so that every observation has the same censoring structure.
#
#   Results are compared to the manuscript's sigma_det values
#   (29, 24, 19 counts for 96, 21, 6 ms) which were fitted from the
#   excess-noise discrepancy at 50 ng/mL (Equation 3, Figure 3).
#
# Author: Randall Julian
# Date: 2026-05-15
rm(list = ls())
library(NADA)

dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

for (file in list.files(funcPath, full.names = TRUE)) {
  source(file = file, echo = FALSE)
}
dataFiles <- list.files(dataPath)

endTime <- 2.2
Pause <- "4ms"
Trans <- "any"

configs <- list(
  list(Chan = "2",    Dwell = "96ms", numDT = 1,  td = 0.096, label = "1x96"),
  list(Chan = "2x4",  Dwell = "21ms", numDT = 4,  td = 0.021, label = "4x21"),
  list(Chan = "2x10", Dwell = "6ms",  numDT = 10, td = 0.006, label = "10x6")
)
concentrations <- c("1500ng", "1000ng", "500ng", "250ng", "50ng")
energies <- c("40eV", "36eV")
compoundNames <- c("Naloxone", "Naloxone-d5")

DL <- 0.5  # detection limit in counts

results <- data.frame(
  Compound = character(),
  Config = character(),
  Conc = character(),
  NaiveMean = numeric(),
  NaiveSD = numeric(),
  ROSMean = numeric(),
  ROSSD = numeric(),
  ZeroPct = numeric(),
  NTotal = integer(),
  stringsAsFactors = FALSE
)

for (eIdx in 1:2) {
  Energy <- energies[eIdx]
  compound <- compoundNames[eIdx]

  for (cfg in configs) {
    for (Conc in concentrations) {
      allCounts <- NULL

      for (rep in 1:3) {
        Replicate <- as.character(rep)
        params <- list(
          Conc = Conc, Chan = cfg$Chan, Dwell = cfg$Dwell,
          Pause = Pause, Replicate = Replicate,
          Trans = Trans, Energy = Energy
        )

        indSelected <- tryCatch(
          getIndexSelectedFilenameFlow(dataFiles, params),
          error = function(e) NULL
        )
        if (is.null(indSelected)) next

        filenames <- dataFiles[indSelected]
        firstFile <- read.delim(paste0(dataPath, filenames[1]), header = FALSE)
        keepInd <- which(firstFile$V1 <= endTime)
        if (length(keepInd) == 0) next

        # read each individual dwell-time trace separately
        for (i in seq_along(filenames)) {
          data <- read.delim(paste0(dataPath, filenames[i]), header = FALSE)
          rates <- data$V2[keepInd]
          counts <- rates * cfg$td  # individual dwell-time counts
          allCounts <- c(allCounts, counts)
        }
      }

      if (is.null(allCounts) || length(allCounts) == 0) next

      nZero <- sum(allCounts == 0)
      nTotal <- length(allCounts)
      pctZero <- 100 * nZero / nTotal

      naiveMean <- mean(allCounts)
      naiveSD <- sd(allCounts)

      # ROS: treat zeros as left-censored at DL
      censored <- (allCounts == 0)
      obs <- ifelse(allCounts == 0, DL, allCounts)

      rosMean <- NA
      rosSD <- NA
      if (sum(!censored) >= 3 && pctZero < 95) {
        fit <- tryCatch({
          ros(obs, censored, forwardT = "log")
        }, error = function(e) NULL)
        if (!is.null(fit)) {
          rosMean <- mean(fit)
          rosSD <- sd(fit)
        }
      }

      results <- rbind(results, data.frame(
        Compound = compound,
        Config = cfg$label,
        Conc = Conc,
        NaiveMean = naiveMean,
        NaiveSD = naiveSD,
        ROSMean = rosMean,
        ROSSD = rosSD,
        ZeroPct = pctZero,
        NTotal = nTotal,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- Output ---

cat("\n======================================================================\n")
cat("Pre-Injection Baseline: Individual Dwell-Time Counts (t <= 2.2 min)\n")
cat("Naive sd() vs NADA::ros() with DL =", DL, "counts\n")
cat("======================================================================\n")

for (compound in compoundNames) {
  sub <- results[results$Compound == compound, ]
  cat("\n---", compound, "---\n\n")
  cat(sprintf("%-7s  %-8s  %9s  %9s  %9s  %9s  %7s  %6s\n",
              "Config", "Conc", "NaiveMn", "NaiveSD", "ROS_Mn", "ROS_SD",
              "Zero%", "N"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  for (i in 1:nrow(sub)) {
    r <- sub[i, ]
    rosm <- ifelse(is.na(r$ROSMean), "   --", sprintf("%9.2f", r$ROSMean))
    ross <- ifelse(is.na(r$ROSSD),   "   --", sprintf("%9.2f", r$ROSSD))
    cat(sprintf("%-7s  %-8s  %9.2f  %9.2f  %s  %s  %6.1f%%  %6d\n",
                r$Config, r$Conc, r$NaiveMean, r$NaiveSD,
                rosm, ross, r$ZeroPct, r$NTotal))
  }
}

# --- Comparison with manuscript sigma_det ---

cat("\n\n======================================================================\n")
cat("Comparison with Manuscript sigma_det (fitted from Figure 3)\n")
cat("======================================================================\n\n")
cat("Manuscript sigma_det:  1x96 = 29 counts,  4x21 = 24 counts,  10x6 = 19 counts\n")
cat("Per individual dwell:  96ms = 29 counts,  21ms = 6.0 counts,  6ms = 1.9 counts\n\n")

cat("Baseline ROS SD (individual dwell-time counts), naloxone only:\n\n")
cat(sprintf("%-7s  %-8s  %9s  %12s\n",
            "Config", "Conc", "ROS_SD", "Manu sigma"))
cat(paste(rep("-", 42), collapse = ""), "\n")

manu_sigma <- c("1x96" = 29, "4x21" = 6.0, "10x6" = 1.9)
for (cfg in configs) {
  nalox <- results[results$Compound == "Naloxone" & results$Config == cfg$label, ]
  for (i in 1:nrow(nalox)) {
    r <- nalox[i, ]
    ross <- ifelse(is.na(r$ROSSD), "   --", sprintf("%9.2f", r$ROSSD))
    cat(sprintf("%-7s  %-8s  %s  %12.1f\n",
                r$Config, r$Conc, ross, manu_sigma[cfg$label]))
  }
}

cat("\nNote: Manuscript sigma_det values are per-sum (n*td), not per-dwell.\n")
cat("      For 4x21: 24 counts over 4 dwells = 6 counts/dwell\n")
cat("      For 10x6: 19 counts over 10 dwells = 1.9 counts/dwell\n")
cat("      These per-dwell values assume uncorrelated detector noise.\n")
