# Table7_Raw_vs_Detrended_Comparison.R
#   Compares Table 7 correlation values computed on raw plateau rates
#   vs. linearly detrended residuals. Runs all 10 combinations
#   (5 concentrations x 2 dwell-time configs) and produces a
#   side-by-side comparison table.
#
# Author: Randall Julian
# Date: 2026-05-15
rm(list = ls())

dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

for (file in list.files(funcPath, full.names = TRUE)) {
  source(file = file, echo = FALSE)
}
dataFiles <- list.files(dataPath)

startTime <- 3
Pause <- "4ms"
Energy <- "40eV"
Trans <- "any"

configs <- list(
  list(Chan = "2x4",  Dwell = "21ms", maxInd = 4),
  list(Chan = "2x10", Dwell = "6ms",  maxInd = 10)
)
concentrations <- c("1500ng", "1000ng", "500ng", "250ng", "50ng")

results <- data.frame(
  DwellConfig = character(),
  Conc = character(),
  AvgCorr_Raw = numeric(),
  AvgCorr_Detrended = numeric(),
  PctVar_Raw = numeric(),
  PctVar_Detrended = numeric(),
  stringsAsFactors = FALSE
)

for (cfg in configs) {
  Chan <- cfg$Chan
  Dwell <- cfg$Dwell
  maxInd <- cfg$maxInd

  for (Conc in concentrations) {

    covar_raw <- array(0, c(maxInd, maxInd))
    corr_raw  <- array(0, c(maxInd, maxInd))
    covar_det <- array(0, c(maxInd, maxInd))
    corr_det  <- array(0, c(maxInd, maxInd))

    for (rep in 1:3) {
      Replicate <- as.character(rep)
      params <- list(
        Conc = Conc, Chan = Chan, Dwell = Dwell,
        Pause = Pause, Replicate = Replicate,
        Trans = Trans, Energy = Energy
      )
      indSelected <- getIndexSelectedFilenameFlow(dataFiles, params)
      filenames <- dataFiles[indSelected]

      firstFile <- read.delim(paste0(dataPath, filenames[1]), header = FALSE)
      keepInd <- which(firstFile$V1 >= startTime)
      timeL <- length(keepInd)
      time <- firstFile$V1[keepInd]

      rates <- array(0, c(timeL, maxInd))
      resid <- array(0, c(timeL, maxInd))

      for (ind in 1:maxInd) {
        allRates <- read.delim(paste0(dataPath, filenames[ind]), header = FALSE)
        amp <- allRates$V2[keepInd]
        rates[, ind] <- amp
        fit <- lm(amp ~ time)
        resid[, ind] <- fit$residuals
      }

      covar_raw <- covar_raw + cov(rates) / 3
      corr_raw  <- corr_raw  + cor(rates) / 3
      covar_det <- covar_det + cov(resid) / 3
      corr_det  <- corr_det  + cor(resid) / 3
    }

    # adjacent-dwell (super-diagonal) average correlation
    cordiag_raw <- sapply(1:(maxInd - 1), function(i) corr_raw[i, i + 1])
    cordiag_det <- sapply(1:(maxInd - 1), function(i) corr_det[i, i + 1])

    pctvar_raw <- 100 * sum(diag(covar_raw)) / sum(covar_raw)
    pctvar_det <- 100 * sum(diag(covar_det)) / sum(covar_det)

    label <- ifelse(Dwell == "21ms", "4x21", "10x6")
    results <- rbind(results, data.frame(
      DwellConfig = label,
      Conc = Conc,
      AvgCorr_Raw = round(mean(cordiag_raw), 3),
      AvgCorr_Detrended = round(mean(cordiag_det), 3),
      PctVar_Raw = round(pctvar_raw, 1),
      PctVar_Detrended = round(pctvar_det, 1),
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n====================================================================\n")
cat("Table 7 Comparison: Raw Plateau Rates vs Linearly Detrended Residuals\n")
cat("====================================================================\n\n")

cat("Adjacent-Dwell Average Correlation:\n\n")
cat(sprintf("%-10s  %-8s  %10s  %10s  %8s\n",
            "Config", "Conc", "Raw", "Detrended", "Delta"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (i in 1:nrow(results)) {
  r <- results[i, ]
  delta <- r$AvgCorr_Raw - r$AvgCorr_Detrended
  cat(sprintf("%-10s  %-8s  %10.3f  %10.3f  %+8.3f\n",
              r$DwellConfig, r$Conc, r$AvgCorr_Raw, r$AvgCorr_Detrended, delta))
}

cat("\n\nVariance as % of Total (Var + Cov):\n\n")
cat(sprintf("%-10s  %-8s  %10s  %10s  %8s\n",
            "Config", "Conc", "Raw", "Detrended", "Delta"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (i in 1:nrow(results)) {
  r <- results[i, ]
  delta <- r$PctVar_Raw - r$PctVar_Detrended
  cat(sprintf("%-10s  %-8s  %9.1f%%  %9.1f%%  %+7.1f%%\n",
              r$DwellConfig, r$Conc, r$PctVar_Raw, r$PctVar_Detrended, delta))
}
cat("\n")
