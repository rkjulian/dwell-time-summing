# chromatogram_figure_05.R
#
#   Figure 5: chromatographic peak for 250 ng/mL, 4x21 configuration, Rep 1.
#   Dual y-axis scatter of summed dwell-time counts (left) and averaged
#   dwell-time rates (right) over the raw trace. No EMG fit curve is drawn;
#   the area is obtained by rectangular-rule numeric integration over the
#   marked window (65.4 s to 78.0 s) after subtracting the pre-peak baseline.
#
#   User functions called: getIndexSelectedFilenameChrom()
rm(list = ls())
library(plotrix)  # twoord.plot for dual y-axis

# directory structure
dataPath <- file.path("data", "chromatogram_data", "")
funcPath <- file.path("code", "functions", "")
for (file in list.files(funcPath, full.names = TRUE)) {
  source(file = file, echo = FALSE)
}
dataFiles <- list.files(dataPath)

# the panel shown in the manuscript
Conc <- "250ng"
Chan <- "2x4"
Dwell <- "21ms"
Replicate <- "1"
Energy <- "40eV"
Pause <- "4ms"
Trans <- "any"
dwellTime <- 0.021   # seconds per dwell for the 4x21 config
liveTime  <- 0.084   # n * td = 4 * 0.021

# fixed parameters (data window and integration window)
startTime     <- 1.0    # data window start (minutes)
endTime       <- 1.3    # data window end   (minutes)
peakStart_min <- 1.09   # integration window start (minutes) = 65.4 s

# select the matching files
params <- list(Conc = Conc, Chan = Chan, Dwell = Dwell, Pause = Pause,
               Replicate = Replicate, Trans = Trans, Energy = Energy)
indSelected <- getIndexSelectedFilenameChrom(dataFiles, params)
filenames <- dataFiles[indSelected]
numDT <- length(filenames)

allTime <- read.delim(paste0(dataPath, filenames[1]), header = FALSE)$V1
keepInd <- which(allTime >= startTime & allTime <= endTime)
time <- allTime[keepInd] * 60   # minutes -> seconds

rates <- array(0, c(numDT, length(time)))
for (i in 1:numDT) {
  data <- read.delim(paste0(dataPath, filenames[i]), header = FALSE)
  rates[i, ] <- data$V2[keepInd]
}

sumDTCounts <- colSums(rates * dwellTime)
avgDTRates  <- colMeans(rates)

# pre-peak baseline (mean of the region before the integration window)
peakStart_s    <- peakStart_min * 60   # 65.4 s
peakEnd_s      <- max(time)            # 78.0 s
baseline_idx   <- which(time < peakStart_s)
baseline_count <- mean(sumDTCounts[baseline_idx])

# publication figure
pdf(file.path("manuscript", "figures", "figure_05.pdf"), width = 6.5, height = 5)
twoord.plot(time, sumDTCounts, time, avgDTRates, xlab = "Time (sec)",
            ylab = "        Sum of Dwell Time Counts",
            rylab = "        Average of Dwell Time Rates",
            type = "p",
            mar = c(5, 4, 2, 4) + 0.1,
            lytickpos = c(0, 200, 600, 1000, 1400),
            rytickpos = c(round(c(0, 200, 600, 1000, 1400) / liveTime)),
            lpch = 1, rpch = 1, lcol = "black", rcol = "black")
grid(col = "gray10")

# mark the integration window (65.4 s to 78.0 s) and the pre-peak baseline
abline(v = peakStart_s, lty = 2, col = "gray40")
abline(v = peakEnd_s,   lty = 2, col = "gray40")
abline(h = baseline_count, lty = 3, col = "gray40")
text(mean(c(peakStart_s, peakEnd_s)), 1300,
     "integration window", cex = 0.8, col = "gray30")
text(mean(c(min(time), peakStart_s)), baseline_count + 60,
     "baseline", cex = 0.8, col = "gray30")
dev.off()

cat(sprintf("Figure 5 written: 250 ng/mL 4x21 Rep 1, baseline = %.1f counts,\n",
            baseline_count))
cat(sprintf("  integration window %.1f-%.1f s (%d points)\n",
            peakStart_s, peakEnd_s, length(which(time >= peakStart_s & time <= peakEnd_s))))
