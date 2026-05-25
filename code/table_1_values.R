# Table 1 Values.R
#   Computes all Table 1 values from raw flow injection data using
#   consistent pooling: pool sigma and mu across replicates first,
#   then derive sigma_pro and CV_pro from pooled statistics.
#   User functions called: getIndexSelectedFilenameFlow()
#   Author: Randall Julian
rm(list=ls())

dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")
for(file in list.files(funcPath, full.names=TRUE)) source(file, echo=FALSE)
dataFiles <- list.files(dataPath)

configs <- list(
  list(conc="1500ng", chan96="2", chan21="2x4", chan06="2x10"),
  list(conc="1000ng", chan96="2", chan21="2x4", chan06="2x10"),
  list(conc="500ng",  chan96="2", chan21="2x4", chan06="2x10"),
  list(conc="250ng",  chan96="2", chan21="2x4", chan06="2x10"),
  list(conc="50ng",   chan96="2", chan21="2x4", chan06="2x10")
)

dwells <- list(
  list(label="96ms", dwell="96ms", td=0.096, chan_key="chan96", ntd=0.096),
  list(label="21ms", dwell="21ms", td=0.021, chan_key="chan21", ntd=0.084),
  list(label="6ms",  dwell="6ms",  td=0.006, chan_key="chan06", ntd=0.06)
)

get_plateau_stats <- function(dataFiles, dataPath, conc, chan, dwell, td, rep) {
  params <- list(Conc=conc, Chan=chan, Dwell=dwell, Pause="4ms",
                 Replicate=as.character(rep), Trans="any", Energy="40eV")
  indSelected <- getIndexSelectedFilenameFlow(dataFiles, params)
  filenames <- dataFiles[indSelected]
  numDT <- length(filenames)

  time <- read.delim(paste0(dataPath, filenames[1]), header=FALSE)$V1
  rates <- array(0, c(numDT, length(time)))
  for(i in 1:numDT) {
    data <- read.delim(paste0(dataPath, filenames[i]), header=FALSE)
    rates[i,] <- data$V2
  }

  plateau <- which(time >= 3)

  # counts: individual dwell time * rate, then sum across dwells
  counts <- rates * td
  if(numDT > 1) {
    sumCounts <- colSums(counts)
  } else {
    sumCounts <- drop(counts)
  }
  fitC <- lm(sumCounts[plateau] ~ time[plateau])
  mu_count <- mean(sumCounts[plateau])
  sd_count <- sd(fitC$residuals)

  # rates: average across dwells
  if(numDT > 1) {
    avgRates <- colMeans(rates)
  } else {
    avgRates <- drop(rates)
  }
  fitR <- lm(avgRates[plateau] ~ time[plateau])
  mu_rate <- mean(avgRates[plateau])
  sd_rate <- sd(fitR$residuals)

  list(mu_count=mu_count, sd_count=sd_count, mu_rate=mu_rate, sd_rate=sd_rate)
}

cat(sprintf("%-6s %-8s %8s %8s %8s %8s %8s %8s %8s\n",
    "Conc", "DT", "mu_sum", "sig_sum", "sig_pro", "CV_pro",
    "mu_avg", "sig_avg", "CV"))
cat(paste0(rep("-", 80), collapse=""), "\n")

for(cfg in configs) {
  for(dw in dwells) {
    chan <- cfg[[dw$chan_key]]
    sds_c <- numeric(3)
    mus_c <- numeric(3)
    sds_r <- numeric(3)
    mus_r <- numeric(3)

    for(rep in 1:3) {
      s <- get_plateau_stats(dataFiles, dataPath, cfg$conc, chan, dw$dwell, dw$td, rep)
      sds_c[rep] <- s$sd_count
      mus_c[rep] <- s$mu_count
      sds_r[rep] <- s$sd_rate
      mus_r[rep] <- s$mu_rate
    }

    mu_sum  <- mean(mus_c)
    sig_sum <- sqrt(sum(sds_c^2)) / 3
    sig_pro <- sqrt(sig_sum^2 - mu_sum)
    cv_pro  <- sig_pro / mu_sum

    mu_avg  <- mean(mus_r)
    sig_avg <- sqrt(sum(sds_r^2)) / 3
    cv      <- sig_sum / mu_sum

    conc_num <- gsub("ng", "", cfg$conc)
    cat(sprintf("%-6s %-8s %8.1f %8.1f %8.2f %8.4f %8.0f %8.1f %8.3f\n",
        conc_num, dw$label, mu_sum, sig_sum, sig_pro, cv_pro,
        mu_avg, sig_avg, cv))
  }
}
