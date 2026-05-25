# Table 2 Slope Corrections.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/10/2022
rm(list=ls())

#   Delta t was empirically determined to be 0.3 ms (see text)
#   Slope values from Table 2 a1 column
ntd <- c(0.096,0.084,0.060)
ntdTrue <- ntd - 0.0003*c(1,4,10)
ntdCor <- ntd/ntdTrue
ntdCorNorm <- ntdCor/ntdCor[1]

rate <- c(115.5,112.8,110.3)
rateRatios <- rate/rate[1]
count <- c(11.09,9.472,6.618)
countRatios <- count/count[1]

rateRatiosCor <- rateRatios*ntdCorNorm
countRatiosCor <- countRatios*ntdCorNorm

cat('\nCorrected ntd for delta time = 0.3 ms',
    '\n\tntd =',ntd,
    '\n\tntd True =',ntdTrue,
    '\n\tntd Corrected =',ntdCor,
    '\n\tntd Corrected and Normalised =',ntdCorNorm,
    '\nSlope Values',
    '\n\tRates =',rate,
    '\n\tRates Normalized =',rateRatios,
    '\n\tCounts =',count,
    '\n\tCounts Normalized =',countRatios,
    '\nCorrected Normalizations',
    '\n\tRate Ratios =',rateRatiosCor,
    '\n\tCount Ratios =',countRatiosCor)