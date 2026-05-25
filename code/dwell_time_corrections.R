# Dwell Time Corrections.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 6/24/2022
rm(list=ls())

delT <- 0.00030
rate <- c(115.5,112.8,110.3)
rateAdj <- c(rate[1]*0.096/(0.096-delT),rate[2]*0.021/(0.021-delT),
              rate[3]*0.006/(0.006-delT))
rateNorm <- (rateAdj/rateAdj[1])
count <- c(11.09,9.472,6.618)
countAdj <- rateAdj*c(0.096,0.084,0.060)
countNorm <- (countAdj/countAdj[1])

# 96 ms
true96 <- 0.096 - delT
ratio96 <- 0.096/true96
true21 <- 0.021 - delT
ratio21 <- 0.021/true21
true06 <- 0.006 - delT
ratio06 <- 0.006/true06
cat('\nDelta Time =',delT)
cat('\n96 Rate:',
    '\n\tTrue =',true96,
    '\n\tRatio =',ratio96,
    '\n\tCorrected Slope =',rateAdj[1],
    '\n21 Rate:',
    '\n\tTrue =',true21,
    '\n\tRatio =',ratio21,
    '\n\tCorrected Slopes =',rateAdj[2],
    '\n06 Rate:',
    '\n\tTrue =',true06,
    '\n\tRatio =',ratio06,
    '\n\tCorrected Slopes =',rateAdj[3],
    '\nSlope Ratios =',rateNorm
)

cat('\n96 Count:',
    '\n\tTrue =',true96,
    '\n\tRatio =',ratio96,
    '\n\tCorrected Slope =',countAdj[1],
    '\n21 Count:',
    '\n\tTrue =',true21,
    '\n\tRatio =',ratio21,
    '\n\tCorrected Slopes =',countAdj[2],
    '\n06 Count:',
    '\n\tTrue =',true06,
    '\n\tRatio =',ratio06,
    '\n\tCorrected Slopes =',countAdj[3],
    '\nSlope Ratios =',countNorm
)
