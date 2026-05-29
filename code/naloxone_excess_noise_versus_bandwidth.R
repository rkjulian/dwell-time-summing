# Naloxone Excess Noise versus Bandwidth.R
# LSQ coefficients from exploratory excess noise analysis (see script-provenance.md)
# Author: Fred Lytle
# Last edit: 9/10/2021
rm(list=ls())

# x-axis values
DT <- c(0.096,0.021,0.006)

# excess counting noise: use figure LSQ coefficients at mean = 1000 counts
meanCount <- 1100
exCount96 <- 58.2 + 0.0527*meanCount
exCount21 <- 24.9 + 0.0958*meanCount
exCount06 <- 19.7 + 0.117*meanCount
excessCountNoise <- c(exCount96,exCount21,exCount06)

plot(sqrt(1/DT),excessCountNoise,xlab='sqrt(1/DwellTime)',
     ylab='Excess Count Noise',pch=1,lwd=2,cex=1.5)
grid(col='darkgray')
fitCount <- lm(excessCountNoise ~ sqrt(1/DT))
coefsCount <- signif(summary(fitCount)$coef,3)
abline(fitCount,lwd=2)
legend('topleft',c(paste('Mean Count:',meanCount),
                   paste('Inter:',coefsCount[1,1]),
                   paste('Slope:',coefsCount[2,1])),bty='n')

cat('\nCount Excess Noise',
    '\nMean Count',meanCount,
    '\n\tCoeffs:',coefsCount[1,1],coefsCount[2,1],
    '\n\tErrors:',coefsCount[1,2],coefsCount[2,2])

# excess rate noise: use figure LSQ coefficients at mean = 10000 counts/sec
meanRate <- 1000
exRate96 <- 668 + 0.0573*meanRate
exRate21 <- 1320 + 0.0965*meanRate
exRate06 <- 2810 + 0.124*meanRate
excessRateNoise <- c(exRate96,exRate21,exRate06)

plot(sqrt(1/DT),excessRateNoise,xlab='sqrt(1/DwellTime)',
     ylab='Excess Rate Noise',pch=1,lwd=2,cex=1.5)
grid(col='darkgray')
fitRate <- lm(excessRateNoise ~ sqrt(1/DT))
coefsRate <- signif(summary(fitRate)$coef,3)
abline(fitRate,lwd=2)
legend('topleft',paste('Mean Rate:',meanRate),
       c(paste('Inter:',coefsRate[1,1]),
                   paste('Slope:',coefsRate[2,1])),bty='n')

cat('\nRate Excess Noise',
    '\nMean Rate:',meanRate,
    '\n\tCoeffs:',coefsRate[1,1],coefsRate[2,1],
    '\n\tErrors:',coefsRate[1,2],coefsRate[2,2])
