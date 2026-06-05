# Chromatogram Reverse Calibration Normalized Simulation.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/06/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function

# user supplied input
N <- 1e8L

# rate computed by rate = inter + slope*50 using 1/x^2 calibration output values
area <- c(-0.00548607+0.00226033*50,-0.02609411+0.002302241*50,-0.04285641+0.002341541*50)

# printed output title
cat('\nChromatography Reverse Rate Normalized\nCalculated at 50 ng/mL Rate')

# 96
mu <- c(-0.00548607,0.00226033)
covar <- matrix(data=c(0.009387198^2,-7.225797e-07,-7.225797e-07,9.800522e-05^2),nrow=2,ncol=2)
randomCoef <- mvrnorm(n=N,mu=mu,Sigma=covar)
C <- (area[1] - randomCoef[,1])/randomCoef[,2]
muC <- mean(C)
sdC <- sd(C)
cv96 <- sdC/muC
cat('\n96',
    '\tArea =',area[1],
    '\n\tConc Mean =',muC,
    '\n\tConc Std =',sdC,
    '\n\tConc CV =',cv96)

# 21
mu <- c(-0.02609411,0.002302241)
covar <- matrix(data=c(0.008506465^2,-5.933516e-07,-5.933516e-07,8.881011e-05^2),nrow=2,ncol=2)
randomCoef <- mvrnorm(n=N,mu=mu,Sigma=covar)
C <- (area[2] - randomCoef[,1])/randomCoef[,2]
muC <- mean(C)
sdC <- sd(C)
cv21 <- sdC/muC
cat('\n21',
    '\tArea =',area[2],
    '\n\tConc Mean =',muC,
    '\n\tConc Std =',sdC,
    '\n\tConc CV =',cv21)

# 06
mu <- c(-0.04285641,0.002341541)
covar <- matrix(data=c(0.009638636^2,-7.61807e-07,-7.61807e-07,0.0001006303^2),nrow=2,ncol=2)
randomCoef <- mvrnorm(n=N,mu=mu,Sigma=covar)
C <- (area[3] - randomCoef[,1])/randomCoef[,2]
muC <- mean(C)
sdC <- sd(C)
cv06 <- sdC/muC
cat('\n06',
    '\tArea =',area[3],
    '\n\tConc Mean =',muC,
    '\n\tConc Std =',sdC,
    '\n\tConc CV =',cv06)

cat('\n\nCV Ratios:',c(cv96,cv21,cv06)/cv96)

