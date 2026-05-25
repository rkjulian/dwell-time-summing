# Chromatogram Reverse Calibration Count Simulation.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/10/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function

# user supplied input
useA50 <- FALSE
target <- 50
N <- 1e8

# rate computed by rate = inter + slope*50 using table 5 1/x^2 calibration values
A50 <- c(-61.06+18.65*target,-182.7+16.20*target,-189.6+10.83*target)
Atable <- c(837.8,601.0,338.7)

# printed output title
cat('\nChromatography Reverse Count',
    ifelse(useA50,'\nCalculated 50 ng/mL Rate','\nCalculated at at Table Rate'))

# 96
mu <- c(-61.06,18.65) # table 5 coefs
# uses coef c(variance,covariance,covariance,variance)
covar <- matrix(data=c(76.63^2,-48.15,-48.15,0.8001^2),nrow=2,ncol=2)
a01 <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[1],Atable[1])
C <- (A - a01[,1])/a01[,2]
muC <- mean(C)
sdC <- sd(C)
cv96 <- sdC/muC
cat('\n96',
    '\tA =',A,
    '\n\tMean =',muC,
    '\n\tStd =',sdC,
    '\n\tCV =',cv96)

# 21
mu <- c(-182.7,16.20) # table 5 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(71.84^2,-42.31,-42.31,0.7500^2),nrow=2,ncol=2)
a01 <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[2],Atable[2])
C <- (A - a01[,1])/a01[,2]
muC <- mean(C)
sdC <- sd(C)
cv21 <- sdC/muC
cat('\n21',
    '\tA =',A,
    '\n\tMean =',muC,
    '\n\tStd =',sdC,
    '\n\tCV =',cv21)

# 06
mu <- c(-189.6,10.83) # table 5 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(45.08^2,-16.67,-16.67,0.4707^2),nrow=2,ncol=2)
a01 <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[3],Atable[3])
C <- (A - a01[,1])/a01[,2]
muC <- mean(C)
sdC <- sd(C)
cv06 <- sdC/muC
cat('\n06',
    '\n\tA =',A,
    '\n\tMean =',muC,
    '\n\tStd =',sdC,
    '\n\tCV =',cv06)

cat('\n\nCV Ratios:',c(cv96,cv21,cv06)/cv96)

