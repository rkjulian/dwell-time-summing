# Chromatogram Reverse Calibration Rate Simulation.R
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
A50 <- c(-635.1+194.2*target,-2175+192.9*target,-3161+180.5*target)
Atable <- c(8728,7155,5645) # values from manuscript Table 4

# printed output title
cat('\nChromatography Reverse Rate',
    ifelse(useA50,'\nCalculated 50 ng/mL Rate','\nCalculated at Table Rate'))

# 96
mu <- c(-635.1,194.2) # table 5 coefs 
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(798.0^2,-5222,-5222,8.332^2),nrow=2,ncol=2)
randomCoef <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[1],Atable[1])
C <- (A - randomCoef[,1])/randomCoef[,2]
muC <- mean(C)
sdC <- sd(C)
cv96 <- sdC/muC
cat('\n96',
    '\tArea =',A,
    '\n\tConc Mean =',muC,
    '\n\tConc Std =',sdC,
    '\n\tConc CV =',cv96)

# 21
mu <- c(-2175,192.9) # table 5 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(853.9^2,-5979,-5979,8.915^2),nrow=2,ncol=2)
randomCoef <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[2],Atable[2])
C <- (A - randomCoef[,1])/randomCoef[,2]
muC <- mean(C)
sdC <- sd(C)
cv21 <- sdC/muC
cat('\n21',
    '\tArea =',A,
    '\n\tConc Mean =',muC,
    '\n\tConc Std =',sdC,
    '\n\tCond CV =',cv21)

# 06
mu <- c(-3161,180.5) # table 5 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(752.3^2,-4640,-4640,7.853^2),nrow=2,ncol=2)
randomCoef <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[3],Atable[3])
C <- (A - randomCoef[,1])/randomCoef[,2]
muC <- mean(C)
sdC <- sd(C)
cv06 <- sdC/muC
cat('\n06',
    '\tArea =',A,
    '\n\tConc Mean =',muC,
    '\n\tConc Std =',sdC,
    '\n\tConc CV =',cv06)

cat('\n\nCV Ratios:',c(cv96,cv21,cv06)/cv96)

