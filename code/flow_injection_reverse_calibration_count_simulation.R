# Flow Injection Reverse Calibration Count Simulation.R
#   User functions called: none
# Author: Fred Lytle
# Last Edit: 10/10/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function

# user supplied input
useA50 <- FALSE
target <- 50
N <- 1e8

# rate computed by rate = inter + slope*50 using unweighted calibration values
A50 <- c(45.98 + 11.86*target,-49.70 + 10.21*target,-81.67 + 7.155*target)
Atable <- c(624.3,448.8,267.6)

cat('\nFlow Injection Reverse Count:',
    ifelse(useA50,'\nCalculated at 50 ng/mL Count','\nCalculated at Table Count'))

# 96
mu <- c(45.98,11.86) # Table 2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(76.53^2,-32.41,-32.41,0.7026^2),nrow=2,ncol=2)
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
mu <- c(-49.70,10.21) # Table 2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(62.20^2,-21.41,-21.41,0.5710^2),nrow=2,ncol=2)
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
mu <- c(-81.67,7.155) # Table 2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(42.97^2,-10.21,-10.21,0.3945^2),nrow=2,ncol=2)
a01 <- mvrnorm(n=N,mu=mu,Sigma=covar)
A <- ifelse(useA50,A50[3],Atable[3])
C <- (A - a01[,1])/a01[,2]
muC <- mean(C)
sdC <- sd(C)
cv06 <- sdC/muC
cat('\n06',
    '\tA =',A,
    '\n\tMean =',muC,
    '\n\tStd =',sdC,
    '\n\tCV =',cv06)

cat('\n\nCV Ratios:',c(cv96,cv21,cv06)/cv96)

