# Flow Injection Reverse Calibration Rate Simulation.R
#   User functions called: none
# Author: Fred Lytle
# Last Edit: 10/10/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function

# user supplied input
useA50 <- TRUE
target <- 50
N <- 1e8

# rate computed by rate = inter + slope*50 using Table 2 1/x^2 calibration values
A50 <- c(478.9 + 123.5*target,-591.2+121.6*target,-1360+119.3*target)
Atable <- c(6503,5343,4461) # values from manuscript Table 1

# printed output title
cat('\nFlow Injection Reverse Rate:',
    ifelse(useA50,'\nCalculated at 50 ng/mL Rate','\nCalculated at Table Rate'))

# 96
mu <- c(478.9,123.5) # Table 2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(796.8^2,-3513,-3513,7.315^2),nrow=2,ncol=2)
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
mu <- c(-591.2,121.6) # Table 2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(739.8^2,-3029,-3029,6.792^2),nrow=2,ncol=2)
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
mu <- c(-1360,119.3) # Table 2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(716.9^2,-2844,-2844,6.582^2),nrow=2,ncol=2)
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

