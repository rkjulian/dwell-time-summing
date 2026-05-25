# Flow Injection Reverse Calibration Rate Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last Edit: 10/10/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function

# user supplied input
useA50 <- TRUE # FALSE is not implemented as of 10/10/2022
target <- 50
N <- 1e8

# rate computed by rate = inter + slope*50 using values from the manuscript text
A50 <- c(0.0240 + 0.00243*50,-0.0015 + 0.00248 *50,-0.0195 + 0.00255*50)
Atable <- c(0.1430,0.1198,0.1053)

cat('\nFlow Injection Reverse Normalized:',
    ifelse(useA50,'\nCalculated at 50 ng/mL Rate','\nCalculated at Table Rate'))

# 96
mu <- c(0.0240,0.00243) # A50 coef values
# uses coef c(variance,covariance,variance) from manuscript text
covar <- matrix(data=c(0.0187^2,-1.936e-06,-1.936e-06,0.000171^2),nrow=2,ncol=2)
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
mu <- c(-0.0015,0.00248) # A50 coef values
# uses coef c(variance,covariance,variance) from manuscript text
covar <- matrix(data=c(0.0188^2,-1.964e-06,-1.964e-06,0.000173^2),nrow=2,ncol=2)
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
mu <- c(-0.0195,0.00255) # A50 coef values
# uses coef c(variance,covariance,variance) from manuscript text
covar <- matrix(data=c(0.0200^2,-2.204e-06,-2.204e-06,0.000183^2),nrow=2,ncol=2)
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


