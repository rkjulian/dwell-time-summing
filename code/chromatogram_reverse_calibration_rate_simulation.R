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

# rate computed by rate = inter + slope*50 using Table S3 1/x^2 calibration values
A50 <- c(-286.5+179.5*target,-1785+183.2*target,-2921+172.7*target)
Atable <- c(8380,7104,5476) # values from manuscript Table 3

# printed output title
cat('\nChromatography Reverse Rate',
    ifelse(useA50,'\nCalculated 50 ng/mL Rate','\nCalculated at Table Rate'))

# 96
mu <- c(-286.5,179.5) # Table S3 1/x^2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(728.1^2,-4348,-4348,7.602^2),nrow=2,ncol=2)
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
mu <- c(-1785,183.2) # Table S3 1/x^2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(792.2^2,-5146,-5146,8.271^2),nrow=2,ncol=2)
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
mu <- c(-2921,172.7) # Table S3 1/x^2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(749.8^2,-4609,-4609,7.828^2),nrow=2,ncol=2)
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

