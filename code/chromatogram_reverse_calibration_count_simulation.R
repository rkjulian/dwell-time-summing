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

# count computed by count = inter + slope*50 using Table S3 1/x^2 calibration values
A50 <- c(-27.50+17.24*target,-149.9+15.39*target,-175.3+10.36*target)
Atable <- c(804.4,596.7,328.6) # values from manuscript Table 3

# printed output title
cat('\nChromatography Reverse Count',
    ifelse(useA50,'\nCalculated 50 ng/mL Rate','\nCalculated at at Table Rate'))

# 96
mu <- c(-27.50,17.24) # Table S3 1/x^2 coefs
# uses coef c(variance,covariance,covariance,variance)
covar <- matrix(data=c(69.90^2,-40.07,-40.07,0.7298^2),nrow=2,ncol=2)
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
mu <- c(-149.9,15.39) # Table S3 1/x^2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(66.54^2,-36.31,-36.31,0.6947^2),nrow=2,ncol=2)
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
mu <- c(-175.3,10.36) # Table S3 1/x^2 coefs
# uses coef c(variance,covariance,variance)
covar <- matrix(data=c(44.99^2,-16.59,-16.59,0.4697^2),nrow=2,ncol=2)
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

