# Chromatogram Reverse Calibration Normalized Simulation.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/06/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function

# user supplied input
N <- 1e8L

# rate computed by rate = inter + slope*50 using R calibration output values
area <- c(-0.008193414+0.002150751*50,-0.02694965+0.002191745*50,-0.0410818+0.002238781*50)

# printed output title
cat('\nChromatography Reverse Rate Normalized\nCalculated at 50 ng/mL Rate')

# 96
mu <- c(-0.008193414,0.002150751)
covar <- matrix(data=c(0.008601558^2,-6.066918e-07,-6.066918e-07,8.980291e-05^2),nrow=2,ncol=2)
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
mu <- c(-0.02694965,0.002191745)
covar <- matrix(data=c(0.008989008^2,-6.625786e-07,-6.625786e-07,9.3848e-05^2),nrow=2,ncol=2)
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
mu <- c(-0.0410818,0.002238781)
covar <- matrix(data=c(0.008205447^2,-5.521008e-07,-5.521008e-07,8.566739e-05^2),nrow=2,ncol=2)
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

