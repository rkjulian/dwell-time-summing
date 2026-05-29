# Flow Injection Calibration Rate Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/07/2022
rm(list=ls())

# user supplied input
# the choice 'var' was not provided due to the complex propagation of variance
weighted <- 'none' # choice: 'none','1/x^2'

# data from flow_injection_full_trace_rate.R (analyte and d5 transitions)
# Note: the avg-values are the average of three replicate normalized areas
conc <- c(50,250,500,1000,1500)
rt96 <- array(0,dim=c(5,3))
rt96[1,] <- c(6307.973,6545.178,6654.84)
rt96[2,] <- c(35645.29,37122,36286.24)
rt96[3,] <- c(60509.3,60420.24,60193.1)
rt96[4,] <- c(115741.5,116766.9,116220.1)
rt96[5,] <- c(175564,178903.5,178558.9)
rt21 <- array(0,dim=c(5,3))
rt21[1,] <- c(5039.117,5277.944,5711.935)
rt21[2,] <- c(33833.14,34394.75,35041.03)
rt21[3,] <- c(59612.09,58478.45,58572.45)
rt21[4,] <- c(115268.5,114373.9,114779.9)
rt21[5,] <- c(171939.1,171224.4,171727.7)
rt06 <- array(0,dim=c(5,3))
rt06[1,] <- c(4456.322,4500.709,4425.052)
rt06[2,] <- c(32946.86,33146.68,32628.33)
rt06[3,] <- c(58044.18,57027.62,56628.3)
rt06[4,] <- c(112052.3,111417.4,112418)
rt06[5,] <- c(166650.7,167501.8,166920.8)
d5rt96 <- array(0,dim=c(5,3))
d5rt96[1,] <- c(45787.06,45191.27,45458.75)
d5rt96[2,] <- c(48633.37,49992.8,49577.19)
d5rt96[3,] <- c(55267.77,55615.47,55350.59)
d5rt96[4,] <- c(46208.36,46544.78,46366.29)
d5rt96[5,] <- c(49895.2,51225.62,51229.97)
d5rt21 <- array(0,dim=c(5,3))
d5rt21[1,] <- c(44854.28,44314.9,44588.63)
d5rt21[2,] <- c(47209.49,47490.28,48038.19)
d5rt21[3,] <- c(54761.75,53431.72,53764.6)
d5rt21[4,] <- c(45522.28,45587.77,44800.99)
d5rt21[5,] <- c(48755.14,48583.06,48624.83)
d5rt06 <- array(0,dim=c(5,3))
d5rt06[1,] <- c(42351.2,42303.8,42420.1)
d5rt06[2,] <- c(45265.8,45091.72,44864.05)
d5rt06[3,] <- c(52774.89,51724.48,51857.88)
d5rt06[4,] <- c(43696.85,43156.45,43113.14)
d5rt06[5,] <- c(46149.32,45999.91,46167.24)
norm96 <- rt96/d5rt96
norm21 <- rt21/d5rt21
norm06 <- rt06/d5rt06
avg96 <- rowMeans(norm96)
avg21 <- rowMeans(norm21)
avg06 <- rowMeans(norm06)

# calculate the least-squares weights
if(weighted == 'none'){weights <- rep(1,length(conc))}else{weights <- 1/conc^2}

# generate the main graph
title <- paste0('Flow Injection Calibration Rate Normalized\nRepetition Avg')
if(weighted == '1/x^2'){title <- paste0(title,' (1/x)^2 Weighted')
}else{title <- paste0(title,' Unweighted')}
plot(conc,avg96,xlim=c(0,1500),ylim=c(0,1.1*max(avg96)),
     xlab='Conc (ng/mL)',ylab='Normalized Rates',cex=1.25,
     main=title,font.main=1,cex.main=1)
grid(col='gray10')
points(conc,avg21,pch=16,cex=0.75)
points(conc,avg06,pch=2,cex=1.25)
legend('topleft',c('96','21','06'),pch=c(1,16,2),bty='n')

# perform the least-squares fits, errors, and covariance
fit96 <- lm(avg96 ~ conc,weights=weights)
coeff96 <- summary(fit96)$coef
covar96 <- summary(fit96)$cov.unscaled[1,2]*summary(fit96)$sigma^2
fit21 <- lm(avg21 ~ conc,weights=weights)
coeff21 <- summary(fit21)$coef
covar21 <- summary(fit21)$cov.unscaled[1,2]*summary(fit21)$sigma^2
fit06 <- lm(avg06 ~ conc,weights=weights)
coeff06 <- summary(fit06)$coef
covar06 <- summary(fit06)$cov.unscaled[1,2]*summary(fit06)$sigma^2
abline(fit96)
abline(fit21)
abline(fit06)

# print results
{
  cat('\n',title,
      '\nDwell Time: \t\t  Value      Error       CV',   'Covar',
      '\n\t96Intercept:\t',coeff96[1,1],coeff96[1,2],coeff96[1,2]/coeff96[1,1],covar96,
      '\n\t96Slope:\t',coeff96[2,1],coeff96[2,2],coeff96[2,2]/coeff96[2,1],
      '\n\t21Intercept:\t',coeff21[1,1],coeff21[1,2],coeff21[1,2]/coeff21[1,1],covar21,
      '\n\t21Slope:\t',coeff21[2,1],coeff21[2,2],coeff21[2,2]/coeff21[2,1],
      '\n\t06Intercept:\t',coeff06[1,1],coeff06[1,2],coeff06[1,2]/coeff06[1,1],covar06,
      '\n\t06Slope:\t',coeff06[2,1],coeff06[2,2],coeff06[2,2]/coeff06[2,1])
  cat('\n\nRate Slope Ratios:\t',
      c(coeff96[2,1],coeff21[2,1],coeff06[2,1])/coeff96[2,1])
  cat('\n\nAmplitude at 50ng/mL',
      '\n\t96:',avg96[1],
      '\n\t21:',avg21[1],
      '\n\t06:',avg06[1])
}


