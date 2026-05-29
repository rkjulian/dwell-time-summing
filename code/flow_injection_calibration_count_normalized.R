# Flow Injection Calibration Count Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/07/2022
rm(list=ls())

# user supplied input
# the choice 'var' was not provided due to the complex propagation of variance
weighted <- 'none' # choice: 'none','1/x^2'

# data from flow_injection_full_trace_count.R (analyte and d5 transitions)
# Note: the avg-values are the average of three replicate normalized areas
conc <- c(50,250,500,1000,1500)
cnt96 <- array(0,dim=c(5,3))
cnt96[1,] <- c(605.5654,628.3371,638.8646)
cnt96[2,] <- c(3421.948,3563.712,3483.479)
cnt96[3,] <- c(5808.893,5800.343,5778.538)
cnt96[4,] <- c(11111.19,11209.63,11157.13)
cnt96[5,] <- c(16854.14,17174.74,17141.65)
cnt21 <- array(0,dim=c(5,3))
cnt21[1,] <- c(423.2859,443.3473,479.8026)
cnt21[2,] <- c(2841.984,2889.159,2943.447)
cnt21[3,] <- c(5007.415,4912.19,4920.086)
cnt21[4,] <- c(9682.554,9607.409,9641.513)
cnt21[5,] <- c(14442.89,14382.85,14425.12)
cnt06 <- array(0,dim=c(5,3))
cnt06[1,] <- c(267.3793,270.0425,265.5031)
cnt06[2,] <- c(1976.812,1988.801,1957.7)
cnt06[3,] <- c(3482.651,3421.657,3397.698)
cnt06[4,] <- c(6723.141,6685.046,6745.082)
cnt06[5,] <- c(9999.039,10050.11,10015.25)
d5cnt96 <- array(0,dim=c(5,3))
d5cnt96[1,] <- c(4395.558,4338.362,4364.04)
d5cnt96[2,] <- c(4668.804,4799.309,4932.795)
d5cnt96[3,] <- c(5305.706,5339.085,5313.657)
d5cnt96[4,] <- c(4436.002,4468.299,4451.164)
d5cnt96[5,] <- c(4789.939,4917.66,4918.077)
d5cnt21 <- array(0,dim=c(5,3))
d5cnt21[1,] <- c(3767.759,3722.452,3745.445)
d5cnt21[2,] <- c(3965.597,3989.183,4035.208)
d5cnt21[3,] <- c(4599.987,4488.264,4516.226)
d5cnt21[4,] <- c(3823.872,3829.372,3763.283)
d5cnt21[5,] <- c(4095.432,4080.977,4084.486)
d5cnt06 <- array(0,dim=c(5,3))
d5cnt06[1,] <- c(2541.072,2538.228,2545.206)
d5cnt06[2,] <- c(2715.948,2705.503,2691.843)
d5cnt06[3,] <- c(3166.493,3103.469,3111.473)
d5cnt06[4,] <- c(2621.811,2589.387,2586.788)
d5cnt06[5,] <- c(2768.959,2759.994,2770.034)
norm96 <- cnt96/d5cnt96
norm21 <- cnt21/d5cnt21
norm06 <- cnt06/d5cnt06
avg96 <- rowMeans(norm96)
avg21 <- rowMeans(norm21)
avg06 <- rowMeans(norm06)

# calculate the least-squares weights
if(weighted == 'none'){weights <- rep(1,length(conc))}else{weights <- 1/conc^2}

# generate the main graph
title <- paste0('Flow Injection Calibration Count Normalized\nRepetition Avg')
if(weighted=='1/x^2'){title <- paste0(title,' (1/x)^2 Weighted')
}else{title <- paste0(title,' Unweighted')}
plot(conc,avg96,xlim=c(0,1500),ylim=c(0,1.1*max(avg96)),
     xlab='Conc (ng/mL)',ylab='Normalized Counts',cex=1.25,
     main=title,font.main=1,cex.main=1)
grid(col='gray10')
points(conc,avg21,pch=16,cex=0.75)
points(conc,avg06,pch=2,cex=1.25)
legend('topleft',c('96','21','06'),pch=c(1,16,2),bty='n')

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

{
  cat('\n',title,
      '\nDwell Time:\t\tValue;\t  Error;\tCV;\tCovar',
      '\n\t96 Intercept:\t',coeff96[1,1],coeff96[1,2],coeff96[1,2]/coeff96[1,1],covar96,
      '\n\t96 Slope:\t',coeff96[2,1],coeff96[2,2],coeff96[2,2]/coeff96[2,1],
      '\n\t21Intercept:\t',coeff21[1,1],coeff21[1,2],coeff21[1,2]/coeff21[1,1],covar21,
      '\n\t21Slope:\t',coeff21[2,1],coeff21[2,2],coeff21[2,2]/coeff21[2,1],
      '\n\t06Intercept:\t',coeff06[1,1],coeff06[1,2],coeff06[1,2]/coeff06[1,1],covar06,
      '\n\t06Slope:\t',coeff06[2,1],coeff06[2,2],coeff06[2,2]/coeff06[2,1]
  )
  
  cat('\n\nCount Slope Ratios:\t',
      c(coeff96[2,1],coeff21[2,1],coeff06[2,1])/coeff96[2,1])
  
}
