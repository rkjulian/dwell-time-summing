# Chromatogram Calibration Peak Counts Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/06/2022
rm(list=ls())

# user supplied input
# the choice 'var' was not provided due to the complex propagation of variance
weight <- '1/x^2' # choices: 'none','1/x^2'

# data from chromatogram_numeric_areas.R (analyte and d5 transitions)
# Note: the avg-values are the average of four replicate normalized areas
conc <- c(50,100,200,250,500)
cnt96 <- array(0,dim=c(5,4))
cnt96[1,] <- c(762.4192,828.0896,736.0960,891.1936)
cnt96[2,] <- c(1893.3632,1684.1856,1875.0656,1804.0768)
cnt96[3,] <- c(3481.0560,3254.9184,3585.6896,3845.7152)
cnt96[4,] <- c(3895.7312,4292.6272,4221.2864,4456.6720)
cnt96[5,] <- c(8106.1952,7786.2833,8375.4240,8228.8832)
cnt21 <- array(0,dim=c(5,4))
cnt21[1,] <- c(619.7562,527.6362,545.9146,693.5614)
cnt21[2,] <- c(1565.5864,1360.2596,1516.7796,1505.6496)
cnt21[3,] <- c(2733.5238,3204.3956,3097.0982,3296.8040)
cnt21[4,] <- c(3525.8062,3465.3430,3448.9644,3410.1032)
cnt21[5,] <- c(7393.4266,7341.4964,7227.7884,7521.8066)
cnt06 <- array(0,dim=c(5,4))
cnt06[1,] <- c(300.7840,355.8224,342.8100,314.8176)
cnt06[2,] <- c(975.6564,773.7772,954.1020,909.4184)
cnt06[3,] <- c(1978.3104,1970.9824,1952.8316,2273.7396)
cnt06[4,] <- c(2411.6000,2405.4984,2267.4832,2349.5032)
cnt06[5,] <- c(5045.0640,4601.7236,4620.8908,4464.2360)
d5cnt96 <- array(0,dim=c(5,4))
d5cnt96[1,] <- c(7267.5264,7652.7488,7752.6720,8367.4944)
d5cnt96[2,] <- c(8524.2240,7485.2736,7066.3616,7855.8464)
d5cnt96[3,] <- c(7425.2352,6939.6480,7777.7088,9183.1168)
d5cnt96[4,] <- c(7172.2240,7378.8288,7401.9520,7618.2080)
d5cnt96[5,] <- c(8029.8944,7807.0716,7307.3344,7915.3152)
d5cnt21 <- array(0,dim=c(5,4))
d5cnt21[1,] <- c(6832.1596,6796.0242,7152.5888,7227.3894)
d5cnt21[2,] <- c(7017.9158,6454.7798,6275.9984,7052.6932)
d5cnt21[3,] <- c(6686.1564,7189.2968,6640.7894,7835.9890)
d5cnt21[4,] <- c(6098.9418,6506.2900,6413.9474,6668.7222)
d5cnt21[5,] <- c(6980.0178,7156.4948,6351.3632,6629.5656)
d5cnt06 <- array(0,dim=c(5,4))
d5cnt06[1,] <- c(4246.0912,4569.9036,4704.3028,4931.3140)
d5cnt06[2,] <- c(4589.4884,4561.3900,4290.0864,4797.9116)
d5cnt06[3,] <- c(4373.8764,4551.5592,4585.7996,4685.4556)
d5cnt06[4,] <- c(4239.9968,4355.6732,4247.1668,4164.4268)
d5cnt06[5,] <- c(4575.2516,4657.4220,4265.8004,4588.4320)
norm96 <- cnt96/d5cnt96
norm21 <- cnt21/d5cnt21
norm06 <- cnt06/d5cnt06
avg96 <- rowMeans(norm96)
avg21 <- rowMeans(norm21)
avg06 <- rowMeans(norm06)

# calculate the least-squares weights
if(weight == 'none'){
  weights96 <- weights21 <- weights06 <- rep(1,length(conc))
}else{
  weights96 <- weights21 <- weights06 <- 1/conc^2
}

# generate the main graph
title <- 'Chromatography Calibration Normalized Count'
if(weight=='1/x^2'){title <- paste0(title,'\nWeighted (1/x)^2')
}else{
  title <- paste0(title,'\nUnweighted')
}
plot(conc,avg96,xlim=c(0,500),ylim=c(0,1.1),
     main=title,font.main=1,cex.main=1,
     xlab='Conc (ng/mL)',ylab='Normalized Peak Area  (count x sec)',cex=1.25)
grid(col='gray10')
points(conc,avg21,pch=16,cex=0.75)
points(conc,avg06,pch=4,cex=1.25)

# perform the least-squares fits, errors, and covariance
fit96 <- lm(avg96 ~ conc,weights=weights96)
coef96 <- summary(fit96)$coef[,1]
error96 <- summary(fit96)$coef[,2]
cov96 <- summary(fit96)$cov[1,2]*summary(fit96)$sigma^2
fit21 <- lm(avg21 ~ conc,weights=weights21)
coef21 <- summary(fit21)$coef[,1]
error21 <- summary(fit21)$coef[,2]
cov21 <- summary(fit21)$cov[1,2]*summary(fit21)$sigma^2
fit06 <- lm(avg06 ~ conc,weights=weights06)
coef06 <- summary(fit06)$coef[,1]
error06 <- summary(fit06)$coef[,2]
cov06 <- summary(fit06)$cov[1,2]*summary(fit06)$sigma^2

# add regression lines to the main graph
abline(fit96)
abline(fit21)
abline(fit06)

# print the results
{
  cat(title,
      '\nDwell Time:\t\t  Value       Error         Covar',
      '\n\t96 Intercept:\t',coef96[1],'; ',error96[1],'; ',cov96,
      '\n\t96 Slope:\t',coef96[2],'; ',error96[2],
      '\n\t21 Intercept:\t',coef21[1],'; ',error21[1],'; ',cov21,
      '\n\t21 Slope:\t',coef21[2],'; ',error21[2],
      '\n\t06 Intercept:\t',coef06[1],'; ',error06[1],'; ',cov06,
      '\n\t06 Slope:\t',coef06[2],';  ',error06[2],'\n',sep='')
  
  cat('\nCount Slope Ratios:\t',
      c(coef96[2],coef21[2],coef06[2])/coef96[2])
}