# Chromatogram Calibration Peak Rates Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/05/2022
rm(list=ls())

# user supplied input
# the choice 'var' was not provided due to the complex propagation of variance
weighted <- 'none' # choice: 'none','1/x^2'

# data from chromatogram_peak_rates.R (analyte and d5 transitions)
# Note: the avg-values are the average of four replicate normalized areas
conc <- c(50,100,200,250,500)
rt96 <- array(0,dim=c(5,4))
rt96[1,] <- c(8355.273,9169.481,8407.73,8977.846)
rt96[2,] <- c(21454.14,18679.04,20152.7,20625.71)
rt96[3,] <- c(38602.73,36674.27,39840.1,42479.37)
rt96[4,] <- c(43676.66,46874.77,47864.88,49746.92)
rt96[5,] <- c(92393.42,88087.61,91157.42,94763.41)
rt21 <- array(0,dim=c(5,4))
rt21[1,] <- c(7531.032,6230.405,6907.982,7951.981)
rt21[2,] <- c(19627.87,16659.93,18726.65,18774.93)
rt21[3,] <- c(35471.93,39182.34,38563.4,39976.3)
rt21[4,] <- c(44674.63,41978.58,43435.12,43242.93)
rt21[5,] <- c(91689.77,91968.16,90067.12,92077.86)
rt06 <- array(0,dim=c(5,4))
rt06[1,] <- c(5400.47,5865.36,5780.342,5533.517)
rt06[2,] <- c(16518.69,13619.18,16459.95,15380.95)
rt06[3,] <- c(34573.27,34487.68,34210.47,39023.74)
rt06[4,] <- c(41447.54,41758.53,39920.78,39742.96)
rt06[5,] <- c(86573.97,81025.41,81659.93,79084.91)
d5rt96 <- array(0,dim=c(5,4))
d5rt96[1,] <- c(87631.71,89990.39,88856.1,99384.03)
d5rt96[2,] <- c(97155.76,88638.97,84234.05,93212.8)
d5rt96[3,] <- c(89863.05,85058.3,90519.37,107529.7)
d5rt96[4,] <- c(85919.64,86161.59,87753.2,91678.37)
d5rt96[5,] <- c(93273.66,92104.23,86211.57,92154.88)
d5rt21 <- array(0,dim=c(5,4))
d5rt21[1,] <- c(88809.52,89640.09,91787.09,93846.48)
d5rt21[2,] <- c(92775.95,83910.93,82145.75,92800.53)
d5rt21[3,] <- c(86196.56,91527.69,87698.15,100829.5)
d5rt21[4,] <- c(82064.27,85379.64,84840.64,87241.02)
d5rt21[5,] <- c(91602.84,93933.39,84529.79,88520.79)
d5rt06 <- array(0,dim=c(5,4))
d5rt06[1,] <- c(77486.76,81300.86,83654.98,87347.94)
d5rt06[2,] <- c(82754.96,84107.72,77999.09,85704.37)
d5rt06[3,] <- c(78277.48,82358.7,84088.82,85107.73)
d5rt06[4,] <- c(77857.17,78437.02,77815.34,77340.78)
d5rt06[5,] <- c(82295.69,82932.47,77685.42,84291.87)
norm96 <- rt96/d5rt96
norm21 <- rt21/d5rt21
norm06 <- rt06/d5rt06
avg96 <- rowMeans(norm96)
avg21 <- rowMeans(norm21)
avg06 <- rowMeans(norm06)

# calculate the least-squares weights
if(weighted == 'none'){
  weights96 <- weights21 <- weights06 <- rep(1,length(conc))
}else{
  weights96 <- weights21 <- weights06 <- 1/conc^2
}

# generate the main graph
title <- 'Chromatography Calibration Normalized Rate'
if(weighted=='1/x^2'){title <- paste0(title,'\nWeighted (1/x)^2')
}else{
  title <- paste0(title,'\nUnweighted')
}
plot(conc,avg96,xlim=c(0,500),ylim=c(0,1.1),
     main=title,font.main=1,cex.main=1,
     xlab='ng/mL',ylab='Normalized Peak Area (count)',cex=1.25)
grid(col='gray10')
points(conc,avg21,pch=16,cex=0.75)
points(conc,avg06,pch=4,cex=1.25)

# perform the least-squares fits and obtain errors and covariance
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
      '\nDwell Time:\t\t  Value       Error        Covar',
      '\n\t96 Intercept:\t',coef96[1],'; ',error96[1],'; ',cov96,
      '\n\t96 Slope:\t',coef96[2],'; ',error96[2],
      '\n\t21 Intercept:\t',coef21[1],'; ',error21[1],'; ',cov21,
      '\n\t21 Slope:\t',coef21[2],'; ',error21[2],
      '\n\t06 Intercept:\t',coef06[1],'; ',error06[1],'; ',cov06,
      '\n\t06 Slope:\t',coef06[2],';  ',error06[2],'\n',sep='')
  
  cat('\nCount Slope Ratios:\t',
      c(coef96[2],coef21[2],coef06[2])/coef96[2])
}