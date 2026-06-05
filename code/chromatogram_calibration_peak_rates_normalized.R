# Chromatogram Calibration Peak Rates Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/05/2022
rm(list=ls())

# user supplied input
# the choice 'var' was not provided due to the complex propagation of variance
weighted <- '1/x^2' # choice: 'none','1/x^2'

# data from chromatogram_numeric_areas.R (analyte and d5 transitions)
# Note: the avg-values are the average of four replicate normalized areas
conc <- c(50,100,200,250,500)
rt96 <- array(0,dim=c(5,4))
rt96[1,] <- c(7941.8667,8625.9333,7667.6667,9283.2667)
rt96[2,] <- c(19722.5333,17543.6000,19531.9333,18792.4667)
rt96[3,] <- c(36261.0000,33905.4000,37350.9333,40059.5333)
rt96[4,] <- c(40580.5333,44714.8667,43971.7333,46423.6667)
rt96[5,] <- c(84439.5333,81107.1180,87244.0000,85717.5333)
rt21 <- array(0,dim=c(5,4))
rt21[1,] <- c(7378.0500,6281.3833,6498.9833,8256.6833)
rt21[2,] <- c(18637.9333,16193.5667,18056.9000,17924.4000)
rt21[3,] <- c(32541.9500,38147.5667,36870.2167,39247.6667)
rt21[4,] <- c(41973.8833,41254.0833,41059.1000,40596.4667)
rt21[5,] <- c(88016.9833,87398.7667,86045.1000,89545.3167)
rt06 <- array(0,dim=c(5,4))
rt06[1,] <- c(5013.0667,5930.3733,5713.5000,5246.9600)
rt06[2,] <- c(16260.9400,12896.2867,15901.7000,15156.9733)
rt06[3,] <- c(32971.8400,32849.7067,32547.1933,37895.6600)
rt06[4,] <- c(40193.3333,40091.6400,37791.3867,39158.3867)
rt06[5,] <- c(84084.4000,76695.3933,77014.8467,74403.9333)
d5rt96 <- array(0,dim=c(5,4))
d5rt96[1,] <- c(75703.4000,79716.1333,80757.0000,87161.4000)
d5rt96[2,] <- c(88794.0000,77971.6000,73607.9333,81831.7333)
d5rt96[3,] <- c(77346.2000,72288.0000,81017.8000,95657.4667)
d5rt96[4,] <- c(74710.6667,76862.8000,77103.6667,79356.3333)
d5rt96[5,] <- c(83644.7333,81323.6620,76118.0667,82451.2000)
d5rt21 <- array(0,dim=c(5,4))
d5rt21[1,] <- c(81335.2333,80905.0500,85149.8667,86040.3500)
d5rt21[2,] <- c(83546.6167,76842.6167,74714.2667,83960.6333)
d5rt21[3,] <- c(79597.1000,85586.8667,79057.0167,93285.5833)
d5rt21[4,] <- c(72606.4500,77455.8333,76356.5167,79389.5500)
d5rt21[5,] <- c(83095.4500,85196.3667,75611.4667,78923.4000)
d5rt06 <- array(0,dim=c(5,4))
d5rt06[1,] <- c(70768.1867,76165.0600,78405.0467,82188.5667)
d5rt06[2,] <- c(76491.4733,76023.1667,71501.4400,79965.1933)
d5rt06[3,] <- c(72897.9400,75859.3200,76429.9933,78090.9267)
d5rt06[4,] <- c(70666.6133,72594.5533,70786.1133,69407.1133)
d5rt06[5,] <- c(76254.1933,77623.7000,71096.6733,76473.8667)
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