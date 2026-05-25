# Chromatogram Calibration Peak Counts Normalized.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 10/06/2022
rm(list=ls())

# user supplied input
# the choice 'var' was not provided due to the complex propagation of variance
weight <- '1/x^2' # choices: 'none','1/x^2'

# data obtained from "Chromatogram Peak Counts.R"
# summarized in "Chromatogram Peak Counts.docx"
#   and "Chromatogram Peak Counts d5.docx"
# Note: the avg-values are the average of four replicate normalized areas
conc <- c(50,100,200,250,500)
cnt96 <- array(0,dim=c(5,4))
cnt96[1,] <- c(802.1061,880.2702,807.1422,861.8732)
cnt96[2,] <- c(2059.597,1793.188,1934.66,1980.068)
cnt96[3,] <- c(3705.862,3520.73,3824.649,4078.019)
cnt96[4,] <- c(4192.96,4499.978,4595.029,4775.704)
cnt96[5,] <- c(8869.768,8456.41,8751.112,9097.287)
cnt21 <- array(0,dim=c(5,4))
cnt21[1,] <- c(632.6068,523.3541,580.2705,667.9665)
cnt21[2,] <- c(1648.741,1399.434,1573.039,1577.094)
cnt21[3,] <- c(2979.642,3291.317,3239.325,3358.009)
cnt21[4,] <- c(3752.669,3526.201,3648.55,3632.406)
cnt21[5,] <- c(7701.941,7725.326,7565.638,7734.54)
cnt06 <- array(0,dim=c(5,4))
cnt06[1,] <- c(324.0281,351.9216,346.8206,332.0111)
cnt06[2,] <- c(991.1216,817.151,987.5972,922.857)
cnt06[3,] <- c(2074.396,2069.261,2052.628,2341.424)
cnt06[4,] <- c(2486.853,2505.512,2395.247,2384.577)
cnt06[5,] <- c(5194.438,4861.525,4899.596,4745.094)
d5cnt96 <- array(0,dim=c(5,4))
d5cnt96[1,] <- c(8412.644,8639.077,8530.186,9540.867)
d5cnt96[2,] <- c(9326.953,8509.341,8086.469,8948.429)
d5cnt96[3,] <- c(8626.853,8165.597,8689.859,10322.86)
d5cnt96[4,] <- c(8248.285,8271.512,8424.308,8801.123)
d5cnt96[5,] <- c(8954.271,8842.006,8276.311,8846.869)
d5cnt21 <- array(0,dim=c(5,4))
d5cnt21[1,] <- c(7460,7529.767,7710.115,7883.104)
d5cnt21[2,] <- c(7793.18,7048.519,6900.243,7795.245)
d5cnt21[3,] <- c(7240.511,7688.326,7366.645,8469.682)
d5cnt21[4,] <- c(6893.398,7171.89,7126.614,7328.246)
d5cnt21[5,] <- c(7694.639,7890.405,7100.502,7435.746)
d5cnt06 <- array(0,dim=c(5,4))
d5cnt06[1,] <- c(4649.206,4878.051,5019.299,5240.877)
d5cnt06[2,] <- c(4965.298,5046.463,4679.946,5142.262)
d5cnt06[3,] <- c(4696.649,4941.522,5045.329,5106.464)
d5cnt06[4,] <- c(4671.43,4706.221,4668.92,4640.447)
d5cnt06[5,] <- c(4937.742,4975.948,4661.125,5057.512)
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