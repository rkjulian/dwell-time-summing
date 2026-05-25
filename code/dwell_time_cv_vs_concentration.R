# Dwell Time CV vs Concentration.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 11/04/2022
rm(list=ls())

#   Data from ..\Flow Injection\Flow Injection Count Sum Performance.docx
#     Table 1 columns 3 & 6
conc <- c(50,250,500,1000,1500)
mean96 <- c(624.3,3490,5796,11159,17057)
mean21 <- c(448.8,2891,4946,9644,14417)
mean06 <- c(267.6,1975,3434,6718,10021)
obsNoise96 <- c(45.33,168.1,247.0,476.4,661.2)
obsNoise21 <- c(36.42,140.5,208.4,402.2,600.8)
obsNoise06 <- c(27.10,100.4,153.4,276.0,392.0)
conNoise <- 300
# these values were optimized to agree with the drop in covariance
#   in manuscript Table 7.  
conNoise96 <- conNoise*0.096
conNoise21 <- conNoise*0.079
conNoise06 <- conNoise*0.063
proNoise96All <- sqrt(obsNoise96^2-mean96-conNoise96^2)
proNoise21All <- sqrt(obsNoise21^2-mean21-conNoise21^2)
proNoise06All <- sqrt(obsNoise06^2-mean06-conNoise06^2)
proNoise96 <- sqrt(obsNoise96^2-mean96)
proNoise21 <- sqrt(obsNoise21^2-mean21)
proNoise06 <- sqrt(obsNoise06^2-mean06)
cv96All <- proNoise96All/mean96
cv21All <- proNoise21All/mean21
cv06All <- proNoise06All/mean06
cv96 <- proNoise96/mean96
cv21 <- proNoise21/mean21
cv06 <- proNoise06/mean06

plot(conc,cv96,ylim=c(0.03,0.09),cex=1.5,
     xlab='Concentration (ng/mL)',ylab='Proportional CV',
     main='\nCV of Count Sums',font.main=1,cex.main=1)
grid(col='darkgray')
points(conc,cv21,pch=16,col='gray50')
points(conc,cv06,pch=4)
points(conc[1:2],cv96All[1:2],pch=1,lwd=2,cex=1.25,col='red')
points(conc[1:2],cv21All[1:2],pch=16,col='red')
points(conc[1:2],cv06All[1:2],pch=4,col='red')
legend('topright',title='Dwell Time',c('96','21','06'),pch=c(1,16,4),bty='n')

# publication figure
if(FALSE){
  pdf(file.path("manuscript", "figures", "figure_03.pdf"),
      width = 6.5, height=5)
  constRatio <- round(c(conNoise96,conNoise21,conNoise06)/conNoise96,3)
  cat('\nConstant Noise Counts =',c(conNoise96,conNoise21,conNoise06))
  cat('\nConst Ratio =',constRatio)
  cat('\nCV Line is at',round(mean(cv21All[3:5]),4))
  cat('\ncv96 circle, cv21 square, cv06 triangle')
  cat('\nCV after subtraction, plus sign')
  par(mar=c(5,4,2,2)+0.1)# reduce upper margin
  {
    plot(conc,cv96,ylim=c(0.03,0.085),cex=1.5,lwd=1,pch=1,
       xlab='Concentration (ng/mL)',ylab='Coefficient of Variation')
    grid(col='gray10')
    points(conc,cv21,pch=16)
    points(conc,cv06,pch=2)
    points(conc[1:2],cv96All[1:2],pch=3,lwd=2,cex=1.25)
    abline(h=mean(cv21[3:5]),lty=1)
  }
  par(mar=c(5,4,4,2)+0.1)# restore upper margin
  dev.off()
}
