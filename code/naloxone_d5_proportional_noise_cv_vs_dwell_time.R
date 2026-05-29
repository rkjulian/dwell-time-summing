# Naloxone d5 Proportional Noise CV vs Dwell Time.R
#   Proportional CV values from flow injection plateau analysis (d5 transition)
# Author: Fred Lytle
# Last edit: 7/19/2022
rm(list=ls())

propCV96 <- c(0.04150306,0.04493080,0.04397789,0.04451547,0.04772897)
propCV21 <- c(0.04312719,0.04343089,0.04443251,0.04644427,0.04955814)
propCV06 <- c(0.04467535,0.04608802,0.04526295,0.04812809,0.04826732)

breaks <- seq(0.040,0.051,0.001)
h96 <- hist(round(propCV96,3),breaks=breaks)
h21 <- hist(round(propCV21,3),breaks=breaks)
h06 <- hist(round(propCV06,3),breaks=breaks)
heightMat <- rbind(h96$counts,h21$counts,h06$counts)
colnames(heightMat) <- c(as.character(h96$mids+0.0005))
barplot(height=heightMat,beside=FALSE,legend.text=c('96','21','06'),
        col=c('blue','green','red'),main='\nProportional CV vs. Dwell Time',
        font.main=1,cex.main=1,args.legend=list(bty='n'),
        xlab='Proportional CV',ylab='Count Out of Five')
grid(col='darkgray')
