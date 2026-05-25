# Flow Injection Full Trace Count.R
#   User functions called: getIndexSelectedFilenameChrom()
# Author: Fred Lytle
# Last edit: 10/07/2022   # fixed baseline CV equation
rm(list=ls())

# directory structure
dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

# user supplied input
Conc <- '250ng' # choices: '1500ng', '1000ng', '500ng', '200ng', '50ng'
Chan <- '2x4'   # choices: '2', '2x4', '2x10', number of cmpds x number dwells
Dwell <- '21ms' # choices: '96ms' '21ms'  '6ms'
Replicate <- '1' # choices: '1','2','3','4'
eVs <- c('40eV','36eV') # don't change
Energy <- eVs[1] # choices: naloxone [1] or d5 [2]

# fixed parameters
Pause <- '4ms'
transitions <- c('328.1104-212.0927','333.1227-258.1961')# naloxone-d5
Trans <- 'any'


# input functions and data filenames
for(file in list.files(funcPath,full.names = TRUE)){
  source(file=file,echo=FALSE)
}
dataFiles <- list.files(dataPath)

# identify files that match requested parameters
params <- list(Conc=Conc,Chan=Chan,Dwell=Dwell,Pause=Pause,Replicate=Replicate,
               Trans=Trans,Energy=Energy)
indSelected <- getIndexSelectedFilenameFlow(dataFiles,params)
filenames <- dataFiles[indSelected]
numDT <- length(filenames)

# use the first file to obtain indices for time >= start time
name <- filenames[1]
file <- paste0(dataPath,name)
time <- read.delim(file,header=FALSE)$V1 # data file time is in minutes

# obtain the rate and count data for each dwell time
rates <- array(0,c(numDT,length(time)))
for(i in 1:numDT){
  name <- filenames[i]
  file <- paste0(dataPath,name)
  data <- read.delim(file,header=FALSE)
  rates[i,] <- data$V2 # instrument-generated data is a rate, cnts/sec
}
dwellTime <- ifelse(Dwell == '96ms',0.096,ifelse (Dwell == '21ms',0.021,0.006))
counts <- rates*dwellTime

# process the data, rows are rate or count vs. time for one dwell time
if(numDT > 1){ #F){# 
  sumDTCounts <- colSums(counts)
  
  # generate the main graph
  if(Energy == '40eV' || Energy == '24eV'){
    title <- paste('\nNaloxone Count:',Conc,Chan,Dwell,'Rep =',Replicate,Energy)
  }else{
    title <- paste('\nNaloxone-d5 Count:',Conc,Chan,Dwell,'Rep =',Replicate,Energy)
  }
  plot(time,sumDTCounts,ylim=c(0,max(sumDTCounts)),xlab='Time (min)',
       ylab='Sum of Dwell Time Counts',
       main=title,font.main=1,cex.main=1)
  grid(col='gray10')
  abline(v=c(2.2,3.0),lty=2,lwd=2)

  # fit counts and compute stdev stats from residuals
  fit <- lm(sumDTCounts[which(time >= 3)] ~ time[which(time >= 3)])
  inter <- fit$coef[1]
  slope <- fit$coef[2]
  resid <- fit$residuals
  abline(fit,lwd=2,lty=2)
  
  # calculate statistical parameters
  meanCount <- mean(sumDTCounts[which(time >= 3)])
  stdevCount <- sd(resid)
  cvCount <- stdevCount/meanCount
  meanBase <- mean(sumDTCounts[which(time <= 2.2)])
  stdevBase <- sd(sumDTCounts[which(time <= 2.2)])
  cvBase <- stdevBase/meanBase
}

if(numDT == 1){
  counts <- drop(counts)
  
  # generate the main graph
  if(Energy == '40eV' || Energy == '24eV'){
    title <- paste('\nNaloxone Count:',Conc,Dwell,'Rep =',Replicate,Energy)
  }else{
      title <- paste('\nNaloxone-d5 Count:',Conc,Dwell,'Rep =',Replicate,Energy)
  }
  # DT rate vs. time
  plot(time,counts,ylim=c(0,max(counts)),xlab='Time (min)',
       ylab='96ms Counts',
       main=paste0(title),
       font.main=1,cex.main=1)
  grid(col='gray10')
  abline(v=c(2.2,3.0),lty=2,lwd=2)

  # fit counts and compute stdev stats from residuals
  fit <- lm(counts[which(time >= 3)] ~ time[which(time >= 3)])
  inter <- fit$coef[1]
  slope <- fit$coef[2]
  resid <- fit$residuals
  abline(fit,lwd=2,lty=2)

  # calculate statistical parameters
  meanCount <- mean(counts[which(time >= 3)])
  stdevCount <- sd(resid)
  cvCount <- stdevCount/meanCount
  meanBase <- mean(counts[which(time <= 2.2)])
  stdevBase <- sd(counts[which(time <= 2.2)])
  cvBase <- stdevBase/meanBase
}

slopeRatio <- slope*(4.5-3.0)/meanCount
cat('\n',title,
    '\n\tPlateau N =',length(time[which(time >= 3.0)]),
    '\n\tIntercept =',inter,
    '\n\tSlope =',slope,
    '\n\tCount Mean (t >= 3 min) =',round(meanCount,1),
    '\n\tSlope Ratio =',100*signif(slopeRatio,4),'%',
    '\n\tResidual Count Stdev  (t >= 3 min) =',stdevCount,
    '\n\tResidual Count CV  (t >= 3 min) =',round(cvCount,4),
    '\n\tPoisson Noise =',sqrt(meanCount),
    '\n\tExcess Noise =',sqrt(stdevCount^2 - meanCount),
    '\n\tExcess Noise CV =',sqrt(stdevCount^2 - meanCount)/meanCount,
    '\n\tBaseline Mean (t <= 2.2 min) =',round(meanBase,1),
    '\n\tBaseline Stdev (t <= 2.2 min) =',round(stdevBase,1),
    '\n\tBaseline CV (t <= 2.2 min) =',round(cvBase,3))

# plot for publication
if(FALSE){
    pdf(file.path("manuscript", "figures", "figure_02.pdf"),
        width = 6.5, height=5)
    library(plotrix)
    sumDTCounts <- colSums(counts)
    avgDTRates <- colMeans(rates)
    twoord.plot(time,sumDTCounts,time,avgDTRates,xlab='Time (min)',
                ylab='        Sum of Dwell Time Counts',
                rylab='        Average of Dwell Time Rates',
                type='p',
                mar=c(5,4,2,4)+0.1,lytickpos=c(0,500,1500,2500,3500),
                rytickpos=c(round(c(0,500,1500,2500,3500)/0.084)),
                lpch=1,rpch=1,lcol='black',rcol='black')
    grid(col='gray10')
    abline(v=3.0,lty=2,lwd=2)
    abline(v=2.2,lty=2,lwd=2)
    
    # fit counts and compute stdev stats from residuals
    fit <- lm(sumDTCounts[which(time >= 3)] ~ time[which(time >= 3)])
    inter <- fit$coef[1]
    slope <- fit$coef[2]
    resid <- fit$residuals
    abline(fit,lwd=2,lty=2)
    dev.off()
}
