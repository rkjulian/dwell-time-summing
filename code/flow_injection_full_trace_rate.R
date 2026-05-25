# Flow Injection Full Trace Rate.R
#   User function called: getIndexSelectedFilenameFlow()
# Author: Fred Lytle
# Last edit: 10/31/2022
rm(list=ls())

# directory structure
dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

# user supplied input
Conc <- '1500ng' # choices: '1500ng', '1000ng', '500ng', '250ng', '50ng'
Chan <- '2x4'   # choices: '2', '2x4', '2x10', number of cmpds x number dwells
Dwell <- '21ms' # choices: '96ms' '21ms'  '6ms'
Replicate <- '1' # choices: '1','2','3'
eVs <- c('40eV','36eV') # don't change
Energy <- eVs[1] # choices: naloxone [1] or d5 [2]

# fixed parameters
Pause <- '4ms'
transitions <- c('328.1104-212.0927','328.1104-253.1204',# naloxone
                 '333.1227-212.1013','333.1227-258.1961')# naloxone-d5
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

# use the first file to obtain time
name <- filenames[1]
file <- paste0(dataPath,name)
time <- read.delim(file,header=FALSE)$V1 # data file time is in minutes

# obtain the rate for each dwell time
rates <- array(0,c(numDT,length(time)))
for(i in 1:numDT){
  name <- filenames[i]
  file <- paste0(dataPath,name)
  data <- read.delim(file,header=FALSE)
  rates[i,] <- data$V2 # instrument-generated data is a rate, cnts/sec
}

# process the data, rows are rate vs. time for one dwell time
dwellTime <- ifelse(Dwell == '96ms',0.096,ifelse (Dwell == '21ms',0.084,0.06))
if(numDT > 1){ #F){# 
  avgDTRates <- colMeans(rates)
  
  # generate the main graph
  if(Energy == '40eV' || Energy == '24eV'){
    title <- paste('\nNaloxone Rate:',Conc,Chan,Dwell,'Rep =',Replicate,Energy)
  }else{
    title <- paste('\nNaloxone-d5 Rate:',Conc,Chan,Dwell,'Rep =',Replicate,Energy)
  }
  plot(time,avgDTRates,ylim=c(0,max(avgDTRates)),xlab='Time (min)',
       ylab='Avg of Dwell Time Rates',
       main=title,font.main=1,cex.main=1)
  grid(col='gray10')
  abline(v=c(2.2,3.0),lty=2,lwd=2)
  
  # fit rates and compute stdev stats from residuals
  fit <- lm(avgDTRates[which(time >= 3)] ~ time[which(time >= 3)])
  inter <- fit$coef[1]
  slope <- fit$coef[2]
  resid <- fit$residuals
  abline(fit,lwd=2,lty=2)
  
  # calculate statistical parameters
  meanRate <- mean(avgDTRates[which(time >= 3)])
  stdevRate <- sd(resid)
  cvRate <- stdevRate/meanRate
  meanBase <- mean(avgDTRates[which(time <= 2.2)])
  stdevBase <- sd(avgDTRates[which(time <= 2.2)])
  cvBase <- stdevBase/meanBase
}

if(numDT == 1){
  # generate main graph
  if(Energy == '40eV' || Energy == '24eV'){
    title <- paste('\nNaloxone Rate:',Conc,Dwell,'Rep =',Replicate,Energy)
  }else{
    title <- paste('\nNaloxone-d5 Rate:',Conc,Dwell,'Rep =',Replicate,Energy)
  }
  # DT rate vs. time
  plot(time,rates,ylim=c(0,max(rates)),xlab='Time (min)',
       ylab='96ms Rates',
       main=paste0(title),
       font.main=1,cex.main=1)
  grid(col='gray10')
  abline(v=c(2.2,3.0),lty=2,lwd=2)
  
  # fit rates and compute stdev stats from residuals
  fit <- lm(rates[which(time >= 3)] ~ time[which(time >= 3)])
  inter <- fit$coef[1]
  slope <- fit$coef[2]
  resid <- fit$residuals
  abline(fit,lwd=2,lty=2)

  # calculate statistical parameters
  meanRate <- mean(rates[which(time >= 3)])
  stdevRate <- sd(resid)
  cvRate <- stdevRate/meanRate
  meanBase <- mean(rates[which(time <= 2.2)])
  stdevBase <- sd(rates[which(time <= 2.2)])
  cvBase <- stdevBase/meanBase
}

# print the results
slopeRatio <- slope*(4.5-3.0)/meanRate
cat('\n',title,
    '\n\t\n\tPlateau N =',length(time[which(time >= 3.0)]),
    '\n\tIntercept =',inter,
    '\n\tSlope =',slope,
    '\n\tRate Mean (t >= 3 min) =',meanRate,
    '\n\tSlope Ratio =',100*signif(slopeRatio,4),'%',
    '\n\tResidual Rate Stdev  (t >= 3 min) =',stdevRate,
    '\n\tResidual Rate CV  (t >= 3 min) =',round(cvRate,4),
    # rate equivalents have to be used to make units the same in subtraction
    '\n\tPoisson Noise =',sqrt(meanRate*0.096)/0.096,
    '\n\tProportional Noise =',sqrt(stdevRate^2 - meanRate/0.096),
    '\n\tProportional Noise CV =',sqrt(stdevRate^2 - meanRate/0.096)/meanRate,
    '\n\tBaseline Mean (t <= 2.2 min) =',round(meanBase,1),
    '\n\tBaseline Stdev (t <= 2.2 min) =',round(stdevBase,1),
    '\n\tBaseline CV (t <= 2.2 min) =',round(cvBase,3))

