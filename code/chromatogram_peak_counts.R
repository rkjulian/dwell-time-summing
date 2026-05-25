# Chromatogram Peak Counts.R
#   User functions called: getIndexSelectedFilenameChrom(),EMG(),marquardtSearch()
#       marquardtSearch() requires getChisqr() and getDeriv()
# Author: Fred Lytle
# Last edit: 10/17/2022
rm(list=ls())
library(pracma) # erfc() required within the EMG function
library(plotrix)  # needed for dual y-axis plots

# directory structure
dataPath <- file.path("data", "chromatogram_data", "")
funcPath <- file.path("code", "functions", "")

# user supplied input
Conc <- '250ng' # choices: '500ng', '250ng', '200ng', '100ng', '50ng', '20mg'
Chan <- '2x4'   # choices: '2', '2x4', '2x10', number of cmpds x number dwells
Dwell <- '21ms' # choices: '96ms' '21ms'  '6ms'
Replicate <- '1' # choices: '1','2','3','4'
eVs <- c('40eV','36eV') # don't change
Energy <- eVs[1] # choices: naloxone [1] or d5 [2]

# fixed parameters
startTime <- 1.0
endTime <- 1.3
peakStart <- 1.09
peakEnd <- 1.75
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
indSelected <- getIndexSelectedFilenameChrom(dataFiles,params)
filenames <- dataFiles[indSelected]
numDT <- length(filenames) # number of dwell times

# use the first file to obtain indices for time >= start time
name <- filenames[1]
file <- paste0(dataPath,name)
allTime <- read.delim(file,header=FALSE)$V1 # data file time is in minutes
keepInd <- which(allTime >= startTime & allTime <= endTime)
time <- allTime[keepInd]*60 # convert minutes to seconds

# obtain the rate and count data for each dwell time
rates <- array(0,c(numDT,length(time)))
for(i in 1:numDT){
  name <- filenames[i]
  file <- paste0(dataPath,name)
  data <- read.delim(file,header=FALSE)
  rates[i,] <- data$V2[keepInd] # instrument-generated data is a rate, cnts/sec
}
dwellTime <- ifelse(Dwell == '96ms',0.096,ifelse (Dwell == '21ms',0.021,0.006))
counts <- rates*dwellTime

# process the data, rows are counts vs. dwell time index
if(numDT > 1){  
  sumDTCounts <- colSums(counts)
  
  # generate the main graph
  if(Energy == '40eV' || Energy == '24eV'){
    title <- paste('\nNaloxone Count:',Conc,Chan,Dwell,'Rep =',Replicate,Energy)
  }else{
    title <- paste('\nNaloxone-d5 Count:',Conc,Chan,Dwell,'Rep =',Replicate,Energy)
  }
  
  plot(time,sumDTCounts,ylim=c(0,max(sumDTCounts)),xlab='Time (sec)',
       ylab='Sum of Dwell Time Counts',
       main=title,font.main=1,cex.main=1)
  grid(col='gray10')
  
  # fit peak to an EMG: guesses are amplitude, mean, stdev and skew
  guess <- c(max(sumDTCounts),time[which.max(sumDTCounts)],0.5,0.5)
  {
    options(warn=-1)
    out <- tryCatch(marquardtSearch(time,sumDTCounts,EMG,guess,verbose=FALSE),
                    error=function(e) out <<- NULL)
    options(warn=0)
  }
  if(!is.null(out)){ # EMG fit is okay
    # draw regression line
    y <- EMG(time,out$coefs)
    lines(time,y,lwd=2) 
  }
}

if(numDT == 1){
  counts <- drop(counts) # converts a one-row matrix into a vector
  
  # generate the main graph
  if(Energy == '40eV' || Energy == '24eV'){
    title <- paste('\nNaloxone Count:',Conc,Dwell,'Rep =',Replicate,Energy)
  }else{
    title <- paste('\nNaloxone-d5 Count:',Conc,Dwell,'Rep =',Replicate,Energy)
  }
  plot(time,counts,ylim=c(0,max(counts)),xlab='Time (sec)',
       ylab='96 ms Counts',
       main=title,font.main=1,cex.main=1)
  grid(col='gray10')
  
  # fit peak to an EMG: guesses are amplitude, mean, stdev and skew
  guess <- c(max(counts),time[which.max(counts)],0.5,0.5)
  {
    options(warn=-1)
    out <- tryCatch(marquardtSearch(time,counts,EMG,guess),
                    error=function(e) out <<- NULL)
    options(warn=0)
  }
  if(!is.null(out)){ # fit to EMG is okay
    # draw regression line
    y <- EMG(time,out$coefs)
    lines(time,y,lwd=2)
  }
}

# print the results
if(!is.null(out)){
  cat('\n',title,
      '\n\tN = ',length(time),
      '\n\tIter = ',out$iter,
      '\n\tApex Amp = ',max(y),
      '\n\tApex Location = ',time[which.max(y)],
      '\n\tNumeric Area = ',ifelse(numDT == 1,
                  sum(counts[which(time >= 1.09*60 & time <= 1.75*60)])*0.2,
                  sum(sumDTCounts[which(time >= 1.09*60 & time <= 1.75*60)])*0.2),
      '\n\tEMG Area = ',out$coefs[1],'+-',out$error[1],
      '\n\tEMG Mean = ',out$coefs[2],
      '\n\tEMG Stdev = ',out$coefs[3],
      '\n\tEMG Tau = ',out$coefs[4],
      '\n\tError of Fit = ',sqrt(out$fitVar),
      '\n\tChisqr = ',out$chi,sep='')
}else{
  peakStart <- 1.09*60
  peakEnd <- 1.75*60
  deltaT <- time[2] - time[1]
  cat('\n',title,
      '\n\tN =',length(time),
      '\n\tNumeric Area =',ifelse(numDT==1,
          sum(counts[which(time >= peakStart & time <= peakEnd)]*deltaT),
          sum(sumDTCounts[which(time >= peakStart & time <= peakEnd)]*deltaT)))
}


# plot for publication
if(FALSE){
  pdf(file.path("manuscript", "figures", "figure_05.pdf"),
      width = 6.5, height=5)
  sumDTCounts <- colSums(counts)
  avgDTRates <- colMeans(rates)
  twoord.plot(time,sumDTCounts,time,avgDTRates,xlab='Time (sec)',
              ylab='        Sum of Dwell Time Counts',
              rylab='        Average of Dwell Time Rates',
              type='p',
              mar=c(5,4,2,4)+0.1,lytickpos=c(0,200,600,1000,1400),
              rytickpos=c(round(c(0,200,600,1000,1400)/0.084)),
              lpch=1,rpch=1,lcol='black',rcol='black')
  grid(col='gray10')
  lines(time,y,lwd=2)
  dev.off()
}
