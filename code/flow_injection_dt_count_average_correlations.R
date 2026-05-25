# Flow Injection DT Count Average Correlations.R
#   User function called: getIndexSelectedFilenameFlow()
# Author: Fred Lytle
# Last edit: 11/14/2022
rm(list=ls())

# directory structure
dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

# user supplied input
Conc <- '1500ng' # choices: '1500ng', '1000ng', '500ng', '250ng', '50ng'
Chan <- '2x4'   # choices: '2', '2x4', '2x10', number of cmpds x number dwells
Dwell <- '21ms' # choices: '96ms' '21ms'  '6ms'

# fixed parameters
startTime <- 3
eVs <- c('40eV','36eV') # don't change
Energy <- eVs[1] # choices: naloxone [1] or d5 [2]
Pause <- '4ms'
transitions <- c('328.1104-212.0927','328.1104-253.1204',# naloxone
                 '333.1227-212.1013','333.1227-258.1961')# naloxone-d5
Trans <- 'any'

# input functions and data filenames
for(file in list.files(funcPath,full.names = TRUE)){
  source(file=file,echo=FALSE)
}
dataFiles <- list.files(dataPath)

maxInd <- ifelse(Dwell == '21ms',4,10)
dwellTime <- ifelse(Dwell == '21ms',0.021,0.006)
covar <- array(0,c(maxInd,maxInd))
corr <- array(0,c(maxInd,maxInd))
var <- array(0,c(3,maxInd))
mu<- array(0,c(3,maxInd))
for(rep in 1:3){
  Replicate <- as.character(rep)
  # obtain the length of the data vector by processing first file
  params <- list(Conc=Conc,Chan=Chan,Dwell=Dwell,Pause=Pause,
                 Replicate=Replicate,Trans=Trans,Energy=Energy)
  indSelected <- getIndexSelectedFilenameFlow(dataFiles,params)
  filenames <- dataFiles[indSelected[1]]
  file <- paste0(dataPath,filenames)
  allRates <- read.delim(file,header=FALSE)
  keepInd <- which(allRates$V1 >= startTime)
  timeL <- length(keepInd)
  time <- allRates$V1[keepInd]
  amp <- allRates$V2[keepInd]
  counts <- array(0,c(timeL,maxInd))
  counts[,1] <- amp*dwellTime
  
  # input data file for remaining transitions
  for(ind in 2:maxInd){
    params <- list(Conc=Conc,Chan=Chan,Dwell=Dwell,Pause=Pause,
                   Replicate=Replicate,Trans=Trans,Energy=Energy)
    indSelected <- getIndexSelectedFilenameFlow(dataFiles,params)
    filenames <- dataFiles[indSelected]
    file <- paste0(dataPath,filenames[ind])
    allRates <- read.delim(file,header=FALSE)
    amp <- allRates$V2[keepInd]
    counts[,ind] <- amp*dwellTime
  }
  
  # compute statistical parameters
  mu[rep,] <- colMeans(counts)
  var[rep,] <- apply(counts,2,var)
  covar <- covar + cov(counts)/3
  corr <- corr + cor(counts)/3
}
cordiag <- array(0,maxInd-1)
for(i in 1:(maxInd-1)){cordiag[i] <- corr[i,i+1]}
avgVar <- colMeans(var)
stdev <- sqrt(avgVar)
avgMean <- colMeans(mu)


{
  pow <- round(log10(max(covar))-5)
  if(pow < 2) pow <- 0
  cat('\nCount Parameters with Replicates Averaged:',
      ' \n\tConc = ',Conc,
      '; Chan =  ',Chan,
      '; Dwell = ',Dwell,
      '; Pause = ',Pause,
      sep='')
  
  colnames(covar) <- as.character(1:maxInd)
  rownames(covar) <- as.character(1:maxInd)
  colnames(corr) <- as.character(1:maxInd)
  rownames(corr) <- as.character(1:maxInd)
  cat('\n\nMeans:',signif(avgMean,4),
      '\nStdevs:',signif(stdev,4),
      '\n')
  cat('\nCovariance/10^',as.character(pow),':\n',sep='')
  print(round(covar/10^pow))
  cat('\nCorrelation:\n')
  print(round(corr,3))
  cat('\n2nd Diagonal:\n')
  print(paste(round(cordiag,3),',',sep=''),quot=FALSE)
  cat('\nCovar Matrix Sum =',formatC(sum(covar),digits=4),
      '\nCovar Diagonal Sum =',formatC(sum(diag(covar)),digits=4),
      '\n% Variance =',100*sum(diag(covar))/sum(covar),
      '\n% Covariance =',100*(1-sum(diag(covar))/sum(covar)) )
}
