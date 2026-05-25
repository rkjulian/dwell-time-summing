# Flow Injection DT96 Correlation Investigation.R
#   User function called: getIndexSelectedFilenameFlow()
# Author: Fred Lytle
# Last edit: 10/31/2022
rm(list=ls())

# directory structure
dataPath <- file.path("data", "flow_injection_data", "")
funcPath <- file.path("code", "functions", "")

# user supplied input
Conc <- '50ng' # choices: '1500ng', '1000ng', '500ng', '250ng', '50ng'

# fixed parameters
Chan <- '2'   # choices: '2', '2x4', '2x10', number of cmpds x number dwells
Dwell <- '96ms' # choices: '96ms' '21ms'  '6ms'
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

residCor <- matrix(0,nrow=4,ncol=4)
residCov <- matrix(0,nrow=4,ncol=4)
for(rep in 1){
  Replicate <- as.character(rep)
  # identify the file that matches requested parameters
  params <- list(Conc=Conc,Chan=Chan,Dwell=Dwell,Pause=Pause,Replicate=Replicate,
                 Trans=Trans,Energy=Energy)
  indSelected <- getIndexSelectedFilenameFlow(dataFiles,params)
  filename <- dataFiles[indSelected]
  data <- read.delim(paste0(dataPath,filename),header=FALSE)
  allTime <- data$V1 # data file time is in minutes
  allRates <- data$V2
  . <- which(allTime >= 3.0)
  time <- allTime[.]
  rates <- allRates[.]
  
  # compute the residuals to get rid of drift artifact
  fit <- lm(rates~time)
  resid <- fit$resid
  
  # break the rates into four groups with each group differing by one index
  indices <- seq(1:length(rates))
  ind1 <- seq(1,449,4)
  ind1 <- ind1[-113]
  ind2 <- seq(2,449,4)
  ind3 <- seq(3,449,4)
  ind4 <- seq(4,449,4)
  resid1 <- resid[ind1]
  resid2 <- resid[ind2]
  resid3 <- resid[ind3]
  resid4 <- resid[ind4]
  residMat <- cbind(resid1,resid2,resid3,resid4)
  
  # compute the correlation and covariance matrices
  residCor <- residCor + cor(residMat)
  residCov <- residCor + cov(residMat)
}

cat('\nDT 96 Correlations for',Conc,'and Rep',rep,
    '\n')
# print(round(residCor/3,3))
print(round(residCor,3))

