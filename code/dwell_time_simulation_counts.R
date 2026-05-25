# Dwell Time Simulation Counts.R
# Author: Fred Lytle
# Last edit: 12/3/2021
rm(list=ls())

setSeed <- TRUE
seed <- 314159
muRate <- 177676#17500
noiseCV <- 0.067
numReps <- 10000
CT <- 0.1
muCntsCT <- muRate*CT
numIons <- ceiling(muCntsCT*2)#1.6

cntsCT <- array(0,numReps)
cntsDT96 <- array(0,numReps)
cntsDT21 <- array(0,c(numReps,4))
cntsDT06 <- array(0,c(numReps,10))
offset <- 0.01
if(setSeed) set.seed(seed)
for(i in 1:numReps){
  if(i %% 1000 == 0) cat('\n',100*i/numReps,'%',sep='')
  ############# original code
  # noisyRates <- rnorm(numIons,muRate,noiseCV*muRate)
  noisyRates <- rnorm(1,muRate,noiseCV*muRate)
  deltaTime <- rexp(numIons,noisyRates)
  #############
  # noisyRates <- rnorm(numIons,muRate,noiseCV*muRate)
  # noiseSpec <- fft(noisyRates)
  # noiseSpec[(14000-13990):(14000+13990)] <- 0
  # filteredRates <- Re(fft(noiseSpec,inverse=T)/numIons)
  # normRates <- (filteredRates-mean(filteredRates))*sd(noisyRates)/
  #                 sd(filteredRates) + mean(noisyRates)
  # deltaTime <- rexp(numIons,normRates)
  #############
  arrivalTime <- cumsum(deltaTime)
  if(max(arrivalTime) < 0.11) stop('Too Few Ions!')
  cntsCT[i] <- length(which(arrivalTime >= 0+offset & 
                              arrivalTime <= 0.1+offset))
  cntsDT96[i] <- length(which(arrivalTime >= 0+offset & 
                                arrivalTime <= 0.096+offset))
  cntsDT21[i,1] <- length(which(arrivalTime >= 0+offset & 
                                  arrivalTime <= 0.021+offset))
  cntsDT21[i,2] <- length(which(arrivalTime >= 0.025+offset & 
                                  arrivalTime <= 0.046+offset))
  cntsDT21[i,3] <- length(which(arrivalTime >= 0.050+offset & 
                                  arrivalTime <= 0.071+offset))
  cntsDT21[i,4] <- length(which(arrivalTime >= 0.075+offset & 
                                  arrivalTime <= 0.096+offset))
  cntsDT06[i,1] <- length(which(arrivalTime >= 0+offset & 
                                  arrivalTime <= 0.006+offset))
  cntsDT06[i,2] <- length(which(arrivalTime >= 0.010+offset & 
                                  arrivalTime <= 0.016+offset))
  cntsDT06[i,3] <- length(which(arrivalTime >= 0.020+offset & 
                                  arrivalTime <= 0.026+offset))
  cntsDT06[i,4] <- length(which(arrivalTime >= 0.030+offset & 
                                  arrivalTime <= 0.036+offset))
  cntsDT06[i,5] <- length(which(arrivalTime >= 0.040+offset & 
                                  arrivalTime <= 0.046+offset))
  cntsDT06[i,6] <- length(which(arrivalTime >= 0.050+offset & 
                                  arrivalTime <= 0.056+offset))
  cntsDT06[i,7] <- length(which(arrivalTime >= 0.060+offset & 
                                  arrivalTime <= 0.066+offset))
  cntsDT06[i,8] <- length(which(arrivalTime >= 0.070+offset & 
                                  arrivalTime <= 0.076+offset))
  cntsDT06[i,9] <- length(which(arrivalTime >= 0.080+offset & 
                                  arrivalTime <= 0.086+offset))
  cntsDT06[i,10] <- length(which(arrivalTime >= 0.090+offset & 
                                  arrivalTime <= 0.096+offset))
}

# 100 ms cycle time results
muCT <- mean(cntsCT)
sdCT <- sd(cntsCT)
cvCT <- sdCT/muCT
# 96 ms dwell time results
muDT96 <- mean(cntsDT96)
sdDT96 <- sd(cntsDT96)
cvDT96 <- sdDT96/muDT96
# 21 ms dwell time results
muDT21 <- apply(cntsDT21,2,mean)
sdDT21 <- apply(cntsDT21,2,sd)
cvDT21 <- sdDT21/muDT21
dt21Avg <- rowMeans(cntsDT21)
muDT21Avg <- mean(dt21Avg)
sdDT21Avg <- sd(dt21Avg)
cvDT21Avg <- sdDT21Avg/muDT21Avg
dt21Sum <- rowSums(cntsDT21)
muDT21Sum <- mean(dt21Sum)
sdDT21Sum <- sd(dt21Sum)
cvDT21Sum <- sdDT21Sum/muDT21Sum
# 6 ms dwell time results
muDT06 <- apply(cntsDT06,2,mean)
sdDT06 <- apply(cntsDT06,2,sd)
cvDT06 <- sdDT06/muDT06
dt06Avg <- rowMeans(cntsDT06)
muDT06Avg <- mean(dt06Avg)
sdDT06Avg <- sd(dt06Avg)
cvDT06Avg <- sdDT06Avg/muDT06Avg
dt06Sum <- rowSums(cntsDT06)
muDT06Sum <- mean(dt06Sum)
sdDT06Sum <- sd(dt06Sum)
cvDT06Sum <- sdDT06Sum/muDT06Sum
# 21 ms correlations
corDT21 <- array(0,3)
for(i in 2:4) corDT21[i-1] <- cor(cntsDT21[,1],cntsDT21[,i])
# 6 ms correlations
corDT06 <- array(0,9)
for(i in 2:10) corDT06[i-1] <- cor(cntsDT06[,1],cntsDT06[,i])

cat('\n\nMeasuring Counts',
    '\nSet Seed =',setSeed,ifelse(setSeed,seed,''),
    '\nMean Rate =',muRate,
    '\nNum Reps =',numReps,
    '\nNoise CV =',noiseCV,
    '\nCT =',CT,
    '\nCT Pop Mean =\t',muCT,'\tStdev =',sdCT,'\tCV =',cvCT,
    '\n\nIndividual Dwell Time Statistics',
    '\nDT96 Pop Mean =\t',muDT96,'\tStdev =',sdDT96,'\tCV =',cvDT96,
    '\nDT21 Pop Mean =\t',muDT21,
    '\nDT21 Pop Stdev =',sdDT21,
    '\nDT21 Pop CV =\t',cvDT21,
    '\nDT06 Pop Mean =\t',muDT06[1:5],'\n\t\t',muDT06[6:10],
    '\nDT06 Pop Stdev =',sdDT06[1:5],'\n\t\t',sdDT06[6:10],
    '\nDT06 Pop CV =\t',cvDT06[1:5],'\n\t\t',cvDT06[6:10],
    '\n\nDwell Time Average Statistics',
    '\nDT21 Average Pop Mean =',muDT21Avg,
    '\nDT21 Average Pop Stdev =',sdDT21Avg,
    '\nDT21 Average Pop CV =',cvDT21Avg,
    '\nDT06 Average Pop Mean =',muDT06Avg,
    '\nDT06 Average Pop Stdev =',sdDT06Avg,
    '\nDT06 Average Pop CV =',cvDT06Avg,
    '\n\nDwell Time Sum Statistics',
    '\nDT21 Sum Pop Mean =',muDT21Sum,
    '\nDT21 Sum Pop Stdev =',sdDT21Sum,
    '\nDT21 Sum Pop CV =',cvDT21Sum,
    '\nDT06 Sum Pop Mean =',muDT06Sum,
    '\nDT06 Sum Pop Stdev =',sdDT06Sum,
    '\nDT06 Sum Pop CV =',cvDT06Sum,
    '\n\nDwell Time Correlation Coefficients',
    '\nDT21 =\t',corDT21,
    '\nDT06 =\t',corDT06[1:5],
    '\n\t',corDT06[6:9],
    '\n')

