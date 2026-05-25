# getIndexSelectedFilenameChrom.R
# Author: Fred Lytle
# Last edit: 10/05/2120

getIndexSelectedFilenameChrom <- function(fileNames,params){
  # obtain parameters from file names
  L <- length(fileNames)
  conc <- array('',L)
  channels <- array('',L)
  DT <- array('',L)
  PT <- array('',L)
  replicate <- array('',L)
  transition <- array('',L)
  energy <- array('',L)
  for(i in 1:L){
    tokens <- unlist(strsplit(fileNames[i],'_'))
    conc[i] <- tokens[2]
    channels[i] <- tokens[7]
    DT[i] <- tokens[10]
    PT[i] <- tokens[12]
    replicate[i] <- tokens[13]
    transition[i] <- tokens[14]
    energy[i] <- tokens[15]
  }

  # get indices of filenames satisfying selected parameters
  indSelected <- NULL
  for(i in 1:L){
    logic1 <- conc[i] == params$Conc || params$Conc == 'any'
    logic2 <- channels[i] == params$Chan || params$Chan == 'any'
    logic3 <- DT[i] == params$Dwell || params$Dwell == 'any'
    logic4 <- PT[i] == params$Pause || params$Pause == 'any'
    logic5 <- replicate[i] == params$Replicate || params$Replicate == 'any'
    logic6 <- transition[i] == params$Trans || params$Trans == 'any'
    logic7 <- energy[i] == params$Energy || params$Energy == 'any'
    if(logic1 && logic2 && logic3 && logic4 && 
       logic5 && logic6 && logic7){
      indSelected <- c(indSelected,i)
    }
  }

  if(is.null(indSelected)) stop('No File Selected')
  return(indSelected)
}