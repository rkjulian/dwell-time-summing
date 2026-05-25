# EMG.R
# Defines the Exponentially Modified Gaussian (EMG) function
# Author: Fred Lytle
# Last edit: 8/30/2021
# params: 1=amplitude,2=mean,3=stdev,4=decay constant
EMG <- function(x,params){
  require('pracma')
  A <- params[1]/(2*params[4])
  EXP <- exp(0.5*(params[3]/params[4])^2 - (x - params[2])/params[4])
  ERFC <- erfc((1/sqrt(2))*((params[3]/params[4])-(x - params[2])/params[3]))
  return(A*EXP*ERFC)
}
