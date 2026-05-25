# getDeriv.R
# Author: Fred Lytle
# Last edit: 11/20/2020
getDeriv <- function(x,FUNC,params,delta){
  deriv <- array(0,c(length(x),length(params)))
  for(r in 1:length(x)){
    for(c in 1:length(params)){
      paramsPlus <- params
      paramsMinus <- params
      paramsPlus[c] <- params[c] + delta/2
      funcPlus <- FUNC(x[r],paramsPlus)
      paramsMinus[c] <- params[c] - delta/2
      funcMinus <- FUNC(x[r],paramsMinus)
      deriv[r,c] <- (funcPlus - funcMinus)/delta
    }
  }
  return(deriv)
}
