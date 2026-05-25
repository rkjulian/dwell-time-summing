# getChisqr.R
# Author: Fred Lytle
# Last edit: 11/19/2020
getChisqr <- function(x,yn,FUNC,params){
  chisqr <- sum((yn - FUNC(x,params))^2)
}
