# marquardtSearch.R
# Author: Fred Lytle
# Last edit: 11/25/2020
marquardtSearch <- function(x,yn,FUNC,trialCoefs,tol=0.0001,maxIter=20,
                            deltaDeriv=0.0001,lambda=0.1,verbose=FALSE,
                            plot=FALSE,lty=1){
  if(verbose) cat('\nMarquard Search: Iter/Lambda/Chisqr/Coefs')
  mu <- trialCoefs
  L <- length(mu)
  beta <- array(0,L)
  alpha <- array(0,c(L,L))
  iteration <- -1
  repeat{
    delY <- yn - FUNC(x,mu)
    deriv <- getDeriv(x,FUNC,mu,deltaDeriv)
    for(r in 1:L){
      beta[r] <- sum(delY*deriv[,r])
    }
    for(r in 1:L){
      for(c in 1:L){
        alpha[r,c] <- sum(deriv[,r]*deriv[,c])
      }
    }
    alphaSave <- alpha######
    # apply lambda
    for(i in 1:L){alpha[i,i] <- abs(alpha[i,i]*(1+lambda))}
    alphaInv <- solve(alpha)
    delta <- alphaInv%*%beta
    iteration <- iteration + 1
    if(all(abs(delta) <= tol) || iteration == maxIter) break
    chiOld <- getChisqr(x,yn,FUNC,mu)
    chiNew <- getChisqr(x,yn,FUNC,mu+delta)
    if(chiNew >= chiOld){
      lambda <- lambda*10
      iteration <- iteration -1
      next
    }
    if(verbose)cat('\n\t',iteration,'\t',lambda,'\t',chiNew,'\t',mu)
    muOld <- mu
    mu <- mu + delta
    if(plot){
      if(.Devices[[2]] == 'RStudioGD' & length(trialCoefs) == 2){
        arrows(x0=muOld[1],y=muOld[2],x1=mu[1],y1=mu[2],length=0,lwd=2,lty=lty)
      }
    }
    if(chiNew < chiOld){
      lambda <- lambda/10
    }
  }
  yFit <- FUNC(x,mu)
  fitVar <- sum((yFit-yn)^2)/(length(x)-length(trialCoefs))
  alphaInv <- solve(alphaSave)
  coefErrors <- sqrt(fitVar*diag(alphaInv))
  return(list(iter=iteration-1,coefs=as.vector(mu),error=coefErrors,
              fitVar=fitVar,alphaInv=alphaInv,chi=sum((yFit-yn)^2)))
}
