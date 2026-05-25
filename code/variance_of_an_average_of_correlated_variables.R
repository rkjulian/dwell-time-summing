# Variance of an Average of Correlated Variables.R
#   User functions called: none
# Author: Fred Lytle
# Last edit: 11/08/2022
rm(list=ls())
library(MASS) # provides mvrnorm() function
set.seed(118)

# the mu and covar values for the three traces are arbitrary
N <- 1000000
muInput <- c(10,11,12)
covarInput <- matrix(data=c(9,2,1,
                            2,8,3,
                            1,3,7),nrow=3,ncol=3)

# generate the correlated data vectors, Nx3
data <- mvrnorm(n=N,mu=muInput,Sigma=covarInput)

# compute the mean and covariance the incorrect way
mu <- colMeans(data)        # three values
muCol <- mean(mu)
varCol <- apply(data,2,var) # three values
varMuCol <- sum(varCol)/3   # equals sum of matrix diagonal

# compute the mean and covariance the correct way
mu <- rowMeans(data)        # N values
muRow <- mean(mu)
varRow <- cov(data)         # 3x3 matrix
varMuRow <- sum(varRow)/3       # equals sum of the entire matrix

cat('\nInput Values',
    '\n\tMeans =',muInput,
    '\n\tInput Covar =\n')
print(covarInput)
cat('\nComputed Covar =\n')
print(varRow)
cat('\nVariance Computed from Data Columns',
    '\n\tMean =',muCol,
    '\n\tVar of Mean =',varMuCol,
    '\nVariance Computed using entire Data Matrix',
    '\n\tMean =',muRow,
    '\n\tVariance of Mean =',varMuRow)
