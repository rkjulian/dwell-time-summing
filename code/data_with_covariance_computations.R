# Data with Covariance Computations.R
# Author: Fred Lytle
# Last edit: 11/16/2022
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
cat('\nData Dim =',dim(data),
    '\nTrue Means =',muInput,
    '\nTrue Var/Cov Matrix\n')
print(covarInput)
cat('\nExp Var/Cov Matrix\n')
print(cov(data))

# process sums by columns (incorrect - misses covariance)
muCol <- colMeans(data)# note: sum(muCol) = muRow
varCol <- apply(data,2,var)# sum(varCol) gives the diagonal sum
cat('\nProcess Sum by Column (incorrect - misses covariance)',
    '\n\tMean =',sum(muCol),
    '\n\tVar =',varCol,
    '\n\tMatrix Diag Sum =',sum(diag(cov(data))))

# process by rows (correct - includes covariance)
sumRow <- rowSums(data)
muRow <- mean(sumRow)
varRow <- var(sumRow)# gives the matrix sum
cat('\nProcess Sum by Row (correct - includes covariance)',
    '\n\tMean =',muRow,
    '\n\tVar =',sum(varRow),
    '\n\tMatrix Full Sum =',sum(cov(data)))

# process averages by rows looking for 1/3^2 term in Eq. 4
avgRow <- rowMeans(data)
varAvgRow <- var(avgRow)
cat('\nProcess Avg by Row (correct 1/3^2 sum result)',
    '\n\tVar =',varAvgRow,
    '\n\tSumVar/3^2 =',varRow/9)

