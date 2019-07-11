# 
# Author: Alejandro Calle Saldarriaga. 22-01-2019.
# 
# This file has the GenerateCurves function. The purpose of this file is to
# generate the simulations of the curves described in Flores et al. [2018]
# which will be used in other files to test the size and power of our
# homogeinity test.

library(pracma)
library(MASS)

GenerateCurves <- function(){
  # Generate the curves described in Flores et al. [2018] for homogeinity tests.
  #
  # Args:
  # none
  # 
  # Returns:
  #
  # S, a list with each of the samples, where the first one is the reference 
  # sample (S0), the next 5 ones are the other samples, and the last one
  # is another realization of the reference sample.
  #
  # Number of observations for each curve.
  kNs = 30
  # Number of curves.
  kNc = 50
  # 30 equidistant points (ts = timesteps).
  ts <- linspace(0, 1, n = kNs) 
  # Initialize covariance matrix for first e(t) with zeroes.
  cov.e <- matrix(rep(0, len=kNs*kNs), nrow=kNs, ncol=kNs)
  # Initializ covariance matrix for h(t) with zeroes.
  cov.h <- matrix(rep(0, len=kNs*kNs), nrow=kNs, ncol=kNs) 
  E1 <- numeric(kNs) 
  E2 <- numeric(kNs)
  # Fill up covariances and expected values as described in the paper.
  for (i in 1:kNs) {
    E1[i] <- 30*ts[i]^(3/2)*(1 - ts[i])
    E2[i] <- 30*ts[i]*(1 - ts[i])^2
    for (j in 1:kNs) {
      cov.e[i,j] <- 0.3*exp(-(abs(ts[i] - ts[j]))/0.3)
      cov.h[i,j] <- 0.5*exp(-(abs(ts[i] - ts[j]))/0.2)
      } 
  }
  E3 <- E1 + 1
  E4 <- E1 + 0.5
  # Vector of zeros, all gaussians simualations will be centered around 0 so
  # we need this vector for later.
  mu <- numeric(kNs)
  # Sample 0 error simulation with mean 0 and covariance as in cove.
  e.S0 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e) 
  # Sample 1 to 5 error, first three with error e(t), last two with error h(t).
  e.S1 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S2 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S3 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S4 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.h)
  e.S5 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.h)
  
  # This is just another realization of the first population,
  # needed for computing the size of the test
  e.S6 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  
  # Now sweep to sum the errors and the means in the correct manner in order to
  # get 50x30 matrices where each row is a function.
  S0 <- sweep(e.S0, 2, E1, "+")
  S1 <- sweep(e.S1, 2, E3, "+")
  S2 <- sweep(e.S2, 2, E4, "+")
  S3 <- sweep(e.S3, 2, E2, "+")
  S4 <- sweep(e.S4, 2, E2, "+")
  S5 <- sweep(e.S5, 2, E1, "+")
  
  S6 <- sweep(e.S6, 2, E1, "+")
  
  S <- list(S0, S1, S2, S3, S4, S5, S6)
  
  return(S)
  
}

# References
#
# Flores, Ram?n, Lillo, Rosa and Romo, Juan. 2018. Homogeinity test for
# functional data. Journal of Applied Statistics, 45(5), 868-883.
