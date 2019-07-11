rm(list=ls())
#source("generate_curves.R")
library(fda.usc)
library(parallel)
library(car)
library(pracma)


Bootstrapper <- function(J, G, B=1000, depth.function=depth.FM, nc = 4){
  # Computes the bootstrap statistics parallely.
  #
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # B: The amount of bootstrap statistics we're going to use.
  # stat: The statistic we are going to use. You can choose P1-P4.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  # nc: Number of clusters. In my PC I have 4 so this is the default value I
  # will use.
  #
  # Returns:
  # J.fun and G.fun: A list containing this two structures, which are
  # the resampled versions of J and G.
  #
  # Number of functions in J.
  H <- c(J, G)
  kN <- length(J)
  kM <- length(G)
  kH <- kN + kM
  cl <- makeCluster(nc)
  bs <- numeric(B)
  depths.inJ <- depth.FM(H, fdataori=G, trim = 0) #with trim = 0.25
  depths.inG <- depth.FM(H, fdataori=J, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  plot(depy, depx)
  lm.ori <- lm(depy ~ depx)
  summary(lm.ori)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  #Need to export the global enviroment to the cluster environment
  clusterExport(cl=cl, list("J", "G", "nc", "depth.function", "H", "kN", "kM", "kH", "b1.ori", "b0.ori"), envir=environment())
  #Need to export the libraries used to the cluster.
  clusterEvalQ(cl, library(pracma))
  clusterEvalQ(cl, library(fda.usc))
  #parSapply to obtain a vector.
  bs <- parSapply(cl, seq_len(B), function(i){
    #Resample and compute the statistic. 
    Hs <- H[sample(1:kH,size=kH,replace=TRUE),]
    Js <- Hs[1:kN, ]
    aux <- kN + 1
    #Next kH functions in H are the functions in the resampled version of H. 
    Gs <- Hs[aux:kH, ]
    depths.inJs <- depth.FM(Hs, fdataori=Gs, trim = 0) #with trim = 0.25
    depths.inGs <- depth.FM(Hs, fdataori=Js, trim = 0)
    depys <- depths.inJs$dep
    depxs <- depths.inGs$dep
    lm.bs <- lm(depys ~ depxs)
    beta1 <- lm.bs$coefficients[2]
    beta0 <- lm.bs$coefficients[1]
    sum_aux <- sum((depxs - mean(depxs))^2)
    beta1.std <- sqrt((sum(lm.bs$residuals^2))/((kH-2)*sum_aux))
    beta0.std <- sqrt((sum(lm.bs$residuals^2)*sum(depxs^2))/((kH-2)*sum_aux))
    t1 <- (beta1 - 1)/beta1.std
    t0 <- (beta0)/beta0.std
    x <- c(beta0, beta1, t0, t1)
  }
  )
  stopCluster(cl)
  return(bs)
}

Tester <- function(J, G, B=1000, depth.function=depth.FM, nc=4){
  # Tests homogeinity using bootstrap confidence intervals. H0: Homogeinity.
  #
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # B: The amount of bootstrap statistics we're going to use.
  # stat: The statistic we are going to use. You can choose P1-P4.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  # nc: Number of clusters. In my PC I have 4 so this is the default value I
  # will use.
  #
  # Returns:
  # TRUE if we do not reject the H0, meaning that J and G come from the same
  # population. False is we reject the H0, meaning that J and G are not 
  # homogeneous and come from different populations.
  #
  # True P4 between J and G.
  # Compute the bootstrap interval
  H <- c(J, G)
  depths.inJ <- depth.FM(H, fdataori=G, trim = 0) #with trim = 0.25
  depths.inG <- depth.FM(H, fdataori=J, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  lm.ori <- lm(depy ~ depx)
  plot(depx, depy)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  stats <- Bootstrapper(J, G, B=1000, depth.function=depth.FM, nc=nc)
  b0 <- stats[1,]
  b1 <- stats[2,]
  t0 <- stats[3,]
  t1 <- stats[4,]
  cvalb0 <- quantile(t0, probs=c(0.025, 0.975))
  cvalb1 <- quantile(t1, probs=c(0.025, 0.975))
  beta0.l <- b0.ori - cvalb0[2]*sqrt(var(b0))
  beta0.u <- b0.ori - cvalb0[1]*sqrt(var(b0))
  beta1.l <- b1.ori - cvalb1[2]*sqrt(var(b1))
  beta1.u <- b1.ori - cvalb1[1]*sqrt(var(b1))
  print('beta0')
  print(beta0.l)
  print(beta0.u)
  print('beta1')
  print(beta1.l)
  print(beta1.u)
  return(beta1.l <= 1 && beta1.u >= 1) #&& beta0.l <= 0 && beta0.u >= 0)
}

Tester3 <- function(J, G, B=1000, depth.function=depth.FM, nc=4){
  # Tests homogeinity using bootstrap confidence intervals. H0: Homogeinity.
  #
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # B: The amount of bootstrap statistics we're going to use.
  # stat: The statistic we are going to use. You can choose P1-P4.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  # nc: Number of clusters. In my PC I have 4 so this is the default value I
  # will use.
  #
  # Returns:
  # TRUE if we do not reject the H0, meaning that J and G come from the same
  # population. False is we reject the H0, meaning that J and G are not 
  # homogeneous and come from different populations.
  #
  # True P4 between J and G.
  # Compute the bootstrap interval
  H <- c(J, G)
  depths.inJ <- depth.FM(H, fdataori=G, trim = 0) #with trim = 0.25
  depths.inG <- depth.FM(H, fdataori=J, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  stats <- Bootstrapper(J, G, B=1000, depth.function=depth.FM, nc=nc)
  b0 <- stats[1,]
  b1 <- stats[2,]
  t0 <- stats[3,]
  t1 <- stats[4,]
  cval0 <- quantile(t0, probs=c(0.025, 0.975))
  cval1 <- quantile(t1, probs=c(0.025, 0.975))
  beta0.l <- cval0[1]
  beta0.u <- cval0[2]
  beta1.l <- cval1[1]
  beta1.u <- cval1[2]
  print('beta0')
  print(beta0.l)
  print(beta0.u)
  print('beta1')
  print(beta1.l)
  print(beta1.u)
  return(beta1.l <= 1 && beta1.u >= 1 && beta0.l <= 0 && beta0.u >= 0)
}


Tester2 <- function(J, G, B=1000, depth.function=depth.FM, nc=4){
  first_test <- Tester(J, G, nc) #Test one way
  if (first_test){
    other_test <- Tester(G, J, nc)
    return(other_test) 
  } else{
    return(FALSE)
  }
}


kNs <- 30
# Number of curves.
kNc <- 50
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

mu <- numeric(kNs)
# Sample 0 error simulation with mean 0 and covariance as in cove.
tests <- logical(100)
for (k in 1:100){
  print('iter')
  print(k)
  e.S0 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e) 
  # Sample 1 to 5 error, first three with error e(t), last two with error h(t).
  e.S1 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.h)
  S0 <- sweep(e.S0, 2, E1, "+")
  S1 <- sweep(e.S1, 2, E1, "+")
  J <- fdata(S0)
  G <- fdata(S1)
  H <- c(J, G)
  tests[k] <- Tester(J, G, nc=8)
  print('result')
  print(tests[k])
}
1 - sum(tests)/100
