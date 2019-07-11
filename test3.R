rm(list=ls())
#source("generate_curves.R")
library(fda.usc)
library(parallel)
library(car)
library(pracma)
library(graphics)
library(ddalpha)

ReSample <- function(H){
  # Number of functions in J.
  dimH <- dim(S)
  kH <- dimH[1]
  kP <- dimH[2]
  r <- (kH)*runif(kH)
  # Initialize the matrix when we will store the resampled functions.
  new.H <- matrix(rep(0, len=kH*kP), nrow=kH, ncol=kP)
  for (i in 1:kH){
    # R to int to choose function in this position
    kr <- ceil(r[i])
    new.H[i, ] <- as.vector(H[kr, ])
  }
  return(new.H)
}

Bootstrapper <- function(J, G, B=1000, nc = 4){
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
  S <- rbind(J, G)
  kN <- dim(J)[1]
  kM <- dim(G)[2]
  cl <- makeCluster(nc)
  bs <- numeric(B)
  J <- rawfd2dataf(J, c(0,1))
  G <- rawfd2dataf(G, c(0,1))
  H <- rawfd2dataf(S, c(0,1))
  depths.inJ <- depthf.hM2(H, J, range=c(0,1)) #with trim = 0.25
  depths.inG <- depthf.hM2(H, G, range=c(0,1))
  depy <- depths.inJ$hM_norm2
  depx <- depths.inG$hM_norm2
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  #Need to export the global enviroment to the cluster environment
  clusterExport(cl=cl, list("ReSample","J","G", "H", "kN", "kM", "b1.ori", "b0.ori", "S"), envir=environment())
  #Need to export the libraries used to the cluster.
  clusterEvalQ(cl, library(pracma))
  clusterEvalQ(cl, library(fda.usc))
  clusterEvalQ(cl, library(ddalpha))
  #parSapply to obtain a vector.
  #print(b1.ori)
  #print(b0.ori)
  bs <- parSapply(cl, seq_len(B), function(i){
    #Resample and compute the statistic. 
    Hs <- ReSample(S)
    Js <- Hs[1:kN, ]
    aux <- kN + 1
    kH <- kN + kM
    #Next kH functions in H are the functions in the resampled version of H. 
    Gs <- Hs[aux:kH, ]
    Js <- rawfd2dataf(Js, c(0,1))
    Gs <- rawfd2dataf(Gs, c(0,1))
    Hs <- rawfd2dataf(Hs, c(0,1))
    depths.inJs <- depthf.hM2(Hs, Js, range=c(0,1)) #with trim = 0.25
    depths.inGs <- depthf.hM2(Hs, Gs, range=c(0,1))
    depxs <- depths.inGs$hM_norm2
    depys <- depths.inJs$hM_norm2
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

Tester <- function(J, G, B=1000, nc=4){
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
  stats <- Bootstrapper(J, G, B=1000, nc=nc)
  betas0 <- stats[1,]
  betas1 <- stats[2,]
  t0 <- stats[3,]
  t1 <- stats[4,]
  S <- rbind(J, G)
  J <- rawfd2dataf(J, c(0,1))
  G <- rawfd2dataf(G, c(0,1))
  H <- rawfd2dataf(S, c(0,1))
  depths.inJ <- depthf.hM2(H, J, range=c(0,1)) #with trim = 0.25
  depths.inG <- depthf.hM2(H, G, range=c(0,1))
  depx <- depths.inG$hM_norm2
  depy <- depths.inJ$hM_norm2
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  cvalb0 <- quantile(t0, probs=c(0.025,0.975))
  cvalb1 <- quantile(t1, probs=c(0.025,0.975))
  beta0.l <- b0.ori - cvalb0[2]*sqrt(var(betas0))
  beta0.u <- b0.ori - cvalb0[1]*sqrt(var(betas0))
  beta1.l <- b1.ori - cvalb1[2]*sqrt(var(betas1))
  beta1.u <- b1.ori - cvalb1[1]*sqrt(var(betas1))
  #print(b.ori)
  print('beta1')
  #print(b1.ori)
  print(beta1.l)
  print(beta1.u)
  print('beta0')
  #print(b0.ori)
  print(beta0.l)
  print(beta0.u)
  return(beta1.l <= 1 && beta1.u >= 1) #&& beta0.l <= 0 && beta0.u >= 0)
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

#cons[k] <- 100000
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
# Fill up covariances and expected values as described in the paper.
for (i in 1:kNs) {
  E1[i] <- 30*ts[i]^(3/2)*(1 - ts[i])
  #E2[i] <- 30*ts[i]*(1 - ts[i])^2
  for (j in 1:kNs) {
    cov.e[i,j] <- 0.3*exp(-3.33*(abs(ts[i] - ts[j])))
    cov.h[i,j] <- 0.3*exp(-(3^100)*(abs(ts[i] - ts[j])))
  } 
}

E2 <- E1 + 0.4
tests <- logical(25)
for (k in 1:100){
  print('iter')
  print(k)
  mu <- numeric(kNs)
  e.S0 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S1 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  J <- sweep(e.S0, 2, E1, "+")
  G <- sweep(e.S1, 2, E1, "+")
  tests[k] <- Tester(J, G, B=1000, nc = 8)
  print(tests[k])
}
1 - sum(tests)/100