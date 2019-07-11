# 
# Author: Alejandro Calle Saldarriaga. 27-06-2019.
# 
# This file implements an original homogeinity test for functional data, based on
# the ideas in [Liu et al., 1999]. 
# The function here are:
# Bootstrapper, which implements a resampling scheme for our statistic.
# Tester, which uses the boostrapped statistics to reject or not reject the test.

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
  # Note: Not currently using beta0 for anything, but some versions of this code
  # used it. If you want to experiment with beta1 and beta0 simulataneously
  # its quite easy as the t-values for beta0 are computed here too.
  #
  H <- c(J, G)
  kN <- length(J)
  kM <- length(G)
  kH <- kN + kM
  cl <- makeCluster(nc)
  bs <- numeric(B)
  depths.inJ <- depth.function(H, fdataori=J, trim = 0) #with trim = 0.25
  depths.inG <- depth.function(H, fdataori=G, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  #Need to export the global enviroment to the cluster environment
  clusterExport(cl=cl, list("J", "G", "nc", "depth.function", "H", "kN",
                            "kM", "kH", "b1.ori", "b0.ori"), envir=environment())
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
    #Resampled depths
    depths.inJs <- depth.function(Hs, fdataori=Js, trim = 0) 
    depths.inGs <- depth.function(Hs, fdataori=Gs, trim = 0)
    depys <- depths.inJs$dep
    depxs <- depths.inGs$dep
    lm.bs <- lm(depys ~ depxs)
    # Betas
    beta1 <- lm.bs$coefficients[2]
    beta0 <- lm.bs$coefficients[1]
    sum_aux <- sum((depxs - mean(depxs))^2)
    # Formula for the standard deviation of both betas. Standard formula from
    # regression analysis.
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
  H <- c(J, G)
  depths.inJ <- depth.function(H, fdataori=J, trim = 0) #with trim = 0.25
  depths.inG <- depth.function(H, fdataori=G, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  stats <- Bootstrapper(J, G, B=1000, depth.function=depth.function, nc=nc)
  # Extract all statistics from the bootstrapper.
  b0 <- stats[1,]
  b1 <- stats[2,]
  t0 <- stats[3,]
  t1 <- stats[4,]
  # 5% cutoff (2.5% upper and lower).
  cvalb0 <- quantile(t0, probs=c(0.025, 0.975))
  cvalb1 <- quantile(t1, probs=c(0.025, 0.975))
  beta0.l <- b0.ori - cvalb0[2]*sqrt(var(b0))
  beta0.u <- b0.ori - cvalb0[1]*sqrt(var(b0))
  beta1.l <- b1.ori - cvalb1[2]*sqrt(var(b1))
  beta1.u <- b1.ori - cvalb1[1]*sqrt(var(b1))
  return(beta1.l <= 1 && beta1.u >= 1) 
}

