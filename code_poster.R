# modele simple
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

don <- read_xlsx("/Users/negarsoleimani/Documents/phd/paper1/Ball_drops_data.xlsx", sheet = 2)
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop)
don <- don[don$drop == 1, ]

t <- don$time
y <- don$Height
length(t)
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)

##############################################################################################

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

##############################################################################################

GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
}

##############################################################################################

mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000) {
  # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T")
  
  # Total iterations = burn-in + desired samples
  total_iter <- n_burnin + n_iter
  
  theta <- init
  delta <- rep(0, length(y))
  chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
  colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
  chain_zeta <- matrix(NA, nrow = total_iter, ncol = length(y))
  
  loglik_chain <- numeric(total_iter)
  accept_psi <- 0
  
  for (iter in 1:total_iter) {
    g <- theta[1]; h0 <- theta[2]; sigma_sq_err <- theta[3]
    alpha_param <- theta[4]; psi_delta <- theta[5]; k <- theta[6]
    sigma_sq_delta <- sigma_sq_err / k
    
    Sigma_delta <- GP_covariance(t, sigma_sq_delta, psi_delta)
    f_theta <- balldropg(t, c(g, h0))
    mean1 <- f_theta
    mean2 <- f_theta + delta

#-------------------------------- probability of the zeta --------------------------------# 
    
    # Method1 for calculate the probability of the zeta:
    # prob_zeta <- 1 / (1 + (alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err))) /
    #                     ((1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))) #P(M2) = ((1-alpha)*M2)/((alpha*M1)+(1-alpha)*M2)
    # zeta <- 1 + (runif(length(y)) < prob_zeta)
    # chain_zeta[iter, ] = zeta
    
    #I don't use Method1 for computing prob_zeta,because this direct probability 
    #computation is numerically unstable, when the parameter sigma_sq_err is small 
    #(sigma_sq_err = 0.01), it leads to overflow and produces NaN/NA probabilities, 
    #and these NA values then propagate into zeta and later cause singular matrices 
    #and when computing solve(A), leading to errors.
    
    # Method2 for calculate the probability of the zeta:(log_sum_exp: https://rpubs.com/FJRubio/LSE and https://en.wikipedia.org/wiki/LogSumExp)
    log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
    log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    #prob_zeta_1 <- exp(log_w1 - log_den)  #  P(zeta=1|y)
    prob_zeta_2 <- exp(log_w2 - log_den)  #  P(zeta=2|y)
    #zeta <- ifelse(runif(length(y)) < prob_zeta_1, 1, 2)
    zeta <- ifelse(runif(length(y)) < prob_zeta_2, 2, 1) #= zeta <- 1 + (runif(length(y)) < prob_zeta_1) # sample zeta: 2 with prob post_p2, otherwise 1
    chain_zeta[iter, ] = zeta

#-------------------------------- log_likelihood --------------------------------#   
    
    log_likelihood <- sum(log(ifelse(zeta == 1,
                                     alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
                                     (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    loglik_chain[iter] <- log_likelihood

#-------------------------------- delta(discrepancy) --------------------------------#  
    
    zeta_2_indices <- which(zeta == 2)
    if (length(zeta_2_indices) > 0) {
      y_m <- y[zeta_2_indices]
      Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
        Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
      Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
      Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
      mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
      Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
      Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
      delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    }
    else delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
    
#-------------------------------- Gibbs step for theta --------------------------------# 

    # prior Normal distribution
    # zeta_1_indices <- which(zeta == 1)
    # zeta_2_indices <- which(zeta == 2)
    # X <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
    # x1 <- X[zeta_1_indices, , drop = FALSE]
    # x2 <- X[zeta_2_indices, , drop = FALSE]
    # theta_hat     <- matrix(c(46.46, 9.8), ncol = 1)
    # inv_sigma_theta <- solve(Sigma_theta)
    # A <- ((t(x1) %*% x1) / sigma_sq_err) + ((t(x2) %*% x2) / sigma_sq_err) + inv_sigma_theta
    # Sigmapost_theta <- solve(A)
    # y1 <- matrix(y[zeta_1_indices], ncol = 1)
    # y2 <- matrix(y[zeta_2_indices], ncol = 1)
    # d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    # B <- (t(x1) %*% y1) / sigma_sq_err +
    #   (t(x2) %*% y2) / sigma_sq_err -
    #   (t(x2) %*% d2) / sigma_sq_err +
    #   inv_sigma_theta %*% theta_hat
    # Mupost_theta <- Sigmapost_theta %*% B
    # theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
    # h0 <- theta_sample[1];  g <- theta_sample[2]
    # theta[1] <- g
    # theta[2] <- h0

    #-------------------------------- Page 22 - part 8.1 --------------------------------#     
    # flat prior on theta 
    zeta_1_indices <- which(zeta == 1)
    zeta_2_indices <- which(zeta == 2)
    
    X <- cbind(1, -0.5 * (t * t_range)^2)
    x1 <- X[zeta_1_indices, , drop = FALSE]
    x2 <- X[zeta_2_indices, , drop = FALSE]
    
    y1 <- matrix(y[zeta_1_indices], ncol = 1)
    y2 <- matrix(y[zeta_2_indices], ncol = 1)
    d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    
    A <- (t(x1) %*% x1 + t(x2) %*% x2) / sigma_sq_err
    B <- (t(x1) %*% y1 + t(x2) %*% y2 - t(x2) %*% d2) / sigma_sq_err
    
    Sigmapost_theta <- solve(A)             
    Mupost_theta    <- Sigmapost_theta %*% B 
    
    theta_sample <- rmvnorm(1, mean = Mupost_theta,
                                     sigma = Sigmapost_theta)
    
    h0 <- theta_sample[1];  g <- theta_sample[2]
    theta[1] <- g
    theta[2] <- h0
    
    
    if(mcmc_parameters[1] == FALSE){
      g <- init[1]
      h0 <- init[2]
    }
    
    ## g = fixer
    # if (mcmc_parameters[1] == FALSE) {
    #   g <- init[1]
    #   h0 <- theta_sample[1]
    # }
    
    ## h0 = fixed
    #if (mcmc_parameters[1] == FALSE) {
    #  h0 <- init[2] 
    #  g <- theta_sample[2] 
    #}

#-------------------------------- Gibbs step for sigma_sq_err(lambda^2) --------------------------------#   
#-------------------------------- Page 23 - part 8.2 --------------------------------#     
    
    R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
    if(n > 0){
      R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
      
      quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
    } 
    else {
      quad_form_delta <- 0
    }
    
    f_theta <- balldropg(t, c(g, h0))
    idx1 <- which(zeta == 1)
    idx2 <- which(zeta == 2)
    residual1 <- y[idx1] - f_theta[idx1]
    rss1   <- sum(residual1^2)

    residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
    rss2   <- sum(residual2^2)

    # When we consider the Jeffreys prior, I use the following code:
    rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    shape_err <- n + 1/2
    
    # When the prior is the sigma parameter of the "inverse gamma distribution", I use the following code:
    # rate_err <- 1 + (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    # shape_err <- 2 + n
    # sigma_sq_err <- rinvgamma(1, shape = shape_err, scale = rate_err)
    sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    theta[3] <- sigma_sq_err
    # 
    if(mcmc_parameters[2] == FALSE){
      sigma_sq_err <- init[3]
      theta[3] <- sigma_sq_err
    }
 
#-------------------------------- Gibbs step for psi_delta --------------------------------#   
    
    # When the prior is the psi_delta of the "uniform(0.1, 1)", I use the following code:
    psi_prop <- rtruncnorm(1, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    sigma_sq_delta_prop <- sigma_sq_err / k
    Sigma_delta <- GP_covariance(t, sigma_sq_delta_prop, psi_delta)
    Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
    log_prop_current <- log(dtruncnorm(psi_delta, a = 0.1, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
    log_prop_prop <- log(dtruncnorm(psi_prop, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
    log_prior_current <- dunif(psi_delta, min = 0.1, max = 1, log = TRUE)
    log_prior_prop    <- dunif(psi_prop,  min = 0.1, max = 1, log = TRUE)
    log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
    log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + (log_prop_current - log_prop_prop)
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      psi_delta <- psi_prop
      theta[5] <- psi_delta
      accept_psi <- accept_psi + 1
    }

    # When the prior is the psi_delta of the "Beta(8, 4)", I use the following code:
    
    # psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    # sigma_sq_delta_prop <- sigma_sq_err / k
    # Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
    # #log_prior_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = 0.5, sd = 0.2))
    # #log_prior_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = 0.5, sd = 0.2))
    # log_prior_current <- dbeta(psi_delta, shape1 = 8, shape2 = 4, log=TRUE)
    # log_prior_prop <- dbeta(psi_prop, shape1 = 8, shape2 = 4, log=TRUE)
    # log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
    # log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    # log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current)
    # if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
    #   psi_delta <- psi_prop
    #   accept_psi <- accept_psi + 1
    # }

    if(mcmc_parameters[3] == FALSE){
      psi_delta <- init[5]
      theta[5] <- psi_delta
    }

#-------------------------------- Gibbs step for k --------------------------------#   
#-------------------------------- Page 24-part 8.3 --------------------------------#     
    
    alpha_k <- (n / 2) + 1
    beta_k <- (1 / (2 * sigma_sq_err)) * quad_form_delta
    for (try_k in 1:100) {
      k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
      #if (k_prop >= 0.1 && k_prop <= 0.9) {
      if (k_prop > 0 && k_prop < 1) {
        k <- k_prop
        theta[6] <- k
        break
      }
    }
    
    if(mcmc_parameters[4] == FALSE){
      k <- init[6]
      theta[6] <- k
    }
 
#-------------------------------- Gibbs step for alpha --------------------------------#   
#-------------------------------- Page 25 - part 8.4 --------------------------------#     
    
    alpha_param <- rbeta(1, sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
    theta[4] <- alpha_param
    
    if(mcmc_parameters[5] == FALSE){
      alpha_param <- init[4]
      theta[4] <- alpha_param
    }
    
    
    theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
    chain_theta[iter, ] <- theta
    chain_delta[iter, ] <- delta
    
  }
 
#----------------------------------------------------------------#   
  
  # Remove burn-in samples if n_burnin > 0
  if (n_burnin > 0) {
    keep_indices <- (n_burnin + 1):total_iter
    chain_theta <- chain_theta[keep_indices, , drop = FALSE]
    chain_delta <- chain_delta[keep_indices, , drop = FALSE]
    chain_zeta <- chain_zeta[keep_indices, , drop = FALSE]
    loglik_chain <- loglik_chain[keep_indices]
    
    # Adjust acceptance rate calculation for post-burn-in period only
    cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
  } else {
    cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
  }
  return(list(theta = chain_theta, delta = chain_delta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_psi / n_iter))
}
##############################################################
##############################################################
###################### Real Data #############################
##############################################################
##############################################################
set.seed(12345)

init_base       <- c(9.8, 46.46, 0.01, 0.5, 0.3, 0.2)
sigma_proposals <- c(NA,NA,NA,NA,0.4,NA)
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta     <- matrix(c(0.5,0,0,0.5), 2)
n_iter          <- 10000
burn_in         <- 5000

results_real_sh5 <- mcmc_step6(
  y = y, t = t,
  n_iter = n_iter,
  init   = init_base,
  sigma_proposals = sigma_proposals,
  mcmc_parameters = mcmc_parameters,
  Sigma_theta     = Sigma_theta,
  n_burnin        = burn_in
)

#save(results_real_sh5,file = "/Users/negarsoleimani/Documents/PhD/Paper1/Real Data/results_real_sh5.RData")
# load("/Users/negarsoleimani/Documents/PhD/Paper1/Real Data/results_real_sh5.RData")
# load("/Users/negarsoleimani/Documents/PhD/Paper1/Real Data/results_real_sh5.RData")
# load("/Users/negarsoleimani/Documents/PhD/Paper1/Real Data/results_real_sh5.RData")
# load("/Users/negarsoleimani/Documents/PhD/Paper1/Real Data/results_real_sh5.RData")
# load("/Users/negarsoleimani/Documents/PhD/Paper1/Real Data/results_real_sh5.RData")

g_chain <- results_real_sh5$theta[,1]
h0_chain <- results_real_sh5$theta[,2]
sigma_sq_err_chain <- results_real_sh5$theta[,3]
alpha_chain <- results_real_sh5$theta[,4]
psi_delta_chain <- results_real_sh5$theta[,5]
k_chain <- results_real_sh5$theta[,6]
delta_chain <- results_real_sh5$delta
zeta_chain <- results_real_sh5$zeta
loglik_chain <- results_real_sh5$loglik
accept_rate_psi <- results_real_sh5$accept_rate_psi

par(mfrow = c(1, 1))
boxplot(results_real_sh5$delta, col = "#CCDAFF")
title(ylab = expression(delta), xlab = "n", line = 1.75)

prob_zeta_model1 <- colMeans(results_real_sh5$zeta == 1)
print(prob_zeta_model1)

prob_zeta_model2 <- colMeans(results_real_sh5$zeta == 2)
print(prob_zeta_model2)

par(mfrow = c(3, 2))
plot(results_real_sh5$theta[, 1], type = "l", main = "Trace of g")
plot(results_real_sh5$theta[, 2], type = "l", main = "Trace of h0")
plot(results_real_sh5$theta[, 3], type = "l", ylim = c(0, 1), main = "Trace of sigma_sq_err")
plot(results_real_sh5$theta[, 4], type = "l", ylim = c(0, 1), main = "alpha")
plot(results_real_sh5$theta[, 5], type = "l", main = "Trace of psi_delta")
plot(results_real_sh5$theta[, 6], type = "l", ylim = c(0, 1), main = "k")

par(mfrow = c(2, 3))
boxplot(results_real_sh5$theta[, 1])
boxplot(results_real_sh5$theta[, 2])
abline(h = 46.46)
boxplot(results_real_sh5$theta[, 3])
boxplot(results_real_sh5$theta[, 4])
boxplot(results_real_sh5$theta[, 5])
boxplot(results_real_sh5$theta[, 6])

library("vioplot")
par(mfrow = c(1, 6))
vioplot(results_real_sh5$theta[, 1], col = "#C9E2FF", xlab = "", xaxt = "n")
title(ylab = expression(g), line = 1.75)
abline(h = 9.8)
vioplot(results_real_sh5$theta[, 2], col = "#C9E2FF", xlab = "", xaxt = "n")
title(ylab = expression(h[0]), line = 1.75)
abline(h = 46.46)
vioplot(results_real_sh5$theta[, 3], col = "#C9E2FF", xlab = "", xaxt = "n")
title(ylab = expression(sigma[err]^2), line = 1.75)
abline(h = 0.01)
vioplot(results_real_sh5$theta[, 4], col = "#C9E2FF", xlab = "", xaxt = "n")
title(ylab = expression(alpha), line = 1.75)
vioplot(results_real_sh5$theta[, 5], col = "#C9E2FF", xlab = "", xaxt = "n")
title(ylab = expression(psi[delta]), line = 1.75)
vioplot(results_real_sh5$theta[, 6], col = "#C9E2FF", xlab = "", xaxt = "n")
title(ylab = expression(k), line = 1.75)
hist(results_real_sh5$theta[, 4])

##############################################################
##############################################################
##############################################################
# Simulation the data of model 2 with different psi_delta
# Code to verify the result

set.seed(12345)
k = 0.2
sim_psi_delta <- 0.9
sigma_sq_err <- (0.1)^2
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 5
n_iter <- 10000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA)
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.1, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)

  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))

  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  #y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)*(1+1/k))

  y_obs[, v] <- y_1
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)

  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
res <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

g <- res[[1]]
h <- res[[2]]
sig <- res[[3]]
alpha <- res[[4]]
psi <- res[[5]]
k <- res[[6]]

hist(k)
hist(sig)
hist(psi)
hist(alpha)

boxplot(colMeans(g))
abline(h = 9.8)
boxplot(colMeans(h))
abline(h = 46.46)
boxplot(colMeans(sig))
abline(h = 0.01)
boxplot(colMeans(alpha))
boxplot(colMeans(psi))
abline(h = 0.7)
boxplot(colMeans(k))
abline(h = 0.2)
# # 
# plot(y_1)
# zeta <- res[[8]]
# prob_zeta_m2 <- colMeans(res[[8]][[1]] == 2)
# delta <- res[[7]][[1]]

##############################################################
##############################################################
# Simulation the data of model 2 with different psi_delta
# ## Psi1 ####################################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1
sim_psi_delta <- 0.01
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000
# 
sigma_props <- c(NA, NA, NA, NA, 0.02, NA)
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, sim_psi_delta, 0.2)
# 
g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)
# 
y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)

  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))

  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)

  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi1_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

# result_m2_sh2_psi1_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
# save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")


## Psi2 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.1 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi2_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi2_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi2_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi2_simple.RData")
## Psi3 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.2
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi3_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi3_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi3_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi3_simple.RData")


## Psi4 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.3 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi4_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi4_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi4_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi4_simple.RData")


## Psi5 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.4
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi5_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi5_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi5_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi5_simple.RData")


## Psi6 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.5 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi6_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi6_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi6_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi6_simple.RData")


## Psi7 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.6 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi7_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi7_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi7_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi7_simple.RData")


## Psi8 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.7 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi8_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi8_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi8_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi8_simple.RData")


## Psi9 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.8 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
  
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi9_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi9_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi9_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi9_simple.RData")


## Psi10 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1
sim_psi_delta <- 0.9
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000
# 
sigma_props <- c(NA, NA, NA, NA, 0.4, NA)
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.1, 0.5, 0.2, 0.2)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)
# 
y_obs <- matrix(NA, n, n_samples)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)

  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))

  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1

  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in)
#   
  g_chain[, v]     <- res$theta[, 1]
  h0_chain[, v]    <- res$theta[, 2]
  sigma_chain[, v] <- res$theta[, 3]
  alpha_chain[, v] <- res$theta[, 4]
  psi_chain[, v]   <- res$theta[, 5]
  k_chain[, v]     <- res$theta[, 6]
  delta_list[[v]]  <- res$delta
  zeta_list[[v]]   <- res$zeta
  loglik_mat[, v]  <- res$loglik
  accept_rate[v]   <- res$accept_rate_psi
}
result_m2_sh2_psi10_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi10_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi10_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi10_simple.RData")
# #####################################################################################################

load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi1_simple.RData")
#load("/home/nsoleimani/my_project/result_m2_sh2_psi1_simple.RData")

load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi2_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi3_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi4_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi5_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi6_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi7_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi8_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi9_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi10_simple.RData")

boxplot(result_m2_sh2_psi8_simple1[[7]][[1]])

g1 <- result_m2_sh2_psi1_simple1[[1]]
g2 <- result_m2_sh2_psi2_simple1[[1]]
g3 <- result_m2_sh2_psi3_simple1[[1]]
g4 <- result_m2_sh2_psi4_simple1[[1]]
g5 <- result_m2_sh2_psi5_simple1[[1]]
g6 <- result_m2_sh2_psi6_simple1[[1]]
g7 <- result_m2_sh2_psi7_simple1[[1]]
g8 <- result_m2_sh2_psi8_simple1[[1]]
g9 <- result_m2_sh2_psi9_simple1[[1]]
g10 <- result_m2_sh2_psi10_simple1[[1]]

par(mfrow = c(2,3))
boxplot(colMeans(g1), colMeans(g2), colMeans(g3), colMeans(g4), colMeans(g5), 
        colMeans(g6), colMeans(g7), colMeans(g8), colMeans(g9), colMeans(g10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "g",
        col = "#800020")
abline(h=9.8, col = "orange")

par(mgp = c(3, 0.7, 0)) 
library(vioplot)

vioplot(colMeans(g1), colMeans(g2), colMeans(g3), colMeans(g4), colMeans(g5),
        colMeans(g6), colMeans(g7), colMeans(g8), colMeans(g9), colMeans(g10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "",  
        ylab = "",
        col = "#E18E96")
title(xlab = expression(psi[delta]), line = 1.5)
title(ylab = expression(g), line = 1.5)
abline(h = 9.8, col = "#800020")




h01 <- result_m2_sh2_psi1_simple1[[2]]
h02 <- result_m2_sh2_psi2_simple1[[2]]
h03 <- result_m2_sh2_psi3_simple1[[2]]
h04 <- result_m2_sh2_psi4_simple1[[2]]
h05 <- result_m2_sh2_psi5_simple1[[2]]
h06 <- result_m2_sh2_psi6_simple1[[2]]
h07 <- result_m2_sh2_psi7_simple1[[2]]
h08 <- result_m2_sh2_psi8_simple1[[2]]
h09 <- result_m2_sh2_psi9_simple1[[2]]
h010 <- result_m2_sh2_psi10_simple1[[2]]

boxplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05), 
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(psi_delta),
        ylab = "h0",
        col = "#800020")
abline(h = 46.46, col = "orange")

vioplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05), 
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(psi[delta]), line = 1.5)
title(ylab = expression(h), line = 1.5)
abline(h = 46.46, col = "#800020")

sigma1 <- result_m2_sh2_psi1_simple1[[3]]
sigma2 <- result_m2_sh2_psi2_simple1[[3]]
sigma3 <- result_m2_sh2_psi3_simple1[[3]]
sigma4 <- result_m2_sh2_psi4_simple1[[3]]
sigma5 <- result_m2_sh2_psi5_simple1[[3]]
sigma6 <- result_m2_sh2_psi6_simple1[[3]]
sigma7 <- result_m2_sh2_psi7_simple1[[3]]
sigma8 <- result_m2_sh2_psi8_simple1[[3]]
sigma9 <- result_m2_sh2_psi9_simple1[[3]]
sigma10 <- result_m2_sh2_psi10_simple1[[3]]

boxplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "sigma",
        col = "#800020")
abline(h = 0.01, col = "orange")

vioplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(psi[delta]), line = 1.5)
title(ylab = expression(sigma[err]^2), line = 1.5)
abline(h = 0.01, col = "#800020")

alpha1 <- result_m2_sh2_psi1_simple1[[4]]
alpha2 <- result_m2_sh2_psi2_simple1[[4]]
alpha3 <- result_m2_sh2_psi3_simple1[[4]]
alpha4 <- result_m2_sh2_psi4_simple1[[4]]
alpha5 <- result_m2_sh2_psi5_simple1[[4]]
alpha6 <- result_m2_sh2_psi6_simple1[[4]]
alpha7 <- result_m2_sh2_psi7_simple1[[4]]
alpha8 <- result_m2_sh2_psi8_simple1[[4]]
alpha9 <- result_m2_sh2_psi9_simple1[[4]]
alpha10 <- result_m2_sh2_psi10_simple1[[4]]

boxplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8", "0.9"),
        xlab = "psi_delta",
        ylab = "alpha",
        ylim = c(0,1),
        col = "#800020")

vioplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8", "0.9"),
        xlab = "",
        ylab = "",
        ylim = c(0,1),
        col = "#E18E96")
title(xlab = expression(psi[delta]), line = 1.5)
title(ylab = expression(alpha), line = 1.5)

psi1 <- result_m2_sh2_psi1_simple1[[5]]
psi2 <- result_m2_sh2_psi2_simple1[[5]]
psi3 <- result_m2_sh2_psi3_simple1[[5]]
psi4 <- result_m2_sh2_psi4_simple1[[5]]
psi5 <- result_m2_sh2_psi5_simple1[[5]]
psi6 <- result_m2_sh2_psi6_simple1[[5]]
psi7 <- result_m2_sh2_psi7_simple1[[5]]
psi8 <- result_m2_sh2_psi8_simple1[[5]]
psi9 <- result_m2_sh2_psi9_simple1[[5]]
psi10 <- result_m2_sh2_psi10_simple1[[5]]


boxplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "psi",
        col = "#800020")

vioplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(psi[delta]), line = 1.5)
title(ylab = expression(psi[delta]), line = 1.5)

k1 <- result_m2_sh2_psi1_simple1[[6]]
k2 <- result_m2_sh2_psi2_simple1[[6]]
k3 <- result_m2_sh2_psi3_simple1[[6]]
k4 <- result_m2_sh2_psi4_simple1[[6]]
k5 <- result_m2_sh2_psi5_simple1[[6]]
k6 <- result_m2_sh2_psi6_simple1[[6]]
k7 <- result_m2_sh2_psi7_simple1[[6]]
k8 <- result_m2_sh2_psi8_simple1[[6]]
k9 <- result_m2_sh2_psi9_simple1[[6]]
k10 <- result_m2_sh2_psi10_simple1[[6]]


boxplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), #colMeans(k10),
        names = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "k",
        col = "#800020")
abline(h = 0.2, col = "orange")

vioplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), #colMeans(k10),
        names = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(psi[delta]), line = 1.5)
title(ylab = expression(k), line = 1.5)
abline(h = 0.2, col = "#800020")

res_psi8 <- result_m2_sh2_psi8_simple
g_psi8 <- result_m2_sh2_psi8_simple[[1]][,10]
h0_psi8 <- result_m2_sh2_psi8_simple[[2]][,10]
sigma_psi8 <- result_m2_sh2_psi8_simple[[3]][,10]
alpha_psi8 <- result_m2_sh2_psi8_simple[[4]][,10]
psi_psi8 <- result_m2_sh2_psi8_simple[[5]][,10]
k_psi8 <- result_m2_sh2_psi8_simple[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi8, type = "l", ylab = expression(g ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi8, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 46.46, col = "red")

plot(sigma_psi8, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(alpha_psi8, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.7), xlab = "iteration")

plot(psi_psi8, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.7, col = "red")

plot(k_psi8, type = "l", ylab = expression(k ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.2, col = "red")


res_psi9 <- result_m2_sh2_psi9_simple
g_psi9 <- result_m2_sh2_psi9_simple[[1]][,10]
h0_psi9 <- result_m2_sh2_psi9_simple[[2]][,10]
sigma_psi9 <- result_m2_sh2_psi9_simple[[3]][,10]
alpha_psi9 <- result_m2_sh2_psi9_simple[[4]][,10]
psi_psi9 <- result_m2_sh2_psi9_simple[[5]][,10]
k_psi9 <- result_m2_sh2_psi9_simple[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi9, type = "l", ylab = expression(g ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi9, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 46.46, col = "red")

plot(sigma_psi9, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(alpha_psi9, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.8), xlab = "iteration")

plot(psi_psi9, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.8, col = "red")

plot(k_psi9, type = "l", ylab = expression(k ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.2, col = "red")

