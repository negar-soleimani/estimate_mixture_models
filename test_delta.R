# final model 1:
# simulation = classic , algorithm = classic
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

# ----------------------------- data -------------------------------------#
don <- read_xlsx("/Users/negarsoleimani/Documents/phd/paper1/Ball_drops_data.xlsx", sheet = 2)
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop)
don <- don[don$drop == 1, ]

t <- don$time
y <- don$Height
n <- length(t)
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)
a <- t_range
# ----------------------------- Helper functions ------------------------------
GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
}

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

# ----------------------------- MCMC Step -------------------------------------#
mcmc_step_fix_delta <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, delta_true, n_burnin=1000) {
  
  freeze_delta_zeta <- mcmc_parameters[6]
  
  # Total iterations = burn-in + desired samples
  total_iter <- n_burnin + n_iter
  
  theta <- init
  delta <- delta_true
  #delta <- rep(0, length(y))
  chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
  colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  #chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
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
    
    # Method2 for calculate the probability of the zeta:(log_sum_exp: https://rpubs.com/FJRubio/LSE and https://en.wikipedia.org/wiki/LogSumExp)
    log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
    log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    #prob_zeta_1 <- exp(log_w1 - log_den)  #  P(zeta=1|y)
    prob_zeta_2 <- exp(log_w2 - log_den)  #  P(zeta=2|y)
    #zeta <- ifelse(runif(length(y)) < prob_zeta_1, 1, 2)
    zeta <- ifelse(runif(length(y)) < prob_zeta_2, 1, 0) #= zeta <- 1 + (runif(length(y)) < prob_zeta_1) # sample zeta: 2 with prob post_p2, otherwise 1
    
    #--------------------
    if (freeze_delta_zeta) {
      zeta <- rep(0, length(y))
    }
    #--------------------
    chain_zeta[iter, ] = zeta
    
    #-------------------------------- log_likelihood --------------------------------#   
    
    log_likelihood <- sum(log(ifelse(zeta == 0,
                                     alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
                                     (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    loglik_chain[iter] <- log_likelihood
    
    # zeta_1_indices <- which(zeta == 2)
    # if (length(zeta_1_indices) > 0) {
    #   y_m <- y[zeta_1_indices]
    #   Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_1_indices)) +
    #     Sigma_delta[zeta_1_indices, zeta_1_indices, drop = FALSE]
    #   Sigma_delta_ym <- Sigma_delta[, zeta_1_indices, drop = FALSE]
    #   Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
    #   mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_1_indices])
    #   Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
    #   Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
    #   delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    # } else {
    #   delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))}
    # 
    # #--------------------
    # if (freeze_delta_zeta) {
    #   delta <- rep(0, n)
    # }
    # #--------------------
    # -------------------------------- Gibbs step for theta (g, h0) ---------------- #
    # flat prior on theta 
    zeta_0_indices <- which(zeta == 0)
    zeta_1_indices <- which(zeta == 1)
    
    X <- cbind(1, -0.5 * (t * t_range)^2)
    x1 <- X[zeta_0_indices, , drop = FALSE]
    x2 <- X[zeta_1_indices, , drop = FALSE]
    
    y1 <- matrix(y[zeta_0_indices], ncol = 1)
    y2 <- matrix(y[zeta_1_indices], ncol = 1)
    d2 <- matrix(delta[zeta_1_indices], ncol = 1)
    
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
    
    # -------------------------------- Gibbs step for sigma_sq_err ---------------- #
    R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
    if(n > 0){
      R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
      
      quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
    } 
    else {
      quad_form_delta <- 0
    }
    
    if (freeze_delta_zeta){
      f_theta <- balldropg(t, c(g, h0))
      idx0 <- which(zeta == 0)
      idx1 <- which(zeta == 1)
      n1   <- length(idx0)
      n2   <- length(idx1)
      residual1 <- y[idx0] - f_theta[idx0]
      rss1   <- sum(residual1^2)
      shape_err <- (n / 2) + 1
      rate_err  <- 0.5 * rss1
    } else {
      f_theta <- balldropg(t, c(g, h0))
      idx0 <- which(zeta == 0)
      idx1 <- which(zeta == 1)
      n1   <- length(idx0)
      n2   <- length(idx1)
      residual1 <- y[idx0] - f_theta[idx0]
      rss1   <- sum(residual1^2)
      
      residual2 <- y[idx1] - f_theta[idx1] - delta[idx1]
      rss2   <- sum(residual2^2)
      d <- 2
      # When we consider the Jeffreys prior, I use the following code:
      rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
      shape_err <- n + (d/2)
    }
    sigma_sq_err <- 1/rgamma(1, shape = shape_err, rate = rate_err)
    theta[3] <- sigma_sq_err
    
    if(mcmc_parameters[2] == FALSE){
      sigma_sq_err <- init[3]
      theta[3] <- sigma_sq_err
    }
    
    # -------------------------------- Gibbs step for psi_delta ---------------- #
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
    if(mcmc_parameters[3] == FALSE){
      psi_delta <- init[5]
      theta[5] <- psi_delta
    }
    
    # -------------------------------- Gibbs step for k ------------------------- #
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
    
    # -------------------------------- Gibbs step for alpha --------------------- #
    alpha_param <- rbeta(1, sum(zeta == 0) + 0.5, sum(zeta == 1) + 0.5)
    theta[4] <- alpha_param
    
    if(mcmc_parameters[5] == FALSE){
      alpha_param <- init[4]
      theta[4] <- alpha_param
    }
    
    theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
    chain_theta[iter, ] <- theta
    #chain_delta[iter, ] <- delta
  }
  
  #----------------------------------------------------------------#   
  
  # Remove burn-in samples if n_burnin > 0
  if (n_burnin > 0) {
    keep_indices <- (n_burnin + 1):total_iter
    chain_theta <- chain_theta[keep_indices, , drop = FALSE]
    #chain_delta <- chain_delta[keep_indices, , drop = FALSE]
    chain_zeta <- chain_zeta[keep_indices, , drop = FALSE]
    loglik_chain <- loglik_chain[keep_indices]
    
    # Adjust acceptance rate calculation for post-burn-in period only
    cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
  } else {
    cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
  }
  return(list(theta = chain_theta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_psi / n_iter))
}

# ----------------------------- Parameters & Run -------------------------------------#
set.seed(12345)
n_samples <- 50
n_iter <- 20000
burn_in <- 2000
#(g,h), sig, psi,k, alpha, freeze
mcmc_parameters <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)
#g, h, sig, alpha, psi, k
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.2)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)

# Exemple : delta_true simulé à partir d'une GP
k_true <- 0.2
sigma_sq_err_true <- 0.01
psi_delta_true <- 0.5
sigma_sq_delta_true <- sigma_sq_err_true / k_true

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
delta_true <- as.vector(rmvnorm(1, rep(0, n), GP_covariance(t, sigma_sq_delta_true, psi_delta_true)))

y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err_true)) + delta_true
# alpha_true <- 0.5
# z_true <- rbinom(n,1,alpha_true)
# 
# f <- balldropg(t, c(9.8, 46.46))
# eps <- rnorm(n,0,sqrt(sigma_sq_err_true))
# 
# y_1 <- f + eps + z_true * delta_true

v=1
y_obs[, v] <- y_1
res <- mcmc_step_fix_delta(y_1, t, n_iter, init, sigma_props, mcmc_parameters, Sigma_theta, delta_true=delta_true, n_burnin=burn_in)

g_chain[, v]     <- res$theta[, 1]
h0_chain[, v]    <- res$theta[, 2]
sigma_chain[, v] <- res$theta[, 3]
alpha_chain[, v] <- res$theta[, 4]
psi_chain[, v]   <- res$theta[, 5]
k_chain[, v]     <- res$theta[, 6]
zeta_list[[v]]   <- res$zeta
loglik_mat[, v]  <- res$loglik
accept_rate[v]   <- res$accept_rate_psi

results <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, zeta_list, loglik_mat, accept_rate)

# ----------------------------- Résultats -------------------------------------#
theta_chain <- res$theta
apply(theta_chain, 2, mean)  # Moyennes postérieures


boxplot(delta_true, ylab = "delta")
plot(delta_true, type = "l")

g <- res[["theta"]][,1]
h0 <- res[["theta"]][,2]
sig <- res[["theta"]][,3]
alpha <- res[["theta"]][,4]
psi <- res[["theta"]][,5]
k <- res[["theta"]][,6]

par(mfrow = c(2,3))
plot(g, type = "l")
abline(h = 9.8, col = "red")
plot(h0, type = "l")
abline(h = 46.46, col = "red")
plot(sig, type = "l")
abline(h = 0.01, col = "red")
plot(alpha, type = "l")
plot(psi, type = "l")
abline(h = 0.5, col = "red")
plot(k, type = "l")
abline(h = 0.2, col = "red")
