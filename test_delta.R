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
  x_vec <- cbind(1, -0.5 * (t * t_range)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

# ----------------------------- MCMC Step -------------------------------------#
mcmc_step_fix_delta <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters,
                                Sigma_theta, delta_true, n_burnin=1000) {
  
  # mcmc_parameters = c(theta(g,h0), sigma_sq_err, psi_delta, k, alpha, freeze_delta_zeta)
  freeze_delta_zeta <- mcmc_parameters[6]
  
  total_iter <- n_burnin + n_iter
  
  theta <- init
  delta <- delta_true
  
  chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
  colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  
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
    
    # --- sample zeta (mixture allocation) ---
    log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
    log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    prob_zeta_2 <- exp(log_w2 - log_den)  # P(zeta=1 | y) in your coding below
    
    zeta <- ifelse(runif(length(y)) < prob_zeta_2, 1, 0)  # 1 = discrepancy component, 0 = no-discrepancy
    
    # -------------------- IMPORTANT CHANGE --------------------
    # For identifiability diagnostic conditional on delta_true:
    # force all points to belong to discrepancy component.
    if (freeze_delta_zeta) {
      zeta <- rep(1, length(y))
    }
    # ----------------------------------------------------------
    chain_zeta[iter, ] <- zeta
    
    # --- log-likelihood (conditional on zeta) ---
    log_likelihood <- sum(log(ifelse(zeta == 0,
                                     alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
                                     (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    loglik_chain[iter] <- log_likelihood
    
    # --- Gibbs step for theta (g, h0): flat prior ---
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
    Mupost_theta <- Sigmapost_theta %*% B
    
    theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
    h0 <- theta_sample[1]; g <- theta_sample[2]
    theta[1] <- g
    theta[2] <- h0
    
    if (mcmc_parameters[1] == FALSE) {
      theta[1] <- init[1]
      theta[2] <- init[2]
      g <- theta[1]; h0 <- theta[2]
    }
    
    # --- Gibbs step for sigma_sq_err ---
    R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
    R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
    quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
    
    f_theta <- balldropg(t, c(g, h0))
    idx0 <- which(zeta == 0)
    idx1 <- which(zeta == 1)
    
    residual1 <- y[idx0] - f_theta[idx0]
    rss1 <- sum(residual1^2)
    
    residual2 <- y[idx1] - f_theta[idx1] - delta[idx1]
    rss2 <- sum(residual2^2)
    
    d_dim <- 2
    rate_err <- 0.5 * (rss1 + rss2 + (k * quad_form_delta))
    shape_err <- n + (d_dim / 2)
    
    sigma_sq_err <- 1 / rgamma(1, shape = shape_err, rate = rate_err)
    theta[3] <- sigma_sq_err
    
    if (mcmc_parameters[2] == FALSE) {
      theta[3] <- init[3]
      sigma_sq_err <- theta[3]
    }
    
    # --- MH step for psi_delta ---
    psi_prop <- rtruncnorm(1, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    
    sigma_sq_delta_prop <- sigma_sq_err / k
    Sigma_delta <- GP_covariance(t, sigma_sq_delta_prop, psi_delta)
    Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
    
    log_prop_current <- log(dtruncnorm(psi_delta, a = 0.1, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
    log_prop_prop <- log(dtruncnorm(psi_prop, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
    
    log_prior_current <- dunif(psi_delta, min = 0.1, max = 1, log = TRUE)
    log_prior_prop <- dunif(psi_prop, min = 0.1, max = 1, log = TRUE)
    
    log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
    log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    
    log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) +
      (log_prop_current - log_prop_prop)
    
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      psi_delta <- psi_prop
      theta[5] <- psi_delta
      accept_psi <- accept_psi + 1
    }
    
    if (mcmc_parameters[3] == FALSE) {
      theta[5] <- init[5]
      psi_delta <- theta[5]
    }
    
    # --- Gibbs step for k (truncated gamma on (0,1)) ---
    alpha_k <- (n / 2) + 1
    beta_k <- (1 / (2 * sigma_sq_err)) * quad_form_delta
    
    for (try_k in 1:100) {
      k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
      if (k_prop > 0 && k_prop < 1) {
        k <- k_prop
        theta[6] <- k
        break
      }
    }
    
    if (mcmc_parameters[4] == FALSE) {
      theta[6] <- init[6]
      k <- theta[6]
    }
    
    # --- Gibbs step for alpha ---
    alpha_param <- rbeta(1, sum(zeta == 0) + 0.5, sum(zeta == 1) + 0.5)
    theta[4] <- alpha_param
    
    if (mcmc_parameters[5] == FALSE) {
      theta[4] <- init[4]
      alpha_param <- theta[4]
    }
    
    theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
    chain_theta[iter, ] <- theta
  }
  
  # Remove burn-in
  if (n_burnin > 0) {
    keep_indices <- (n_burnin + 1):total_iter
    chain_theta <- chain_theta[keep_indices, , drop = FALSE]
    chain_zeta <- chain_zeta[keep_indices, , drop = FALSE]
    loglik_chain <- loglik_chain[keep_indices]
    cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
  } else {
    cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
  }
  
  return(list(theta = chain_theta, zeta = chain_zeta, loglik = loglik_chain,
              accept_rate_psi = accept_psi / n_iter))
}

# -----------------------------------------------------------------#
set.seed(12345)

n_samples <- 50
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA)

# (g,h0), sig2err, psi, k, alpha, freeze_delta_zeta
# In this diagnostic: freeze_delta_zeta = TRUE (force zeta=1), and fix alpha (set alpha update OFF)
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)

# g, h0, sig, alpha, psi, k
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)

# True values for simulation
k_true <- 0.1
sigma_sq_err_true <- 0.01
psi_delta_true <- 0.5
sigma_sq_delta_true <- sigma_sq_err_true / k_true

g_mean <- numeric(n_samples)
h0_mean <- numeric(n_samples)
sig_mean <- numeric(n_samples)
psi_mean <- numeric(n_samples)
k_mean <- numeric(n_samples)

accept_rate <- numeric(n_samples)

for (v in 1:n_samples) {
  
  delta_true <- as.vector(rmvnorm(1, rep(0, n),
                                  GP_covariance(t, sigma_sq_delta_true, psi_delta_true)))
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err_true)) + delta_true
  
  res <- mcmc_step_fix_delta(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                             Sigma_theta, delta_true = delta_true, n_burnin = burn_in)
  
  theta_chain <- res$theta
  
  g_mean[v]   <- mean(theta_chain[,1])
  h0_mean[v]  <- mean(theta_chain[,2])
  sig_mean[v] <- mean(theta_chain[,3])
  psi_mean[v] <- mean(theta_chain[,5])
  k_mean[v]   <- mean(theta_chain[,6])
  
  accept_rate[v] <- res$accept_rate_psi
}

cat("Mean acceptance rate psi_delta:", mean(accept_rate), "\n")

par(mfrow=c(2,3))
boxplot(g_mean, ylab="k", col = "lightseagreen")
abline(h=9.8, col="orange")
boxplot(h0_mean, ylab="k", col = "lightseagreen")
abline(h=46.45, col="orange")
boxplot(sig_mean, ylab="k", col = "lightseagreen")
abline(h=0.01, col="orange")
boxplot(psi_mean, ylab=expression(gamma[delta]), col = "lightseagreen")
abline(h=psi_delta_true, col="orange")
boxplot(k_mean, ylab="k", col = "lightseagreen")
abline(h=k_true, col="orange")


# theta_chain <- res$theta
# 
# psi_delta_true <- 0.5 
# k_true <- 0.1     

par(mfrow = c(1,2))

# Trace plot for psi_delta
plot(theta_chain[,5], type = "l",
     xlab = "iteration", ylab = expression(psi_delta),
     main = expression("Trace plot: " * psi_delta))
abline(h = psi_delta_true, col = "orange", lwd = 2)

# Trace plot for k
plot(theta_chain[,6], type = "l",
     xlab = "iteration", ylab = "k",
     main = "Trace plot: k")
abline(h = k_true, col = "orange", lwd = 2)