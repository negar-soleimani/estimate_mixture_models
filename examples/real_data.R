# mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000) {
#   # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T")
#   
#   # Total iterations = burn-in + desired samples
#   total_iter <- n_burnin + n_iter
#   
#   theta <- init
#   delta <- rep(0, length(y))
#   chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
#   colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
#   chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
#   chain_zeta <- matrix(NA, nrow = total_iter, ncol = length(y))
#   
#   loglik_chain <- numeric(total_iter)
#   accept_psi <- 0
#   
#   for (iter in 1:total_iter) {
#     g <- theta[1]; h0 <- theta[2]; sigma_sq_err <- theta[3]
#     alpha_param <- theta[4]; psi_delta <- theta[5]; k <- theta[6]
#     sigma_sq_delta <- sigma_sq_err / k
#     
#     Sigma_delta <- GP_covariance(t, sigma_sq_delta, psi_delta)
#     f_theta <- balldropg(t, c(g, h0))
#     mean1 <- f_theta
#     mean2 <- f_theta + delta
#     s <- 0.3
#     
#     #prob_zeta <- 1 / (1 + (alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err))) /
#     #                    ((1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err))))*(abs(delta)>s)
#     #zeta <- 1 + (runif(length(y)) < prob_zeta)
#     
#     # log_sum_exp <- function(x, y) {
#     #   m <- pmax(x, y)
#     #   m + log(exp(x - m) + exp(y - m))
#     # }
#     # log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
#     # log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
#     # 
#     # log_denominator <- log_sum_exp(log_w1, log_w2)
#     # 
#     # prob_zeta_base <- exp(log_w2 - log_denominator)
#     # 
#     # prob_zeta <- prob_zeta_base * (abs(delta) > s)
#     # zeta <- 1 + (runif(length(y)) < prob_zeta)
#     # #zeta <- ifelse(runif(length(y)) < prob_zeta, 1, 2)
#     # chain_zeta[iter, ] <- zeta
#     log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE) # comp 1
#     log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE) # comp 2
#     log_max <- pmax(log_w1, log_w2)
#     log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
#     
#     if (seuil) {
#       prob_zeta_base <- exp(log_w2 - log_den)
#       prob_zeta <- prob_zeta_base * (abs(delta) > s)
#       zeta <- ifelse(runif(length(y)) < prob_zeta, 2, 1)
#       #zeta <- 1 + (runif(length(y)) < prob_zeta)
#     } else {
#       prob_zeta <- exp(log_w2 - log_den)
#       #zeta <- 1 + (runif(length(y)) < prob_zeta)
#       zeta <- ifelse(runif(length(y)) < prob_zeta, 2, 1)
#     }
#     
#     chain_zeta[iter, ] <- zeta
#     
#     
#     log_likelihood <- sum(log(ifelse(zeta == 1,
#                                      alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
#                                      (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
#     loglik_chain[iter] <- log_likelihood
#     
#     zeta_2_indices <- which(zeta == 2)
#     if (length(zeta_2_indices) > 0) {
#       y_m <- y[zeta_2_indices]
#       Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
#         Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
#       Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
#       Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
#       mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
#       Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
#       delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
#     }
#     else delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
#     
#     # Gibbs for theta
#     zeta_1_indices <- which(zeta == 1)
#     zeta_2_indices <- which(zeta == 2)
#     X <- cbind(1, -0.5 * (t * t_range + t_min)^2)
#     x1 <- X[zeta_1_indices, , drop = FALSE]  
#     x2 <- X[zeta_2_indices, , drop = FALSE]
#     theta_hat     <- matrix(c(46.45, 9.8), ncol = 1)
#     inv_sigma_theta <- solve(Sigma_theta)
#     A <- ((t(x1) %*% x1) / theta[3]) + ((t(x2) %*% x2) / theta[3]) + inv_sigma_theta
#     Sigmapost_theta <- solve(A)
#     y1 <- matrix(y[zeta_1_indices], ncol = 1)
#     y2 <- matrix(y[zeta_2_indices], ncol = 1)
#     d2 <- matrix(delta[zeta_2_indices], ncol = 1)
#     B <- (t(x1) %*% y1) / theta[3] + 
#       (t(x2) %*% y2) / theta[3] -
#       (t(x2) %*% d2) / theta[3] +
#       inv_sigma_theta %*% theta_hat
#     Mupost_theta <- Sigmapost_theta %*% B
#     theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
#     h0 <- theta_sample[1];  g <- theta_sample[2]
#     
#     #if(mcmc_parameters[1] == FALSE){
#     #  g <- init[1]
#     #  h0 <- init[2]
#     #}
#     
#     ## g = fixer
#     if (mcmc_parameters[1] == FALSE) {
#       g <- init[1]
#       h0 <- theta_sample[1]
#     }
#     
#     ## h0 = fixed
#     #if (mcmc_parameters[1] == FALSE) {
#     #  h0 <- init[2] 
#     #  g <- theta_sample[2] 
#     #}
#     
#     # Gibbs for sigma_sq_err
#     
#     R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
#     if(n > 0){
#       R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
#       
#       quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
#     } 
#     else {
#       quad_form_delta <- 0
#     }
#     
#     f_theta <- balldropg(t, c(g, h0))
#     idx1 <- which(zeta == 1)
#     idx2 <- which(zeta == 2)
#     residual1 <- y[idx1] - f_theta[idx1]
#     rss1   <- sum(residual1^2)
#     
#     residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
#     rss2   <- sum(residual2^2)
#     
#     rate_err <- 0.5 + (0.5 * ( rss1 + rss2 + (theta[6] * quad_form_delta)))
#     
#     shape_err <- 4 + n
#     sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
#     
#     if(mcmc_parameters[2] == FALSE){
#       sigma_sq_err <- init[3]
#     }
#     
#     # MH for psi_delta
#     psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
#     sigma_sq_delta_prop <- sigma_sq_err / k
#     Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
#     log_prop_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
#     log_prop_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
#     log_prior_current <- dbeta(psi_delta, shape1 = 4, shape2 = 6, log=TRUE)
#     log_prior_prop <- dbeta(psi_prop, shape1 = 5, shape2 = 6, log=TRUE)
#     log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
#     log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
#     log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + log_prop_current - log_prop_prop
#     if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
#       psi_delta <- psi_prop
#       #theta[5] <- psi_delta
#       accept_psi <- accept_psi + 1
#     }
#     
#     if(mcmc_parameters[3] == FALSE){
#       psi_delta <- init[5]
#     }
#     
#     # Gibbs for k
#     R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
#     alpha_k <- (n / 2) + 1
#     beta_k <- (1 / (2 * sigma_sq_err)) * sum(delta * (solve(R, delta)))
#     for (try_k in 1:100) {
#       k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
#       #if (k_prop >= 0.1 && k_prop <= 0.9) {
#       if (k_prop > 0 && k_prop < 1) {
#         k <- k_prop
#         break
#       }
#     }
#     
#     
#     if(mcmc_parameters[4] == FALSE){
#       k <- init[6]
#     }
#     
#     # Gibbs step for alpha
#     
#     zetabis <- zeta
#     zetabis[abs(delta)<s] <- 0
#     alpha_prime <- rbeta(1, sum(zetabis == 1) + 0.5, sum(zetabis == 2) + 0.5)
#     alpha_param <- rbeta(1,sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
#     
#     if(mcmc_parameters[5] == FALSE){
#       alpha_param <- init[4]
#     }
#     
#     theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
#     chain_theta[iter, ] <- theta
#     chain_delta[iter, ] <- delta
#     loglik_chain[iter] <- sum(log(
#       ifelse(zeta==1,
#              alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
#              (1-alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))
#     ))
#   }
#   
#   # Remove burn-in samples if n_burnin > 0
#   if (n_burnin > 0) {
#     keep_indices <- (n_burnin + 1):total_iter
#     chain_theta <- chain_theta[keep_indices, , drop = FALSE]
#     chain_delta <- chain_delta[keep_indices, , drop = FALSE]
#     chain_zeta <- chain_zeta[keep_indices, , drop = FALSE]
#     loglik_chain <- loglik_chain[keep_indices]
#     
#     # Adjust acceptance rate calculation for post-burn-in period only
#     cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
#   } else {
#     cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
#   }
#   return(list(theta = chain_theta, delta = chain_delta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_psi / n_iter))
# }
# sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
# mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
# Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2)
# ress <- mcmc_step6(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000)


source("data/prepare_data.R")
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

don <- read_xlsx("/Users/negar/Documents/phd/estimate_mixture_models-main/data/Ball_drops_data.xlsx", sheet = 2)
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
a <- t_range
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

set.seed(12345)
n_iter <- 10000
burn_in <- 2000
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)

# init = c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

res <- mcmc_step6(
  y = y, t = t, n_iter = n_iter, init = init, sigma_proposals = sigma_props,
  g_init = TRUE, 
  h0_init = TRUE,
  sig2er_init = FALSE,
  alpha_init = FALSE,
  psi_init = FALSE,
  k_init = FALSE,
  Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
  n_burnin = burn_in,
  seuil = FALSE,  
  s = 0.3       
)

