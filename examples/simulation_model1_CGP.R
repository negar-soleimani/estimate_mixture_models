rm(list = ls())
source("data/prepare_data.R")
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

# Simulation the data of model 2 with different psi_delta
# ## Psi1 ####################################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1
sim_psi_delta <- 0.01 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000
# 
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)
#(g,h0), sigma, psi, k, alpha
# mcmc parameter (g,h), sig2err, psidelta, k, alpha, freeze delta-zeta
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)
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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi1_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi1_simple.RData")

## Psi2 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.1 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi2_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi2_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi2_simple.RData")

# library(ggplot2)
# library(tidyr)

# res_all <- result_m2_sh2_psi1_simple_ex
# 
# g_mat       <- res_all[[1]]
# h0_mat      <- res_all[[2]]
# sigma_mat   <- res_all[[3]]
# alpha_mat   <- res_all[[4]]
# psi_mat     <- res_all[[5]]
# k_mat       <- res_all[[6]]
# delta_list  <- res_all[[7]]
# zeta_list   <- res_all[[8]]
# loglik_mat  <- res_all[[9]]
# accept_rate <- res_all[[10]]

# n_iter    <- nrow(g_mat)
# n_samples <- ncol(g_mat)
# n_obs     <- ncol(delta_list[[1]])
# 
# cat("MCMC diagnostics\n")
# cat("  iterations (post-burn-in): ", n_iter, "\n")
# cat("  replicate chains:          ", n_samples, "\n")
# cat("  observations:              ", n_obs, "\n")
# cat("  psi acceptance rate(s):    ",
#     paste(round(accept_rate, 3), collapse = ", "), "\n\n")
# 
#
# true_g     <- 9.8
# true_h0    <- 46.46
# true_sigma <- 0.01
# true_psi   <- 0.5
# true_k     <- 0.1

# =============================================================================
# 1) TRACE PLOTS  - Markov chains for each parameter
# =============================================================================

# op <- par(mfrow = c(2, 3),
#           mar   = c(3.5, 4, 2.5, 1),
#           mgp   = c(2.2, 0.8, 0))
# plot(g_mat[, 1],     type = "l", col = "steelblue",
#      xlab = "iteration", ylab = "g",
#      main = "Trace: g")
# abline(h = true_g, col = "red", lty = 2, lwd = 1.5)
# 
# 
# plot(h0_mat[, 1],    type = "l", col = "steelblue",
#      xlab = "iteration", ylab = expression(h[0]),
#      main = expression(paste("Trace: ", h[0])))
# abline(h = true_h0, col = "red", lty = 2, lwd = 1.5)
# 
# plot(sigma_mat[, 1], type = "l", col = "steelblue",
#      xlab = "iteration", ylab = expression(lambda^2),
#      main = expression(paste("Trace: ", lambda^2)))
# abline(h = true_sigma, col = "red", lty = 2, lwd = 1.5)
# 
# plot(alpha_mat[, 1], type = "l", col = "steelblue",
#      xlab = "iteration", ylab = expression(alpha),
#      main = expression(paste("Trace: ", alpha)))
# 
# plot(psi_mat[, 1],   type = "l", col = "steelblue",
#      xlab = "iteration", ylab = expression(gamma[delta]),
#      main = expression(paste("Trace: ", gamma[delta])))
# abline(h = true_psi, col = "red", lty = 2, lwd = 1.5)
# 
# plot(k_mat[, 1],     type = "l", col = "steelblue",
#      xlab = "iteration", ylab = "k",
#      main = "Trace: k")
# abline(h = true_k, col = "red", lty = 2, lwd = 1.5)
# 
# par(op)

# plot_overlay <- function(M, ylab_expr, true_val = NA, main_txt = "") {
#   matplot(M, type = "l", lty = 1,
#           col  = rainbow(ncol(M), alpha = 0.7),
#           xlab = "iteration", ylab = ylab_expr,
#           main = main_txt)
#   if (!is.na(true_val)) abline(h = true_val, col = "black", lty = 2, lwd = 2)
# }
# 
# op <- par(mfrow = c(2, 3),
#           mar   = c(3.5, 4, 2.5, 1),
#           mgp   = c(2.2, 0.8, 0))
# plot_overlay(g_mat,     "g",                 true_g,     "All chains: g")
# plot_overlay(h0_mat,    expression(h[0]),    true_h0,    expression(paste("All chains: ", h[0])))
# plot_overlay(sigma_mat, expression(lambda^2), true_sigma, expression(paste("All chains: ", lambda^2)))
# plot_overlay(alpha_mat, expression(alpha),   NA,         expression(paste("All chains: ", alpha)))
# plot_overlay(psi_mat,   expression(gamma[delta]), true_psi, expression(paste("All chains: ", gamma[delta])))
# plot_overlay(k_mat,     "k",                 true_k,     "All chains: k")
# par(op)

## Psi3 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.2
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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

result_m2_sh2_psi3_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi3_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi3_simple.RData")


## Psi4 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.3 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi4_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi4_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi4_simple.RData")


## Psi5 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.4
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi5_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi5_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi5_simple.RData")


## Psi6 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.5 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi6_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi6_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi6_simple.RData")


## Psi7 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.6 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi7_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi7_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi7_simple.RData")


## Psi8 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.7 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi8_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi8_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi8_simple.RData")


## Psi9 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.8 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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
  y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))  + delta
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

result_m2_sh2_psi9_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi9_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi9_simple.RData")


## Psi10 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1
sim_psi_delta <- 0.9
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 50
n_iter <- 10000
burn_in <- 2000
#
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.5, 0.1)

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

result_m2_sh2_psi10_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi10_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi10_simple.RData")