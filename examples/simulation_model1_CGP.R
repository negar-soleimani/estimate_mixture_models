source("data/prepare_data.R")
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_CGP.R")

# Simulation the data of model 2 with different psi_delta
# ## Psi1 ####################################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1
sim_psi_delta <- 0.1 # psi_delta -> [0.1, 1]
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
init <- c(9.8, 46.46, 0.01, 0.5, 0.1, 0.1)
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
#result_m2_sh2_psi1_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi1_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
# save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")


## Psi2 ############################################
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
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.1)

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
#result_m2_sh2_psi2_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi2_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi2_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi2_simple.RData")
## Psi3 ############################################
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
#result_m2_sh2_psi3_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi3_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi3_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi3_simple.RData")


## Psi4 ############################################
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
#result_m2_sh2_psi4_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi4_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi4_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi4_simple.RData")


## Psi5 ############################################
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
#result_m2_sh2_psi5_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi5_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi5_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi5_simple.RData")


## Psi6 ############################################
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
#result_m2_sh2_psi6_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi6_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi6_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi6_simple.RData")


## Psi7 ############################################
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
#result_m2_sh2_psi7_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi7_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi7_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi7_simple.RData")


## Psi8 ############################################
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
#result_m2_sh2_psi8_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi8_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi8_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi8_simple.RData")


## Psi9 ############################################
set.seed(12345)
k = 0.1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.9 
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
#result_m2_sh2_psi9_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi9_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi9_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi9_simple.RData")


## Psi10 ############################################
# set.seed(12345)
# k = 0.2
# #sigma_sq_delta <- 0.1
# sim_psi_delta <- 0.9
# sigma_sq_err <- 0.01
# sigma_sq_delta <- sigma_sq_err / k
# n_samples <- 50
# n_iter <- 10000
# burn_in <- 2000
# # 
# sigma_props <- c(NA, NA, NA, NA, 0.5, NA)
# #(g,h0), sigma, psi, k, alpha
# mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
# Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# init <- c(9.8, 46.46, 0.1, 0.5, 0.2, 0.2)
# 
# g_chain     <- matrix(NA, n_iter, n_samples)
# h0_chain    <- matrix(NA, n_iter, n_samples)
# sigma_chain <- matrix(NA, n_iter, n_samples)
# alpha_chain <- matrix(NA, n_iter, n_samples)
# psi_chain   <- matrix(NA, n_iter, n_samples)
# k_chain     <- matrix(NA, n_iter, n_samples)
# loglik_mat  <- matrix(NA, n_iter, n_samples)
# delta_list  <- vector("list", n_samples)
# zeta_list   <- vector("list", n_samples)
# accept_rate <- numeric(n_samples)
# # 
# y_obs <- matrix(NA, n, n_samples)
# for (v in 1:n_samples) {
#   
#   y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err))
#   y_obs[, v] <- y_1
#   
#   res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
#                     Sigma_theta, n_burnin = burn_in)
#   #   
#   g_chain[, v]     <- res$theta[, 1]
#   h0_chain[, v]    <- res$theta[, 2]
#   sigma_chain[, v] <- res$theta[, 3]
#   alpha_chain[, v] <- res$theta[, 4]
#   psi_chain[, v]   <- res$theta[, 5]
#   k_chain[, v]     <- res$theta[, 6]
#   delta_list[[v]]  <- res$delta
#   zeta_list[[v]]   <- res$zeta
#   loglik_mat[, v]  <- res$loglik
#   accept_rate[v]   <- res$accept_rate_psi
# }
# result_m2_sh2_psi10_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi10_simple <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi10_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/m2_simplegp/result_m2_sh2_psi10_simple.RData")
