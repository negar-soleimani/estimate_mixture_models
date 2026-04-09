rm(list = ls())

source("data/prepare_data.R")
source("scripts/physics_model.R")
source("scripts/helper_function_OGP.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_OGP.R")

set.seed(12345)
n_iter <- 10000
burn_in <- 2500
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
a <- t_range
# init = c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)
# (g, h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta     <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)

# res <- mcmc_step6(
#   y = y, t = t, n_iter = n_iter, init = init, sigma_proposals = sigma_props,
#   g_init = TRUE, 
#   h0_init = FALSE,
#   sig2er_init = FALSE,
#   alpha_init = FALSE,
#   psi_init = FALSE,
#   k_init = FALSE,
#   Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
#   n_burnin = burn_in,
#   seuil = TRUE,  
#   s = 0.3       
# )
result <- mcmc_step6(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000)
#View(result)

g <- result[["theta"]][, 1]
h0 <- result[["theta"]][, 2]
sig <- result[["theta"]][, 3]
alpha <- result[["theta"]][, 4]
psi <- result[["theta"]][, 5]
k <- result[["theta"]][, 6]

delta <- result[["delta"]]
zeta <- result[["zeta"]]

boxplot(delta)

par(mfrow = c(2,3))
boxplot(g)
boxplot(h0)
boxplot(sig)
boxplot(alpha)
boxplot(psi)
boxplot(k)


g1 <- rep(1, n)
g2 <- -0.5 * (t * t_range)^2

w <- rep(1/n, n)

proj1 <- apply(delta, 1, function(d) sum(w * g1 * d))
proj2 <- apply(delta, 1, function(d) sum(w * g2 * d))

summary(proj1)
summary(proj2)
#########################################
set.seed(12345)

k <- 0.2
sigma_sq_err   <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2000

# (g, h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta     <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)

psi_vec <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

results_list <- vector("list", length(psi_vec))
names(results_list) <- paste0("psi", 1:length(psi_vec))

for (i in seq_along(psi_vec)) {
  sim_psi_delta <- psi_vec[i]
  cat("Running psi_delta =", sim_psi_delta, "\n")
  
  sigma_props <- c(NA, NA, NA, NA, 0.4, NA)
  init        <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.1)
  # if (i == 1) {
  #   # Psi1
  #   sigma_props <- c(NA, NA, NA, NA, 0.02, NA)
  # } else if (i == 10) {
  #   # Psi10
  #   sigma_props <- c(NA, NA, NA, NA, 0.3, NA)
  # } else {
  #   # Psi2 ... Psi9
  #   sigma_props <- c(NA, NA, NA, NA, 0.4, NA)
  # }
  
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
    
    res <- mcmc_step6(
      y_1, t, n_iter, init,
      sigma_props, mcmc_parameters,
      Sigma_theta, n_burnin = burn_in
    )
    
    g_chain[, v]     <- res$theta[, 1]
    h0_chain[, v]    <- res$theta[, 2]
    sigma_chain[, v] <- res$theta[, 3]
    alpha_chain[, v] <- res$theta[, 4]
    psi_chain[, v]   <- res$theta[, 5]
    k_chain[, v]     <- res$theta[, 6]
    
    delta_list[[v]] <- res$delta
    zeta_list[[v]]  <- res$zeta
    loglik_mat[, v] <- res$loglik
    accept_rate[v]  <- res$accept_rate_psi
  }
  
  results_list[[i]] <- list(
    g_chain     = g_chain,
    h0_chain    = h0_chain,
    sigma_chain = sigma_chain,
    alpha_chain = alpha_chain,
    psi_chain   = psi_chain,
    k_chain     = k_chain,
    delta_list  = delta_list,
    zeta_list   = zeta_list,
    loglik_mat  = loglik_mat,
    accept_rate = accept_rate,
    y_obs       = y_obs,
    psi_delta   = sim_psi_delta
  )
  
  # save(results_list[[i]], file = paste0("/Users/negarsoleimani/Documents/phd/paper1/", "result_m2_sh2_psi", i, "_simple1.RData"))
}
####---------------------------------------------
#or:
####---------------------------------------------
set.seed(12345)
Sigma_theta <- matrix(c(0.5,0,0,0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.2, 0.1)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
n_samples <- 10
n_iter <- 10000
burn_in <- 2000
# FALSE= fixed parameter
# mcmc parameter (g,h), sig2err, psidelta, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)

g            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
h0           <- matrix(NA, n_iter, n_samples, byrow = FALSE)
sigma_sq_err <- matrix(NA, n_iter, n_samples, byrow = FALSE)
alpha        <- matrix(NA, n_iter, n_samples, byrow = FALSE)
psi_delta    <- matrix(NA, n_iter, n_samples, byrow = FALSE)
k            <- matrix(NA, n_iter, n_samples, byrow = FALSE)


y_obs <- matrix(NA, length(t), n_samples, byrow = FALSE)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta)
  delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  y_1 = balldropg(t,c(9.8, 46.45)) + rnorm(n, 0, sqrt(0.01)) + delta
  y_obs[,v] <- y_1
  
  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, g_init=FALSE, h0_init= FALSE, sig2er_init = FALSE,
                        alpha_init = FALSE, psi_init = FALSE, k_init = FALSE, Sigma_theta, n_burnin = burn_in, seuil = FALSE, s = 0.3)
  g[,v] = results$theta[,1]
  h0[,v]=results$theta[,2]
  sigma_sq_err[,v]=results$theta[,3]
  alpha[,v] <- results$theta[,4]
  psi_delta[,v] <- results$theta[,5]
  k[,v] <- results$theta[,6]
}

result_m1_sh2_presentation <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs)

#############################################################
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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.1)
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
#result_m2_sh2_psi1_simple1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi1_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi1_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi1_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi2_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi2_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi2_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi2_ortho.RData")
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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi3_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi3_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi3_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi3_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi4_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi4_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi4_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi4_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi5_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi5_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi5_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi5_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi6_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi6_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi6_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi6_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi7_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi7_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi7_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi7_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi8_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi8_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi8_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi8_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
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
#result_m2_sh2_psi9_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi9_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi9_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality/result_m2_sh2_psi9_ortho.RData")


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
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.46, 0.1, 0.5, 0.2, 0.1)

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
#result_m2_sh2_psi10_ortho1 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

result_m2_sh2_psi10_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
save(result_m2_sh2_psi10_ortho,file = "/Users/negar/Documents/phd/Result/Model1/Orthogonality//result_m2_sh2_psi10_ortho.RData")

