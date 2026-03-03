source("data/prepare_data.R")
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.01 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500
#n_flat <- 10

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  #Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  #delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  
  #delta[1:n_flat] <- 0
  
  #y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  #for fig 1 rep 50 times
  ndis <- 30
  n <- length(y)
  #y = y + (46-y)*.1# Fig 2:y = y + delta_i ; i = 43, when delta_i ~ GP(0, (sig2err/k)*exp(-(ti-tj)2/psi))
  #y_1 = y + (46-y)*.3 * c(rep(0, ndis), rep(1, length(y) - ndis))
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))

  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
    g_init = FALSE,
    h0_init = FALSE,
    sig2er_init = FALSE,
    alpha_init = FALSE,
    psi_init = FALSE,
    k_init = FALSE,
    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
    n_burnin = burn_in,
    seuil = TRUE,
    s = 0.1 
  )
  
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

result_m2_sh2_psi1_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

# par(mfrow=c(1,1))
# boxplot(result_m2_sh2_psi1_simple_seuil[[7]][[1]])
# abline(h=0)
# g <- result_m2_sh2_psi1_simple_seuil[[1]]
# h0 <- result_m2_sh2_psi1_simple_seuil[[2]]
# sig <- result_m2_sh2_psi1_simple_seuil[[3]]
# alpha <- result_m2_sh2_psi1_simple_seuil[[4]]
# psi <- result_m2_sh2_psi1_simple_seuil[[5]]
# k <- result_m2_sh2_psi1_simple_seuil[[6]]
# 
# par(mfrow=c(2,3))
# boxplot(colMeans(g))
# abline(h=9.8, col = "orange")
# boxplot(colMeans(h0))
# abline(h=46.45, col = "orange")
# boxplot(colMeans(sig))
# abline(h=0.01, col = "orange")
# boxplot(colMeans(alpha))
# boxplot(colMeans(psi))
# boxplot(colMeans(k))
# abline(h=0.1, col = "orange")
# 
# par(mfrow=c(2,3))
# plot(g[,50], type = "l")
# abline(h=9.8, col = "orange")
# plot(h0[,50], type = "l")
# abline(h=46.45, col = "orange")
# plot(sig[,50], type = "l")
# abline(h=0.01, col = "orange")
# plot(alpha[,50], type = "l")
# plot(psi[,50], type = "l")
# plot(k[,50], type = "l")
# abline(h=0.2, col = "orange")
# 
# prob_zeta_model2 <- colMeans(res$zeta == 2)
# print(prob_zeta_model2)


set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.1 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi2_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.2 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi3_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.3 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi4_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.4 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi5_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.5 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi6_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")


set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.6 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi7_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")


set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.7 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi8_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.8 
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi9_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.9
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 50
n_iter <- 10000
burn_in <- 2500

init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)

y_obs <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {
  
  y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  ndis <- 30
  n <- length(y)
  y_1 = y + 2 * c(rep(0, ndis), rep(1, length(y) - ndis))
  
  y_obs[, v] <- y_1
  
  res <- mcmc_step6(y = y_1, t = t, n_iter = n_iter, init = init,
                    sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
                    g_init = FALSE,
                    h0_init = FALSE,
                    sig2er_init = FALSE,
                    alpha_init = FALSE,
                    psi_init = FALSE,
                    k_init = FALSE,
                    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
                    n_burnin = burn_in,
                    seuil = TRUE,
                    s = 0.1 
  )
  
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

result_m2_sh2_psi10_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")
