source("scripts/helper_function.R")
source("scripts/main_function.R")
set.seed(12345)
k <- 0.2
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err/k
psi_delta_vec <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
n_samples <- 1
n_iter <- 1000
burn_in <- 100
init <- c(9.8, 46.45, 0.08, 0.7, 0.5, 0.2)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)

results_all <- list()
for (i in 1:length(psi_delta_vec)) {
  g            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  h0           <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  sigma_sq_err <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  alpha        <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  psi_delta    <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  k            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  loglik1      <- matrix(NA, n_iter, n_samples, byrow = FALSE)
  delta1       <- list()
  zeta1        <- list()
  accept_rate  <- c()
  
  y_obs <- matrix(NA, length(t), n_samples, byrow = FALSE)
  
  for (v in 1:n_samples) {
    Sigma_delta <- GP_covariance(t, sigma_sq_delta, psi_delta_vec[i])
    delta <- rmvnorm(1, rep(3, n), Sigma_delta)
    y_1 = balldropg(t, c(9.8, 46.45)) + rnorm(length(t), 0, sqrt(0.1)) + delta 
    y_obs[, v] <- y_1
    a_psi = i
    b_psi = 10 - i
    results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, g_init = TRUE, h0_init = FALSE, sig2er_init = FALSE,
                          alpha_init = FALSE, psi_init = FALSE, k_init = FALSE, Sigma_theta, n_burnin = burn_in, a_psi, b_psi, seuil = TRUE, s = 0.3)
    
    g[, v] <- results$theta[, 1]
    h0[, v] <- results$theta[, 2]
    sigma_sq_err[, v] <- results$theta[, 3]
    alpha[, v] <- results$theta[, 4]
    psi_delta[, v] <- results$theta[, 5]
    k[, v] <- results$theta[, 6]
    delta1[[v]] <- results$delta
    zeta1[[v]] <- results$zeta
    loglik1[, v] <- results$loglik
    accept_rate[v] <- results$accept_rate_psi
  }
  
  results_all[[i]] <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, delta1, zeta1, loglik1, accept_rate)
}