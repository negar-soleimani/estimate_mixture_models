source("scripts/helper_function_CGP.R")
source("scripts/main_function_CGP.R")

set.seed(12345)
Sigma_theta <- matrix(c(0.5,0,0,0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.2)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
n_samples       <- 50
burn_in         <- 2000
n_iter          <- 10000
# FALSE= fixed parameter
# mcmc parameter (g,h), sig2err, psidelta, k, alpha, freeze delta-zeta
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

g            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
h0           <- matrix(NA, n_iter, n_samples, byrow = FALSE)
sigma_sq_err <- matrix(NA, n_iter, n_samples, byrow = FALSE)
alpha        <- matrix(NA, n_iter, n_samples, byrow = FALSE)
psi_delta    <- matrix(NA, n_iter, n_samples, byrow = FALSE)
k            <- matrix(NA, n_iter, n_samples, byrow = FALSE)


y_obs <- matrix(NA, length(t), n_samples, byrow = FALSE)
for (v in 1:n_samples) {
  y_1 = balldropg(t,c(9.8, 46.45)) + rnorm(n, 0, sqrt(0.01))
  y_obs[,v] <- y_1
  
  #results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, g_init=FALSE, h0_init= FALSE, sig2er_init = FALSE,
  #                      alpha_init = FALSE, psi_init = FALSE, k_init = FALSE, Sigma_theta, n_burnin = burn_in, a_psi, b_psi, seuil = FALSE, s = 0.3)
  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=2000)
  g[,v] = results$theta[,1]
  h0[,v]=results$theta[,2]
  sigma_sq_err[,v]=results$theta[,3]
  alpha[,v] <- results$theta[,4]
  psi_delta[,v] <- results$theta[,5]
  k[,v] <- results$theta[,6]
}

result_m1_sh2 <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs)
y_obs_m1_sh2 <- y_obs
