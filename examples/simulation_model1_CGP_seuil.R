source("data/prepare_data.R")
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

# ------ Scenario (I): Baseline mixture inference without thresholding; seuil == FALSE -------------#
# ------ Scenario (II): Thresholded allocations; seuil == TRUE -------------#
# ------ Scenario (III): Effect of theta and delta confounding under thresholded allocations; seuil == TRUE, g, h0 == FALSE -------------#

set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.5
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

# discrepancy is ~0 for first ndis points, then active afterwards
ndis <- 20      
envelope_type <- "step"

make_envelope <- function(n, ndis, type = c("step","ramp")) {
  type <- match.arg(type)
  if (type == "step") {
    w <- c(rep(0, ndis), rep(1, n - ndis))
  } else {
    # smoothish ramp from 0 to 1 after ndis
    u <- seq(0, 1, length.out = n - ndis)
    w <- c(rep(0, ndis), u)
  }
  return(w)
}

w <- make_envelope(n, ndis, envelope_type)

n_samples <- 10
n_iter <- 10000
burn_in <- 2500
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)

# init = c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)
y_obs       <- matrix(NA, n, n_samples)

for (v in 1:n_samples) {

  y0 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  
    Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  tilde_delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  delta_true  <- w * tilde_delta
  
  y_1 <- y0 + delta_true
  y_obs[, v] <- y_1
  
  # -------- Scenario (I), (II), (III) --------
  res <- mcmc_step6(
    y = y_1, t = t, n_iter = n_iter, init = init, sigma_proposals = sigma_props,
    g_init = TRUE, 
    h0_init = TRUE,
    sig2er_init = FALSE,
    alpha_init = FALSE,
    psi_init = FALSE,
    k_init = FALSE,
    Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
    n_burnin = burn_in,
    seuil = TRUE,  
    s = 0.3       
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

# result_scenario_I <- list(
#   g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain,
#   delta_list, zeta_list, loglik_mat, accept_rate, y_obs = y_obs,
#   envelope = w
# )

# result_scenario_II <- list(
#   g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain,
#   delta_list, zeta_list, loglik_mat, accept_rate, y_obs = y_obs,
#   envelope = w
# )

result_scenario_III <- list(
  g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain,
  delta_list, zeta_list, loglik_mat, accept_rate, y_obs = y_obs,
  envelope = w
)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# set.seed(12345)
# k <- 0.1
# sim_psi_delta <- 0.5 
# sigma_sq_err <- 0.01
# sigma_sq_delta <- sigma_sq_err / k
# 
# n_samples <- 10
# n_iter <- 10000
# burn_in <- 2500
# #n_flat <- 10
# 
# init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)
# 
# g_chain     <- matrix(NA, n_iter, n_samples)
# h0_chain    <- matrix(NA, n_iter, n_samples)
# alpha_chain <- matrix(NA, n_iter, n_samples)
# psi_chain   <- matrix(NA, n_iter, n_samples)
# k_chain     <- matrix(NA, n_iter, n_samples)
# sigma_chain <- matrix(NA, n_iter, n_samples)
# zeta_list   <- vector("list", n_samples)
# delta_list  <- vector("list", n_samples)
# loglik_mat  <- matrix(NA, n_iter, n_samples)
# accept_rate <- numeric(n_samples)
# 
# y_obs <- matrix(NA, n, n_samples)
# 
# for (v in 1:n_samples) {
#   
#   #Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
#   #delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
#   
#   #delta[1:n_flat] <- 0
#   
#   #y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
#   
#   # y = balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) #+ delta
#   # #for fig 1 rep 50 times
#   # ndis <- 10
#   # n <- length(y)
#   # #y = y + (46-y)*.1# Fig 2:y = y + delta_i ; i = 45, when delta_i ~ GP(0, (sig2err/k)*exp(-(ti-tj)2/psi))
#   # y_1 = y + (46-y)*.3 * c(rep(0, ndis), rep(1, length(y) - ndis))
#   # #y_1 = y + 4 * c(rep(0, ndis), rep(1, length(y) - ndis))
#   # 
#   # y_obs[, v] <- y_1
#   y = balldropg(t,c(9.8,46.45)) + rnorm(45,0,sqrt(.01)) # for fig 1 rep 50 times
#   ndis <- 20
#   y_1 = y + (46-y)*.1 * c(rep(0, ndis), rep(1, length(y) - ndis))
#   y_obs[, v] <- y_1
#   
#   res <- mcmc_step6(y = y, t = t, n_iter = n_iter, init = init,
#     sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
#     g_init = FALSE,
#     h0_init = FALSE,
#     sig2er_init = FALSE,
#     alpha_init = FALSE,
#     psi_init = FALSE,
#     k_init = FALSE,
#     Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
#     n_burnin = burn_in,
#     seuil = FALSE,
#     s = 0.3 
#   )
#   
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
# 
# boxplot(res[["delta"]])
# colMeans(res[["zeta"]] == 2)
# colMeans(res[["zeta"]] == 1)
# 
# g <- res[["theta"]][,1]
# h0 <- res[["theta"]][,2]
# sigma <- res[["theta"]][,3]
# alpha <- res[["theta"]][,4]
# psi <- res[["theta"]][,5]
# k <- res[["theta"]][,6]
# 
# par(mfrow=c(2,3))
# boxplot(g)
# boxplot(h0)
# abline(h = 46.45)
# boxplot(sigma)
# abline(h = 0.01)
# boxplot(alpha)
# boxplot(psi)
# boxplot(k)
# abline(h = 0.1)
# 
# result_m2_sh2_psi1_simple_seuil <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
# #save(result_m2_sh2_psi1_simple,file = "/Users/negarsoleimani/Documents/phd/paper1/result_m2_sh2_psi1_simple.RData")
# 
# par(mfrow=c(1,1))
# boxplot(result_m2_sh2_psi1_simple_seuil[[7]][[1]])
# abline(h=0)
# 
# # simulation des données
# 
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
# prob_zeta_model2 <- colMeans(res$zeta == 2)
# print(prob_zeta_model2)
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
# # 
# 
# 
# 
