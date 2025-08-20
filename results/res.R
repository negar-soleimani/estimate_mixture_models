#################################################################################
#################################################################################
###################### Results for model 1 ######################################
#################################################################################
#################################################################################
set.seed(12345)
Sigma_theta <- matrix(c(0.1,0,0,0.1), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 48, 0.08, 0.7, 0.5, 0.2)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
n_samples       <- 3
burn_in         <- 100
n_iter          <- 1000
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
  y_1 = balldropg(t,c(9.8, 46.45)) + rnorm(n,0,.1)
  y_obs[,v] <- y_1
  
  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin = burn_in)
  
  g[,v] = results$theta[,1]
  h0[,v]=results$theta[,2]
  sigma_sq_err[,v]=results$theta[,3]
  alpha[,v] <- results$theta[,4]
  psi_delta[,v] <- results$theta[,5]
  k[,v] <- results$theta[,6]
}

result_m1_sh2_presentation <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs)
y_obs_m1_sh2_presentation <- y_obs

g_mean <- colMeans(g)
boxplot(g_mean)
h0_mean <- colMeans(h0)
boxplot(h0_mean)
sigma_sq_err_mean <- colMeans(sigma_sq_err)
boxplot(sigma_sq_err_mean)
alpha_mean <- colMeans(alpha)
boxplot(alpha_mean)
psi_delta_mean <- colMeans(psi_delta)
boxplot(psi_delta_mean)
k_mean <- colMeans(k)
boxplot(k_mean)


prob_zeta_model1 <- colMeans(results$zeta == 1)
print(prob_zeta_model1)

#################################################################################
#################################################################################
###################### Results for model 2 ######################################
#################################################################################
#################################################################################


#psi5#########################################
set.seed(12345)
k <- 0.2
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err/k
psi_delta_vec <-c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
n_samples <- 10
n_iter <- 1000
burn_in <- 100
#y_1 <- y[,1]
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
# FALSE= fixed parameter
# mcmc parameter (g,h), sig2err, psidelta, k, alpha
#mcmc_parameters <- c(FALSE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
# h0 = 47.7, 46.45, 46.44, 48, 46.13
init <- c(9.8, 46.45, 0.08, 0.7, 0.5, 0.2)

g            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
h0           <- matrix(NA, n_iter, n_samples, byrow = FALSE)
sigma_sq_err    <- matrix(NA, n_iter, n_samples, byrow = FALSE)
alpha        <- matrix(NA, n_iter, n_samples, byrow = FALSE)
psi_delta    <- matrix(NA, n_iter, n_samples, byrow = FALSE)
k            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
loglik1      <- matrix(NA, n_iter, n_samples, byrow = FALSE)
delta1            <- list()
zeta1            <- list()
accept_rate  <- c()

y_obs<- matrix(NA, length(t), n_samples, byrow = FALSE)
for (v in 1:n_samples) {
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, psi_delta_vec[5])
  delta<-rmvnorm(1,rep(3,n),Sigma_delta)
  y_1 = balldropg(t,c(9.8,46.45)) + rnorm(length(t),0,sqrt(0.1)) + delta 
  y_obs[,v]<- y_1
  a_psi = 4
  b_psi = 6
  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, g_init = TRUE, h0_init= FALSE, sig2er_init = FALSE,
                        alpha_init = FALSE, psi_init = FALSE, k_init = FALSE, Sigma_theta, n_burnin=burn_in, a_psi, b_psi, seuil = TRUE, s = 0.3)
  
  g[,v] = results$theta[,1]
  h0[,v]=results$theta[,2]
  sigma_sq_err[,v]=results$theta[,3]
  alpha[,v] <- results$theta[,4]
  psi_delta[,v] <- results$theta[,5]
  k[,v] <- results$theta[,6]
  delta1[[v]] <-results$delta
  zeta1[[v]] <- results$zeta
  loglik1[, v] <- results$loglik
  accept_rate[v] <- results$accept_rate_psi
}

result_m2_sh2_psi5_new <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, delta1, zeta1, loglik1, accept_rate)

g <- result_m2_sh2_psi5_new[[1]]
h <- result_m2_sh2_psi5_new[[2]]
sigma <- result_m2_sh2_psi5_new[[3]]
alpha <- result_m2_sh2_psi5_new[[4]]
psi <- result_m2_sh2_psi5_new[[5]]
k <- result_m2_sh2_psi5_new[[6]]

par(mfrow=c(2,3))
boxplot(g)
abline(h=9.8)
boxplot(h)
abline(h=46.45)
boxplot(sigma)
abline(h=0.1)
boxplot(alpha)
boxplot(psi)
abline(h=0.4)
boxplot(k)
abline(h=0.2)