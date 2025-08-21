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
# Define a function to plot boxplots for each parameter
plot_boxplot <- function(results_all, param_index, ylab, abline_value, ylim_value = NULL) {
  results_list <- lapply(1:10, function(i) results_all[[i]][[param_index]])
  
  par(mfrow = c(1,1), mai = c(1,1,0.5,0.5))
  boxplot(lapply(results_list, colMeans), 
          names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
          xlab = expression(psi[delta]),
          ylab = ylab,
          ylim = ylim_value,
          cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
  abline(h = abline_value)
}

# Plot for g
plot_boxplot(results_all, 1, expression(g), 9.8)

# Plot for h0
plot_boxplot(results_all, 2, expression(h), 46.45, ylim_value = c(45, 50))

# Plot for sigma_sq_err
plot_boxplot(results_all, 3, expression(sigma[err]^2), 0.1)

# Plot for alpha
plot_boxplot(results_all, 4, expression(alpha), 0.4, ylim_value = c(0, 1))

# Plot for psi_delta
plot_boxplot(results_all, 5, expression(psi[delta]), 0.4)

# Plot for k
plot_boxplot(results_all, 6, expression(k), 0.2)

