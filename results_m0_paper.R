# result - model 1:
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

# ----------------------------- data -------------------------------------#
don <- read_xlsx("/Users/negarsoleimani/Documents/phd/paper1/Ball_drops_data.xlsx", sheet = 2)
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop)
don <- don[don$drop == 1, ]

t <- don$time
# dt_est <- median(diff(don$time))
# t <- seq(from = don$time[1],
#              by   = dt_est,
#              length.out = 100)
# 
# t

# t_raw  <- don$time          # 0 to ~2.903
# t_min  <- min(t_raw)
# t_max  <- max(t_raw)
# 
# make_t <- function(n) {
#   seq(from = t_min, to = t_max, length.out = n)
# }
# 
# t <- make_t(100) 

#t <- seq(0, 1, 0.01)

y <- don$Height
n <- length(t)
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)

# ----------------------------- Helper functions ------------------------------
GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
}

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

# =============================================================================
# MCMC
# =============================================================================

mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000) {
  # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T", freeze_delta_zeta = "F")
  #use_discrepancy <- mcmc_parameters[6]
  freeze_delta_zeta <- mcmc_parameters[6]
 
  total_iter <- n_burnin + n_iter
  
  theta <- init
  delta <- rep(0, length(y))
  chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
  colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
  chain_zeta <- matrix(NA, nrow = total_iter, ncol = length(y))
  
  loglik_chain <- numeric(total_iter)
  accept_psi <- 0
  
  for (iter in 1:total_iter) {
    g <- theta[1]; h0 <- theta[2]; sigma_sq_err <- theta[3]
    alpha_param <- theta[4]; psi_delta <- theta[5]; k <- theta[6]
    sigma_sq_delta <- sigma_sq_err / k
    
    Sigma_delta <- GP_covariance(t, sigma_sq_delta, psi_delta)
    f_theta <- balldropg(t, c(g, h0))
    mean1 <- f_theta
    mean2 <- f_theta + delta
    
    #-------------------------------- probability of the zeta --------------------------------# 
    
    # Method1 for calculate the probability of the zeta:
    # prob_zeta <- 1 / (1 + (alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err))) /
    #                     ((1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))) #P(M2) = ((1-alpha)*M2)/((alpha*M1)+(1-alpha)*M2)
    # zeta <- 1 + (runif(length(y)) < prob_zeta)
    # chain_zeta[iter, ] = zeta
    
    #I don't use Method1 for computing prob_zeta,because this direct probability 
    #computation is numerically unstable, when the parameter sigma_sq_err is small 
    #(sigma_sq_err = 0.01), it leads to overflow and produces NaN/NA probabilities, 
    #and these NA values then propagate into zeta and later cause singular matrices 
    #and when computing solve(A), leading to errors.
    
    # Method2 for calculate the probability of the zeta:(log_sum_exp: https://rpubs.com/FJRubio/LSE and https://en.wikipedia.org/wiki/LogSumExp)
    #if (use_discrepancy) {
    log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
    log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    #prob_zeta_1 <- exp(log_w1 - log_den)  #  P(zeta=1|y)
    prob_zeta_2 <- exp(log_w2 - log_den)  #  P(zeta=2|y)
    #zeta <- ifelse(runif(length(y)) < prob_zeta_1, 1, 2)
    zeta <- ifelse(runif(length(y)) < prob_zeta_2, 2, 1) #= zeta <- 1 + (runif(length(y)) < prob_zeta_1) # sample zeta: 2 with prob post_p2, otherwise 1
    #} else {
      
    #  zeta <- rep(1, length(y))   # MODEL WITHOUT DISCREPANCY
      
    #}
    #--------------------
    if (freeze_delta_zeta) {
      zeta <- rep(1, length(y))
    }
    #--------------------
    
    chain_zeta[iter, ] = zeta
    
    #-------------------------------- log_likelihood --------------------------------#   
    
    log_likelihood <- sum(log(ifelse(zeta == 1,
                                     alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
                                     (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    loglik_chain[iter] <- log_likelihood
    
    #-------------------------------- delta(discrepancy) --------------------------------#  
    #if (use_discrepancy) {
    zeta_2_indices <- which(zeta == 2)
    if (length(zeta_2_indices) > 0) {
      y_m <- y[zeta_2_indices]
      Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
        Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
      Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
      Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
      mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
      Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
      Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
      delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    } else {
      delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))}
    #} else {
    #  delta <- rep(0, n)  
    #}
    #--------------------
    if (freeze_delta_zeta) {
      delta <- rep(0, n)
    }
    #--------------------
    #-------------------------------- Gibbs step for theta --------------------------------# 
    
    # # prior Normal distribution
    # zeta_1_indices <- which(zeta == 1)
    # zeta_2_indices <- which(zeta == 2)
    # X <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
    # x1 <- X[zeta_1_indices, , drop = FALSE]
    # x2 <- X[zeta_2_indices, , drop = FALSE]
    # theta_hat     <- matrix(c(46.46, 9.8), ncol = 1)
    # inv_sigma_theta <- solve(Sigma_theta)
    # A <- ((t(x1) %*% x1) / sigma_sq_err) + ((t(x2) %*% x2) / sigma_sq_err) + inv_sigma_theta
    # Sigmapost_theta <- solve(A)
    # y1 <- matrix(y[zeta_1_indices], ncol = 1)
    # y2 <- matrix(y[zeta_2_indices], ncol = 1)
    # d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    # B <- (t(x1) %*% y1) / sigma_sq_err +
    #   (t(x2) %*% y2) / sigma_sq_err -
    #   (t(x2) %*% d2) / sigma_sq_err +
    #   inv_sigma_theta %*% theta_hat
    # Mupost_theta <- Sigmapost_theta %*% B
    # theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
    # h0 <- theta_sample[1];  g <- theta_sample[2]
    # theta[1] <- g
    # theta[2] <- h0
    
    #-------------------------------- Page 22 - part 8.1 --------------------------------#     
    # flat prior on theta
    zeta_1_indices <- which(zeta == 1)
    zeta_2_indices <- which(zeta == 2)

    X <- cbind(1, -0.5 * (t * t_range)^2)
    x1 <- X[zeta_1_indices, , drop = FALSE]
    x2 <- X[zeta_2_indices, , drop = FALSE]

    y1 <- matrix(y[zeta_1_indices], ncol = 1)
    y2 <- matrix(y[zeta_2_indices], ncol = 1)
    d2 <- matrix(delta[zeta_2_indices], ncol = 1)

    A <- (t(x1) %*% x1 + t(x2) %*% x2) / sigma_sq_err
    B <- (t(x1) %*% y1 + t(x2) %*% y2 - t(x2) %*% d2) / sigma_sq_err

    Sigmapost_theta <- solve(A)
    Mupost_theta    <- Sigmapost_theta %*% B

    theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)

    h0 <- theta_sample[1];  g <- theta_sample[2]
    theta[1] <- g
    theta[2] <- h0


    if(mcmc_parameters[1] == FALSE){
      g <- init[1]
      h0 <- init[2]
    }
    
    # # g = fixer
    # if (mcmc_parameters[1] == FALSE) {
    #   g <- init[1]
    #   h0 <- theta_sample[1]
    # }
    
    # h0 = fixed
    # if (mcmc_parameters[1] == FALSE) {
    #  h0 <- init[2]
    #  g <- theta_sample[2]
    # }
    
    #-------------------------------- Gibbs step for sigma_sq_err(lambda^2) --------------------------------#   
    #-------------------------------- Page 23 - part 8.2 --------------------------------#     
    
    R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
    if(n > 0){
      R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
      
      quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
    } 
    else {
      quad_form_delta <- 0
    }
    
    # #if (use_discrepancy) {
    # if (freeze_delta_zeta == FALSE) {
    # f_theta <- balldropg(t, c(g, h0))
    # idx1 <- which(zeta == 1)
    # idx2 <- which(zeta == 2)
    # n1   <- length(idx1)
    # n2   <- length(idx2)
    # residual1 <- y[idx1] - f_theta[idx1]
    # rss1   <- sum(residual1^2)
    # 
    # residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
    # rss2   <- sum(residual2^2)
    # d <- 2
    # # When we consider the Jeffreys prior, I use the following code:
    # rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    # shape_err <- n + (d/2)
    # # When the prior is the sigma parameter of the "inverse gamma distribution", I use the following code:
    # # rate_err <- 0.05 + (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    # # shape_err <- 3 + n
    # # sigma_sq_err <- rinvgamma(1, shape = shape_err, scale = rate_err)
    # } else {
    #   #residual <- y - f_theta
    #   #rss <- sum(residual^2)
    #   f_theta <- balldropg(t, c(g, h0))
    #   idx1 <- which(zeta == 1)
    #   idx2 <- which(zeta == 2)
    #   n1   <- length(idx1)
    #   n2   <- length(idx2)
    #   residual1 <- y[idx1] - f_theta[idx1]
    #   rss1   <- sum(residual1^2)
    #   shape_err <- (n / 2) + 1
    #   rate_err  <- 0.5 * rss1
    # }
    # #sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    # sigma_sq_err <- 1/rgamma(1, shape = shape_err, rate = rate_err)
    # theta[3] <- sigma_sq_err
    
    if (freeze_delta_zeta){
      f_theta <- balldropg(t, c(g, h0))
      idx1 <- which(zeta == 1)
      idx2 <- which(zeta == 2)
      n1   <- length(idx1)
      n2   <- length(idx2)
      residual1 <- y[idx1] - f_theta[idx1]
      rss1   <- sum(residual1^2)
      shape_err <- (n / 2) + 1
      rate_err  <- 0.5 * rss1
    } else {
      f_theta <- balldropg(t, c(g, h0))
      idx1 <- which(zeta == 1)
      idx2 <- which(zeta == 2)
      n1   <- length(idx1)
      n2   <- length(idx2)
      residual1 <- y[idx1] - f_theta[idx1]
      rss1   <- sum(residual1^2)
      
      residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
      rss2   <- sum(residual2^2)
      d <- 2
      # When we consider the Jeffreys prior, I use the following code:
      rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
      shape_err <- n + (d/2)
    }
    sigma_sq_err <- 1/rgamma(1, shape = shape_err, rate = rate_err)
    theta[3] <- sigma_sq_err
    
    if(mcmc_parameters[2] == FALSE){
      sigma_sq_err <- init[3]
      theta[3] <- sigma_sq_err
    }
    
    #-------------------------------- Gibbs step for psi_delta --------------------------------#   
    
    # When the prior is the psi_delta of the "uniform(0.1, 1)", I use the following code:
    psi_prop <- rtruncnorm(1, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    sigma_sq_delta_prop <- sigma_sq_err / k
    Sigma_delta <- GP_covariance(t, sigma_sq_delta_prop, psi_delta)
    Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
    log_prop_current <- log(dtruncnorm(psi_delta, a = 0.1, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
    log_prop_prop <- log(dtruncnorm(psi_prop, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
    log_prior_current <- dunif(psi_delta, min = 0.1, max = 1, log = TRUE)
    log_prior_prop    <- dunif(psi_prop,  min = 0.1, max = 1, log = TRUE)
    log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
    log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + (log_prop_current - log_prop_prop)
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      psi_delta <- psi_prop
      theta[5] <- psi_delta
      accept_psi <- accept_psi + 1
    }
    
    # When the prior is the psi_delta of the "Beta(8, 4)", I use the following code:
    
    # psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    # sigma_sq_delta_prop <- sigma_sq_err / k
    # Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
    # #log_prior_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = 0.5, sd = 0.2))
    # #log_prior_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = 0.5, sd = 0.2))
    # log_prior_current <- dbeta(psi_delta, shape1 = 8, shape2 = 4, log=TRUE)
    # log_prior_prop <- dbeta(psi_prop, shape1 = 8, shape2 = 4, log=TRUE)
    # log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
    # log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    # log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current)
    # if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
    #   psi_delta <- psi_prop
    #   accept_psi <- accept_psi + 1
    # }
    
    if(mcmc_parameters[3] == FALSE){
      psi_delta <- init[5]
      theta[5] <- psi_delta
    }
    
    #-------------------------------- Gibbs step for k --------------------------------#   
    #-------------------------------- Page 24-part 8.3 --------------------------------#     
    
    alpha_k <- (n / 2) + 1
    beta_k <-  (1 / (2 * sigma_sq_err)) * quad_form_delta
    for (try_k in 1:100) {
      k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
      #if (k_prop >= 0.1 && k_prop <= 0.9) {
      if (k_prop > 0 && k_prop < 1) {
        k <- k_prop
        theta[6] <- k
        break
      }
    }
    
    if(mcmc_parameters[4] == FALSE){
      k <- init[6]
      theta[6] <- k
    }
    
    #-------------------------------- Gibbs step for alpha --------------------------------#   
    #-------------------------------- Page 25 - part 8.4 --------------------------------#     
    
    alpha_param <- rbeta(1, sum(zeta == 1) + 7, sum(zeta == 2) + 1)
    theta[4] <- alpha_param
    
    if(mcmc_parameters[5] == FALSE){
      alpha_param <- init[4]
      theta[4] <- alpha_param
    }
    
    #if (!use_discrepancy) {
    #  alpha_param <- 1
    #  k <- 0.2
    #  psi_delta <- 0.5
    #}
    
    theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
    chain_theta[iter, ] <- theta
    chain_delta[iter, ] <- delta
    
  }
  
  #----------------------------------------------------------------#   
  
  # Remove burn-in samples if n_burnin > 0
  if (n_burnin > 0) {
    keep_indices <- (n_burnin + 1):total_iter
    chain_theta <- chain_theta[keep_indices, , drop = FALSE]
    chain_delta <- chain_delta[keep_indices, , drop = FALSE]
    chain_zeta <- chain_zeta[keep_indices, , drop = FALSE]
    loglik_chain <- loglik_chain[keep_indices]
    
    # Adjust acceptance rate calculation for post-burn-in period only
    cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
  } else {
    cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
  }
  return(list(theta = chain_theta, delta = chain_delta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_psi / n_iter))
}


set.seed(123)
Sigma_theta <- matrix(c(0.5,0,0,0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.4, 0.2)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
n_samples       <- 30
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
loglik_mat  <- matrix(NA, n_iter, n_samples)
delta_list  <- vector("list", n_samples)
zeta_list   <- vector("list", n_samples)
accept_rate <- numeric(n_samples)


y_obs <- matrix(NA, length(t), n_samples, byrow = FALSE)
for (v in 1:n_samples) {
  y_1 = balldropg(t,c(9.8, 46.45)) + rnorm(n, 0, sqrt(0.01))
  y_obs[,v] <- y_1
  
  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=2000)
  g[,v] = results$theta[,1]
  h0[,v]=results$theta[,2]
  sigma_sq_err[,v]=results$theta[,3]
  alpha[,v] <- results$theta[,4]
  psi_delta[,v] <- results$theta[,5]
  k[,v] <- results$theta[,6]
  delta_list[[v]]  <- results$delta
  zeta_list[[v]]   <- results$zeta
  loglik_mat[, v]  <- results$loglik
  accept_rate[v]   <- results$accept_rate_psi
}
result_m0_sh2_classic_classic <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs, delta_list, zeta_list, loglik_mat, accept_rate)
# y_obs_m1_sh2_ex <- y_obs
# save(result_m0_sh2_classic_classic,file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m0_sh2_classic_classic.RData")
# load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m0_sh2_classic_classic.RData")
View(result_m0_sh2_classic_classic)

g_sh2 <- result_m0_sh2_classic_classic[[1]]
h0_sh2 <- result_m0_sh2_classic_classic[[2]]
sigma_sq_err_sh2 <- result_m0_sh2_classic_classic[[3]]
alpha_sh2 <- result_m0_sh2_classic_classic[[4]]
psi_delta_sh2 <- result_m0_sh2_classic_classic[[5]]
k_sh2 <- result_m0_sh2_classic_classic[[6]]
zeta <- result_m0_sh2_classic_classic[[9]][[30]]
  
mean(colMeans(sigma_sq_err_sh2))
par(mfrow = c(1, 3),
    mar   = c(3, 4, 1, 1) 
)
boxplot(
  colMeans(g_sh2),
  ylab = "g",
  #xlab = "Bowling Ball",
  col  = "lightseagreen",
  main = ""
)
abline(h = 9.8, lty = 2)

boxplot(
  colMeans(h0_sh2),
  ylab = "h0",
  #xlab = "",
  col  = "lightseagreen",
  main = ""
)
abline(h = 46.45045, lty = 2)

boxplot(
  colMeans(sigma_sq_err_sh2),
  ylab = expression(lambda^2),
  #ylim = c(0.01, 0.0452),
  #xlab = "",
  col  = "lightseagreen",
  main = ""
)
abline(h = 0.01, lty = 2)

boxplot(
  colMeans(psi_delta_sh2),
  ylab = expression(lambda^2),
  #ylim = c(0.01, 0.0452),
  #xlab = "",
  col  = "lightseagreen",
  main = ""
)
abline(h = 0.4, lty = 2)

alpha <- result_m0_sh2_classic_classic[[4]]

library(ggplot2)

# 1) flatten the matrix to a single vector
alpha_vec <- as.vector(alpha)

df_alpha <- data.frame(alpha = alpha_vec)

ggplot(df_alpha, aes(x = alpha)) +
  geom_density(
    fill  = "#4CCDC9",  
    color = "lightseagreen", 
    alpha = 0.6,        
    size  = 1
  ) +
  labs(
    x     = expression(alpha),
    y     = "Density",
    title = "Posterior Density of alpha"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    plot.margin  = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

#result jeffry
# result_m1_sh2 <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs)
# y_obs_m1_sh2 <- y_obs
# save(result_m1_sh2,file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m1_sh2.RData")
# load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m1_sh2.RData")

# result Inverse gamma:
# result_m1_sh2_IG <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs)
# y_obs_m1_sh2_IG <- y_obs
# save(result_m1_sh2_IG,file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m1_sh2_IG.RData")
# load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m1_sh2_IG.RData")

#################################################################################
#################################################################################
###################### Results for model 1 ######################################
#################################################################################
#################################################################################
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh2_presentation.RData")
load("E:Phd_Paris Saclay/resultat/Results_M1_Simulation/y_obs_m1_sh2_presentation.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh1.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh2.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh3.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh4.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh5.RData")

sheet_names <- c("Baseball",
                 "Blue Basketball",
                 "Green Basketball",
                 "Volleyball",
                 "Bowling Ball")

#View(result_m1_sh1)

g_sh1 <- result_m1_sh1[[1]]
h0_sh1 <- result_m1_sh1[[2]]
sigma_sq_err_sh1 <- result_m1_sh1[[3]]
alpha_sh1 <- result_m1_sh1[[4]]
psi_delta_sh1 <- result_m1_sh1[[5]]
k_sh1 <- result_m1_sh1[[6]]

g_sh2 <- result_m1_sh2_IG[[1]]
h0_sh2 <- result_m1_sh2_IG[[2]]
sigma_sq_err_sh2 <- result_m1_sh2_IG[[3]]
alpha_sh2 <- result_m1_sh2_IG[[4]]
psi_delta_sh2 <- result_m1_sh2_IG[[5]]
k_sh2 <- result_m1_sh2_IG[[6]]

g_sh3 <- result_m1_sh3[[1]]
h0_sh3 <- result_m1_sh3[[2]]
sigma_sq_err_sh3 <- result_m1_sh3[[3]]
alpha_sh3 <- result_m1_sh3[[4]]
psi_delta_sh3 <- result_m1_sh3[[5]]
k_sh3 <- result_m1_sh3[[6]]

g_sh4 <- result_m1_sh4[[1]]
h0_sh4 <- result_m1_sh4[[2]]
sigma_sq_err_sh4 <- result_m1_sh4[[3]]
alpha_sh4 <- result_m1_sh4[[4]]
psi_delta_sh4 <- result_m1_sh4[[5]]
k_sh4 <- result_m1_sh4[[6]]

g_sh5 <- result_m1_sh5[[1]]
h0_sh5 <- result_m1_sh5[[2]]
sigma_sq_err_sh5 <- result_m1_sh5[[3]]
alpha_sh5 <- result_m1_sh5[[4]]
psi_delta_sh5 <- result_m1_sh5[[5]]
k_sh5 <- result_m1_sh5[[6]]

g_sh2 <- result_m1_sh2_presentation[[1]]
h0_sh2 <- result_m1_sh2_presentation[[2]]
sigma_sq_err_sh2 <- result_m1_sh2_presentation[[3]]
alpha_sh2 <- result_m1_sh2_presentation[[4]]
psi_delta_sh2 <- result_m1_sh2_presentation[[5]]
k_sh2 <- result_m1_sh2_presentation[[6]]

par(mfrow = c(1, 3),
    mar   = c(3, 4, 1, 1) 
)
boxplot(
  colMeans(g_sh2),
  ylab = "g",
  #xlab = "Bowling Ball",
  col  = "lightseagreen",
  main = NULL
)
abline(h = 9.8, lty = 2)

boxplot(
  colMeans(h0_sh2),
  ylab = "h0",
  #xlab = "",
  col  = "lightseagreen",
  main = NULL
)
abline(h = 46.45045, lty = 2)

boxplot(
  colMeans(sigma_sq_err_sh2),
  ylab = expression(lambda^2),
  #ylim = c(0.01, 0.0452),
  #xlab = "",
  col  = "lightseagreen",
  main = NULL
)
abline(h = 0.01, lty = 2)


boxplot(
  colMeans(alpha_sh2),
  ylab = expression(alpha),
  xlab = "alpha",
  col  = "#800020",
  main = NULL
)
abline(h = 0.7, lty = 2)

boxplot(
  colMeans(psi_delta_sh2),
  ylab = expression(psi[delta]),
  xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 0.3, lty = 2)

boxplot(
  colMeans(k_sh2),
  ylab = "k",
  xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 0.2, lty = 2)

don <- read_xlsx("E:/Phd_Paris Saclay/10Sep2024_Cours Calibration/Ball_drops_data.xls", sheet = 2)
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop)
don <- don[don$drop == 1, ]

t <- don$time
y <- don$Height
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

g_mat <- result_m0_sh2_classic_classic[[1]]   # n_iter × n_sims
h0_mat<- result_m0_sh2_classic_classic[[2]]

g_last <- g_mat[, 3]
h0_last <- h0_mat[, 3]
y_obs <- result_m0_sh2_classic_classic$y_obs
y_obs30 <- y_obs[,3]
y_true <- balldropg(t, c(9.8, 46.45))

y_pred <- matrix(NA, length(t), length(g_last), byrow = FALSE)

for (i in 1:length(g_last)) {
  theta <- c(g_last[i], h0_last[i])
  y_pred[,i] <- balldropg(t, theta)
}

y_pred1  = t(y_pred)
par(mfrow = c(1,1))
boxplot(y_pred1, col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 1)
legend(x=30, y=45, legend=c("Simulated data",
                            "True code", "Predictions"), lwd=rep(2,2), col=c("gold","blue", "orange2"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))

boxplot(y_pred1[, 1:5], col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 3)
legend(x=4.2, y=46.55, legend=c("Simulated data", "Predictions",
                            "True code"), lwd=rep(2,2), col=c("gold", "orange2","blue"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))

y_pred2 <- rowMeans(y_pred) 
y_obs2 <- rowMeans(y_obs)

df_cmp <- data.frame(
  time   = t,
  y_obs2  = y_obs2,
  y_pred2  = y_pred2,
  y_true = y_true
)
library(ggplot2)
p <- ggplot(df_cmp, aes(x = time)) +
  geom_point(aes(y = y_obs2, color = "Observed"),      size = 1) +
  geom_line( aes(y = y_pred2,  color = "Estimate"),     linetype = "dashed", size = 0.5) +
  geom_line( aes(y = y_true, color = "True"),         linetype = "solid",  size = 0.5) +
  scale_color_manual(
    name   = NULL,
    values = c("Observed" = "#800020",
               "Estimate" = "green",
               "True"     = "blue")
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x     = "Time",
    y     = ""
  ) +
  theme(
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.6)),
    legend.key           = element_rect(fill = NA)
  )

print(p)

boxplot(colMeans(g_sh1), colMeans(g_sh2), colMeans(g_sh3), colMeans(g_sh4), colMeans(g_sh5),
        names  = sheet_names,
        ylab = "g",
        col    = rainbow(5))

boxplot(colMeans(h0_sh1), colMeans(h0_sh2), colMeans(h0_sh3), colMeans(h0_sh4), colMeans(h0_sh5),
        names  = sheet_names,
        ylab = "h0",
        col    = rainbow(5))

boxplot(colMeans(sigma_sq_err_sh1), colMeans(sigma_sq_err_sh2), colMeans(sigma_sq_err_sh3), 
        colMeans(sigma_sq_err_sh4), colMeans(sigma_sq_err_sh5),
        names  = sheet_names,
        ylab = "sigma_sq_err",
        col    = rainbow(5))

boxplot(colMeans(alpha_sh1), colMeans(alpha_sh2), colMeans(alpha_sh3), colMeans(alpha_sh4), 
        colMeans(alpha_sh5),
        names  = sheet_names,
        ylab = "alpha",
        ylim = c(0.4,1),
        col    = rainbow(5))

boxplot(colMeans(psi_delta_sh1), colMeans(psi_delta_sh2), colMeans(psi_delta_sh3), 
        colMeans(psi_delta_sh4), colMeans(psi_delta_sh5),
        names  = sheet_names,
        ylab = "psi_delta",
        col    = rainbow(5))

boxplot(colMeans(k_sh1), colMeans(k_sh2), colMeans(k_sh3), colMeans(k_sh4), colMeans(k_sh5),
        names  = sheet_names,
        ylab = "k",
        col    = rainbow(5))

load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh2_presentation.RData")
load("E:Phd_Paris Saclay/resultat/Results_M1_Simulation/y_obs_m1_sh2_presentation.RData")

alpha <- result_m1_sh2_IG[[4]]
hist(alpha)
plot(density(alpha))

library(ggplot2)

# 1) flatten the matrix to a single vector
alpha_vec <- as.vector(alpha)

df_alpha <- data.frame(alpha = alpha_vec)

ggplot(df_alpha, aes(x = alpha)) +
  geom_density(
    fill  = "#4CCDC9",  
    color = "lightseagreen", 
    alpha = 0.6,        
    size  = 1
  ) +
  labs(
    x     = expression(alpha),
    y     = "Density",
    title = "Posterior Density of " ~ alpha
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    plot.margin  = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

t <- seq(0,1,0.01)

#============================================================================#
#model 1 - result - OGP
#============================================================================#
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

don <- read_xlsx("/Users/negarsoleimani/Documents/phd/paper1/Ball_drops_data.xlsx", sheet = 2)
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

GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
}

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}


J_x <- function(x, psi_delta) {
  psi <- psi_delta
  psi * ( 2 - exp(-x/psi) - exp(-(1 - x)/psi) )
}

I_x <- function(x, psi_delta) {
  psi <- psi_delta
  A1 <- function(b) { psi * (1 - exp(-b/psi)) }
  A2 <- function(b) { psi^2 * (1 - exp(-b/psi)) - psi * b * exp(-b/psi) }
  A3 <- function(b) { 2*psi^3 - exp(-b/psi) * (b^2*psi + 2*b*psi^2 + 2*psi^3) }
  
  term1 <- x^2 * ( A1(x) + A1(1 - x) )
  term2 <- 2*x * ( -A2(x) + A2(1 - x) )
  term3 <- A3(x) + A3(1 - x)
  
  term1 + term2 + term3
}

A_scalar <- function(psi_delta) {
  psi <- psi_delta
  2*psi - ((2*psi^2) * (1 - exp(-1/psi)))
}

B_scalar <- function(psi_delta) {
  psi <- psi_delta
  E <- exp(-1/psi)
  ((2/3)*psi) - (psi^2) + (2*psi^3) - (4*psi^4) + E*(psi^2 + 2*psi^3 + 4*psi^4)
}

D_scalar <- function(psi_delta) {
  psi <- psi_delta
  E <- exp(-1/psi)
  ((2/5)*psi) - (psi^2) + ((4/3)*psi^3) - (8*(psi^6)) + E*(4*(psi^4) + 8*(psi^5) + 8*(psi^6))
}

h_vec_x <- function(x, psi_delta) {
  psi <- psi_delta
  a   <- t_range 
  J <- J_x(x, psi)
  I <- I_x(x, psi)
  rbind(J, -0.5 * (a^2) * I)
}

H_matrix <- function(psi_delta) {
  A <- A_scalar(psi_delta); B <- B_scalar(psi_delta); D <- D_scalar(psi_delta)
  H <- matrix(c(A, -((a^2)*B)/2, -((a^2)*B)/2, ((a^4)*D)/4), nrow = 2, byrow = TRUE)
  H <- 0.5 * (H + t(H))
  H + diag(1e-10, 2)
}


GP_correlation <- function(t, psi_delta) {
  R_se <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
  
  h_mat <- t(sapply(t, function(x) h_vec_x(x, psi_delta)))
  H <- H_matrix(psi_delta)
  H_inv <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(H_inv)) return(NULL)
  
  K_star <- R_se - h_mat %*% H_inv %*% t(h_mat)
  #K_star <- 0.5 * (K_star + t(K_star))  # symmetrize
  K_star
}

GP_covariance_star_complete <- function(t, sigma_sq_err, k, psi_delta) {
  K_star <- GP_correlation(t, psi_delta)
  if (is.null(K_star)) {
    return(NULL)
  }
  (sigma_sq_err / k) * K_star
}

mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000, a_psi, b_psi) {
  # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T")
  freeze_delta_zeta <- mcmc_parameters[6]
  
  total_iter <- n_burnin + n_iter
  n <- length(y)
  
  theta <- init
  delta <- rep(0, length(y))
  
  chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
  colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
  chain_zeta <- matrix(NA, nrow = total_iter, ncol = length(y))
  
  loglik_chain <- numeric(total_iter)
  accept_psi <- 0
  
  for (iter in 1:total_iter) {
    g <- theta[1]; h0 <- theta[2]; sigma_sq_err <- theta[3]
    alpha_param <- theta[4]; psi_delta <- theta[5]; k <- theta[6]
    sigma_sq_delta <- sigma_sq_err / k
    
    #-------------------------------- probability of the zeta --------------------------------# 
    
    Sigma_delta <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta) #c*(x,x')
    #Sigma_delta <- 0.5 * (Sigma_delta + t(Sigma_delta)) + diag(1e-10, n)
    
    f_theta <- balldropg(t, c(g, h0))
    mean1 <- f_theta
    mean2 <- f_theta + delta
    
    # ##s <- 0.3
    # prob_zeta <- 1 / (1 + (alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err))) /
    #                     ((1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err))))
    # ##*(abs(delta) > s)
    # zeta <- 1 + (runif(length(y)) < prob_zeta)
    # 
    # chain_zeta[iter, ] = zeta
    
    
    # Method2 for calculate the probability of the zeta:(log_sum_exp: https://rpubs.com/FJRubio/LSE and https://en.wikipedia.org/wiki/LogSumExp)
    log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
    log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    #prob_zeta_1 <- exp(log_w1 - log_den)  #  P(zeta=1|y)
    prob_zeta_2 <- exp(log_w2 - log_den)  #  P(zeta=2|y)
    #zeta <- ifelse(runif(length(y)) < prob_zeta_1, 1, 2)
    zeta <- ifelse(runif(length(y)) < prob_zeta_2, 2, 1) #= zeta <- 1 + (runif(length(y)) < prob_zeta_1) # sample zeta: 2 with prob post_p2, otherwise 1
    
    #--------------------
    if (freeze_delta_zeta) {
      zeta <- rep(1, length(y))
    }
    #--------------------
    
    chain_zeta[iter, ] = zeta
    
    #-------------------------------- log_likelihood --------------------------------#   
    
    log_likelihood <- sum(log(ifelse(zeta == 1,
                                     alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
                                     (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    loglik_chain[iter] <- log_likelihood
    
    #-------------------------------- delta(discrepancy) --------------------------------#  
    
    zeta_2_indices <- which(zeta == 2)
    if (length(zeta_2_indices) > 0) {
      y_m <- y[zeta_2_indices]
      Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
        Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
      Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
      Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
      mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
      Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
      Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
      delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    } else {
      delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
    }
    
    #--------------------
    if (freeze_delta_zeta) {
      delta <- rep(0, n)
    }
    #--------------------
    #-------------------------------- Gibbs step for theta --------------------------------#
    
    # # prior Normal distribution
    # zeta_1_indices <- which(zeta == 1)
    # zeta_2_indices <- which(zeta == 2)
    # X <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
    # x1 <- X[zeta_1_indices, , drop = FALSE]
    # x2 <- X[zeta_2_indices, , drop = FALSE]
    # theta_hat     <- matrix(c(46.46, 9.8), ncol = 1)
    # inv_sigma_theta <- solve(Sigma_theta)
    # A <- ((t(x1) %*% x1) / sigma_sq_err) + ((t(x2) %*% x2) / sigma_sq_err) + inv_sigma_theta
    # Sigmapost_theta <- solve(A)
    # y1 <- matrix(y[zeta_1_indices], ncol = 1)
    # y2 <- matrix(y[zeta_2_indices], ncol = 1)
    # d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    # B <- (t(x1) %*% y1) / sigma_sq_err +
    #   (t(x2) %*% y2) / sigma_sq_err -
    #   (t(x2) %*% d2) / sigma_sq_err +
    #   inv_sigma_theta %*% theta_hat
    # Mupost_theta <- Sigmapost_theta %*% B
    # theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
    # h0 <- theta_sample[1];  g <- theta_sample[2]
    # theta[1] <- g
    # theta[2] <- h0
    # 
    
    # flat prior on theta 
    zeta_1_indices <- which(zeta == 1)
    zeta_2_indices <- which(zeta == 2)
    
    X <- cbind(1, -0.5 * (t * t_range)^2)
    x1 <- X[zeta_1_indices, , drop = FALSE]
    x2 <- X[zeta_2_indices, , drop = FALSE]
    
    y1 <- matrix(y[zeta_1_indices], ncol = 1)
    y2 <- matrix(y[zeta_2_indices], ncol = 1)
    d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    
    A <- (t(x1) %*% x1 + t(x2) %*% x2) / sigma_sq_err
    B <- (t(x1) %*% y1 + t(x2) %*% y2 - t(x2) %*% d2) / sigma_sq_err
    
    Sigmapost_theta <- solve(A)
    Mupost_theta    <- Sigmapost_theta %*% B
    
    theta_sample <- rmvnorm(1, mean = Mupost_theta,
                            sigma = Sigmapost_theta)
    
    h0 <- theta_sample[1];  g <- theta_sample[2]
    theta[1] <- g
    theta[2] <- h0
    
    
    if(mcmc_parameters[1] == FALSE){
      g <- init[1]
      h0 <- init[2]
    }
    ## g = fixer
    # if (mcmc_parameters[1] == FALSE) {
    #  g <- init[1]
    #  h0 <- theta_sample[1]
    # }
    # 
    ## h0 = fixed
    #if (mcmc_parameters[1] == FALSE) {
    #  h0 <- init[2] 
    #  g <- theta_sample[2] 
    #}
    
    #-------------------------------- Gibbs step for sigma_sq_err --------------------------------#   
    
    K_star_psi <- GP_correlation(t, psi_delta)
    inv_Kstar_psi <- tryCatch(chol2inv(chol(K_star_psi)), error = function(e) NULL)
    
    if (is.null(inv_Kstar_psi)) {
      quad_form_delta <- -Inf
    } else {
      quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
    }
    
    if (quad_form_delta == -Inf) next
    
    if (freeze_delta_zeta){
      f_theta <- balldropg(t, c(g, h0))
      idx1 <- which(zeta == 1)
      idx2 <- which(zeta == 2)
      n1   <- length(idx1)
      n2   <- length(idx2)
      residual1 <- y[idx1] - f_theta[idx1]
      rss1   <- sum(residual1^2)
      shape_err <- (n / 2) + 1
      rate_err  <- 0.5 * rss1
    } else {
      f_theta <- balldropg(t, c(g, h0))
      idx1 <- which(zeta == 1)
      idx2 <- which(zeta == 2)
      residual1 <- y[idx1] - f_theta[idx1]
      rss1 <- sum(residual1^2)
      residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
      rss2 <- sum(residual2^2)
      d <- 2
      # When we consider the Jeffreys prior, I use the following code:
      # rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
      # shape_err <- n + 1#/2
      # When the prior is the sigma parameter of the "inverse gamma distribution", I use the following code:
      rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
      shape_err <- n + (d/2)
    }
    # rate_err <- 0.05 + (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    # shape_err <- 3 + n
    # sigma_sq_err <- rinvgamma(1, shape = shape_err, scale = rate_err)
    sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    theta[3] <- sigma_sq_err
    
    if (mcmc_parameters[2] == FALSE) {
      sigma_sq_err <- init[3]
      theta[3] <- sigma_sq_err
    }
    
    #-------------------------------- Gibbs step for psi_delta --------------------------------#   
    
    # # # MH for psi_delta
    # psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    # 
    # K_star_prop <- GP_correlation(t, psi_prop)
    # 
    # # if (is.null(K_star_prop)) {
    # #   log_acc <- -Inf
    # # } else {
    # Sigma_delta_cur <- (sigma_sq_err / k) * K_star_psi
    # Sigma_delta_prop <- (sigma_sq_err / k) * K_star_prop
    # Sigma_delta_prop <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta = psi_prop)
    # Sigma_delta_cur <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta)
    # log_prop_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
    # log_prop_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
    # log_prior_current <- dbeta(psi_delta, shape1 = a_psi, shape2 = b_psi, log=TRUE)
    # log_prior_prop <- dbeta(psi_prop, shape1 = a_psi, shape2 = b_psi, log=TRUE)
    # log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_cur, log = TRUE), error = function(e) -Inf)
    # log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    # log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + log_prop_current - log_prop_prop
    # # }
    # if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
    #   psi_delta <- psi_prop
    #   theta[5] <- psi_delta
    #   accept_psi <- accept_psi + 1
    # }
    
    # When the prior is the psi_delta of the "uniform(0.1, 1)", I use the following code:
    psi_prop <- rtruncnorm(1, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    
    K_star_prop <- GP_correlation(t, psi_prop)
    Sigma_delta_cur  <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta)
    Sigma_delta_prop <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_prop)
    log_prop_current <- log(dtruncnorm(psi_delta, a = 0.1, b = 1, mean = psi_prop,  sd = sigma_proposals[5]))
    log_prop_prop    <- log(dtruncnorm(psi_prop,  a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
    log_prior_current <- dunif(psi_delta, min = 0.1, max = 1, log = TRUE)
    log_prior_prop    <- dunif(psi_prop,  min = 0.1, max = 1, log = TRUE)
    log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_cur,  log = TRUE), error = function(e) -Inf)
    log_like_prop    <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    
    log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) +
      (log_prop_current - log_prop_prop)
    
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      psi_delta <- psi_prop
      theta[5]  <- psi_delta
      accept_psi <- accept_psi + 1
    }
    
    if(mcmc_parameters[3] == FALSE){
      psi_delta <- init[5]
      theta[5] <- psi_delta
    }
    
    #-------------------------------- Gibbs step for k --------------------------------#   
    
    K_star_psi <- GP_correlation(t, psi_delta)
    inv_Kstar <- tryCatch(chol2inv(chol(K_star_psi)), 
                          error = function(e) NULL)
    beta_k <- 1 + (1 / (2 * sigma_sq_err)) * quad_form_delta
    alpha_k <- (n / 2) + 0.1
    for (try_k in 1:100) {
      k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
      if (k_prop > 0 && k_prop < 1) { 
        k <- k_prop
        theta[6] <- k
        break 
      }
    }
    if(mcmc_parameters[4] == FALSE){
      k <- init[6]
      theta[6] <- k
    }
    
    #-------------------------------- Gibbs step for alpha --------------------------------#   
    
    alpha_param <- rbeta(1, sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
    theta[4] <- alpha_param
    
    if(mcmc_parameters[5] == FALSE){
      alpha_param <- init[4]
    }
    
    theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
    chain_theta[iter, ] <- theta
    chain_delta[iter, ] <- delta
    
  }
  ##################################################################
  # Remove burn-in 
  if (n_burnin > 0) {
    
    keep <- (n_burnin + 1):total_iter
    chain_theta <- chain_theta[keep, , drop = FALSE]
    chain_delta <- chain_delta[keep, , drop = FALSE]
    chain_zeta <- chain_zeta[keep, , drop = FALSE]
    loglik_chain <- loglik_chain[keep]
    
    cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
  } else {
    cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
  }
  return(list(theta = chain_theta, delta = chain_delta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_psi / n_iter))
}



set.seed(12345)
Sigma_theta <- matrix(c(0.5,0,0,0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.2, 0.2)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
n_samples       <- 50
burn_in         <- 2000
n_iter          <- 20000
# FALSE= fixed parameter
# mcmc parameter (g,h), sig2err, psidelta, k, alpha
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
  
  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000)
  g[,v] = results$theta[,1]
  h0[,v]=results$theta[,2]
  sigma_sq_err[,v]=results$theta[,3]
  alpha[,v] <- results$theta[,4]
  psi_delta[,v] <- results$theta[,5]
  k[,v] <- results$theta[,6]
}

result_m0_sh2_ortho_classic <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs)

save(result_m0_sh2_ortho_classic,file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m0_sh2_ortho_classic.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m0_sh2_ortho_classic.RData")
#load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m0_sh2_OGP_jef.RData")

g_sh2 <- result_m0_sh2_ortho_classic[[1]]
h0_sh2 <- result_m0_sh2_ortho_classic[[2]]
sigma_sq_err_sh2 <- result_m0_sh2_ortho_classic[[3]]
alpha_sh2 <- result_m0_sh2_ortho_classic[[4]]
psi_delta_sh2 <- result_m0_sh2_ortho_classic[[5]]
k_sh2 <- result_m0_sh2_ortho_classic[[6]]

par(mfrow = c(1, 3),
    mar   = c(3, 4, 1, 1) 
)
boxplot(
  colMeans(g_sh2),
  ylab = "g",
  #xlab = "Bowling Ball",
  col  = "lightseagreen",
  main = NULL
)
abline(h = 9.8, lty = 2)

boxplot(
  colMeans(h0_sh2),
  ylab = "h0",
  #xlab = "",
  col  = "lightseagreen",
  main = NULL
)
abline(h = 46.45045, lty = 2)

boxplot(
  colMeans(sigma_sq_err_sh2),
  ylab = expression(lambda^2),
  #ylim = c(0.01, 0.0452),
  #xlab = "",
  col  = "lightseagreen",
  main = NULL
)
abline(h = 0.01, lty = 2)

g_mat <- result_m0_sh2_ortho_classic[[1]]   # n_iter × n_sims
h0_mat<- result_m0_sh2_ortho_classic[[2]]

g_last <- g_mat[, 50]
h0_last <- h0_mat[, 50]
y_obs <- result_m0_sh2_ortho_classic$y_obs
y_obs30 <- y_obs[,50]
y_true <- balldropg(t, c(9.8, 46.45))

y_pred <- matrix(NA, length(t), length(g_last), byrow = FALSE)

for (i in 1:length(g_last)) {
  theta <- c(g_last[i], h0_last[i])
  y_pred[,i] <- balldropg(t, theta)
}

y_pred1  = t(y_pred)
par(mfrow = c(1,1))
boxplot(y_pred1, col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 1)
legend(x=30, y=45, legend=c("Simulated data",
                            "True code", "Predictions"), lwd=rep(2,2), col=c("gold","blue", "orange2"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))

boxplot(y_pred1[, 1:5], col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 3)
legend(x=4.2, y=46.55, legend=c("Simulated data", "Predictions",
                                "True code"), lwd=rep(2,2), col=c("gold", "orange2","blue"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))


alpha <- result_m0_sh2_ortho_classic[[4]]
hist(alpha)
plot(density(alpha))

library(ggplot2)

# 1) flatten the matrix to a single vector
alpha_vec <- as.vector(alpha)

# 2) put into a data.frame
df_alpha <- data.frame(alpha = alpha_vec)

# 3) ggplot density
ggplot(df_alpha, aes(x = alpha)) +
  geom_density(
    fill  = "#4CCDC9",   # burgundy fill
    color = "lightseagreen",   # darker border
    alpha = 0.6,         # semi‐transparent
    size  = 1
  ) +
  labs(
    x     = expression(alpha),
    y     = "Density",
    title = "Posterior Density of " ~ alpha
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    plot.margin  = unit(c(0.2,0.2,0.2,0.2), "cm")
  )
##############################################################
##-------------------------------------------------
## Empirical coverage of 50% and 95% credible CIs
##-------------------------------------------------

g_true       <- 9.8
h0_true      <- 46.45      
lambda2_true <- 0.01

# Number of simulated datasets
n_sims <- ncol(g_sh2)

# Coverage indicators (0/1) for each dataset and each CI level
cover_g_50       <- logical(n_sims)
cover_h0_50      <- logical(n_sims)
cover_lambda2_50 <- logical(n_sims)

cover_g_95       <- logical(n_sims)
cover_h0_95      <- logical(n_sims)
cover_lambda2_95 <- logical(n_sims)

ci_g_50_mat       <- matrix(NA, nrow = n_sims, ncol = 2)
ci_h0_50_mat      <- matrix(NA, nrow = n_sims, ncol = 2)
ci_lambda2_50_mat <- matrix(NA, nrow = n_sims, ncol = 2)

ci_g_95_mat       <- matrix(NA, nrow = n_sims, ncol = 2)
ci_h0_95_mat      <- matrix(NA, nrow = n_sims, ncol = 2)
ci_lambda2_95_mat <- matrix(NA, nrow = n_sims, ncol = 2)

for (v in 1:n_sims) {
  # --- 50% central credible interval: 25% and 75% quantiles ---
  ci_g_50       <- quantile(g_sh2[, v],            probs = c(0.25, 0.75))
  ci_h0_50      <- quantile(h0_sh2[, v],           probs = c(0.25, 0.75))
  ci_lambda2_50 <- quantile(sigma_sq_err_sh2[, v], probs = c(0.25, 0.75))
  
  ci_g_50_mat[v, ]       <- ci_g_50
  ci_h0_50_mat[v, ]      <- ci_h0_50
  ci_lambda2_50_mat[v, ] <- ci_lambda2_50
  
  cover_g_50[v]       <- (g_true       >= ci_g_50[1]       && g_true       <= ci_g_50[2])
  cover_h0_50[v]      <- (h0_true      >= ci_h0_50[1]      && h0_true      <= ci_h0_50[2])
  cover_lambda2_50[v] <- (lambda2_true >= ci_lambda2_50[1] && lambda2_true <= ci_lambda2_50[2])
  
  # --- 95% central credible interval: 2.5% and 97.5% quantiles ---
  ci_g_95       <- quantile(g_sh2[, v],            probs = c(0.025, 0.975))
  ci_h0_95      <- quantile(h0_sh2[, v],           probs = c(0.025, 0.975))
  ci_lambda2_95 <- quantile(sigma_sq_err_sh2[, v], probs = c(0.025, 0.975))
  
  ci_g_95_mat[v, ]       <- ci_g_95
  ci_h0_95_mat[v, ]      <- ci_h0_95
  ci_lambda2_95_mat[v, ] <- ci_lambda2_95
  
  cover_g_95[v]       <- (g_true       >= ci_g_95[1]       && g_true       <= ci_g_95[2])
  cover_h0_95[v]      <- (h0_true      >= ci_h0_95[1]      && h0_true      <= ci_h0_95[2])
  cover_lambda2_95[v] <- (lambda2_true >= ci_lambda2_95[1] && lambda2_true <= ci_lambda2_95[2])
}

# Empirical coverage across the 50 datasets
cat("50% CI coverage:\n")
cat("  g       :", mean(cover_g_50),       "\n")
cat("  h0      :", mean(cover_h0_50),      "\n")
cat("  lambda^2:", mean(cover_lambda2_50), "\n\n")

cat("95% CI coverage:\n")
cat("  g       :", mean(cover_g_95),       "\n")
cat("  h0      :", mean(cover_h0_95),      "\n")
cat("  lambda^2:", mean(cover_lambda2_95), "\n")

coverage_table <- data.frame(
  dataset        = 1:n_sims,
  cover_g_50     = cover_g_50,
  cover_h0_50    = cover_h0_50,
  cover_lambda2_50 = cover_lambda2_50,
  cover_g_95     = cover_g_95,
  cover_h0_95    = cover_h0_95,
  cover_lambda2_95 = cover_lambda2_95
)

print(coverage_table)

######################################################

get_ci <- function(samples, level = 0.95) {
  alpha <- 1 - level
  q <- quantile(samples, probs = c(alpha/2, 1 - alpha/2))
  c(lower = q[1], upper = q[2])
}

get_central_ci <- function(samples, probs = c(0.25, 0.75)) {
  q <- quantile(samples, probs = probs)
  c(lower = q[1], upper = q[2])
}

# آیا مقدار واقعی داخل بازه هست یا نه؟
is_covered <- function(true_value, ci_vec) {
  (true_value >= ci_vec[1]) && (true_value <= ci_vec[2])
}

# lower / upper
ci_g_50       <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_g_95       <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_h0_50      <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_h0_95      <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_lambda2_50 <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_lambda2_95 <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))

cover_g_50       <- logical(n_sims)
cover_g_95       <- logical(n_sims)
cover_h0_50      <- logical(n_sims)
cover_h0_95      <- logical(n_sims)
cover_lambda2_50 <- logical(n_sims)
cover_lambda2_95 <- logical(n_sims)

for (s in 1:n_sims) {
  g_s       <- g_sh2[, s]
  h0_s      <- h0_sh2[, s]
  lambda2_s <- sigma_sq_err_sh2[, s]
  
  ## ---- CI 50
  ci_g_50[s, ]       <- get_central_ci(g_s,       probs = c(0.25, 0.75))
  ci_h0_50[s, ]      <- get_central_ci(h0_s,      probs = c(0.25, 0.75))
  ci_lambda2_50[s, ] <- get_central_ci(lambda2_s, probs = c(0.25, 0.75))
  
  ## ---- CI 95٪ 
  ci_g_95[s, ]       <- get_ci(g_s,       level = 0.95)
  ci_h0_95[s, ]      <- get_ci(h0_s,      level = 0.95)
  ci_lambda2_95[s, ] <- get_ci(lambda2_s, level = 0.95)
  
  ## ---- coverage (hit / miss) ----
  cover_g_50[s]       <- is_covered(g_true,       ci_g_50[s, ])
  cover_g_95[s]       <- is_covered(g_true,       ci_g_95[s, ])
  cover_h0_50[s]      <- is_covered(h0_true,      ci_h0_50[s, ])
  cover_h0_95[s]      <- is_covered(h0_true,      ci_h0_95[s, ])
  cover_lambda2_50[s] <- is_covered(lambda2_true, ci_lambda2_50[s, ])
  cover_lambda2_95[s] <- is_covered(lambda2_true, ci_lambda2_95[s, ])
}

len_g_50       <- ci_g_50[, "upper"]       - ci_g_50[, "lower"]
len_g_95       <- ci_g_95[, "upper"]       - ci_g_95[, "lower"]
len_h0_50      <- ci_h0_50[, "upper"]      - ci_h0_50[, "lower"]
len_h0_95      <- ci_h0_95[, "upper"]      - ci_h0_95[, "lower"]
len_lambda2_50 <- ci_lambda2_50[, "upper"] - ci_lambda2_50[, "lower"]
len_lambda2_95 <- ci_lambda2_95[, "upper"] - ci_lambda2_95[, "lower"]

g_mean       <- colMeans(g_sh2)
h0_mean      <- colMeans(h0_sh2)
lambda2_mean <- colMeans(sigma_sq_err_sh2)

bias_g       <- g_mean       - g_true
bias_h0      <- h0_mean      - h0_true
bias_lambda2 <- lambda2_mean - lambda2_true

summary_table <- data.frame(
  parameter = rep(c("g", "h0", "lambda^2"), each = 2),
  CI_level  = rep(c("50%", "95%"), times = 3),
  coverage  = c(
    mean(cover_g_50),
    mean(cover_g_95),
    mean(cover_h0_50),
    mean(cover_h0_95),
    mean(cover_lambda2_50),
    mean(cover_lambda2_95)
  ) * 100,  
  mean_CI_length = c(
    mean(len_g_50),
    mean(len_g_95),
    mean(len_h0_50),
    mean(len_h0_95),
    mean(len_lambda2_50),
    mean(len_lambda2_95)
  ),
  mean_abs_bias = c(
    mean(abs(bias_g)),
    mean(abs(bias_g)),
    mean(abs(bias_h0)),
    mean(abs(bias_h0)),
    mean(abs(bias_lambda2)),
    mean(abs(bias_lambda2))
  )
)

summary_table

library(ggplot2)

## -------- g: 50% credible intervals --------
df_ci_g <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_g_95[, "lower"],
  upper   = ci_g_95[, "upper"],
  covered = factor(cover_g_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)


# مرکز هر بازه (فقط برای رسم نقطه)
df_ci_g$center <- (df_ci_g$lower + df_ci_g$upper) / 2

par(mfrow = c(1,3))
p_g_95 <- ggplot(df_ci_g,
                 aes(y = dataset,
                     x = center,
                     xmin = lower,
                     xmax = upper,
                     color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = g_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = "g",
    y = "dataset index",
    color = "coverage",
    title = "95% credible intervals for g across 50 datasets"
  ) +
  theme_minimal()

p_g_95

## -------- h0: 50% credible intervals --------
df_ci_h0 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_h0_50[, "lower"],
  upper   = ci_h0_50[, "upper"],
  covered = factor(cover_h0_50,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_h0$center <- (df_ci_h0$lower + df_ci_h0$upper) / 2

p_h0_50 <- ggplot(df_ci_h0,
                  aes(y = dataset,
                      x = center,
                      xmin = lower,
                      xmax = upper,
                      color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = h0_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(h[0]),
    y = "dataset index",
    color = "coverage",
    title = "50% credible intervals for h0 across 50 datasets"
  ) +
  theme_minimal()

p_h0_50

## -------- h0: 95% credible intervals --------
df_ci_h0 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_h0_95[, "lower"],
  upper   = ci_h0_95[, "upper"],
  covered = factor(cover_h0_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_h0$center <- (df_ci_h0$lower + df_ci_h0$upper) / 2

p_h0_95 <- ggplot(df_ci_h0,
                  aes(y = dataset,
                      x = center,
                      xmin = lower,
                      xmax = upper,
                      color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = h0_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(h[0]),
    y = "dataset index",
    color = "coverage",
    title = "95% credible intervals for h0 across 50 datasets"
  ) +
  theme_minimal()

p_h0_95

## -------- lambda^2: 50% credible intervals --------
df_ci_lambda2 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_lambda2_50[, "lower"],
  upper   = ci_lambda2_50[, "upper"],
  covered = factor(cover_lambda2_50,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_lambda2$center <- (df_ci_lambda2$lower + df_ci_lambda2$upper) / 2

p_lambda2_50 <- ggplot(df_ci_lambda2,
                       aes(y = dataset,
                           x = center,
                           xmin = lower,
                           xmax = upper,
                           color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = lambda2_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(lambda^2),
    y = "dataset index",
    color = "coverage",
    title = expression("50% credible intervals for " * lambda^2 * " across 50 datasets")
  ) +
  theme_minimal()

p_lambda2_50

## -------- lambda^2: 95% credible intervals --------
df_ci_lambda2 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_lambda2_95[, "lower"],
  upper   = ci_lambda2_95[, "upper"],
  covered = factor(cover_lambda2_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_lambda2$center <- (df_ci_lambda2$lower + df_ci_lambda2$upper) / 2

p_lambda2_95 <- ggplot(df_ci_lambda2,
                       aes(y = dataset,
                           x = center,
                           xmin = lower,
                           xmax = upper,
                           color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = lambda2_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(lambda^2),
    y = "dataset index",
    color = "coverage",
    title = expression("95% credible intervals for " * lambda^2 * " across 50 datasets")
  ) +
  theme_minimal()

p_lambda2_95
#########################################################



