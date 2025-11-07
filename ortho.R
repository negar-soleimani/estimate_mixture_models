
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

don <- read_xlsx("/Users/negarsoleimani/Documents/phd/paper1/Ball_drops_data.xlsx", sheet = 5)
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop); don <- don[don$drop == 1, ]

GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
}

t <- don$time
y <- don$Height
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
a <- t_range
n <- length(y)


balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  #h[h < 0] <- 0
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

# I_x <- function(x, psi_delta){
#   psi <- psi_delta
#   p1 <- (psi * (1 - exp(-x/psi)) + psi * (1 - exp(-(1-x)/psi))) 
#   p2 <- (((psi^2) * (1 - exp(-x/psi))) + (psi * x * exp(-x/psi)) + ((psi^2) * (1 - exp((1-x)/psi))) - (psi * (1-x) * exp(-(1-x)/psi)))
#   p3 <- ((2 * psi^3) - exp(-x/psi) * (((x^2) * psi) + (2 * x * psi^2)) + (2 * psi^3)) + (2 * psi^3) - (exp((1-x)/psi) * ((((1-x)^2) * psi) + (2 * (1-x) * psi^2) + (2 * psi^3)))
#   term1 <- ((x^2) * p1)
#   term2 <- ((2*x) * p2)
#   term1 + term2 + p3
# }

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

# H_matrix <- function(psi_delta) {
#   A <- A_scalar(psi_delta)
#   B <- B_scalar(psi_delta)
#   D <- D_scalar(psi_delta)
#   matrix(c(A, -((a^2)*B)/2, -((a^2)*B)/2, ((a^4)*D)/4), nrow = 2, byrow = TRUE)
# }

H_matrix <- function(psi_delta) {
  A <- A_scalar(psi_delta); B <- B_scalar(psi_delta); D <- D_scalar(psi_delta)
  H <- matrix(c(A, -((a^2)*B)/2, -((a^2)*B)/2, ((a^4)*D)/4), nrow = 2, byrow = TRUE)
  H <- 0.5 * (H + t(H))
  H + diag(1e-10, 2)
}


# GP_correlation <- function(t, psi_delta) {
#   R_se <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
#   
#   h_mat <- t(sapply(t, function(x) h_vec_x(x, psi_delta)))
#   H <- H_matrix(psi_delta)
#   H_inv <- tryCatch(solve(H), error = function(e) NULL)
#   if (is.null(H_inv)) return(NULL)
#   
#   K_star <- R_se - h_mat %*% H_inv %*% t(h_mat)
#   #K_star <- 0.5 * (K_star + t(K_star))  # symmetrize
#   K_star
# }

GP_correlation <- function(t, psi_delta) {
  R_se  <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
  h_mat <- t(sapply(t, function(x) h_vec_x(x, psi_delta)))
  H     <- H_matrix(psi_delta)
  H_inv <- tryCatch(chol2inv(chol(H)), error = function(e) NULL)
  if (is.null(H_inv)) return(NULL)
  K_star <- R_se - h_mat %*% H_inv %*% t(h_mat)
  K_star <- 0.5 * (K_star + t(K_star))   
  eig <- eigen(K_star, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < 1e-10) {                 
    K_star <- K_star + diag(1e-8 - min(0, min(eig)) + 1e-12, nrow(K_star))
  }
  K_star
}


GP_covariance_star_complete <- function(t, sigma_sq_err, k, psi_delta) {
  K_star <- GP_correlation(t, psi_delta)
  if (is.null(K_star)) {
    return(NULL)
  }
  (sigma_sq_err / k) * K_star
}

# psi_delta <- 0.1
# sigma_sq_err <- 0.1
# k <- 0.2
# sigma_sq_delta <- sigma_sq_err / k
# t <- seq(0,1,0.25)
# t_range <- max(t) - min(t)
# a <- t_range
# 
# J_x(t, psi_delta)
# I_x(t, psi_delta)
# #
# H_matrix(psi_delta)
# GP_covariance(t, sigma_sq_delta, psi_delta)
# GP_correlation(t, psi_delta)
# GP_covariance(t, 1, psi_delta)-t(h_vec_x(t, psi_delta))%*%solve(H_matrix(psi_delta))%*%h_vec_x(t, psi_delta)
# GP_covariance_star_complete(t, 1, 1, psi_delta)
# GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta)


mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000, a_psi, b_psi) {
  # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T")
  
  # Total iterations = burn-in + desired samples
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
    
    Sigma_delta <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta) #c*(x,x')
    Sigma_delta <- 0.5 * (Sigma_delta + t(Sigma_delta)) + diag(1e-10, n)
    
    
    # if (is.null(Sigma_delta)) {
    #   log_likelihood <- -Inf
    # } else {
    f_theta <- balldropg(t, c(g, h0))
    mean1 <- f_theta
    mean2 <- f_theta + delta
    
    ###s <- 0.3
    prob_zeta <- 1 / (1 + (alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err))) /
                        ((1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err))))
    ##*(abs(delta) > s)
    zeta <- 1 + (runif(length(y)) < prob_zeta)
    
    chain_zeta[iter, ] = zeta
    
    
    log_likelihood <- sum(log(ifelse(zeta == 1,
                                     alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
                                     (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    loglik_chain[iter] <- log_likelihood
    
    zeta_2_indices <- which(zeta == 2)
    if (length(zeta_2_indices) > 0) {
      y_m <- y[zeta_2_indices]
      Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
        Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
      Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
      #Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
      Sigma_inv <- tryCatch(chol2inv(chol(Sigma_delta_ymym)),
                            error = function(e) diag(1, nrow(Sigma_delta_ymym)))
      #Sigma_inv <- 0.5 * (Sigma_inv + t(Sigma_inv))
      mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
      Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
      # jitter <- 1e-8
      # eig <- eigen(Sigma_delta_hat, symmetric = TRUE)
      # min_eig <- min(eig$values)
      # if (min_eig < 0) {
      #   Sigma_delta_hat <- Sigma_delta_hat + diag(abs(min_eig) + jitter, nrow(Sigma_delta_hat))
      # }
      # Sigma_delta_hat <- Sigma_delta_hat + diag(jitter, nrow(Sigma_delta_hat))
      # Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
      Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
      delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    } else {
      delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
    }
    #}
    
    K_star_psi <- GP_correlation(t, psi_delta)
    #inv_Kstar_psi <- tryCatch(solve(K_star_psi), error = function(e) NULL)
    inv_Kstar_psi <- tryCatch(chol2inv(chol(K_star_psi)), error = function(e) NULL)
    
    if (is.null(inv_Kstar_psi)) {
      quad_form_delta <- -Inf
    } else {
      #quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
      quad_form_delta <- sum(delta * (inv_Kstar_psi %*% delta))
    }
    
    #if (is.null(K_star_psi)) {
    #  quad_form_delta <- -Inf
    #} else {
    #  quad_form_delta <- as.numeric(t(delta) %*% solve(K_star_psi) %*% delta)
    #}
    if (quad_form_delta == -Inf) next
    
    # # Gibbs for theta
    zeta_1_indices <- which(zeta == 1)
    zeta_2_indices <- which(zeta == 2)
    X <- cbind(1, -0.5 * (t * t_range)^2) #cbind(1, -0.5 * (t * t_range + t_min)^2)
    x1 <- X[zeta_1_indices, , drop = FALSE]
    x2 <- X[zeta_2_indices, , drop = FALSE]
    theta_hat     <- matrix(c(46.14, 9.8), ncol = 1)
    inv_sigma_theta <- solve(Sigma_theta)
    A <- ((t(x1) %*% x1) / theta[3]) + ((t(x2) %*% x2) / theta[3]) + inv_sigma_theta
    Sigmapost_theta <- solve(A)
    y1 <- matrix(y[zeta_1_indices], ncol = 1)
    y2 <- matrix(y[zeta_2_indices], ncol = 1)
    d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    B <- (t(x1) %*% y1) / theta[3] +
         (t(x2) %*% y2) / theta[3] -
         (t(x2) %*% d2) / theta[3] +
      inv_sigma_theta %*% theta_hat
    Mupost_theta <- Sigmapost_theta %*% B
    theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
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
    
    # Gibbs for sigma_sq_err
    f_theta <- balldropg(t, c(g, h0))
    idx1 <- which(zeta == 1)
    idx2 <- which(zeta == 2)
    residual1 <- y[idx1] - f_theta[idx1]
    rss1 <- sum(residual1^2)
    residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
    rss2 <- sum(residual2^2)
    
    # rate_err <- 2 + (0.5 * ( rss1 + rss2 + (theta[6] * quad_form_delta)))
    # shape_err <- 1 + n
    rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    shape_err <- n + 1/2
    ##sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    #sigma_sq_err <- 1/rgamma(1, shape = shape_err, rate = rate_err)
    sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)

    theta[3] <- sigma_sq_err
    if (mcmc_parameters[2] == FALSE) {
      sigma_sq_err <- init[3]
    }
    
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
    
    # ## ---- MH for psi_delta (Uniform(0.1, 1) prior) ----
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
    }
    
    # # Gibbs for k
    # R_se <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
    # alpha_k <- (n / 2) + 1
    # #beta_k <- (1 / (2 * sigma_sq_err)) * as.numeric(t(delta) %*% solve(R_se) %*% delta)
    # beta_k <- (1 / (2 * sigma_sq_err)) * sum(delta * (solve(R_se, delta)))
    # for (try_k in 1:100) {
    #   k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
    #   if (k_prop > 0 && k_prop < 1) {
    #     k <- k_prop
    #     theta[6] <- k
    #     break
    #   }
    # }
    # if(mcmc_parameters[4] == FALSE){
    #   k <- init[6]
    # }
    
    K_star_psi <- GP_correlation(t, psi_delta)
    #K_star_psi <- 0.5 * (K_star_psi + t(K_star_psi))
    
    inv_Kstar <- tryCatch(chol2inv(chol(K_star_psi)), 
                          error = function(e) NULL)
    #if (is.null(inv_Kstar)) {
    #  beta_k <- Inf
    #} else {
    beta_k <- (1 / (2 * sigma_sq_err)) * sum(t(delta) %*% inv_Kstar %*% delta)
    #}
    
    alpha_k <- (n / 2) + 1
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
    }
    
    # Gibbs step for alpha
    alpha_param <- rbeta(1, sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
    theta[4] <- alpha_param
    
    if(mcmc_parameters[5] == FALSE){
      alpha_param <- init[4]
    }
    
    theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
    chain_theta[iter, ] <- theta
    chain_delta[iter, ] <- delta
    
  }
  
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

## Psi1 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.01 
sigma_sq_err <- (0.1)^2
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 1
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, FALSE, FALSE, FALSE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.14, 0.01, 0.5, sim_psi_delta, 0.2)

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
  
  #y_1 <- balldropg(t, c(9.8, 46.14)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_1 <- balldropg(t, c(9.8, 46.14)) + rnorm(n, 0, sqrt(sigma_sq_err)*(1+1/k))
  
  y_obs[, v] <- y_1
  a_psi = 1
  b_psi = 1#99
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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

res <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

g <- res[[1]]
h <- res[[2]]
sig <- res[[3]]
alpha <- res[[4]]
psi <- res[[5]]
k <- res[[6]]

par(mfrow = c(2,3))
plot(g, type = "l")
abline(h = 9.8, col = "red")
plot(h, type = "l")
abline(h = 46.14, col = "red")
plot(sig, type = "l")
abline(h = 0.01, col = "red")
plot(alpha, type = "l")
plot(psi, type = "l")
plot(k, type = "l")
abline(h = 0.2, col = "red")
hist(k)
hist(sig)
hist(psi)
hist(alpha)


GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
}
GP_covariance(t, 0.05, 0.3)

boxplot(colMeans(g))
abline(h = 9.8)
boxplot(colMeans(h))
abline(h = 46.14)
boxplot(colMeans(sig))
abline(h = 0.1)
boxplot(colMeans(alpha))
boxplot(colMeans(psi))
abline(h = 0.7)
boxplot(colMeans(k))
abline(h = 0.2)

plot(y_1)
zeta <- res[[8]]
prob_zeta_m2 <- colMeans(res[[8]][[1]] == 2)
delta <- res[[7]][[1]]
boxplot(delta)
abline(h = 0, col = "red")


delta_mean <- apply(do.call(rbind, delta_list), 2, mean)
plot(t, delta_mean, type="l", main="Posterior mean of delta")


#result_m2_sh2_psi1_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi1_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi1_ortho.RData")

########################
set.seed(12345)

init_base       <- c(9.8, 46.14, 0.01, 0.7, 0.5, 0.2)
sigma_proposals <- c(NA,NA,NA,NA,0.5,NA)
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
#Sigma_theta     <- matrix(c(0.1,0,0,0.1),2)
Sigma_theta     <- matrix(c(0.5,0,0,0.5),2)
n_iter          <- 20000
burn_in         <- 1000

results_real_sh2 <- mcmc_step6(
  y = y, t = t,
  n_iter = n_iter,
  init   = init_base,
  sigma_proposals = sigma_proposals,
  mcmc_parameters = mcmc_parameters,
  Sigma_theta     = Sigma_theta,
  n_burnin        = burn_in
)
g_chain <- results_real_sh2$theta[,1]
h0_chain <- results_real_sh2$theta[,2]
sigma_sq_err_chain <- results_real_sh2$theta[,3]
alpha_chain <- results_real_sh2$theta[,4]
psi_delta_chain <- results_real_sh2$theta[,5]
k_chain <- results_real_sh2$theta[,6]
delta_chain <- results_real_sh2$delta
zeta_chain <- results_real_sh2$zeta
loglik_chain <- results_real_sh2$loglik
accept_rate_psi <- results_real_sh2$accept_rate_psi

par(mfrow = c(1, 1))
boxplot(results_real_sh2$delta, ylab = expression(delta))

par(mfrow = c(2, 3))
plot(results_real_sh2$theta[,1], type = "l",
     xlab = "Blue Basketball", ylab = "g")
abline(h = 9.8, col = "red")
plot(results_real_sh2$theta[,2], type = "l",
     xlab = "Blue Basketball", ylab = "h0")
abline(h = 46.14, col = "red")
plot(results_real_sh2$theta[,3], type = "l",
     xlab = "Blue Basketball", ylab = expression(sigma[err]^2))
abline(h = 0.01, col = "red")
plot(results_real_sh2$theta[,4], type = "l",
     xlab = "Blue Basketball", ylab = expression(alpha))
plot(results_real_sh2$theta[,5], type = "l",
     xlab = "Blue Basketball", ylab = expression(psi[delta]))
plot(results_real_sh2$theta[,6], type = "l",
     xlab = "Blue Basketball", ylab = "k")
######################


## Psi2 ############################################
set.seed(12345)
k = 1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.1 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 1
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 1)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1
  b_psi = 1#9
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res2 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi2_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi2_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi2_ortho.RData")
## Psi3 ############################################
set.seed(12345)
k = 1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.2 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 1
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 1)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#2
  b_psi = 1#8
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res3 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi3_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi3_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi3_ortho.RData")


## Psi4 ############################################
set.seed(12345)
k = 1
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.3 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 1
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 1)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#3
  b_psi = 1#7
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res4 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi4_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi4_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi4_ortho.RData")
View(res2)
g1 <- res[["theta"]][, "g"]
g2 <- res2[[1]]
g3 <- res3[[1]]
g4 <- res4[[1]]

boxplot(g1, g2, g3, g4,
        names = c("0.01","0.1","0.2","0.3"),
        xlab = "psi_delta",
        ylab = "g")
abline(h=9.8)

h1 <- res[["theta"]][, "h0"]
h2 <- res2[[2]]
h3 <- res3[[2]]
h4 <- res4[[2]]

boxplot(h1, h2, h3, h4,
        names = c("0.01","0.1","0.2","0.3"),
        xlab = "psi_delta",
        ylab = "h")
abline(h=46.45)

sig1 <- res[["theta"]][, 3]
sig2 <- res2[[3]]
sig3 <- res3[[3]]
sig4 <- res4[[3]]
boxplot(sig1, sig2, sig3, sig4,
        names = c("0.01","0.1","0.2","0.3"),
        xlab = "psi_delta",
        ylab = "sig")
abline(h=0.1)

a1 <- res[["theta"]][, 4]
a2 <- res2[[4]]
a3 <- res3[[4]]
a4 <- res4[[4]]
boxplot(a1, a2, a3, a4,
        names = c("0.01","0.1","0.2","0.3"),
        xlab = "psi_delta",
        ylab = "alpha")

p1 <- res[["theta"]][, 5]
p2 <- res2[[5]]
p3 <- res3[[5]]
p4 <- res4[[5]]

boxplot(p1, p2, p3, p4,
        names = c("0.01","0.1","0.2","0.3"),
        xlab = "psi_delta",
        ylab = "psi")

k1 <- res[["theta"]][, 6]
k2 <- res2[[6]]
k3 <- res3[[6]]
k4 <- res4[[6]]

boxplot(k1, k2, k3, k4,
        names = c("0.01","0.1","0.2","0.3"),
        xlab = "psi_delta",
        ylab = "k")
abline(h=0.2)

## Psi5 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.4 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.3, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 0.2)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#4
  b_psi = 1#6
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res5 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi5_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi5_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi5_ortho.RData")


## Psi6 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.5 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 0.2)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#5
  b_psi = 1#5
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res6 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi6_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi6_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi6_ortho.RData")


## Psi7 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.6 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 0.2)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#6
  b_psi = 1#4
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res7 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi7_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi7_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi7_ortho.RData")


## Psi8 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.7 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 0.2)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#7
  b_psi = 1#3
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res8 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi8_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi8_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi8_ortho.RData")


## Psi9 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.8 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.5, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 0.2)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#8
  b_psi = 1#2
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res9 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi9_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi9_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi9_ortho.RData")


## Psi10 ############################################
set.seed(12345)
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.2 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 3
n_iter <- 20000
burn_in <- 5000

sigma_props <- c(NA, NA, NA, NA, 0.4, NA) 
#(g,h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
init <- c(9.8, 46.45, 0.1, 0.5, sim_psi_delta, 0.2)

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
  
  y_1 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
  y_obs[, v] <- y_1
  a_psi = 1#9
  b_psi = 1
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi, b_psi)
  
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
res10 <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m2_sh2_psi10_ortho <- list(g_chain, h0_chain, sigma_chain, alpha_chain, psi_chain, k_chain, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m2_sh2_psi10_ortho,file = "/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi10_ortho.RData")
###############################################################################################

load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi1_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi2_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi3_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi4_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi5_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi6_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi7_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi8_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi9_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/Orthogonality/result_m2_sh2_psi10_ortho.RData")

g1 <- result_m2_sh2_psi1_ortho[[1]]
g2 <- result_m2_sh2_psi2_ortho[[1]]
g3 <- result_m2_sh2_psi3_ortho[[1]]
g4 <- result_m2_sh2_psi4_ortho[[1]]
g5 <- result_m2_sh2_psi5_ortho[[1]]
g6 <- result_m2_sh2_psi6_ortho[[1]]
g7 <- result_m2_sh2_psi7_ortho[[1]]
g8 <- result_m2_sh2_psi8_ortho[[1]]
g9 <- result_m2_sh2_psi9_ortho[[1]]
g10 <- result_m2_sh2_psi10_ortho[[1]]

par(mfrow = c(1,1))
boxplot(colMeans(g1), colMeans(g2), colMeans(g3), colMeans(g4), colMeans(g5), 
        colMeans(g6), colMeans(g7), colMeans(g8), colMeans(g9), colMeans(g10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "g")
abline(h=9.8)

h01 <- result_m2_sh2_psi1_ortho[[2]]
h02 <- result_m2_sh2_psi2_ortho[[2]]
h03 <- result_m2_sh2_psi3_ortho[[2]]
h04 <- result_m2_sh2_psi4_ortho[[2]]
h05 <- result_m2_sh2_psi5_ortho[[2]]
h06 <- result_m2_sh2_psi6_ortho[[2]]
h07 <- result_m2_sh2_psi7_ortho[[2]]
h08 <- result_m2_sh2_psi8_ortho[[2]]
h09 <- result_m2_sh2_psi9_ortho[[2]]
h010 <- result_m2_sh2_psi10_ortho[[2]]

par(mfrow = c(1,1))
boxplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05), 
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "h0")
abline(h = 46.45)

sigma1 <- result_m2_sh2_psi1_ortho[[3]]
sigma2 <- result_m2_sh2_psi2_ortho[[3]]
sigma3 <- result_m2_sh2_psi3_ortho[[3]]
sigma4 <- result_m2_sh2_psi4_ortho[[3]]
sigma5 <- result_m2_sh2_psi5_ortho[[3]]
sigma6 <- result_m2_sh2_psi6_ortho[[3]]
sigma7 <- result_m2_sh2_psi7_ortho[[3]]
sigma8 <- result_m2_sh2_psi8_ortho[[3]]
sigma9 <- result_m2_sh2_psi9_ortho[[3]]
sigma10 <- result_m2_sh2_psi10_ortho[[3]]

par(mfrow = c(1,1))
boxplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "sigma")
abline(h = 0.1)

alpha1 <- result_m2_sh2_psi1_ortho[[4]]
alpha2 <- result_m2_sh2_psi2_ortho[[4]]
alpha3 <- result_m2_sh2_psi3_ortho[[4]]
alpha4 <- result_m2_sh2_psi4_ortho[[4]]
alpha5 <- result_m2_sh2_psi5_ortho[[4]]
alpha6 <- result_m2_sh2_psi6_ortho[[4]]
alpha7 <- result_m2_sh2_psi7_ortho[[4]]
alpha8 <- result_m2_sh2_psi8_ortho[[4]]
alpha9 <- result_m2_sh2_psi9_ortho[[4]]
alpha10 <- result_m2_sh2_psi10_ortho[[4]]

par(mfrow = c(1,1))
boxplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "alpha")

psi1 <- result_m2_sh2_psi1_ortho[[5]]
psi2 <- result_m2_sh2_psi2_ortho[[5]]
psi3 <- result_m2_sh2_psi3_ortho[[5]]
psi4 <- result_m2_sh2_psi4_ortho[[5]]
psi5 <- result_m2_sh2_psi5_ortho[[5]]
psi6 <- result_m2_sh2_psi6_ortho[[5]]
psi7 <- result_m2_sh2_psi7_ortho[[5]]
psi8 <- result_m2_sh2_psi8_ortho[[5]]
psi9 <- result_m2_sh2_psi9_ortho[[5]]
psi10 <- result_m2_sh2_psi10_ortho[[5]]

par(mfrow = c(1,1))
boxplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "psi")

k1 <- result_m2_sh2_psi1_ortho[[6]]
k2 <- result_m2_sh2_psi2_ortho[[6]]
k3 <- result_m2_sh2_psi3_ortho[[6]]
k4 <- result_m2_sh2_psi4_ortho[[6]]
k5 <- result_m2_sh2_psi5_ortho[[6]]
k6 <- result_m2_sh2_psi6_ortho[[6]]
k7 <- result_m2_sh2_psi7_ortho[[6]]
k8 <- result_m2_sh2_psi8_ortho[[6]]
k9 <- result_m2_sh2_psi9_ortho[[6]]
k10 <- result_m2_sh2_psi10_ortho[[6]]

par(mfrow = c(1,1))
boxplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), colMeans(k10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "psi_delta",
        ylab = "k")
abline(h = 0.2)


res_psi2 <- result_m2_sh2_psi2_ortho
g_psi2 <- result_m2_sh2_psi2_ortho[[1]][,10]
h0_psi2 <- result_m2_sh2_psi2_ortho[[2]][,10]
sigma_psi2 <- result_m2_sh2_psi2_ortho[[3]][,10]
alpha_psi2 <- result_m2_sh2_psi2_ortho[[4]][,10]
psi_psi2 <- result_m2_sh2_psi2_ortho[[5]][,10]
k_psi2 <- result_m2_sh2_psi2_ortho[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi2, type = "l", ylab = expression(g ~ " with " ~ psi == 0.1), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi2, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.1), xlab = "iteration")
abline(h = 46.45, col = "red")

plot(sigma_psi2, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.1), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(alpha_psi2, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.1), xlab = "iteration")

plot(psi_psi2, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.1), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(k_psi2, type = "l", ylab = expression(k ~ " with " ~ psi == 0.1), xlab = "iteration")
abline(h = 0.2, col = "red")

res_psi8 <- result_m2_sh2_psi8_ortho
g_psi8 <- result_m2_sh2_psi8_ortho[[1]][,10]
h0_psi8 <- result_m2_sh2_psi8_ortho[[2]][,10]
sigma_psi8 <- result_m2_sh2_psi8_ortho[[3]][,10]
alpha_psi8 <- result_m2_sh2_psi8_ortho[[4]][,10]
psi_psi8 <- result_m2_sh2_psi8_ortho[[5]][,10]
k_psi8 <- result_m2_sh2_psi8_ortho[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi8, type = "l", ylab = expression(g ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi8, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 46.45, col = "red")

plot(sigma_psi8, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(alpha_psi8, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.7), xlab = "iteration")

plot(psi_psi8, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.7, col = "red")

plot(k_psi8, type = "l", ylab = expression(k ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.2, col = "red")

res_psi9 <- result_m2_sh2_psi9_ortho
g_psi9 <- result_m2_sh2_psi9_ortho[[1]][,10]
h0_psi9 <- result_m2_sh2_psi9_ortho[[2]][,10]
sigma_psi9 <- result_m2_sh2_psi9_ortho[[3]][,10]
alpha_psi9 <- result_m2_sh2_psi9_ortho[[4]][,10]
psi_psi9 <- result_m2_sh2_psi9_ortho[[5]][,10]
k_psi9 <- result_m2_sh2_psi9_ortho[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi9, type = "l", ylab = expression(g ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi9, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 46.45, col = "red")

plot(sigma_psi9, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(alpha_psi9, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.8), xlab = "iteration")

plot(psi_psi9, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.8, col = "red")

plot(k_psi9, type = "l", ylab = expression(k ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.2, col = "red")
