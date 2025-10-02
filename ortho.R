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
  x_vec <- cbind(1, -0.5 * (t * t_range + t_min)^2)
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

H_matrix <- function(psi_delta) {
  A <- A_scalar(psi_delta)
  B <- B_scalar(psi_delta)
  D <- D_scalar(psi_delta)
  matrix(c(A, -((a^2)*B)/2, -((a^2)*B)/2, ((a^4)*D)/4), nrow = 2, byrow = TRUE)
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
    
    # Gibbs for theta
    zeta_1_indices <- which(zeta == 1)
    zeta_2_indices <- which(zeta == 2)
    X <- cbind(1, -0.5 * (t * t_range + t_min)^2)
    x1 <- X[zeta_1_indices, , drop = FALSE]  
    x2 <- X[zeta_2_indices, , drop = FALSE]
    theta_hat     <- matrix(c(46.45, 9.8), ncol = 1)
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
    
    #rate_err <- 0.5 * (rss1 + rss2 + (k * quad_form_delta))
    rate_err <- 0.5 + (0.5 * ( rss1 + rss2 + (theta[6] * quad_form_delta)))
    shape_err <- 4 + n
    sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    theta[3] <- sigma_sq_err
    if (mcmc_parameters[2] == FALSE) {
      sigma_sq_err <- init[3]
    }
    
    # MH for psi_delta
    psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
    
    K_star_prop <- GP_correlation(t, psi_prop)
    
    # if (is.null(K_star_prop)) {
    #   log_acc <- -Inf
    # } else {
      #Sigma_delta_cur <- (sigma_sq_err / k) * K_star_psi
      #Sigma_delta_prop <- (sigma_sq_err / k) * K_star_prop
      Sigma_delta_prop <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta = psi_prop)
      Sigma_delta_cur <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta)
      log_prop_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
      log_prop_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
      log_prior_current <- dbeta(psi_delta, shape1 = 1, shape2 = 9, log=TRUE)
      log_prior_prop <- dbeta(psi_prop, shape1 = 1, shape2 = 9, log=TRUE)
      log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_cur, log = TRUE), error = function(e) -Inf)
      log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
      log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + log_prop_current - log_prop_prop
   # }
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      psi_delta <- psi_prop
      theta[5] <- psi_delta
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
    
    # Gibbs for k  (use K_* not R)
    K_star_psi <- GP_correlation(t, psi_delta)
    #K_star_psi <- 0.5 * (K_star_psi + t(K_star_psi))       # numerical hygiene
    
    inv_Kstar <- tryCatch(chol2inv(chol(K_star_psi)),       # more stable than solve()
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
  
  # Remove burn-in samples if n_burnin > 0
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
k = 0.2
#sigma_sq_delta <- 0.1 
sim_psi_delta <- 0.1 
sigma_sq_err <- 0.1
sigma_sq_delta <- sigma_sq_err / k
n_samples <- 1
n_iter <- 10000
burn_in <- 1000

sigma_props <- c(NA, NA, NA, NA, 0.1, NA) 
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
  res <- mcmc_step6(y_1, t, n_iter, init, sigma_props, mcmc_parameters,
                    Sigma_theta, n_burnin = burn_in, a_psi = 1, b_psi = 9)
  
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

par(mfrow=c(2,3))
boxplot(g_chain, main="g")
abline(h=9.8)
boxplot(h0_chain, main="h0")
abline(h=46.45)
boxplot(sigma_chain, main=expression(sigma[err]^2))
abline(h=0.1)
boxplot(alpha_chain, main=expression(alpha))
boxplot(psi_chain,   main=expression(psi[delta]))
abline(h=0.1)
boxplot(k_chain,main="k")
abline(h=0.2)

par(mfrow=c(2,3))
plot(g_chain, type = "l", main="g")
abline(h=9.8)
plot(h0_chain, type = "l", main="h0")
abline(h=46.45)
plot(sigma_chain, type = "l",main=expression(sigma[err]^2))
abline(h=0.1)
plot(alpha_chain, type = "l",main=expression(alpha))
plot(psi_chain, type = "l",main=expression(psi[delta]))
abline(h=0.1)
plot(k_chain,type = "l",main="k")
abline(h=0.2)
