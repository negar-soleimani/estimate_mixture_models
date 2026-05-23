mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000) {
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
    
    #-------------------------------- Gibbs step for theta --------------------------------#
    
    # prior Normal distribution
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
    
    # K_star_psi <- GP_correlation(t, psi_delta)
    # inv_Kstar_psi <- tryCatch(chol2inv(chol(K_star_psi)), error = function(e) NULL)
    # 
    # if (is.null(inv_Kstar_psi)) {
    #   quad_form_delta <- -Inf
    # } else {
    #   quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
    # }
    # 
    # if (quad_form_delta == -Inf) next
    
    K_star_psi <- GP_correlation(t, psi_delta)
    
    Kpd <- chol_adapt(K_star_psi, jitter0 = 1e-10, jitter_max = 1e-2)
    inv_Kstar_psi <- inv_from_chol(Kpd$chol)
    
    quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
    
    quad_form_delta <- max(quad_form_delta, 0)
    
    f_theta <- balldropg(t, c(g, h0))
    idx1 <- which(zeta == 1)
    idx2 <- which(zeta == 2)
    residual1 <- y[idx1] - f_theta[idx1]
    rss1 <- sum(residual1^2)
    residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
    rss2 <- sum(residual2^2)
    
    d <- 2
    # When we consider the Jeffreys prior, I use the following code:
    rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    shape_err <- n + (d/2)
    # When the prior is the sigma parameter of the "inverse gamma distribution", I use the following code:
    # rate_err <- 1 + (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    # shape_err <- 2 + n
    # sigma_sq_err <- rinvgamma(1, shape = shape_err, scale = rate_err)
    # sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    sigma_sq_err <- 1/rgamma(1, shape = shape_err, rate = rate_err)
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
    Sigma_delta_cur <- 0.5 * ( Sigma_delta_cur + t( Sigma_delta_cur))
    Sigma_delta_cur <-  Sigma_delta_cur + diag(1e-10, nrow(Sigma_delta_cur))
    Sigma_delta_prop <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_prop)
    Sigma_delta_prop <- 0.5 * (Sigma_delta_prop + t(Sigma_delta_prop))
    Sigma_delta_prop <- Sigma_delta_prop + diag(1e-10, nrow(Sigma_delta_prop))
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
    alpha_k <- (n / 2) + 1
    beta_k <- (1 / (2 * sigma_sq_err)) * quad_form_delta
    # beta_k <- (1 / (2 * sigma_sq_err)) * quad_form_delta
    # alpha_k <- (n / 2) + 1
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
