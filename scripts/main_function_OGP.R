
mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals,
                       g_init=TRUE, h0_init= TRUE, sig2er_init = TRUE,
                       alpha_init = TRUE, psi_init = TRUE, k_init = TRUE, Sigma_theta, n_burnin=1000,
                       continue_chain = FALSE, last_delta = NULL) {
  # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T")
  
  # Total iterations = burn-in + desired samples
  total_iter <- n_burnin + n_iter
  n <- length(y)
  
  theta <- init
  #delta <- rep(0, length(y))
  names(theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  
  if (continue_chain) {
    if (is.null(last_delta)) {
      stop("continue_chain = TRUE but last_delta is NULL.")
    }
    if (length(last_delta) != n) {
      stop("last_delta has wrong length.")
    }
    delta <- as.numeric(last_delta)
  } else {
    delta <- rep(0, n)
  }
  
  chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
  colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
  chain_zeta <- matrix(NA, nrow = total_iter, ncol = length(y))
  
  loglik_chain <- numeric(total_iter)
  accept_psi <- 0
  
  safe_solve <- function(M, jitter = 1e-8) {
    M <- 0.5 * (M + t(M))
    solve(M + jitter * diag(nrow(M)))
  }
  
  d_free <- as.integer(!g_init) + as.integer(!h0_init)
  
  for (iter in 1:total_iter) {
    g <- theta[1]; h0 <- theta[2]; sigma_sq_err <- theta[3]
    alpha_param <- theta[4]; psi_delta <- theta[5]; k <- theta[6]
    sigma_sq_delta <- sigma_sq_err / k
    
    Sigma_delta <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta) #c*(x,x')
    #Sigma_delta <- 0.5 * (Sigma_delta + t(Sigma_delta)) + diag(1e-10, n)
    
    f_theta <- balldropg(t, c(g, h0))
    mean1 <- f_theta
    mean2 <- f_theta + delta
    
    #-------------------------------- probability of the zeta --------------------------------# 
    
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

    #log_likelihood <- sum(log(ifelse(zeta == 1,
    #                                 alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
    #                                 (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
    #loglik_chain[iter] <- log_likelihood
    
    log_likelihood <- sum(ifelse(
      zeta == 1,
      log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE),
      log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
    ))
    loglik_chain[iter] <- log_likelihood
    
    #-------------------------------- delta(discrepancy) --------------------------------#  
    
    zeta_2_indices <- which(zeta == 2)
    if (length(zeta_2_indices) > 0) {
      y_m <- y[zeta_2_indices]
      Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
        Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
      Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
      #Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
      Sigma_inv <- safe_solve(Sigma_delta_ymym)
      mu_delta_hat <- Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
      Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
      #Sigma_delta_hat <- 0.5 * (Sigma_delta_hat + t(Sigma_delta_hat))
      delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    } else {
      delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
    }
    
    #-------------------------------- Gibbs step for theta --------------------------------#
    # new theta
    zeta_1_indices <- which(zeta == 1)
    zeta_2_indices <- which(zeta == 2)
    
    X <- cbind(1, -0.5 * (t * t_range + t_min)^2)
    x1 <- X[zeta_1_indices, , drop = FALSE]
    x2 <- X[zeta_2_indices, , drop = FALSE]
    
    y1 <- matrix(y[zeta_1_indices], ncol = 1)
    y2 <- matrix(y[zeta_2_indices], ncol = 1)
    d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    
    if (g_init && !h0_init) {
      # CASE: g fixed at init[1], sample h0 alone from its conditional
      g_fixed <- init[1]
      # adjust y by subtracting the fixed g term
      y1_adj <- y1 - (x1[, 2, drop = FALSE] * g_fixed)
      y2_adj <- y2 - d2 - (x2[, 2, drop = FALSE] * g_fixed)
      
      sigma_post_h0 <- sigma_sq_err / n
      mu_post_h0    <- (sum(y1_adj) + sum(y2_adj)) / n
      sd_post_h0 <- sqrt(sigma_post_h0)
      
      h0 <- rnorm(1, mu_post_h0, sd = sd_post_h0)
      g  <- g_fixed
      #theta[1] <- g
      #theta[2] <- h0
      
    } else if (!g_init && h0_init) {
      # CASE: h0 fixed at init[2], sample g alone from its conditional
      h0_fixed <- init[2]
      y1_adj <- y1 - (x1[, 1, drop = FALSE] * h0_fixed)
      mean_adj1 <- x1[, 2, drop = FALSE] * y1_adj
      y2_adj <- y2 - d2 - (x2[, 1, drop = FALSE] * h0_fixed)
      mean_adj2 <- x2[, 2, drop = FALSE] * y2_adj
      mean_adj <- sum(mean_adj1) + sum(mean_adj2)
      
      mean_x_sqr <- sum((x1[, 2, drop = FALSE])^2) + sum((x2[, 2, drop = FALSE])^2)
      
      sigma_post_g <- sigma_sq_err / mean_x_sqr
      sd_post_g <- sqrt(sigma_post_g)
      mu_post_g    <- mean_adj / mean_x_sqr
      
      g  <- rnorm(1, mu_post_g, sd = sd_post_g)
      h0 <- h0_fixed
      #theta[1] <- g
      #theta[2] <- h0
      
    } else if (!g_init && !h0_init) {
      # CASE: both free — current bivariate Gibbs sampling
      
      A <- (t(x1) %*% x1 + t(x2) %*% x2) #/ sigma_sq_err
      B <- (t(x1) %*% y1 + t(x2) %*% y2 - t(x2) %*% d2) #/ sigma_sq_err
      #Sigmapost_theta <- safe_solve(A) %*% B
      #Mupost_theta    <- Sigmapost_theta %*% B
      
      A_inv <- safe_solve(A)
      
      Sigmapost_theta <- sigma_sq_err * A_inv
      Mupost_theta    <- A_inv %*% B
      
      theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
      h0 <- theta_sample[1]; g <- theta_sample[2]
      #theta[1] <- g
      #theta[2] <- h0
      
    } else {
      # CASE: both fixed
      g  <- init[1]
      h0 <- init[2]
      #theta[1] <- g
      #theta[2] <- h0
    }

    #-------------------------------- Gibbs step for sigma_sq_err --------------------------------#   
    
    K_star_psi <- GP_correlation(t, psi_delta)
    inv_Kstar_psi <- tryCatch(chol2inv(chol(K_star_psi)), error = function(e) NULL)
    
    if (is.null(inv_Kstar_psi)) {
      quad_form_delta <- -Inf
    } else {
      quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
    }
    
    if (quad_form_delta == -Inf) next
    
    if(sig2er_init){
      sigma_sq_err <- init[3]
    }else{
    
    f_theta <- balldropg(t, c(g, h0))
    idx1 <- which(zeta == 1)
    idx2 <- which(zeta == 2)
    
    residual1 <- y[idx1] - f_theta[idx1]
    rss1 <- sum(residual1^2)
    
    residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
    rss2 <- sum(residual2^2)
    
    #d <- 2
    # When we consider the Jeffreys prior, I use the following code:
    rate_err <- (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    shape_err <- n + (d_free / 2)
    # When the prior is the sigma parameter of the "inverse gamma distribution", I use the following code:
    # rate_err <- 1 + (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
    # shape_err <- 2 + n
    # sigma_sq_err <- rinvgamma(1, shape = shape_err, scale = rate_err)
    # sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
    #sigma_sq_err <- 1/rgamma(1, shape = shape_err, rate = rate_err)
    #theta[3] <- sigma_sq_err
    sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)}
    
    #if (mcmc_parameters[2] == FALSE) {
    #  sigma_sq_err <- init[3]
    #  theta[3] <- sigma_sq_err
    #}

    #-------------------------------- Gibbs step for psi_delta --------------------------------#   

    if(psi_init){
      psi_delta <- init[5]
    }else{
    # When the prior is the psi_delta of the "uniform(0.1, 1)", I use the following code:
    psi_prop <- rtruncnorm(1, a = 0.1, b = 0.5, mean = psi_delta, sd = sigma_proposals[5])
    
    K_star_prop <- GP_correlation(t, psi_prop)
    Sigma_delta_cur  <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_delta)
    #Sigma_delta_cur <- 0.5 * ( Sigma_delta_cur + t( Sigma_delta_cur))
    #Sigma_delta_cur <-  Sigma_delta_cur + diag(1e-10, nrow(Sigma_delta_cur))
    Sigma_delta_prop <- GP_covariance_star_complete(t, sigma_sq_err, k, psi_prop)
    #Sigma_delta_prop <- 0.5 * (Sigma_delta_prop + t(Sigma_delta_prop))
    #Sigma_delta_prop <- Sigma_delta_prop + diag(1e-10, nrow(Sigma_delta_prop))
    log_prop_current <- log(dtruncnorm(psi_delta, a = 0.1, b = 0.5, mean = psi_prop,  sd = sigma_proposals[5]))
    log_prop_prop    <- log(dtruncnorm(psi_prop,  a = 0.1, b = 0.5, mean = psi_delta, sd = sigma_proposals[5]))
    #log_prior_current <- dunif(psi_delta, min = 0.1, max = 0.5, log = TRUE)
    #log_prior_prop    <- dunif(psi_prop,  min = 0.1, max = 0.5, log = TRUE)
    log_prior_current <- dbeta(psi_delta, shape1 = 7, shape2 = 13, log = TRUE)
    log_prior_prop    <- dbeta(psi_prop,  shape1 = 7, shape2 = 13, log = TRUE)
    log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_cur,  log = TRUE), error = function(e) -Inf)
    log_like_prop    <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
    
    log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + (log_prop_current - log_prop_prop)
    
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      psi_delta <- psi_prop
      #theta[5]  <- psi_delta
      accept_psi <- accept_psi + 1
    }}
    
    #if(mcmc_parameters[3] == FALSE){
    #  psi_delta <- init[5]
    #  theta[5] <- psi_delta
    #}

    #-------------------------------- Gibbs step for k --------------------------------#  
    
    # Kpd <- chol_adapt(K_star_psi, jitter0 = 1e-10, jitter_max = 1e-2)
    # inv_Kstar_psi <- inv_from_chol(Kpd$chol)
    # 
    # quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
    # 
    # quad_form_delta <- max(quad_form_delta, 0)
    # 
    # inv_Kstar <- tryCatch(chol2inv(chol(K_star_psi)), 
    #                       error = function(e) NULL)
    
    K_star_psi <- GP_correlation(t, psi_delta)
    inv_Kstar_psi <- tryCatch(chol2inv(chol(K_star_psi)), error = function(e) NULL)
    
    if (is.null(inv_Kstar_psi)) {
      quad_form_delta <- -Inf
    } else {
      quad_form_delta <- as.numeric(t(delta) %*% inv_Kstar_psi %*% delta)
    }
    
    if (quad_form_delta == -Inf) next
    
    if (k_init) {
      k <- init[6]
    } else {
      alpha_k <- (n / 2) + 1
      beta_k  <- (1 / (2 * sigma_sq_err)) * quad_form_delta
      
      F0 <- pgamma(0.02, shape = alpha_k, rate = beta_k)
      F1 <- pgamma(1, shape = alpha_k, rate = beta_k)
      
      if (F1 <= F0 || !is.finite(F1)) {
        warning("Truncated Gamma CDF for k is numerically degenerate; keeping previous k")
      } else {
        u <- runif(1, F0, F1)
        k <- qgamma(u, shape = alpha_k, rate = beta_k)
      }
    }
    
    #-------------------------------- Gibbs step for alpha --------------------------------#   
    if (alpha_init) {
      alpha_param <- init[4]
    } else {
    alpha_param <- rbeta(1, sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
    #theta[4] <- alpha_param
    }
    
    #if(mcmc_parameters[5] == FALSE){
    #  alpha_param <- init[4]
    #}
    
    theta <- c(
      g = g,
      h0 = h0,
      sigma_sq_err = sigma_sq_err,
      alpha = alpha_param,
      psi_delta = psi_delta,
      k = k
    )
    
    #theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
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
  accept_rate_psi <- accept_psi / total_iter
  return(list(theta = chain_theta, delta = chain_delta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_rate_psi))
}
