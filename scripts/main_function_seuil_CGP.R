mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals,
                       g_init=TRUE, h0_init= TRUE, sig2er_init = TRUE,
                       alpha_init = TRUE, psi_init = TRUE, k_init = TRUE, Sigma_theta, n_burnin=1000, seuil = FALSE, s = 0.3) {
  
  # Total iterations = burn-in + desired samples
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
    
    log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE) # comp 1
    log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE) # comp 2
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    
    if (seuil) {
      prob_zeta_base <- exp(log_w2 - log_den)
      prob_zeta <- prob_zeta_base * (abs(delta) > s)
      zeta <- ifelse(runif(length(y)) < prob_zeta, 2, 1)
      #zeta <- 1 + (runif(length(y)) < prob_zeta)
    } else {
      prob_zeta <- exp(log_w2 - log_den)
      #zeta <- 1 + (runif(length(y)) < prob_zeta)
      zeta <- ifelse(runif(length(y)) < prob_zeta, 2, 1)
    }
    
    chain_zeta[iter, ] <- zeta
    
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
      # for scenario II
      Sigma_delta_hat <- Sigma_delta_hat + 1e-8 * diag(n)
      delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
    }
    else delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
    
    #-------------------------------- Gibbs step for theta --------------------------------# 
    
    # zeta_1_indices <- which(zeta == 1)
    # zeta_2_indices <- which(zeta == 2)
    # X <- cbind(1, -0.5 * (t * t_range + t_min)^2)
    # x1 <- X[zeta_1_indices, , drop = FALSE]
    # x2 <- X[zeta_2_indices, , drop = FALSE]
    # theta_hat     <- matrix(c(46.45, 9.8), ncol = 1)
    # inv_sigma_theta <- solve(Sigma_theta)
    # A <- ((t(x1) %*% x1) / theta[3]) + ((t(x2) %*% x2) / theta[3]) + inv_sigma_theta
    # Sigmapost_theta <- solve(A)
    # y1 <- matrix(y[zeta_1_indices], ncol = 1)
    # y2 <- matrix(y[zeta_2_indices], ncol = 1)
    # d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    # B <- (t(x1) %*% y1) / theta[3] +
    #   (t(x2) %*% y2) / theta[3] -
    #   (t(x2) %*% d2) / theta[3] +
    #   inv_sigma_theta %*% theta_hat
    # Mupost_theta <- Sigmapost_theta %*% B
    # theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
    # h0 <- theta_sample[1];  g <- theta_sample[2]

    #-------------------------------- Page 22 - part 8.1 --------------------------------#     
    # flat prior on theta 
    
    #zeta_1_indices <- which(zeta == 1)
    #zeta_2_indices <- which(zeta == 2)
    
    ## X <- cbind(1, -0.5 * (t * t_range)^2)
    #X <- cbind(1, -0.5 * (t * t_range + t_min)^2)
    #x1 <- X[zeta_1_indices, , drop = FALSE]
    #x2 <- X[zeta_2_indices, , drop = FALSE]
    
    #y1 <- matrix(y[zeta_1_indices], ncol = 1)
    #y2 <- matrix(y[zeta_2_indices], ncol = 1)
    #d2 <- matrix(delta[zeta_2_indices], ncol = 1)
    
    #A <- (t(x1) %*% x1 + t(x2) %*% x2) / sigma_sq_err
    #B <- (t(x1) %*% y1 + t(x2) %*% y2 - t(x2) %*% d2) / sigma_sq_err
    
    #Sigmapost_theta <- solve(A)
    #Mupost_theta    <- Sigmapost_theta %*% B
    
    #theta_sample <- rmvnorm(1, mean = Mupost_theta,
    #                        sigma = Sigmapost_theta)
    
    #h0 <- theta_sample[1];  g <- theta_sample[2]
    #theta[1] <- g
    #theta[2] <- h0
    
    #if(g_init){
    #  g <- init[1]
    #}else{g <- theta_sample[2]}
    
    #if(h0_init){
    #  h0 <- init[2]
    #}else{h0 <- theta_sample[1]}
    
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
      y1_adj <- y1 - x1[, 2, drop = FALSE] * g_fixed
      y2_adj <- y2 - d2 - x2[, 2, drop = FALSE] * g_fixed
      # x_intercept columns
      x1_int <- x1[, 1, drop = FALSE]   # all ones
      x2_int <- x2[, 1, drop = FALSE]
      
      A_h0 <- as.numeric(t(x1_int) %*% x1_int + t(x2_int) %*% x2_int) / sigma_sq_err
      B_h0 <- as.numeric(t(x1_int) %*% y1_adj + t(x2_int) %*% y2_adj) / sigma_sq_err
      
      sigma_post_h0 <- 1 / A_h0
      mu_post_h0    <- sigma_post_h0 * B_h0
      
      h0 <- rnorm(1, mu_post_h0, sqrt(sigma_post_h0))
      g  <- g_fixed
      theta[1] <- g
      theta[2] <- h0
      
    } else if (!g_init && h0_init) {
      # CASE: h0 fixed at init[2], sample g alone from its conditional
      h0_fixed <- init[2]
      y1_adj <- y1 - x1[, 1, drop = FALSE] * h0_fixed
      y2_adj <- y2 - d2 - x2[, 1, drop = FALSE] * h0_fixed
      x1_slope <- x1[, 2, drop = FALSE]
      x2_slope <- x2[, 2, drop = FALSE]
      
      A_g <- as.numeric(t(x1_slope) %*% x1_slope + t(x2_slope) %*% x2_slope) / sigma_sq_err
      B_g <- as.numeric(t(x1_slope) %*% y1_adj + t(x2_slope) %*% y2_adj) / sigma_sq_err
      
      sigma_post_g <- 1 / A_g
      mu_post_g    <- sigma_post_g * B_g
      
      g  <- rnorm(1, mu_post_g, sqrt(sigma_post_g))
      h0 <- h0_fixed
      theta[1] <- g
      theta[2] <- h0
      
    } else if (!g_init && !h0_init) {
      # CASE: both free — current bivariate Gibbs sampling
      A <- (t(x1) %*% x1 + t(x2) %*% x2) / sigma_sq_err
      B <- (t(x1) %*% y1 + t(x2) %*% y2 - t(x2) %*% d2) / sigma_sq_err
      Sigmapost_theta <- solve(A)
      Mupost_theta    <- Sigmapost_theta %*% B
      theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
      h0 <- theta_sample[1]; g <- theta_sample[2]
      theta[1] <- g
      theta[2] <- h0
      
    } else {
      # CASE: both fixed
      g  <- init[1]
      h0 <- init[2]
      theta[1] <- g
      theta[2] <- h0
    }


    #-------------------------------- Gibbs step for sigma_sq_err(lambda^2) --------------------------------#   
    #-------------------------------- Page 23 - part 8.2 --------------------------------#     
    
    if(sig2er_init){
      sigma_sq_err <- init[3]
    }else{
      R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
      if(n > 0){
        R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
        
        quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
      } 
      else {
        quad_form_delta <- 0
      }
      f_theta <- balldropg(t, c(g, h0))
      idx1 <- which(zeta == 1)
      idx2 <- which(zeta == 2)
      residual1 <- y[idx1] - f_theta[idx1]
      rss1   <- sum(residual1^2)
      
      residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
      rss2   <- sum(residual2^2)
      
      # lambda^2 ~ IG(a_lambda, b_lambda)
      #rate_err  <- b_lambda + 0.5 * (rss1 + rss2 + k * quad_form_delta)
      #shape_err <- n + a_lambda
      #sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
      
      d <- 2
      # When we consider the Jeffreys prior, I use the following code:
      rate_err <- (0.5 * ( rss1 + rss2 + (theta[6] * quad_form_delta)))
      shape_err <- n + (d/2) #n + 1#/2 #4 + n
      # When the prior is the sigma parameter of the "inverse gamma distribution", I use the following code:
      # rate_err <- 1 + (0.5 * ( rss1 + rss2 + (k * quad_form_delta)))
      # shape_err <- 2 + n
      sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)}
    
    
    #-------------------------------- Gibbs step for psi_delta --------------------------------#   
    
    # When the prior is the psi_delta of the "uniform(0.1, 1)", I use the following code:
    
    if(psi_init){
      psi_delta <- init[5]
    }else{
      #psi_prop <- rtruncnorm(1, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5])
      psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
      sigma_sq_delta_prop <- sigma_sq_err / k
      Sigma_delta <- GP_covariance(t, sigma_sq_delta_prop, psi_delta)
      Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
      #log_prop_current <- log(dtruncnorm(psi_delta, a = 0.1, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
      #log_prop_prop <- log(dtruncnorm(psi_prop, a = 0.1, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
      #log_prior_current <- dunif(psi_delta, min = 0.1, max = 1, log = TRUE)
      #log_prior_prop    <- dunif(psi_prop,  min = 0.1, max = 1, log = TRUE)
      log_prop_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
      log_prop_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
      log_prior_current <- dbeta(psi_delta, shape1 = 7, shape2 = 13, log = TRUE)
      log_prior_prop    <- dbeta(psi_prop,  shape1 = 7, shape2 = 13, log = TRUE)
      log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
      log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
      log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + (log_prop_current - log_prop_prop)
      if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
        psi_delta <- psi_prop
        theta[5] <- psi_delta
        accept_psi <- accept_psi + 1
      }}
    
    #-------------------------------- Gibbs step for k --------------------------------#   
    #-------------------------------- Page 24-part 8.3 --------------------------------#     
    R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
    if(n > 0){
      R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
      
      quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
    } 
    else {
      quad_form_delta <- 0
    }
    
    if(k_init){
      k <- init[6]
    }else{
      alpha_k <- (n / 2) + 1
      beta_k <- (1 / (2 * sigma_sq_err)) * quad_form_delta
      for (try_k in 1:100) {
        k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
        #if (k_prop >= 0.1 && k_prop <= 0.9) {
        if (k_prop > 0 && k_prop < 1) {
          k <- k_prop
          theta[6] <- k
          break
        }
      }}
    
    #-------------------------------- Gibbs step for alpha --------------------------------#   
    #-------------------------------- Page 25 - part 8.4 --------------------------------#     
    
    if (alpha_init) {
      alpha_param <- init[4]
    } else {
      if (seuil) {
        zetabis <- zeta
        zetabis[abs(delta) < s] <- 0
        #alpha_prime <- rbeta(1, sum(zetabis == 1) + 0.5, sum(zetabis == 2) + 0.5)
        alpha_param <- rbeta(1, sum(zetabis == 1) + 0.5, sum(zetabis == 2) + 0.5)
      } else {
        alpha_param <- rbeta(1, sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
      }
      theta[4] <- alpha_param
    }
    
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