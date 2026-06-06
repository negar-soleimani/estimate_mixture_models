# mcmc_step6 <- function(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000) {
#   # mcmc_parameters : c(theta(g, h0) = "TRUE", sigma_sq_err = "T", psi_delta = "T", k = "T", alpha = "T")
#   
#   # Total iterations = burn-in + desired samples
#   total_iter <- n_burnin + n_iter
#   
#   theta <- init
#   delta <- rep(0, length(y))
#   chain_theta <- matrix(NA, nrow = total_iter, ncol = length(init))
#   colnames(chain_theta) <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
#   chain_delta <- matrix(NA, nrow = total_iter, ncol = length(y))
#   chain_zeta <- matrix(NA, nrow = total_iter, ncol = length(y))
#   
#   loglik_chain <- numeric(total_iter)
#   accept_psi <- 0
#   
#   for (iter in 1:total_iter) {
#     g <- theta[1]; h0 <- theta[2]; sigma_sq_err <- theta[3]
#     alpha_param <- theta[4]; psi_delta <- theta[5]; k <- theta[6]
#     sigma_sq_delta <- sigma_sq_err / k
#     
#     Sigma_delta <- GP_covariance(t, sigma_sq_delta, psi_delta)
#     f_theta <- balldropg(t, c(g, h0))
#     mean1 <- f_theta
#     mean2 <- f_theta + delta
#     s <- 0.3
#     
#     #prob_zeta <- 1 / (1 + (alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err))) /
#     #                    ((1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err))))*(abs(delta)>s)
#     #zeta <- 1 + (runif(length(y)) < prob_zeta)
#     
#     # log_sum_exp <- function(x, y) {
#     #   m <- pmax(x, y)
#     #   m + log(exp(x - m) + exp(y - m))
#     # }
#     # log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE)
#     # log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE)
#     # 
#     # log_denominator <- log_sum_exp(log_w1, log_w2)
#     # 
#     # prob_zeta_base <- exp(log_w2 - log_denominator)
#     # 
#     # prob_zeta <- prob_zeta_base * (abs(delta) > s)
#     # zeta <- 1 + (runif(length(y)) < prob_zeta)
#     # #zeta <- ifelse(runif(length(y)) < prob_zeta, 1, 2)
#     # chain_zeta[iter, ] <- zeta
#     log_w1 <- log(alpha_param) + dnorm(y, mean1, sqrt(sigma_sq_err), log = TRUE) # comp 1
#     log_w2 <- log(1 - alpha_param) + dnorm(y, mean2, sqrt(sigma_sq_err), log = TRUE) # comp 2
#     log_max <- pmax(log_w1, log_w2)
#     log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
#     
#     if (seuil) {
#       prob_zeta_base <- exp(log_w2 - log_den)
#       prob_zeta <- prob_zeta_base * (abs(delta) > s)
#       zeta <- ifelse(runif(length(y)) < prob_zeta, 2, 1)
#       #zeta <- 1 + (runif(length(y)) < prob_zeta)
#     } else {
#       prob_zeta <- exp(log_w2 - log_den)
#       #zeta <- 1 + (runif(length(y)) < prob_zeta)
#       zeta <- ifelse(runif(length(y)) < prob_zeta, 2, 1)
#     }
#     
#     chain_zeta[iter, ] <- zeta
#     
#     
#     log_likelihood <- sum(log(ifelse(zeta == 1,
#                                      alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
#                                      (1 - alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))))
#     loglik_chain[iter] <- log_likelihood
#     
#     zeta_2_indices <- which(zeta == 2)
#     if (length(zeta_2_indices) > 0) {
#       y_m <- y[zeta_2_indices]
#       Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) +
#         Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
#       Sigma_delta_ym <- Sigma_delta[, zeta_2_indices, drop = FALSE]
#       Sigma_inv <- tryCatch(solve(Sigma_delta_ymym), error = function(e) diag(1, nrow(Sigma_delta_ymym)))
#       mu_delta_hat <- rep(0, n) + Sigma_delta_ym %*% Sigma_inv %*% (y_m - f_theta[zeta_2_indices])
#       Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
#       delta <- as.vector(rmvnorm(1, mean = mu_delta_hat, sigma = Sigma_delta_hat))
#     }
#     else delta = as.vector(rmvnorm(1, rep(0, n), sigma = Sigma_delta))
#     
#     # Gibbs for theta
#     zeta_1_indices <- which(zeta == 1)
#     zeta_2_indices <- which(zeta == 2)
#     X <- cbind(1, -0.5 * (t * t_range + t_min)^2)
#     x1 <- X[zeta_1_indices, , drop = FALSE]  
#     x2 <- X[zeta_2_indices, , drop = FALSE]
#     theta_hat     <- matrix(c(46.45, 9.8), ncol = 1)
#     inv_sigma_theta <- solve(Sigma_theta)
#     A <- ((t(x1) %*% x1) / theta[3]) + ((t(x2) %*% x2) / theta[3]) + inv_sigma_theta
#     Sigmapost_theta <- solve(A)
#     y1 <- matrix(y[zeta_1_indices], ncol = 1)
#     y2 <- matrix(y[zeta_2_indices], ncol = 1)
#     d2 <- matrix(delta[zeta_2_indices], ncol = 1)
#     B <- (t(x1) %*% y1) / theta[3] + 
#       (t(x2) %*% y2) / theta[3] -
#       (t(x2) %*% d2) / theta[3] +
#       inv_sigma_theta %*% theta_hat
#     Mupost_theta <- Sigmapost_theta %*% B
#     theta_sample <- rmvnorm(1, mean = Mupost_theta, sigma = Sigmapost_theta)
#     h0 <- theta_sample[1];  g <- theta_sample[2]
#     
#     #if(mcmc_parameters[1] == FALSE){
#     #  g <- init[1]
#     #  h0 <- init[2]
#     #}
#     
#     ## g = fixer
#     if (mcmc_parameters[1] == FALSE) {
#       g <- init[1]
#       h0 <- theta_sample[1]
#     }
#     
#     ## h0 = fixed
#     #if (mcmc_parameters[1] == FALSE) {
#     #  h0 <- init[2] 
#     #  g <- theta_sample[2] 
#     #}
#     
#     # Gibbs for sigma_sq_err
#     
#     R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
#     if(n > 0){
#       R_inv <- tryCatch(solve(R), error = function(e) diag(1, n))
#       
#       quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
#     } 
#     else {
#       quad_form_delta <- 0
#     }
#     
#     f_theta <- balldropg(t, c(g, h0))
#     idx1 <- which(zeta == 1)
#     idx2 <- which(zeta == 2)
#     residual1 <- y[idx1] - f_theta[idx1]
#     rss1   <- sum(residual1^2)
#     
#     residual2 <- y[idx2] - f_theta[idx2] - delta[idx2]
#     rss2   <- sum(residual2^2)
#     
#     rate_err <- 0.5 + (0.5 * ( rss1 + rss2 + (theta[6] * quad_form_delta)))
#     
#     shape_err <- 4 + n
#     sigma_sq_err <- rinvgamma(1, shape = shape_err, rate = rate_err)
#     
#     if(mcmc_parameters[2] == FALSE){
#       sigma_sq_err <- init[3]
#     }
#     
#     # MH for psi_delta
#     psi_prop <- rtruncnorm(1, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5])
#     sigma_sq_delta_prop <- sigma_sq_err / k
#     Sigma_delta_prop <- GP_covariance(t, sigma_sq_delta_prop, psi_prop)
#     log_prop_current <- log(dtruncnorm(psi_delta, a = 0, b = 1, mean = psi_prop, sd = sigma_proposals[5]))
#     log_prop_prop <- log(dtruncnorm(psi_prop, a = 0, b = 1, mean = psi_delta, sd = sigma_proposals[5]))
#     log_prior_current <- dbeta(psi_delta, shape1 = 4, shape2 = 6, log=TRUE)
#     log_prior_prop <- dbeta(psi_prop, shape1 = 5, shape2 = 6, log=TRUE)
#     log_like_current <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta, log = TRUE), error = function(e) -Inf)
#     log_like_prop <- tryCatch(dmvnorm(delta, rep(0, n), Sigma_delta_prop, log = TRUE), error = function(e) -Inf)
#     log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) + log_prop_current - log_prop_prop
#     if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
#       psi_delta <- psi_prop
#       #theta[5] <- psi_delta
#       accept_psi <- accept_psi + 1
#     }
#     
#     if(mcmc_parameters[3] == FALSE){
#       psi_delta <- init[5]
#     }
#     
#     # Gibbs for k
#     R <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
#     alpha_k <- (n / 2) + 1
#     beta_k <- (1 / (2 * sigma_sq_err)) * sum(delta * (solve(R, delta)))
#     for (try_k in 1:100) {
#       k_prop <- rgamma(1, shape = alpha_k, rate = beta_k)
#       #if (k_prop >= 0.1 && k_prop <= 0.9) {
#       if (k_prop > 0 && k_prop < 1) {
#         k <- k_prop
#         break
#       }
#     }
#     
#     
#     if(mcmc_parameters[4] == FALSE){
#       k <- init[6]
#     }
#     
#     # Gibbs step for alpha
#     
#     zetabis <- zeta
#     zetabis[abs(delta)<s] <- 0
#     alpha_prime <- rbeta(1, sum(zetabis == 1) + 0.5, sum(zetabis == 2) + 0.5)
#     alpha_param <- rbeta(1,sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
#     
#     if(mcmc_parameters[5] == FALSE){
#       alpha_param <- init[4]
#     }
#     
#     theta <- c(g, h0, sigma_sq_err, alpha_param, psi_delta, k)
#     chain_theta[iter, ] <- theta
#     chain_delta[iter, ] <- delta
#     loglik_chain[iter] <- sum(log(
#       ifelse(zeta==1,
#              alpha_param * dnorm(y, mean1, sqrt(sigma_sq_err)),
#              (1-alpha_param) * dnorm(y, mean2, sqrt(sigma_sq_err)))
#     ))
#   }
#   
#   # Remove burn-in samples if n_burnin > 0
#   if (n_burnin > 0) {
#     keep_indices <- (n_burnin + 1):total_iter
#     chain_theta <- chain_theta[keep_indices, , drop = FALSE]
#     chain_delta <- chain_delta[keep_indices, , drop = FALSE]
#     chain_zeta <- chain_zeta[keep_indices, , drop = FALSE]
#     loglik_chain <- loglik_chain[keep_indices]
#     
#     # Adjust acceptance rate calculation for post-burn-in period only
#     cat("Acceptance rate for psi_delta:", round(accept_psi / total_iter, 4), "\n")
#   } else {
#     cat("Acceptance rate for psi_delta:", round(accept_psi / n_iter, 4), "\n")
#   }
#   return(list(theta = chain_theta, delta = chain_delta, zeta = chain_zeta, loglik = loglik_chain, accept_rate_psi = accept_psi / n_iter))
# }
# sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
# mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
# Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2)
# ress <- mcmc_step6(y, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=1000)



# rm(list = ls())
# library(truncnorm)
# library(readxl)
# library(ggplot2)
# library(tidyr)
# library(MASS)
# library(mvtnorm)
# library(invgamma)
# 
# don <- read_xlsx("/Users/negar/Documents/phd/estimate_mixture_models-main/data/Ball_drops_data.xlsx", sheet = 2)
# names(don) <- c("drop", "time", "Height", "Velocity")
# don$drop <- as.factor(don$drop)
# don <- don[don$drop == 1, ]
# 
# t <- don$time
# y <- don$Height
# length(t)
# t_min <- min(t)
# t_range <- max(t) - min(t)
# t <- (t - t_min) / t_range
# n <- length(y)
# a <- t_range
rm(list = ls())
source("data/prepare_data.R")
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

set.seed(12345)
n_iter <- 10000
burn_in <- 2000
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)

# init = c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

res <- mcmc_step6(
  y = y, t = t, n_iter = n_iter, init = init, sigma_proposals = sigma_props,
  g_init = TRUE, 
  h0_init = FALSE,
  sig2er_init = FALSE,
  alpha_init = FALSE,
  psi_init = FALSE,
  k_init = FALSE,
  Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
  n_burnin = burn_in,
  seuil = TRUE,  
  s = 0.3       
)


###############################################################################
rm(list = ls())
source("data/prepare_data.R")        # provides y, t, n
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

library(coda)
library(mvtnorm)

# =============================================================================
#  GELMAN-RUBIN CONVERGENCE DIAGNOSTIC  --  Real data
#
#  Same workflow as simulation_model0_CGP_gelman.R but on the real (y, t)
#  loaded by prepare_data.R: one MCMC analysis, three parallel chains,
#  run 10,000 iterations first, then check Gelman-Rubin, then extend in
#  blocks of step_iter until R-hat < 1.1 for every monitored parameter or
#  until max_iter is reached.
# =============================================================================

# ---- Gelman-Rubin convergence check -----------------------------------------
check_gelman <- function(chain_list,
                         pars = c("g", "h0", "sigma_sq_err",
                                  "alpha", "psi_delta", "k"),
                         rhat_lim = 1.1) {
  
  mcmc_list <- mcmc.list(
    lapply(chain_list, function(M) mcmc(M[, pars, drop = FALSE]))
  )
  
  gd <- tryCatch(gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE),
                 error = function(e) NULL)
  
  if (is.null(gd)) {
    return(list(ok = FALSE, rhat = rep(NA, length(pars))))
  }
  
  rhat_upper <- gd$psrf[, "Upper C.I."]
  ok <- all(is.finite(rhat_upper)) && all(rhat_upper < rhat_lim)
  
  list(ok = ok, rhat = rhat_upper)
}

# ---- Multi-chain MCMC with Gelman-Rubin extension ---------------------------
run_mcmc_gelman <- function(y, t,
                            init_list,
                            sigma_props,
                            Sigma_theta,
                            min_iter   = 10000,
                            step_iter  = 2000,
                            max_iter   = 120000,
                            rhat_lim   = 1.1,
                            seuil      = TRUE,
                            s          = 0.3,
                            pars_check = c("g", "h0", "sigma_sq_err",
                                           "alpha", "psi_delta", "k")) {
  
  n_chains <- length(init_list)
  
  all_theta  <- vector("list", n_chains)
  all_delta  <- vector("list", n_chains)
  all_zeta   <- vector("list", n_chains)
  all_loglik <- vector("list", n_chains)
  
  current_init  <- init_list
  current_delta <- vector("list", n_chains)
  accept_psi    <- numeric(n_chains)
  
  total_done <- 0
  converged  <- FALSE
  gelman_res <- list(ok = FALSE, rhat = rep(NA, length(pars_check)))
  
  while (!converged && total_done < max_iter) {
    
    if (total_done < min_iter) {
      block_iter <- min_iter - total_done
    } else {
      block_iter <- step_iter
    }
    
    cat("  Running", block_iter, "more iterations on each of",
        n_chains, "chains...\n")
    
    for (c in seq_len(n_chains)) {
      
      continue_now <- total_done > 0
      last_delta_c <- if (continue_now) current_delta[[c]] else NULL
      
      res <- mcmc_step6(
        y               = y,
        t               = t,
        n_iter          = block_iter,
        init            = current_init[[c]],
        sigma_proposals = sigma_props,
        g_init      = TRUE,
        h0_init     = FALSE,
        sig2er_init = FALSE,
        alpha_init  = FALSE,
        psi_init    = FALSE,
        k_init      = FALSE,
        Sigma_theta     = Sigma_theta,
        n_burnin        = 0,
        seuil           = seuil,
        s               = s,
        continue_chain  = continue_now,
        last_delta      = last_delta_c
      )
      
      all_theta[[c]]  <- rbind(all_theta[[c]],  res$theta)
      all_delta[[c]]  <- rbind(all_delta[[c]],  res$delta)
      all_zeta[[c]]   <- rbind(all_zeta[[c]],   res$zeta)
      all_loglik[[c]] <- c(all_loglik[[c]],     res$loglik)
      
      current_init[[c]] <- as.numeric(res$theta[nrow(res$theta), ])
      names(current_init[[c]]) <- colnames(res$theta)
      current_delta[[c]] <- res$delta[nrow(res$delta), ]
      
      accept_psi[c] <- accept_psi[c] + res$accept_rate_psi * block_iter
    }
    
    total_done <- nrow(all_theta[[1]])
    cat("  Total iterations so far:", total_done, "\n")
    
    if (total_done >= min_iter) {
      gelman_res <- check_gelman(all_theta, pars = pars_check, rhat_lim = rhat_lim)
      cat("  Gelman-Rubin R-hat (upper C.I.) on the parameters:\n")
      print(round(gelman_res$rhat, 4))
      converged <- gelman_res$ok
    }
  }
  
  list(
    theta_list   = all_theta,
    delta_list   = all_delta,
    zeta_list    = all_zeta,
    loglik_list  = all_loglik,
    accept_rate  = accept_psi / total_done,
    n_iter_used  = total_done,
    rhat         = gelman_res$rhat,
    converged    = converged
  )
}

# ---- Over-dispersed starting points for the 3 chains ------------------------
#  Same parameter ordering as in the simulation files:
#  c(g, h0, sig2err, alpha, psidelta, k)
make_inits <- function() {
  list(
    chain1 = c(9.8,  46.45, 0.01, 0.5, 0.5, 0.1),
    chain2 = c(9.8,  46.45, 0.01, 0.5, 0.5, 0.1),
    chain3 = c(9.8,  46.45, 0.01, 0.5, 0.5, 0.1)
  )
}

# ---- Trace + density PDF for the real-data run ------------------------------
#  One page, 6 parameters x 2 columns (trace, density) = 6x2 grid.
save_traceplots_pdf <- function(res, filename) {
  
  pdf(filename, width = 12, height = 14)
  
  param_names <- c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
  param_labels <- list(
    g            = expression(g),
    h0           = expression(h[0]),
    sigma_sq_err = expression(lambda^2),
    alpha        = expression(alpha),
    psi_delta    = expression(gamma[delta]),
    k            = expression(k)
  )
  
  chain_mats <- res$theta_list
  n_chains   <- length(chain_mats)
  cols       <- c("steelblue", "darkorange", "forestgreen")
  
  op <- par(mfrow = c(6, 2),
            mar   = c(3.2, 4, 2.0, 1),
            mgp   = c(2.1, 0.7, 0),
            oma   = c(0, 0, 3.0, 0))
  
  for (p in param_names) {
    
    n_iter_p   <- nrow(chain_mats[[1]])
    vals_chain <- lapply(chain_mats, function(M) M[, p])
    ylim_p     <- range(unlist(vals_chain), na.rm = TRUE)
    
    # ---- (a) trace plot ----
    plot(NA, xlim = c(1, n_iter_p), ylim = ylim_p,
         xlab = "iteration",
         ylab = param_labels[[p]],
         main = bquote("Trace: " * .(param_labels[[p]])))
    
    for (c in seq_len(n_chains)) {
      lines(vals_chain[[c]], col = cols[c], lwd = 0.5)
    }
    
    # ---- (b) density plot ----
    dens_list <- lapply(vals_chain, density)
    xlim_d <- range(sapply(dens_list, function(d) range(d$x)))
    ylim_d <- range(sapply(dens_list, function(d) range(d$y)))
    
    plot(NA, xlim = xlim_d, ylim = ylim_d,
         xlab = param_labels[[p]],
         ylab = "Density",
         main = bquote("Density: " * .(param_labels[[p]])))
    
    for (c in seq_len(n_chains)) {
      lines(dens_list[[c]], col = cols[c], lwd = 1.4)
    }
  }
  
  mtext(paste0("Real data  |  iterations = ", res$n_iter_used,
               "  |  R-hat max = ",
               round(max(res$rhat, na.rm = TRUE), 3),
               "  |  converged = ", res$converged),
        outer = TRUE, cex = 1.1, font = 2)
  
  par(op)
  dev.off()
}

# =============================================================================
#  RUN THE REAL-DATA MCMC
# =============================================================================

set.seed(12345)

sigma_props <- c(NA, NA, NA, NA, 0.5, NA)
Sigma_theta <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)

init_list <- make_inits()

res <- run_mcmc_gelman(
  y           = y,
  t           = t,
  init_list   = init_list,
  sigma_props = sigma_props,
  Sigma_theta = Sigma_theta,
  min_iter    = 10000,
  step_iter   = 2000,
  max_iter    = 120000,
  rhat_lim    = 1.1,
  seuil       = TRUE,
  s           = 0.3,
  pars_check  = c("g", "h0", "sigma_sq_err",
                  "alpha", "psi_delta", "k")
)

cat("\n========================\n")
cat("Real-data MCMC finished.\n")
cat("Total iterations per chain:", res$n_iter_used, "\n")
cat("Converged                  :", res$converged, "\n")
cat("R-hat (upper C.I.):\n")
print(round(res$rhat, 4))

# ---- Save chains and trace+density PDF --------------------------------------
out_dir <- "/Users/negarsoleimani/Documents/phd/paper1/github/real_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

result_real_data <- res
save(result_real_data,
     file = file.path(out_dir, "result_real_data.RData"))

save_traceplots_pdf(res,
                    filename = file.path(out_dir, "traceplots_real_data.pdf"))

message("Done. Files written to: ", out_dir)




###################### plot ############################



# =============================================================================
# Posterior plots for the real-data analysis (Gelman-Rubin version).
#
# Input: object `res` returned by run_mcmc_gelman() in real_data_CGP_gelman.R
#   - res$theta_list[[c]] : (n_iter x 6) matrix, columns
#       g, h0, sigma_sq_err, alpha, psi_delta, k          (c = 1..3 chains)
#   - res$delta_list[[c]] : (n_iter x n) matrix
#   - res$zeta_list[[c]]  : (n_iter x n) matrix
#
# We pool the three chains (row-bind) to recover the same flat matrices
# the original single-chain script used. Then every plotting block below is
# UNCHANGED from your original code.
# =============================================================================

library(ggplot2)
library(gridExtra)
library(patchwork)

# Uncomment if running from a fresh R session:
# load("/Users/negarsoleimani/Documents/phd/paper1/github/real_data/result_real_data.RData")
# res <- result_real_data

# =============================================================================
#  POOL THE THREE CHAINS  (the only new block)
# =============================================================================
theta_mat <- do.call(rbind, res$theta_list)   # (3 * n_iter) x 6
delta_mat <- do.call(rbind, res$delta_list)   # (3 * n_iter) x n
zeta_mat  <- do.call(rbind, res$zeta_list)    # (3 * n_iter) x n

g     <- theta_mat[, "g"]
h0    <- theta_mat[, "h0"]
sigma <- theta_mat[, "sigma_sq_err"]
alpha <- theta_mat[, "alpha"]
psi   <- theta_mat[, "psi_delta"]
k     <- theta_mat[, "k"]

n <- ncol(delta_mat)

# =============================================================================
#  EVERYTHING BELOW IS YOUR ORIGINAL CODE, UNCHANGED
# =============================================================================

## ---------- Trace plots ----------
par(mfrow = c(2, 3))
plot(g,     type = "l", main = "trace: g",  xlab = "iter", ylab = "g"); abline(h = 9.8, lty = 2)
plot(h0,    type = "l", main = "trace: h0", xlab = "iter", ylab = "h0")
plot(sigma, type = "l", main = expression("trace: " * lambda^2),
     xlab = "iter", ylab = expression(lambda^2))
plot(alpha, type = "l", main = expression("trace: " * alpha),
     xlab = "iter", ylab = expression(alpha))
plot(psi,   type = "l", main = expression("trace: " * psi[delta]),
     xlab = "iter", ylab = expression(psi[delta]))
plot(k,     type = "l", main = "trace: k",  xlab = "iter", ylab = "k")

# ratio <- sigma / k
# plot(ratio, type = "l")
# boxplot(delta_mat)

n_iter_post <- nrow(theta_mat)

df_alpha <- data.frame(alpha = alpha)
p_alpha <- ggplot(df_alpha, aes(x = alpha)) +
  geom_density(fill = "#4CCDC9", color = "lightseagreen",
               alpha = 0.6, linewidth = 1) +
  labs(x = expression(alpha), y = "Density",
       title = "Posterior density of alpha (real data, Scenario IV)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

## ---------- p_hat_i = P(zeta_i == 1 | Y) ----------
p_hat <- colMeans(zeta_mat == 2)
df_p  <- data.frame(i = 1:n, x = t, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(x = "Index i",
       y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
       title = "Pointwise inclusion probabilities (Scenario IV)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

## ---------- delta(x_i) ----------
delta_mean <- colMeans(delta_mat)
boxplot(delta_mat)
delta_ci   <- t(apply(delta_mat, 2, quantile, probs = c(0.025, 0.975)))

df_delta <- data.frame(
  i      = 1:n,
  x      = t,
  d_mean = delta_mean,
  d_low  = delta_ci[, 1],
  d_high = delta_ci[, 2]
)

p_delta <- ggplot(df_delta, aes(x = i, y = d_mean)) +
  geom_ribbon(aes(ymin = d_low, ymax = d_high),
              fill = "lightseagreen", alpha = 0.25) +
  geom_line(color = "lightseagreen", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Index i", y = expression(delta(x[i])),
       title = expression("Posterior summary of " *
                            delta(x[i]) ~ "(Scenario IV)")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

## ---------- Per-parameter boxplots ----------
make_box <- function(x, ylab, ref = NULL) {
  df <- data.frame(val = x)
  p <- ggplot(df, aes(x = 1, y = val)) +
    geom_boxplot(fill = "lightseagreen", width = 0.25, outlier.size = 1) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(3, 6, 3, 6)
    )
  if (!is.null(ref)) {
    p <- p + geom_hline(yintercept = ref, linetype = "dashed",
                        color = "orange", linewidth = 0.7)
  }
  p
}

p_g     <- make_box(g,     ylab = "g")
p_h0    <- make_box(h0,    ylab = "h0")
p_sig   <- make_box(sigma, ylab = expression(lambda^2))
p_a_box <- make_box(alpha, ylab = expression(alpha))
p_psi   <- make_box(psi,   ylab = expression(gamma[delta]))
p_k     <- make_box(k,     ylab = "k")

left_grid <- (p_g   | p_h0) /
  (p_sig | p_a_box) /
  (p_psi | p_k)

right_col <- p_alpha / p_delta + plot_layout(heights = c(1, 1.25))

fig_all  <- left_grid | right_col
fig_all  <- fig_all + plot_layout(widths = c(1.05, 1.35))

right_col <- p_alpha / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

fig_all2 <- left_grid | right_col
fig_all2 <- fig_all2 + patchwork::plot_layout(widths = c(1.05, 1.45))

fig_all2

