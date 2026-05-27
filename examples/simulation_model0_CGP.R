source("data/prepare_data.R")
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

don <- read_xlsx("/Users/negar/Documents/phd/estimate_mixture_models-main/data/Ball_drops_data.xlsx", sheet = 2)
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

#t_obs <- don$time
#y_obs_real <- don$Height

#t_min <- min(t_obs)
#t_range <- max(t_obs) - min(t_obs)

#t <- seq(0, 1, length.out = 100)

#n <- length(t)

source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

set.seed(12345)
Sigma_theta <- matrix(c(0.5,0,0,0.5), nrow = 2)
# c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)
sigma_proposals <- c(NA, NA, NA, NA, 0.5, NA)
n_samples       <- 50
burn_in         <- 2000
n_iter          <- 10000
# FALSE= fixed parameter
# mcmc parameter (g,h), sig2err, psidelta, k, alpha, freeze delta-zeta
#mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)

#g            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
#h0           <- matrix(NA, n_iter, n_samples, byrow = FALSE)
#sigma_sq_err <- matrix(NA, n_iter, n_samples, byrow = FALSE)
#alpha        <- matrix(NA, n_iter, n_samples, byrow = FALSE)
#psi_delta    <- matrix(NA, n_iter, n_samples, byrow = FALSE)
#k            <- matrix(NA, n_iter, n_samples, byrow = FALSE)
#loglik_mat   <- matrix(NA, n_iter, n_samples)
#delta_list   <- vector("list", n_samples)
#zeta_list    <- vector("list", n_samples)
#accept_rate  <- numeric(n_samples)

#y_obs <- matrix(NA, length(t), n_samples, byrow = FALSE)
#for (v in 1:n_samples) {
#  y_1 = balldropg(t,c(9.8, 46.45)) + rnorm(n, 0, sqrt(0.01))
#  y_obs[,v] <- y_1
  
#  results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals,
#                        g_init=FALSE, h0_init= FALSE, sig2er_init = FALSE,
#                        alpha_init = FALSE, psi_init = FALSE, k_init = FALSE, Sigma_theta, n_burnin=1000, seuil = FALSE, s = 0.3)
#  #results <- mcmc_step6(y_1, t, n_iter, init, sigma_proposals, mcmc_parameters, Sigma_theta, n_burnin=2000)
#  g[,v] = results$theta[,1]
#  h0[,v]=results$theta[,2]
#  sigma_sq_err[,v]=results$theta[,3]
#  alpha[,v] <- results$theta[,4]
#  psi_delta[,v] <- results$theta[,5]
#  k[,v] <- results$theta[,6]
#  delta_list[[v]]  <- results$delta
#  zeta_list[[v]]   <- results$zeta
#  loglik_mat[, v]  <- results$loglik
#  accept_rate[v]   <- results$accept_rate_psi
#}

#result_m0_sh2_classic_classic_ex <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs, delta_list, zeta_list, loglik_mat, accept_rate)

#result_m0_sh2_classic_classic_200 <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m0_sh2_classic_classic_200,file = "/Users/negar/Documents/phd/Result/Model1/Classic/result_m0_sh2_classic_classic_200.RData")

#result_m0_sh2_classic_classic_100 <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m0_sh2_classic_classic_100,file = "/Users/negar/Documents/phd/Result/Model1/Classic/result_m0_sh2_classic_classic_100.RData")

#result_m0_sh2_classic_classic <- list(g, h0, sigma_sq_err, alpha, psi_delta, k, y_obs = y_obs, delta_list, zeta_list, loglik_mat, accept_rate)
#save(result_m0_sh2_classic_classic,file = "/Users/negar/Documents/phd/Result/Model1/Classic/result_m0_sh2_classic_classic.RData")
# y_obs_m1_sh2_ex <- y_obs
# load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Model_1/final_results/result_m0_sh2_classic_classic.RData")
#View(result_m0_sh2_classic_classic_100)

library(coda)

min_iter    <- 10000   # minimum iterations for every sample
step_iter   <- 2000    # extra iterations if not converged
max_iter    <- 50000
diag_window <- 1000    # check Geweke on last 1000 iterations
keep_last   <- 2000    # posterior summaries from last 2000 iterations

check_geweke <- function(theta_mat,
                         pars = c("g", "h0", "sigma_sq_err"),
                         diag_window = 1000,
                         zlim = 2) {
  
  if (nrow(theta_mat) < diag_window) {
    return(list(ok = FALSE, z = rep(NA, length(pars))))
  }
  
  recent_theta <- tail(theta_mat[, pars, drop = FALSE], diag_window)
  
  z_values <- sapply(pars, function(p) {
    geweke.diag(mcmc(recent_theta[, p]))$z
  })
  
  ok <- all(is.finite(z_values)) && all(abs(z_values) < zlim)
  
  return(list(ok = ok, z = z_values))
}

g_list            <- vector("list", n_samples)
h0_list           <- vector("list", n_samples)
sigma_sq_err_list <- vector("list", n_samples)
alpha_list        <- vector("list", n_samples)
psi_delta_list    <- vector("list", n_samples)
k_list            <- vector("list", n_samples)

delta_list        <- vector("list", n_samples)
zeta_list         <- vector("list", n_samples)
loglik_list       <- vector("list", n_samples)

full_theta_list   <- vector("list", n_samples)
full_delta_list   <- vector("list", n_samples)
full_zeta_list    <- vector("list", n_samples)

n_iter_used       <- numeric(n_samples)
geweke_z_list     <- vector("list", n_samples)
accept_rate_list  <- numeric(n_samples)

y_obs <- matrix(NA, length(t), n_samples, byrow = FALSE)

for (v in 1:n_samples) {
  
  cat("\n========================\n")
  cat("Sample:", v, "\n")
  
  y_1 <- balldropg(t, c(9.8, 46.45)) +
    rnorm(n, 0, sqrt(0.01))
  
  y_obs[, v] <- y_1
  
  converged <- FALSE
  
  current_init  <- init
  current_delta <- NULL
  
  all_theta  <- NULL
  all_delta  <- NULL
  all_zeta   <- NULL
  all_loglik <- NULL
  
  total_done <- 0
  
  while (!converged && total_done < max_iter) {
    
    if (total_done < min_iter) {
      block_iter <- min_iter - total_done
    } else {
      block_iter <- step_iter
    }
    
    cat("Running", block_iter, "more iterations...\n")
    
    results <- mcmc_step6(
      y = y_1,
      t = t,
      n_iter = block_iter,
      init = current_init,
      sigma_proposals = sigma_proposals,
      
      g_init = FALSE,
      h0_init = FALSE,
      sig2er_init = FALSE,
      alpha_init = FALSE,
      psi_init = FALSE,
      k_init = FALSE,
      
      Sigma_theta = Sigma_theta,
      n_burnin = 0,
      
      seuil = FALSE,
      s = 0.3,
      
      continue_chain = TRUE,
      last_delta = current_delta
    )
    
    all_theta  <- rbind(all_theta, results$theta)
    all_delta  <- rbind(all_delta, results$delta)
    all_zeta   <- rbind(all_zeta, results$zeta)
    all_loglik <- c(all_loglik, results$loglik)
    
    current_init <- as.numeric(results$theta[nrow(results$theta), ])
    names(current_init) <- colnames(results$theta)
    
    current_delta <- results$delta[nrow(results$delta), ]
    
    total_done <- nrow(all_theta)
    
    cat("Total iterations so far:", total_done, "\n")
    
    if (total_done >= min_iter) {
      
      geweke_res <- check_geweke(
        theta_mat = all_theta,
        pars = c("g", "h0", "sigma_sq_err"),
        diag_window = diag_window,
        zlim = 2
      )
      
      cat("Geweke z-scores on last", diag_window, "iterations:\n")
      print(round(geweke_res$z, 3))
      
      converged <- geweke_res$ok
    }
  }
  
  if (!converged) {
    warning(paste("Sample", v, "did not converge before max_iter"))
  } else {
    cat("Converged for sample", v, "\n")
  }
  
  n_iter_used[v] <- nrow(all_theta)
  geweke_z_list[[v]] <- geweke_res$z
  accept_rate_list[v] <- results$accept_rate_psi
  
  full_theta_list[[v]] <- all_theta
  full_delta_list[[v]] <- all_delta
  full_zeta_list[[v]]  <- all_zeta
  
  n_keep <- min(keep_last, nrow(all_theta))
  keep_idx <- (nrow(all_theta) - n_keep + 1):nrow(all_theta)
  
  theta_keep <- all_theta[keep_idx, , drop = FALSE]
  delta_keep <- all_delta[keep_idx, , drop = FALSE]
  zeta_keep  <- all_zeta[keep_idx, , drop = FALSE]
  loglik_keep <- all_loglik[keep_idx]
  
  g_list[[v]]            <- theta_keep[, "g"]
  h0_list[[v]]           <- theta_keep[, "h0"]
  sigma_sq_err_list[[v]] <- theta_keep[, "sigma_sq_err"]
  alpha_list[[v]]        <- theta_keep[, "alpha"]
  psi_delta_list[[v]]    <- theta_keep[, "psi_delta"]
  k_list[[v]]            <- theta_keep[, "k"]
  
  delta_list[[v]]  <- delta_keep
  zeta_list[[v]]   <- zeta_keep
  loglik_list[[v]] <- loglik_keep
}

result_m0_sh2_classic_classic_ex <- list(
  g = g_list,
  h0 = h0_list,
  sigma_sq_err = sigma_sq_err_list,
  alpha = alpha_list,
  psi_delta = psi_delta_list,
  k = k_list,
  
  y_obs = y_obs,
  
  delta = delta_list,
  zeta = zeta_list,
  loglik = loglik_list,
  
  full_theta = full_theta_list,
  full_delta = full_delta_list,
  full_zeta = full_zeta_list,
  
  n_iter_used = n_iter_used,
  geweke_z = geweke_z_list,
  accept_rate = accept_rate_list
)

par(mfrow = c(1, 3))
boxplot(sapply(g_list, mean), main = "Posterior mean of g")
abline(h = 9.8, col = "red", lty = 2)

boxplot(sapply(h0_list, mean), main = "Posterior mean of h0")
abline(h = 46.45, col = "red", lty = 2)

boxplot(sapply(sigma_sq_err_list, mean), main = expression("Posterior mean of " * lambda^2))
abline(h = 0.01, col = "red", lty = 2)

plot(n_iter_used, type = "b",
     xlab = "Sample",
     ylab = "Number of iterations used",
     main = "Iterations needed until convergence")
abline(h = 10000, col = "red", lty = 2)

v <- 41

par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))

plot(full_theta_list[[v]][, "g"], type = "l",
     main = "trace: g", xlab = "iteration", ylab = "g")
abline(h = 9.8, col = "red", lty = 2)

plot(full_theta_list[[v]][, "h0"], type = "l",
     main = "trace: h0", xlab = "iteration", ylab = "h0")
abline(h = 46.45, col = "red", lty = 2)

plot(full_theta_list[[v]][, "sigma_sq_err"], type = "l",
     main = expression("trace: " * lambda^2),
     xlab = "iteration", ylab = expression(lambda^2))
abline(h = 0.01, col = "red", lty = 2)

plot(full_theta_list[[v]][, "alpha"], type = "l",
     main = expression("trace: " * alpha),
     xlab = "iteration", ylab = expression(alpha))

plot(full_theta_list[[v]][, "psi_delta"], type = "l",
     main = expression("trace: " * psi[delta]),
     xlab = "iteration", ylab = expression(psi[delta]))

plot(full_theta_list[[v]][, "k"], type = "l",
     main = "trace: k", xlab = "iteration", ylab = "k")

pdf("traceplots_sigma_sq_err_final_45.pdf", width = 10, height = 6)

for (i in 1:n_samples) {
  
  plot(
    full_theta_list[[i]][, "sigma_sq_err"],
    type = "l",
    main = paste("Trace plot - sigma_sq_err - Sample", i),
    xlab = "Iteration",
    ylab = expression(lambda^2)
  )
  abline(h = 0.01, lty = 2, col = "red")
  
}

dev.off()

pdf("g_final_45.pdf", width = 10, height = 6)

for (i in 1:n_samples) {
  
  plot(
    full_theta_list[[i]][, "g"],
    type = "l",
    main = paste("Trace plot - g - Sample", i),
    xlab = "Iteration",
    ylab = expression(g)
  )
  abline(h = 9.8, lty = 2, col = "red")
  
}

dev.off()
pdf("h0_final_45.pdf", width = 10, height = 6)

for (i in 1:n_samples) {
  
  plot(
    full_theta_list[[i]][, "h0"],
    type = "l",
    main = paste("Trace plot - h0 - Sample", i),
    xlab = "Iteration",
    ylab = expression(h0)
  )
  abline(h = 46.45, lty = 2, col = "red")
  
}

dev.off()
