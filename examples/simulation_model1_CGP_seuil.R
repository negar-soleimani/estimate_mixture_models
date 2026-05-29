source("data/prepare_data.R")

t <- don$time
y <- don$Height
length(t)
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)
a <- t_range
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")
n <- length(y)

# ------ Scenario (I):  seuil == FALSE, g == FALSE -------------#
# ------ Scenario (II):  g == TRUE, seuil == FALSE -------------#
# ------ Scenario (III): g == FALSE, seuil == TRUE -------------#
# ------ Scenario (IIII): seuil == TRUE, g == TRUE -------------#

# ----------------------- Simulation B ------------------------ #

set.seed(12345)

# -------------------------------------------------------------
# inference parameters (the MCMC is not changed)
# -------------------------------------------------------------
k <- 0.1
sigma_sq_err <- 0.01

# Parameters used only to simulate delta_true
sim_k <- 0.01
sim_psi_delta <- 0.2

# discrepancy is ~0 for first ndis points, then active afterwards
ndis <- 20
envelope_type <- "step"

make_envelope <- function(n, ndis, type = c("step","ramp")) {
  type <- match.arg(type)
  if (type == "step") {
    w <- c(rep(0, ndis), rep(1, n - ndis))
  } else {
    u <- seq(0, 1, length.out = n - ndis)
    w <- c(rep(0, ndis), u)
  }
  return(w)
}

w <- make_envelope(n, ndis, envelope_type)

n_samples <- 50
n_iter <- 10000
burn_in <- 2000
sigma_props <- c(NA, NA, NA, NA, 0.5, NA)

# init = c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

g_chain     <- matrix(NA, n_iter, n_samples)
h0_chain    <- matrix(NA, n_iter, n_samples)
sigma_chain <- matrix(NA, n_iter, n_samples)
alpha_chain <- matrix(NA, n_iter, n_samples)
psi_chain   <- matrix(NA, n_iter, n_samples)
k_chain     <- matrix(NA, n_iter, n_samples)
zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, n_iter, n_samples)
accept_rate <- numeric(n_samples)
y_obs       <- matrix(NA, n, n_samples)
delta_true_mat <- matrix(NA, n, n_samples)

draw_delta_true <- function(t, w, sigma_sq_err, sim_k, sim_psi_delta,
                            keep_same_sign = TRUE,
                            force_positive = FALSE,
                            s = 0.3,
                            require_above_threshold = FALSE,
                            max_try = 5000) {
  
  n <- length(t)
  sigma_sq_delta <- sigma_sq_err / sim_k
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  active_idx <- which(w > 0)
  
  for (iter in 1:max_try) {
    tilde_delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
    delta_true  <- w * tilde_delta
    d_act <- delta_true[active_idx]
    
    ok_sign <- TRUE
    if (keep_same_sign) {
      if (force_positive) {
        ok_sign <- all(d_act >= 0)
      } else {
        ok_sign <- all(d_act >= 0) || all(d_act <= 0)
      }
    }
    
    ok_threshold <- TRUE
    if (require_above_threshold) {
      ok_threshold <- all(abs(d_act) > s)
    }
    
    if (ok_sign && ok_threshold) {
      return(delta_true)
    }
  }
  
  stop("No acceptable delta_true found after max_try draws. Try smaller sim_k or larger sim_psi_delta.")
}

for (v in 1:n_samples) {
  
  y0 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  
  # -----------------------------------------------------------
  # Here we replace the old raw GP simulation
  # with a version better suited to the threshold case
  # -----------------------------------------------------------
  delta_true <- draw_delta_true(
    t = t,
    w = w,
    sigma_sq_err = sigma_sq_err,
    sim_k = sim_k,
    sim_psi_delta = sim_psi_delta,
    keep_same_sign = TRUE,        # delta_true remains a single sign in the active zone
    force_positive = TRUE,        # TRUE if impose delta_true > 0
    s = 0.3,
    require_above_threshold = FALSE,  # Set to TRUE if the entire active area is above the threshold
    max_try = 5000
  )
  
  y_1 <- y0 + delta_true
  y_obs[, v] <- y_1
  delta_true_mat[, v] <- delta_true
  
  # -------- Scenario (I), (II), (III), (IV) --------
  res <- mcmc_step6(
    y = y_1, t = t, n_iter = n_iter, init = init, sigma_proposals = sigma_props,
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

result_scenario_IV <- list(
  g_chain = g_chain,
  h0_chain = h0_chain,
  sigma_chain = sigma_chain,
  alpha_chain = alpha_chain,
  psi_chain = psi_chain,
  k_chain = k_chain,
  delta_list = delta_list,
  zeta_list = zeta_list,
  loglik_mat = loglik_mat,
  accept_rate = accept_rate,
  y_obs = y_obs,
  envelope = w,
  delta_true_mat = delta_true_mat,
  sim_k = sim_k,
  sim_psi_delta = sim_psi_delta
)

#save(result_scenario_I, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_I.RData")
#save(result_scenario_II, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_II.RData")
#save(result_scenario_II_psi0.2, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_II_psi0.2.RData")
#save(result_scenario_III, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_III.RData")
#save(result_scenario_IV, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_IV.RData")
#save(result_scenario_IV_psi0.2, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_IV_psi0.2.RData")

source("data/prepare_data.R")

t <- don$time
y <- don$Height
length(t)
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)
a <- t_range

source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

n <- length(y)

# ------ Scenario (I):  seuil == FALSE, g == FALSE -------------#
# ------ Scenario (II):  g == TRUE, seuil == FALSE -------------#
# ------ Scenario (III): g == FALSE, seuil == TRUE -------------#
# ------ Scenario (IIII): seuil == TRUE, g == TRUE -------------#

# ----------------------- Simulation B ------------------------ #

set.seed(12345)

# -------------------------------------------------------------
# inference parameters (the MCMC is not changed)
# -------------------------------------------------------------
k <- 0.1
sigma_sq_err <- 0.01

# Parameters used only to simulate delta_true
sim_k <- 0.01
sim_psi_delta <- 0.2

# discrepancy is ~0 for first ndis points, then active afterwards
ndis <- 20
envelope_type <- "step"

make_envelope <- function(n, ndis, type = c("step","ramp")) {
  type <- match.arg(type)
  if (type == "step") {
    w <- c(rep(0, ndis), rep(1, n - ndis))
  } else {
    u <- seq(0, 1, length.out = n - ndis)
    w <- c(rep(0, ndis), u)
  }
  return(w)
}

w <- make_envelope(n, ndis, envelope_type)

n_samples <- 50

# -------------------------------------------------------------
# convergence control
# -------------------------------------------------------------
library(coda)

min_iter    <- 10000   # minimum iterations for every sample
step_iter   <- 2000    # extra iterations if not converged
max_iter    <- 50000
diag_window <- 1000    # check Geweke on last 1000 iterations
keep_last   <- 2000    # posterior summaries from last 2000 iterations

burn_in <- 0           # burn-in is not used inside blocks
n_iter  <- keep_last   # stored posterior chain length

sigma_props <- c(NA, NA, NA, NA, 0.5, NA)

# init = c(g, h0, sig2err, alpha, psidelta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)

check_geweke <- function(theta_mat,
                         diag_window = 1000,
                         zlim = 2,
                         psi_lower = 0.1,
                         psi_upper = 0.5,
                         k_lower = 0.005,
                         k_upper = 1,
                         eps = 1e-10,
                         fixed_pars = character(0)) {
  
  z_values <- rep(NA_real_, 6)
  names(z_values) <- c("g", "h0", "log_sigma_sq_err",
                       "logit_alpha", "logit_psi_delta", "logit_k")
  
  if (nrow(theta_mat) < diag_window) {
    return(list(ok = FALSE, z = z_values))
  }
  
  recent_theta <- tail(theta_mat, diag_window)
  
  g_chain  <- recent_theta[, "g"]
  h0_chain <- recent_theta[, "h0"]
  
  sigma_chain <- recent_theta[, "sigma_sq_err"]
  alpha_chain <- recent_theta[, "alpha"]
  psi_chain   <- recent_theta[, "psi_delta"]
  k_chain     <- recent_theta[, "k"]
  
  # transformations for constrained parameters
  log_sigma_chain <- log(pmax(sigma_chain, eps))
  
  alpha_safe <- pmin(pmax(alpha_chain, eps), 1 - eps)
  logit_alpha_chain <- log(alpha_safe / (1 - alpha_safe))
  
  psi_scaled <- (psi_chain - psi_lower) / (psi_upper - psi_lower)
  psi_scaled <- pmin(pmax(psi_scaled, eps), 1 - eps)
  logit_psi_chain <- log(psi_scaled / (1 - psi_scaled))
  
  k_scaled <- (k_chain - k_lower) / (k_upper - k_lower)
  k_scaled <- pmin(pmax(k_scaled, eps), 1 - eps)
  logit_k_chain <- log(k_scaled / (1 - k_scaled))
  
  chains_to_check <- list(
    g = g_chain,
    h0 = h0_chain,
    log_sigma_sq_err = log_sigma_chain,
    logit_alpha = logit_alpha_chain,
    logit_psi_delta = logit_psi_chain,
    logit_k = logit_k_chain
  )
  
  for (p in names(chains_to_check)) {
    
    # In Scenario II, g is fixed. It is not a stochastic MCMC chain.
    if (p %in% fixed_pars) {
      z_values[p] <- 0
      next
    }
    
    x <- chains_to_check[[p]]
    
    if (any(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) {
      z_values[p] <- NA_real_
      next
    }
    
    z_values[p] <- tryCatch(
      as.numeric(geweke.diag(mcmc(x))$z),
      error = function(e) NA_real_
    )
  }
  
  ok <- all(is.finite(z_values)) && all(abs(z_values) < zlim)
  
  return(list(ok = ok, z = z_values))
}

g_chain     <- matrix(NA, keep_last, n_samples)
h0_chain    <- matrix(NA, keep_last, n_samples)
sigma_chain <- matrix(NA, keep_last, n_samples)
alpha_chain <- matrix(NA, keep_last, n_samples)
psi_chain   <- matrix(NA, keep_last, n_samples)
k_chain     <- matrix(NA, keep_last, n_samples)

zeta_list   <- vector("list", n_samples)
delta_list  <- vector("list", n_samples)
loglik_mat  <- matrix(NA, keep_last, n_samples)
accept_rate <- numeric(n_samples)

y_obs          <- matrix(NA, n, n_samples)
delta_true_mat <- matrix(NA, n, n_samples)

full_theta_list <- vector("list", n_samples)
full_delta_list <- vector("list", n_samples)
full_zeta_list  <- vector("list", n_samples)
full_loglik_list <- vector("list", n_samples)

n_iter_used   <- numeric(n_samples)
geweke_z_list <- vector("list", n_samples)

draw_delta_true <- function(t, w, sigma_sq_err, sim_k, sim_psi_delta,
                            keep_same_sign = TRUE,
                            force_positive = FALSE,
                            s = 0.3,
                            require_above_threshold = FALSE,
                            max_try = 5000) {
  
  n <- length(t)
  sigma_sq_delta <- sigma_sq_err / sim_k
  Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  
  active_idx <- which(w > 0)
  
  for (iter in 1:max_try) {
    tilde_delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
    delta_true  <- w * tilde_delta
    d_act <- delta_true[active_idx]
    
    ok_sign <- TRUE
    if (keep_same_sign) {
      if (force_positive) {
        ok_sign <- all(d_act >= 0)
      } else {
        ok_sign <- all(d_act >= 0) || all(d_act <= 0)
      }
    }
    
    ok_threshold <- TRUE
    if (require_above_threshold) {
      ok_threshold <- all(abs(d_act) > s)
    }
    
    if (ok_sign && ok_threshold) {
      return(delta_true)
    }
  }
  
  stop("No acceptable delta_true found after max_try draws. Try smaller sim_k or larger sim_psi_delta.")
}

for (v in 1:n_samples) {
  
  cat("\n========================\n")
  cat("Sample:", v, "\n")
  
  y0 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  
  # -----------------------------------------------------------
  # Here we replace the old raw GP simulation
  # with a version better suited to the threshold case
  # -----------------------------------------------------------
  delta_true <- draw_delta_true(
    t = t,
    w = w,
    sigma_sq_err = sigma_sq_err,
    sim_k = sim_k,
    sim_psi_delta = sim_psi_delta,
    keep_same_sign = TRUE,              # delta_true remains a single sign in the active zone
    force_positive = TRUE,              # TRUE if impose delta_true > 0
    s = 0.3,
    require_above_threshold = FALSE,    # Set to TRUE if the entire active area is above the threshold
    max_try = 5000
  )
  
  y_1 <- y0 + delta_true
  
  y_obs[, v] <- y_1
  delta_true_mat[, v] <- delta_true
  
  converged <- FALSE
  
  current_init  <- init
  current_delta <- NULL
  
  all_theta  <- NULL
  all_delta  <- NULL
  all_zeta   <- NULL
  all_loglik <- NULL
  
  total_done <- 0
  geweke_res <- NULL
  
  while (!converged && total_done < max_iter) {
    
    if (total_done < min_iter) {
      block_iter <- min_iter - total_done
    } else {
      block_iter <- step_iter
    }
    
    cat("Running", block_iter, "more iterations...\n")
    
    # -------- Scenario (II): g == TRUE, seuil == FALSE --------
    res <- mcmc_step6(
      y = y_1,
      t = t,
      n_iter = block_iter,
      init = current_init,
      sigma_proposals = sigma_props,
      
      g_init = TRUE,
      h0_init = FALSE,
      sig2er_init = FALSE,
      alpha_init = FALSE,
      psi_init = FALSE,
      k_init = FALSE,
      
      Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
      n_burnin = 0,
      
      seuil = FALSE,
      s = 0.3,
      
      continue_chain = !is.null(current_delta),
      last_delta = current_delta
    )
    
    all_theta  <- rbind(all_theta, res$theta)
    all_delta  <- rbind(all_delta, res$delta)
    all_zeta   <- rbind(all_zeta, res$zeta)
    all_loglik <- c(all_loglik, res$loglik)
    
    current_init <- as.numeric(res$theta[nrow(res$theta), ])
    names(current_init) <- colnames(res$theta)
    
    current_delta <- res$delta[nrow(res$delta), ]
    
    total_done <- nrow(all_theta)
    
    cat("Total iterations so far:", total_done, "\n")
    
    if (total_done >= min_iter) {
      
      geweke_res <- check_geweke(
        theta_mat = all_theta,
        diag_window = diag_window,
        zlim = 2,
        psi_lower = 0.1,
        psi_upper = 0.5,
        k_lower = 0.005,
        k_upper = 1,
        fixed_pars = c("g")
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
  accept_rate[v] <- res$accept_rate_psi
  
  full_theta_list[[v]] <- all_theta
  full_delta_list[[v]] <- all_delta
  full_zeta_list[[v]]  <- all_zeta
  full_loglik_list[[v]] <- all_loglik
  
  n_keep <- min(keep_last, nrow(all_theta))
  keep_idx <- (nrow(all_theta) - n_keep + 1):nrow(all_theta)
  
  theta_keep  <- all_theta[keep_idx, , drop = FALSE]
  delta_keep  <- all_delta[keep_idx, , drop = FALSE]
  zeta_keep   <- all_zeta[keep_idx, , drop = FALSE]
  loglik_keep <- all_loglik[keep_idx]
  
  g_chain[1:n_keep, v]     <- theta_keep[, "g"]
  h0_chain[1:n_keep, v]    <- theta_keep[, "h0"]
  sigma_chain[1:n_keep, v] <- theta_keep[, "sigma_sq_err"]
  alpha_chain[1:n_keep, v] <- theta_keep[, "alpha"]
  psi_chain[1:n_keep, v]   <- theta_keep[, "psi_delta"]
  k_chain[1:n_keep, v]     <- theta_keep[, "k"]
  
  delta_list[[v]] <- delta_keep
  zeta_list[[v]]  <- zeta_keep
  loglik_mat[1:n_keep, v] <- loglik_keep
}

result_scenario_II <- list(
  g_chain = g_chain,
  h0_chain = h0_chain,
  sigma_chain = sigma_chain,
  alpha_chain = alpha_chain,
  psi_chain = psi_chain,
  k_chain = k_chain,
  delta_list = delta_list,
  zeta_list = zeta_list,
  loglik_mat = loglik_mat,
  accept_rate = accept_rate,
  y_obs = y_obs,
  envelope = w,
  delta_true_mat = delta_true_mat,
  sim_k = sim_k,
  sim_psi_delta = sim_psi_delta,
  
  full_theta = full_theta_list,
  full_delta = full_delta_list,
  full_zeta = full_zeta_list,
  full_loglik = full_loglik_list,
  n_iter_used = n_iter_used,
  geweke_z = geweke_z_list
)

# View(result_scenario_II)
# h0 <- result_scenario_II[["h0_chain"]]
# plot(h0[, 50], type = "l")

#save(result_scenario_II, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_II.RData")
#save(result_scenario_II_psi0.2, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_II_psi0.2.RData")

