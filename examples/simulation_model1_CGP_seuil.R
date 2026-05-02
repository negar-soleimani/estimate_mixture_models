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
sim_psi_delta <- 0.5

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
    require_above_threshold = TRUE,  # Set to TRUE if the entire active area is above the threshold
    max_try = 5000
  )
  
  y_1 <- y0 + delta_true
  y_obs[, v] <- y_1
  delta_true_mat[, v] <- delta_true
  
  # -------- Scenario (I), (II), (III), (IIII) --------
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
    seuil = FALSE,
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
  sim_psi_delta = sim_psi_delta
)

save(result_scenario_II, file = "/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_II.RData")


# ----------------------- Simulation A ------------------------ #
set.seed(12345)
k <- 0.1
sim_psi_delta <- 0.5
sigma_sq_err <- 0.01
sigma_sq_delta <- sigma_sq_err / k

# discrepancy is ~0 for first ndis points, then active afterwards
ndis <- 20      
envelope_type <- "step"

make_envelope <- function(n, ndis, type = c("step","ramp")) {
  type <- match.arg(type)
  if (type == "step") {
    w <- c(rep(0, ndis), rep(1, n - ndis))
  } else {
    # smoothish ramp from 0 to 1 after ndis
    u <- seq(0, 1, length.out = n - ndis)
    w <- c(rep(0, ndis), u)
  }
  return(w)
}

w <- make_envelope(n, ndis, envelope_type)

n_samples <- 10
n_iter <- 10000
burn_in <- 2500
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


for (v in 1:n_samples) {

  y0 <- balldropg(t, c(9.8, 46.45)) + rnorm(n, 0, sqrt(sigma_sq_err))
  
    Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
  tilde_delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
  delta_true  <- w * tilde_delta
  
  y_1 <- y0 + delta_true
  y_obs[, v] <- y_1
  
  # -------- Scenario (I), (II), (III), (IIII) --------
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

# ------ Scenario (I):  seuil == FALSE, g == FALSE -------------#
# result_scenario_I <- list(
#   g_chain = g_chain,
#   h0_chain = h0_chain,
#   sigma_chain = sigma_chain,
#   alpha_chain = alpha_chain,
#   psi_chain = psi_chain,
#   k_chain = k_chain,
#   delta_list = delta_list,
#   zeta_list = zeta_list,
#   loglik_mat = loglik_mat,
#   accept_rate = accept_rate,
#   y_obs = y_obs,
#   envelope = w
# )

# ------ Scenario (II):  g == TRUE, seuil == FALSE -------------#
# result_scenario_II <- list(
#   g_chain = g_chain,
#   h0_chain = h0_chain,
#   sigma_chain = sigma_chain,
#   alpha_chain = alpha_chain,
#   psi_chain = psi_chain,
#   k_chain = k_chain,
#   delta_list = delta_list,
#   zeta_list = zeta_list,
#   loglik_mat = loglik_mat,
#   accept_rate = accept_rate,
#   y_obs = y_obs,
#   envelope = w
# )

# ------ Scenario (III): g == FALSE, seuil == TRUE -------------#
# result_scenario_III <- list(
#   g_chain = g_chain,
#   h0_chain = h0_chain,
#   sigma_chain = sigma_chain,
#   alpha_chain = alpha_chain,
#   psi_chain = psi_chain,
#   k_chain = k_chain,
#   delta_list = delta_list,
#   zeta_list = zeta_list,
#   loglik_mat = loglik_mat,
#   accept_rate = accept_rate,
#   y_obs = y_obs,
#   envelope = w
# )

# ------ Scenario (IIII): seuil == TRUE, g == TRUE -------------#
result_scenario_IIII <- list(
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
  envelope = w
)

