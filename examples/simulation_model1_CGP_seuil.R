source("data/prepare_data.R")
rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)

#don <- read_xlsx("/Users/negar/Documents/phd/estimate_mixture_models-main/data/Ball_drops_data.xlsx", sheet = 2)
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
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")
n <- length(y)

# ------ Scenario (I):  seuil == FALSE, g == FALSE -------------#
# ------ Scenario (II):  g == TRUE, seuil == FALSE -------------#
# ------ Scenario (III): g == FALSE, seuil == TRUE -------------#
# ------ Scenario (IIII): seuil == TRUE, g == TRUE -------------#

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
#save(result_scenario_I, file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_I.RData")
#save(result_scenario_II, file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_II.RData")
#save(result_scenario_III, file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_III.RData")
#save(result_scenario_IIII, file = "/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_IIII.RData")


#load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_I.RData")
#load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_II.RData")
#load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_III.RData")
#load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_IIII.RData")

################################################################################
################################################################################
################################################################################
################################################################################

## =========================================================
## Scenario IIII : seuil == TRUE , g == TRUE
## =========================================================
plot(y_1)
lines(y0)

library(ggplot2)
library(gridExtra)
library(patchwork)

res_obj <- result_scenario_IIII

v <- 10

g     <- res_obj$g_chain[, v]
h0    <- res_obj$h0_chain[, v]
sigma <- res_obj$sigma_chain[, v]
alpha <- res_obj$alpha_chain[, v]
psi   <- res_obj$psi_chain[, v]
k     <- res_obj$k_chain[, v]

delta_mat <- res_obj$delta_list[[v]]   # (n_iter x n)
zeta_mat  <- res_obj$zeta_list[[v]]    # (n_iter x n)
y         <- res_obj$y_obs[, v]

n_iter_post <- length(g)
n <- ncol(delta_mat)

## numerical summaries for 1 dataset

summ_param <- function(x){
  c(mean   = mean(x),
    median = median(x),
    sd     = sd(x),
    q025   = unname(quantile(x, 0.025)),
    q975   = unname(quantile(x, 0.975)))
}

tab_sum <- rbind(
  g     = summ_param(g),
  h0    = summ_param(h0),
  sigma = summ_param(sigma),
  alpha = summ_param(alpha),
  psi   = summ_param(psi),
  k     = summ_param(k)
)

print(round(tab_sum, 4))

hist(alpha)

par(mfrow = c(2, 3), mar = c(3,3,2,1))
plot(g,     type = "l", main = "trace: g", xlab = "iter", ylab = "g", ylim = c(9.30, 9.85))
abline(h = 9.8, lty = 2)

plot(h0,    type = "l", main = "trace: h0", xlab = "iter", ylab = "h0")
abline(h = 46.45)

plot(sigma, type = "l", main = expression("trace: " * lambda^2),
     xlab = "iter", ylab = expression(lambda^2))
plot(alpha, type = "l", main = expression("trace: " * alpha),
     xlab = "iter", ylab = expression(alpha))
plot(psi,   type = "l", main = expression("trace: " * psi[delta]),
     xlab = "iter", ylab = expression(psi[delta]))
plot(k,     type = "l", main = "trace: k", xlab = "iter", ylab = "k")


## posterior densities for one dataset

df_long <- rbind(
  data.frame(value = g,     par = "g"),
  data.frame(value = h0,    par = "h0"),
  data.frame(value = sigma, par = "lambda2"),
  data.frame(value = alpha, par = "alpha"),
  data.frame(value = psi,   par = "psi_delta"),
  data.frame(value = k,     par = "k")
)

p_dens_all <- ggplot(df_long, aes(x = value)) +
  geom_density(fill = "#4CCDC9", color = "lightseagreen", alpha = 0.6, linewidth = 0.9) +
  facet_wrap(~par, scales = "free", ncol = 3) +
  labs(x = "", y = "Density",
       title = "Posterior densities (Scenario IV, one dataset)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_dens_all)


## pooled posterior density of alpha across all datasets

alpha_vec <- as.vector(res_obj$alpha_chain)

df_alpha_pool <- data.frame(alpha = alpha_vec)

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(fill = "#E8C1CF", color = "#8E244D",
               alpha = 0.6, linewidth = 1) +
  labs(x = expression(alpha), y = "Density",
       title = "Pooled posterior density of alpha (Scenario IV)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_alpha_pool)

## zeta==2 means model with discrepancy

p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(i = 1:n, x = t, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(x = "Index i",
       y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
       title = "Pointwise inclusion probabilities (Scenario IV)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_zeta)


## delta(x_i) 

delta_mean <- colMeans(delta_mat)
delta_ci   <- t(apply(delta_mat, 2, quantile, probs = c(0.025, 0.975)))

df_delta <- data.frame(
  i = 1:n,
  x = t,
  d_mean = delta_mean,
  d_low  = delta_ci[, 1],
  d_high = delta_ci[, 2]
)

p_delta <- ggplot(df_delta, aes(x = i, y = d_mean)) +
  geom_ribbon(aes(ymin = d_low, ymax = d_high),
              fill = "#6C63A8", alpha = 0.25) +
  geom_line(color = "lightseagreen", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed") +
  labs(x = "Index i", y = expression(delta(x[i])),
       title = expression("Posterior summary of " * delta(x[i]) ~ "(Scenario IV)")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_delta)


## posterior predictive for one dataset

mu_mat <- matrix(NA, nrow = n_iter_post, ncol = n)

for (tt in 1:n_iter_post) {
  f_tt <- balldropg(t, c(g[tt], h0[tt]))
  mu_mat[tt, ] <- f_tt + (1 - alpha[tt]) * delta_mat[tt, ]
}

pp_mean <- colMeans(mu_mat)
pp_ci   <- t(apply(mu_mat, 2, quantile, probs = c(0.025, 0.975)))

df_pp <- data.frame(
  i  = 1:n,
  x  = t,
  y  = y,
  m  = pp_mean,
  lo = pp_ci[, 1],
  hi = pp_ci[, 2]
)

p_pp <- ggplot(df_pp, aes(x = x, y = y)) +
  geom_point(size = 1.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = 0.25, fill = "lightseagreen") +
  geom_line(aes(y = m), linewidth = 0.9, color = "lightseagreen") +
  labs(x = "Rescaled time x", y = "Height",
       title = "Posterior predictive (Scenario IV)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_pp)

# =========================================================

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

left_grid <- (p_g | p_h0) /
  (p_sig | p_a_box) /
  (p_psi | p_k)

right_col <- p_alpha_pool / p_delta + plot_layout(heights = c(1, 1.25))

fig_all <- left_grid | right_col
fig_all <- fig_all + plot_layout(widths = c(1.05, 1.35))

print(fig_all)


## zeta==1 means model without discrepancy

prob_zeta_model0 <- colMeans(zeta_mat == 1)

df_p0 <- data.frame(
  i = 1:length(prob_zeta_model0),
  p0 = prob_zeta_model0
)

p_zeta0 <- ggplot(df_p0, aes(x = i, y = p0)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(
    x = "Index i",
    y = expression(hat(p)[i] == P(zeta[i] == 0 ~ "|" ~ Y)),
    title = expression("Pointwise allocation probabilities to " * M[0])
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_zeta0)


right_col2 <- p_alpha_pool / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

fig_all2 <- left_grid | right_col2
fig_all2 <- fig_all2 + patchwork::plot_layout(widths = c(1.05, 1.45))

print(fig_all2)
print(right_col2)

right_col2 <- p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1))

fig_all2 <- left_grid | right_col2
fig_all2 <- fig_all2 + patchwork::plot_layout(widths = c(1.05, 1.45))

print(fig_all2)

################################################################################
################################################################################
################################################################################
################################################################################

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

result_scenario_III <- list(
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

