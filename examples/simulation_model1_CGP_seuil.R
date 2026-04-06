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
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")
n <- length(y)

# ------ Scenario (I):  seuil == FALSE, g == FALSE -------------#
# ------ Scenario (II):  g == TRUE, seuil == FALSE -------------#
# ------ Scenario (III): g == FALSE, seuil == TRUE -------------#
# ------ Scenario (IIII): seuil == TRUE, g == TRUE -------------#

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

res_obj <- result_scenario_III

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
plot(g,     type = "l", main = "trace: g", xlab = "iter", ylab = "g")
abline(h = 9.8, lty = 2)

plot(h0,    type = "l", main = "trace: h0", xlab = "iter", ylab = "h0")
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
       title = "Posterior densities (Scenario IIII, one dataset)") +
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
       title = "Pooled posterior density of alpha (Scenario IIII)") +
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
       title = "Pointwise inclusion probabilities (Scenario IIII)") +
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
       title = expression("Posterior summary of " * delta(x[i]) ~ "(Scenario IIII)")) +
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
       title = "Posterior predictive (Scenario IIII)") +
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


right_col2 <- p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1))

fig_all2 <- left_grid | right_col2
fig_all2 <- fig_all2 + patchwork::plot_layout(widths = c(1.05, 1.45))

print(fig_all2)


################################################
################################################
library(ggplot2)
library(gridExtra)

res_obj <- result_scenario_I 
V <- length(res_obj$zeta_list)
n <- ncol(res_obj$zeta_list[[1]])

## p_hat_i^(v) و delta_mean_i^(v)
p_hat_mat <- matrix(NA, nrow = n, ncol = V)
delta_mean_mat <- matrix(NA, nrow = n, ncol = V)

for (v in 1:V) {
  zeta_mat  <- res_obj$zeta_list[[v]]   # iter x n
  delta_mat <- res_obj$delta_list[[v]]  # iter x n
  
  # zeta==2 = with discrepancy
  p_hat_mat[, v]      <- colMeans(zeta_mat == 2)
  delta_mean_mat[, v] <- colMeans(delta_mat)
}


df_agg <- data.frame(
  i = 1:n,
  
  p_mean = rowMeans(p_hat_mat),
  p_med  = apply(p_hat_mat, 1, median),
  p_q25  = apply(p_hat_mat, 1, quantile, probs = 0.25),
  p_q75  = apply(p_hat_mat, 1, quantile, probs = 0.75),
  p_q05  = apply(p_hat_mat, 1, quantile, probs = 0.05),
  p_q95  = apply(p_hat_mat, 1, quantile, probs = 0.95),
  
  d_mean = rowMeans(delta_mean_mat),
  d_med  = apply(delta_mean_mat, 1, median),
  d_q25  = apply(delta_mean_mat, 1, quantile, probs = 0.25),
  d_q75  = apply(delta_mean_mat, 1, quantile, probs = 0.75),
  d_q05  = apply(delta_mean_mat, 1, quantile, probs = 0.05),
  d_q95  = apply(delta_mean_mat, 1, quantile, probs = 0.95)
)

p_hat_agg <- ggplot(df_agg, aes(x = i, y = p_med)) +
  geom_ribbon(aes(ymin = p_q25, ymax = p_q75), alpha = 0.25) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed") +
  ylim(0, 1) +
  labs(
    x = "Index i",
    y = expression(hat(p)[i]),
    title = expression("Aggregated pointwise inclusion probabilities across 50 datasets")
  ) +
  theme_minimal(base_size = 13)

print(p_hat_agg)


p_delta_agg <- ggplot(df_agg, aes(x = i, y = d_med)) +
  geom_ribbon(aes(ymin = d_q25, ymax = d_q75), alpha = 0.25) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed") +
  labs(
    x = "Index i",
    y = expression(delta(x[i])),
    title = expression("Aggregated posterior mean of " * delta(x[i]) * " across 50 datasets")
  ) +
  theme_minimal(base_size = 13)

print(p_delta_agg)



grid.arrange(
  p_delta_agg,
  p_hat_agg,
  ncol = 1,
  heights = c(1, 1)
)

df_heat <- data.frame(
  i = rep(1:n, times = V),
  replication = rep(1:V, each = n),
  p_hat = as.vector(t(p_hat_mat))
)

p_heat <- ggplot(df_heat, aes(x = i, y = replication, fill = p_hat)) +
  geom_tile() +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed") +
  labs(
    x = "Index i",
    y = "Replication",
    fill = expression(hat(p)[i]),
    title = "Pointwise inclusion probabilities across replications"
  ) +
  theme_minimal(base_size = 13)

print(p_heat)


param_df <- data.frame(
  replication = 1:V,
  g     = colMeans(res_obj$g_chain),
  h0    = colMeans(res_obj$h0_chain),
  sigma = colMeans(res_obj$sigma_chain),
  alpha = colMeans(res_obj$alpha_chain),
  psi   = colMeans(res_obj$psi_chain),
  k     = colMeans(res_obj$k_chain)
)

print(param_df)

#######################################
#######################################
col_param   <- "#1F9D94"   # teal اصلی برای پارامترها و alpha
col_param_fill <- "#A8DDD8"

col_delta   <- "#7E5A9B"   # purple برای delta
col_delta_fill <- "#D8C6E6"

col_pihat   <- "#E68A2E"   # orange برای p_hat
col_pihat_fill <- "#F4C38A"

col_ref     <- "#C73E1D"   # red/orange برای مرز واقعی
col_zero    <- "grey30"    # خط صفر


make_box <- function(x, ylab, ref = NULL) {
  df <- data.frame(val = x)
  p <- ggplot(df, aes(x = 1, y = val)) +
    geom_boxplot(
      fill = col_param_fill,
      color = col_param,
      width = 0.25,
      outlier.size = 1,
      outlier.color = col_param
    ) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(3, 6, 3, 6),
      panel.grid.minor = element_line(color = "grey90")
    )
  if (!is.null(ref)) {
    p <- p + geom_hline(
      yintercept = ref,
      linetype = "dashed",
      color = col_ref,
      linewidth = 0.8
    )
  }
  p
}

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(
    fill = col_param_fill,
    color = col_param,
    alpha = 0.8,
    linewidth = 1.2
  ) +
  labs(
    x = expression(alpha),
    y = "Density",
    title = "Pooled posterior density of alpha (Scenario II)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_line(color = "grey90")
  ) + theme_paper


p_delta <- ggplot(df_delta, aes(x = i, y = d_mean)) +
  geom_ribbon(
    aes(ymin = d_low, ymax = d_high),
    fill = col_delta_fill,
    alpha = 0.55
  ) +
  geom_line(
    color = col_delta,
    linewidth = 1.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = col_zero,
    linewidth = 0.9
  ) +
  geom_vline(
    xintercept = ndis + 0.5,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1
  ) +
  labs(
    x = "Index i",
    y = expression(delta(x[i])),
    title = expression("Posterior summary of " * delta(x[i]) ~ "(Scenario II)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_line(color = "grey90")
  ) + theme_paper

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(
    fill = col_pihat,
    color = col_pihat,
    width = 0.85
  ) +
  geom_vline(
    xintercept = ndis + 0.5,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1
  ) +
  ylim(0, 1) +
  labs(
    x = "Index i",
    y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
    title = "Pointwise inclusion probabilities (Scenario II)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_line(color = "grey90")
  ) + theme_paper


p_zeta0 <- ggplot(df_p0, aes(x = i, y = p0)) +
  geom_col(
    fill = col_param,
    color = col_param,
    width = 0.85
  ) +
  geom_vline(
    xintercept = ndis + 0.5,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1
  ) +
  ylim(0, 1) +
  labs(
    x = "Index i",
    y = expression(hat(p)[i] == P(zeta[i] == 0 ~ "|" ~ Y)),
    title = expression("Pointwise allocation probabilities to " * M[0])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_line(color = "grey90")
  )+ theme_paper



#############################
############################
############################
## =========================================================
## FULL PLOTTING SCRIPT
## Compatible with your current saved object structure
## =========================================================


library(ggplot2)
library(gridExtra)
library(patchwork)

## ---------------------------------------------------------
## 1) Load your saved results
## ---------------------------------------------------------
## Choose the file you want
#load("/Users/negarsoleimani/Documents/phd/paper1/Simulation/Seuil_Simulation/result_scenario_II.RData")

## Choose scenario object
res_obj <- result_scenario_II
scenario_name <- "Scenario II"

## ---------------------------------------------------------
## 2) Basic information
## ---------------------------------------------------------
## representative dataset index
v <- 50

## number of replications
V <- ncol(res_obj$g_chain)

## number of points
n <- nrow(res_obj$y_obs)

## activation boundary from the envelope
ndis <- sum(res_obj$envelope == 0)

## ---------------------------------------------------------
## 3) Color palette close to your presentation theme
## ---------------------------------------------------------
col_teal        <- "#2A9D8F"
col_teal_fill   <- "#BFE8E3"

col_maroon      <- "#8E244D"
col_maroon_fill <- "#E8C1CF"

col_purple      <- "#6C63A8"
col_purple_fill <- "#D8D3F0"

col_orange      <- "#E76F51"
col_gray        <- "#4F5D75"

theme_perso <- function(base_size = 13){
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_line(linewidth = 0.2),
      panel.grid.major = element_line(linewidth = 0.35),
      strip.text = element_text(face = "bold")
    )
}

## ---------------------------------------------------------
## 4) Extract one representative dataset
## ---------------------------------------------------------
g     <- res_obj$g_chain[, v]
h0    <- res_obj$h0_chain[, v]
sigma <- res_obj$sigma_chain[, v]
alpha <- res_obj$alpha_chain[, v]
psi   <- res_obj$psi_chain[, v]
k     <- res_obj$k_chain[, v]

delta_mat <- res_obj$delta_list[[v]]   # iterations x n
zeta_mat  <- res_obj$zeta_list[[v]]    # iterations x n
y         <- res_obj$y_obs[, v]

n_iter_post <- length(g)

## ---------------------------------------------------------
## 5) Numerical summaries for one dataset
## ---------------------------------------------------------
summ_param <- function(x){
  c(
    mean   = mean(x),
    median = median(x),
    sd     = sd(x),
    q025   = unname(quantile(x, 0.025)),
    q975   = unname(quantile(x, 0.975))
  )
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

## ---------------------------------------------------------
## 6) Trace plots for one dataset
## ---------------------------------------------------------
par(mfrow = c(2, 3), mar = c(3,3,2,1))

plot(g, type = "l", main = "trace: g", xlab = "iter", ylab = "g")
abline(h = 9.8, lty = 2, col = col_orange, lwd = 2)

plot(h0, type = "l", main = "trace: h0", xlab = "iter", ylab = "h0")
plot(sigma, type = "l", main = expression("trace: " * lambda^2),
     xlab = "iter", ylab = expression(lambda^2))
plot(alpha, type = "l", main = expression("trace: " * alpha),
     xlab = "iter", ylab = expression(alpha))
plot(psi, type = "l", main = expression("trace: " * gamma[delta]),
     xlab = "iter", ylab = expression(gamma[delta]))
plot(k, type = "l", main = "trace: k", xlab = "iter", ylab = "k")

## ---------------------------------------------------------
## 7) Posterior densities for one dataset
## ---------------------------------------------------------
df_long <- rbind(
  data.frame(value = g,     par = "g"),
  data.frame(value = h0,    par = "h0"),
  data.frame(value = sigma, par = "lambda2"),
  data.frame(value = alpha, par = "alpha"),
  data.frame(value = psi,   par = "gamma_delta"),
  data.frame(value = k,     par = "k")
)

p_dens_all <- ggplot(df_long, aes(x = value)) +
  geom_density(fill = col_teal_fill, color = col_teal, alpha = 0.7, linewidth = 1) +
  facet_wrap(~par, scales = "free", ncol = 3) +
  labs(
    x = "",
    y = "Density",
    title = paste("Posterior densities (", scenario_name, ", one dataset)", sep = "")
  ) +
  theme_perso(12)

print(p_dens_all)

## ---------------------------------------------------------
## 8) Pooled posterior density of alpha across all datasets
## ---------------------------------------------------------
alpha_vec <- as.vector(res_obj$alpha_chain)

df_alpha_pool <- data.frame(alpha = alpha_vec)

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(fill = col_maroon_fill, color = col_maroon, alpha = 0.75, linewidth = 1.1) +
  labs(
    x = expression(alpha),
    y = "Density",
    title = paste("Pooled posterior density of alpha (", scenario_name, ")", sep = "")
  ) +
  theme_perso(13)

print(p_alpha_pool)

## ---------------------------------------------------------
## 9) Local diagnostic for one dataset
## NOTE:
## In your code:
## zeta == 2  -> model with discrepancy
## zeta == 1  -> model without discrepancy
##
## In the figure label below, we display it as the probability
## of selecting the discrepancy component.
## ---------------------------------------------------------
p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(i = 1:n, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = col_teal, alpha = 0.95) +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed", color = col_orange, linewidth = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Index i",
    y = expression(hat(p)[i]),
    title = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y))
  ) +
  theme_perso(12)

print(p_zeta)

## ---------------------------------------------------------
## 10) Probability of allocation to M0 for one dataset
## ---------------------------------------------------------
prob_zeta_model0 <- colMeans(zeta_mat == 1)

df_p0 <- data.frame(i = 1:n, p0 = prob_zeta_model0)

p_zeta0 <- ggplot(df_p0, aes(x = i, y = p0)) +
  geom_col(fill = col_gray, alpha = 0.9) +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed", color = col_orange, linewidth = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Index i",
    y = expression(P(zeta[i] == 0 ~ "|" ~ Y)),
    title = expression("Pointwise allocation probabilities to " * M[0])
  ) +
  theme_perso(12)

print(p_zeta0)

## ---------------------------------------------------------
## 11) Posterior summary of delta(x_i) for one dataset
## ---------------------------------------------------------
delta_mean <- colMeans(delta_mat)
delta_ci   <- t(apply(delta_mat, 2, quantile, probs = c(0.025, 0.975)))

df_delta <- data.frame(
  i = 1:n,
  d_mean = delta_mean,
  d_low  = delta_ci[, 1],
  d_high = delta_ci[, 2]
)

p_delta <- ggplot(df_delta, aes(x = i, y = d_mean)) +
  geom_ribbon(aes(ymin = d_low, ymax = d_high),
              fill = col_purple_fill, alpha = 0.7) +
  geom_line(color = col_purple, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = col_gray, linewidth = 0.8) +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed", color = col_orange, linewidth = 0.9) +
  labs(
    x = "Index i",
    y = expression(delta(x[i])),
    title = paste("Posterior summary of delta(x_i) (", scenario_name, ")", sep = "")
  ) +
  theme_perso(12)

print(p_delta)

## ---------------------------------------------------------
## 12) Posterior predictive for one dataset
## ---------------------------------------------------------
## If g is fixed in the scenario, the chain will be constant and this still works
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
  geom_point(size = 1.5, color = col_gray) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = col_teal_fill, alpha = 0.6) +
  geom_line(aes(y = m), color = col_teal, linewidth = 1) +
  labs(
    x = "Rescaled time x",
    y = "Height",
    title = paste("Posterior predictive (", scenario_name, ")", sep = "")
  ) +
  theme_perso(12)

print(p_pp)

## ---------------------------------------------------------
## 13) Boxplots for one dataset
## ---------------------------------------------------------
make_box <- function(x, ylab, ref = NULL) {
  df <- data.frame(val = x)
  
  p <- ggplot(df, aes(x = 1, y = val)) +
    geom_boxplot(fill = col_teal, alpha = 0.9, width = 0.25, outlier.size = 1) +
    labs(x = NULL, y = ylab) +
    theme_perso(11) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(3, 6, 3, 6)
    )
  
  if (!is.null(ref)) {
    p <- p + geom_hline(yintercept = ref, linetype = "dashed",
                        color = col_orange, linewidth = 0.8)
  }
  
  p
}

p_g     <- make_box(g,     ylab = "g", ref = 9.8)
p_h0    <- make_box(h0,    ylab = "h0", ref = 46.45)
p_sig   <- make_box(sigma, ylab = expression(lambda^2), ref = 0.01)
p_a_box <- make_box(alpha, ylab = expression(alpha))
p_psi   <- make_box(psi,   ylab = expression(gamma[delta]))
p_k     <- make_box(k,     ylab = "k", ref = 0.1)

left_grid <- (p_g | p_h0) /
  (p_sig | p_a_box) /
  (p_psi | p_k)

## ---------------------------------------------------------
## 14) Combined figure: parameters + alpha + delta + p_hat
## ---------------------------------------------------------
right_col <- p_alpha_pool / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

fig_all <- left_grid | right_col
fig_all <- fig_all + patchwork::plot_layout(widths = c(1.05, 1.45))

print(fig_all)

## ---------------------------------------------------------
## 15) Alternative combined figure with M0 probabilities
## ---------------------------------------------------------
right_col2 <- p_alpha_pool / p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

fig_all2 <- left_grid | right_col2
fig_all2 <- fig_all2 + patchwork::plot_layout(widths = c(1.05, 1.45))

print(fig_all2)

## ---------------------------------------------------------
## 16) Aggregated summaries across all datasets
## ---------------------------------------------------------
p_hat_mat      <- matrix(NA, nrow = n, ncol = V)
delta_mean_mat <- matrix(NA, nrow = n, ncol = V)

for (vv in 1:V) {
  zeta_tmp  <- res_obj$zeta_list[[vv]]
  delta_tmp <- res_obj$delta_list[[vv]]
  
  p_hat_mat[, vv]      <- colMeans(zeta_tmp == 2)
  delta_mean_mat[, vv] <- colMeans(delta_tmp)
}

df_agg <- data.frame(
  i = 1:n,
  
  p_mean = rowMeans(p_hat_mat),
  p_med  = apply(p_hat_mat, 1, median),
  p_q25  = apply(p_hat_mat, 1, quantile, probs = 0.25),
  p_q75  = apply(p_hat_mat, 1, quantile, probs = 0.75),
  
  d_mean = rowMeans(delta_mean_mat),
  d_med  = apply(delta_mean_mat, 1, median),
  d_q25  = apply(delta_mean_mat, 1, quantile, probs = 0.25),
  d_q75  = apply(delta_mean_mat, 1, quantile, probs = 0.75)
)

## aggregated p_hat
p_hat_agg <- ggplot(df_agg, aes(x = i, y = p_med)) +
  geom_ribbon(aes(ymin = p_q25, ymax = p_q75),
              fill = col_teal_fill, alpha = 0.7) +
  geom_line(color = col_teal, linewidth = 1) +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed", color = col_orange, linewidth = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Index i",
    y = expression(hat(p)[i]),
    title = "Aggregated pointwise inclusion probabilities across replications"
  ) +
  theme_perso(13)

print(p_hat_agg)

## aggregated delta
p_delta_agg <- ggplot(df_agg, aes(x = i, y = d_med)) +
  geom_ribbon(aes(ymin = d_q25, ymax = d_q75),
              fill = col_purple_fill, alpha = 0.7) +
  geom_line(color = col_purple, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = col_gray, linewidth = 0.8) +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed", color = col_orange, linewidth = 0.9) +
  labs(
    x = "Index i",
    y = expression(delta(x[i])),
    title = expression("Aggregated posterior summary of " * delta(x[i]))
  ) +
  theme_perso(13)

print(p_delta_agg)

grid.arrange(
  p_delta_agg,
  p_hat_agg,
  ncol = 1,
  heights = c(1, 1)
)

## ---------------------------------------------------------
## 17) Heatmap across replications
## ---------------------------------------------------------
df_heat <- data.frame(
  i = rep(1:n, times = V),
  replication = rep(1:V, each = n),
  p_hat = as.vector(t(p_hat_mat))
)

p_heat <- ggplot(df_heat, aes(x = i, y = replication, fill = p_hat)) +
  geom_tile() +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed", color = col_orange, linewidth = 0.8) +
  scale_fill_gradient(low = "#F7FBFF", high = col_teal, limits = c(0, 1)) +
  labs(
    x = "Index i",
    y = "Replication",
    fill = expression(hat(p)[i]),
    title = "Pointwise inclusion probabilities across replications"
  ) +
  theme_perso(13)

print(p_heat)

## ---------------------------------------------------------
## 18) Posterior mean of parameters across replications
## ---------------------------------------------------------
param_df <- data.frame(
  replication = 1:V,
  g     = colMeans(res_obj$g_chain),
  h0    = colMeans(res_obj$h0_chain),
  sigma = colMeans(res_obj$sigma_chain),
  alpha = colMeans(res_obj$alpha_chain),
  psi   = colMeans(res_obj$psi_chain),
  k     = colMeans(res_obj$k_chain)
)

print(head(param_df, 10))

## ---------------------------------------------------------
## 19) Optional: save figures
## ---------------------------------------------------------
# ggsave("p_alpha_pool.png", p_alpha_pool, width = 8, height = 5, dpi = 300)
# ggsave("p_delta.png", p_delta, width = 8, height = 5, dpi = 300)
# ggsave("p_zeta.png", p_zeta, width = 8, height = 4.5, dpi = 300)
# ggsave("fig_all.png", fig_all, width = 14, height = 9, dpi = 300)
# ggsave("p_heat.png", p_heat, width = 10, height = 6, dpi = 300)

