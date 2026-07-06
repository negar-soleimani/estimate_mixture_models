# 
# rm(list = ls())
# source("data/prepare_data.R")
# source("scripts/physics_model.R")
# source("scripts/helper_function_CGP.R")
# source("scripts/main_function_seuil_CGP.R")
# 
# #set.seed(12345)
# n_iter <- 20000
# burn_in <- 2000
# sigma_props <- c(NA, NA, NA, NA, 0.5, NA)
# 
# # init = c(g, h0, sig2err, alpha, psidelta, k)
# init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)
# 
# res <- mcmc_step6(
#   y = y, t = t, n_iter = n_iter, init = init, sigma_proposals = sigma_props,
#   g_init = FALSE, 
#   h0_init = FALSE,
#   sig2er_init = FALSE,
#   alpha_init = FALSE,
#   psi_init = FALSE,
#   k_init = FALSE,
#   Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
#   n_burnin = burn_in,
#   seuil = TRUE,  
#   s = 0.3       
# )


###############################################################################
rm(list = ls())
source("data/prepare_data.R")        
source("scripts/physics_model.R")
source("scripts/helper_function_CGP.R")
source("scripts/main_function_seuil_CGP.R")

library(coda)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(patchwork)

# =============================================================================
#  REAL-DATA ANALYSIS  --  Scenario (I): seuil = FALSE, g sampled
#
#  seuil   = FALSE
#  g_init  = FALSE      (g is sampled)
#
#  Convergence is checked on the parameters that actually move; constant
#  parameters (whose R-hat is undefined) are excluded from pars_check.
# =============================================================================

# ---- Gelman-Rubin convergence check -----------------------------------------
check_gelman <- function(chain_list,
                         pars     = c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k"),
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
  
  is_constant <- !is.finite(rhat_upper)
  if (any(is_constant)) {
    cat("  Note: dropping from R-hat check (constant chains):",
        paste(names(rhat_upper)[is_constant], collapse = ", "), "\n")
  }
  rhat_for_check <- rhat_upper[!is_constant]
  ok <- length(rhat_for_check) > 0 && all(rhat_for_check < rhat_lim)
  
  list(ok = ok, rhat = rhat_upper)
}

# ---- Multi-chain MCMC with Gelman-Rubin extension ---------------------------
run_mcmc_gelman <- function(y, t,
                            init_list,
                            sigma_props,
                            Sigma_theta,
                            min_iter   = 10000,
                            step_iter  = 2000,
                            max_iter   = 200000,
                            rhat_lim   = 1.1,
                            seuil      = FALSE,
                            s          = 0.3,
                            pars_check = c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")) {
  
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
        g_init      = FALSE,
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
#  c(g, h0, sig2err, alpha, psidelta, k)
make_inits <- function() {
  list(
    chain1 = c(9.8,  46.45, 0.01, 0.5, 0.5, 0.1),
    chain2 = c(9.8,  46.45, 0.01, 0.5, 0.5, 0.1),
    chain3 = c(9.8,  46.45, 0.01, 0.5, 0.5, 0.1)
  )
}

# =============================================================================
#  RUN THE MCMC
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
  max_iter    = 200000,
  rhat_lim    = 1.1,
  seuil       = FALSE,
  s           = 0.3,
  pars_check  = c("g", "h0", "sigma_sq_err", "alpha", "psi_delta", "k")
)

cat("\n========================\n")
cat("Scenario (I): seuil = FALSE, g sampled finished.\n")
cat("Total iterations per chain:", res$n_iter_used, "\n")
cat("Converged                  :", res$converged, "\n")
cat("R-hat (upper C.I.):\n")
print(round(res$rhat, 4))

# =============================================================================
#  POOL THE THREE CHAINS FOR POSTERIOR SUMMARIES
# =============================================================================

burn_per_chain <- 2000

theta_mat <- do.call(rbind, lapply(res$theta_list, function(M) {
  M[-(1:burn_per_chain), , drop = FALSE]
}))
delta_mat <- do.call(rbind, lapply(res$delta_list, function(M) {
  M[-(1:burn_per_chain), , drop = FALSE]
}))
zeta_mat  <- do.call(rbind, lapply(res$zeta_list, function(M) {
  M[-(1:burn_per_chain), , drop = FALSE]
}))

g     <- theta_mat[, "g"]
h0    <- theta_mat[, "h0"]
sigma <- theta_mat[, "sigma_sq_err"]
alpha <- theta_mat[, "alpha"]
psi   <- theta_mat[, "psi_delta"]
k     <- theta_mat[, "k"]

# =============================================================================
#  TRACE PLOTS (three chains overlaid)
# =============================================================================

plot_trace_chains <- function(param, label, true_val = NULL) {
  chains <- lapply(res$theta_list, function(M) M[, param])
  n_it   <- length(chains[[1]])
  ylim   <- range(unlist(chains))
  cols   <- c("steelblue", "darkorange", "forestgreen")
  plot(NA, xlim = c(1, n_it), ylim = ylim,
       xlab = "iter", ylab = label,
       main = bquote("trace: " * .(label)))
  for (c in seq_along(chains)) lines(chains[[c]], col = cols[c], lwd = 0.4)
  if (!is.null(true_val)) abline(h = true_val, lty = 2)
}

par(
  mfrow = c(2, 3),
  mar   = c(4, 5, 2, 1),
  
  # axis-title size (equivalent to axis.title)
  cex.lab = 1.0,
  
  # axis-number size (equivalent to axis.text)
  cex.axis = 1.2,
  
  # plot title size
  cex.main = 1.0
)
par(mfrow = c(2, 3), mar = c(3, 4, 2, 1))
plot_trace_chains("g",            "g",                 true_val = 9.8)
plot_trace_chains("h0",           "h0")
plot_trace_chains("sigma_sq_err", expression(lambda^2))
plot_trace_chains("alpha",        expression(alpha))
plot_trace_chains("psi_delta",    expression(psi[delta]))
plot_trace_chains("k",            "k")

library(patchwork)
# =============================================================================
#  POSTERIOR DENSITY OF ALPHA
# =============================================================================

df_alpha <- data.frame(alpha = alpha)

p_alpha <- ggplot(df_alpha, aes(x = alpha)) +
  geom_density(
    fill = "#4CCDC9",
    color = "lightseagreen",
    alpha = 0.6,
    linewidth = 1.2
  ) +
  labs(
    title = NULL,
    x = expression(alpha),
    y = "Density"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    
    # =========================
    # AXIS TITLES (x and y)
    # =========================
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 12),
    
    # =========================
    # AXIS TEXT (tick numbers)
    # =========================
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 12),
    
    # =========================
    # CLEAN TITLE
    # =========================
    plot.title = element_blank()
  )
# =============================================================================
#  POINTWISE INCLUSION PROBABILITIES
# =============================================================================

p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(
  i     = 1:n,
  x     = t,
  p_hat = p_hat
)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = NULL,
    x = "Index i",
    y = expression(hat(p)[i])
  ) +
  theme_minimal(base_size = 20) +
  theme(
    
    # =========================
    # AXIS TITLES (x and y)
    # =========================
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 12),
    
    # =========================
    # AXIS TEXT (tick numbers)
    # =========================
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 12),
    
    # =========================
    # CLEAN TITLE
    # =========================
    plot.title = element_blank()
  )
# =============================================================================
#  POSTERIOR SUMMARY OF delta(x_i)
# =============================================================================

delta_mean <- colMeans(delta_mat)

delta_ci <- t(
  apply(
    delta_mat,
    2,
    quantile,
    probs = c(0.025, 0.975)
  )
)

df_delta <- data.frame(
  i      = 1:n,
  x      = t,
  d_mean = delta_mean,
  d_low  = delta_ci[, 1],
  d_high = delta_ci[, 2]
)

p_delta <- ggplot(df_delta, aes(x = i, y = d_mean)) +
  geom_ribbon(
    aes(ymin = d_low, ymax = d_high),
    fill = "lightseagreen",
    alpha = 0.25
  ) +
  geom_line(
    color = "lightseagreen",
    linewidth = 1.1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = NULL,
    x = "Index i",
    y = expression(delta(x[i]))
  ) +
  theme_minimal(base_size = 20) +
  theme(
    
    # =========================
    # AXIS TITLES (x and y)
    # =========================
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 12),
    
    # =========================
    # AXIS TEXT (tick numbers)
    # =========================
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 12),
    
    # =========================
    # CLEAN TITLE
    # =========================
    plot.title = element_blank()
  )


# =============================================================================
#  PER-PARAMETER BOXPLOTS  +  PATCHWORK COMPOSITES
# =============================================================================

make_box <- function(x, ylab, ref = NULL) {
  
  df <- data.frame(val = x)
  
  p <- ggplot(df, aes(x = 1, y = val)) +
    geom_boxplot(
      fill = "lightseagreen",
      width = 0.25,
      outlier.size = 1.5
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    theme_minimal(base_size = 20) +
    theme(
      
      # =========================
      # AXIS TITLES (x and y)
      # =========================
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 12),
      
      # =========================
      # AXIS TEXT (tick numbers)
      # =========================
      axis.text.x  = element_blank(),
      axis.text.y  = element_text(size = 12),
      
      axis.ticks.x = element_blank(),
      
      # =========================
      # CLEAN TITLE
      # =========================
      plot.title = element_blank(),
      
      # spacing
      plot.margin = margin(3, 6, 3, 6)
    )
  
  if (!is.null(ref)) {
    p <- p +
      geom_hline(
        yintercept = ref,
        linetype = "dashed",
        color = "orange",
        linewidth = 0.8
      )
  }
  
  p
}

p_g     <- make_box(g,     ylab = expression(g[e]))
p_h0    <- make_box(h0,    ylab = "h0")
p_sig   <- make_box(sigma, ylab = expression(lambda^2))
p_a_box <- make_box(alpha, ylab = expression(alpha))
p_psi   <- make_box(psi,   ylab = expression(gamma[delta]))
p_k     <- make_box(k,     ylab = "k")
left_grid <- (p_g   | p_h0) /
  (p_sig | p_a_box) /
  (p_psi | p_k)

right_col_1 <- p_alpha / p_delta + plot_layout(heights = c(1, 1.25))
fig_all     <- left_grid | right_col_1
fig_all     <- fig_all + plot_layout(widths = c(1.05, 1.35))

right_col_2 <- p_alpha / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))
fig_all2    <- left_grid | right_col_2
fig_all2    <- fig_all2 + patchwork::plot_layout(widths = c(1.05, 1.45))

fig_all2


