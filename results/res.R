#################################################################################
#################################################################################
###################### Results for model 1 ######################################
#################################################################################
#################################################################################
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh2_presentation.RData")
load("E:Phd_Paris Saclay/resultat/Results_M1_Simulation/y_obs_m1_sh2_presentation.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh1.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh2.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh3.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh4.RData")
load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh5.RData")

sheet_names <- c("Baseball",
                 "Blue Basketball",
                 "Green Basketball",
                 "Volleyball",
                 "Bowling Ball")

#View(result_m1_sh1)

g_sh1 <- result_m1_sh1[[1]]
h0_sh1 <- result_m1_sh1[[2]]
sigma_sq_err_sh1 <- result_m1_sh1[[3]]
alpha_sh1 <- result_m1_sh1[[4]]
psi_delta_sh1 <- result_m1_sh1[[5]]
k_sh1 <- result_m1_sh1[[6]]

g_sh2 <- result_m1_sh2[[1]]
h0_sh2 <- result_m1_sh2[[2]]
sigma_sq_err_sh2 <- result_m1_sh2[[3]]
alpha_sh2 <- result_m1_sh2[[4]]
psi_delta_sh2 <- result_m1_sh2[[5]]
k_sh2 <- result_m1_sh2[[6]]

g_sh3 <- result_m1_sh3[[1]]
h0_sh3 <- result_m1_sh3[[2]]
sigma_sq_err_sh3 <- result_m1_sh3[[3]]
alpha_sh3 <- result_m1_sh3[[4]]
psi_delta_sh3 <- result_m1_sh3[[5]]
k_sh3 <- result_m1_sh3[[6]]

g_sh4 <- result_m1_sh4[[1]]
h0_sh4 <- result_m1_sh4[[2]]
sigma_sq_err_sh4 <- result_m1_sh4[[3]]
alpha_sh4 <- result_m1_sh4[[4]]
psi_delta_sh4 <- result_m1_sh4[[5]]
k_sh4 <- result_m1_sh4[[6]]

g_sh5 <- result_m1_sh5[[1]]
h0_sh5 <- result_m1_sh5[[2]]
sigma_sq_err_sh5 <- result_m1_sh5[[3]]
alpha_sh5 <- result_m1_sh5[[4]]
psi_delta_sh5 <- result_m1_sh5[[5]]
k_sh5 <- result_m1_sh5[[6]]

g_sh2 <- result_m1_sh2_presentation[[1]]
h0_sh2 <- result_m1_sh2_presentation[[2]]
sigma_sq_err_sh2 <- result_m1_sh2_presentation[[3]]
alpha_sh2 <- result_m1_sh2_presentation[[4]]
psi_delta_sh2 <- result_m1_sh2_presentation[[5]]
k_sh2 <- result_m1_sh2_presentation[[6]]

par(mfrow = c(1, 3),
    mar   = c(3, 4, 1, 1) 
)
boxplot(
  colMeans(g_sh2),
  ylab = "g",
  #xlab = "Bowling Ball",
  col  = "#800020",
  main = NULL
)
abline(h = 9.8, lty = 2)

boxplot(
  colMeans(h0_sh2),
  ylab = "h0",
  #xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 46.452, lty = 2)

boxplot(
  colMeans(sigma_sq_err_sh2),
  ylab = expression(sigma[err]^2),
  #xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 0.0537, lty = 2)


boxplot(
  colMeans(alpha_sh5),
  ylab = expression(alpha),
  xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 0.7, lty = 2)

boxplot(
  colMeans(psi_delta_sh5),
  ylab = expression(psi[delta]),
  xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 0.3, lty = 2)

boxplot(
  colMeans(k_sh5),
  ylab = "k",
  xlab = "",
  col  = "#800020",
  main = NULL
)
abline(h = 0.2, lty = 2)

don <- read_xlsx("E:/Phd_Paris Saclay/10Sep2024_Cours Calibration/Ball_drops_data.xls", sheet = 2)
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop)
don <- don[don$drop == 1, ]

t <- don$time
y <- don$Height
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
n <- length(y)

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

g_mat <- result_m1_sh2_presentation[[1]]   # n_iter × n_sims
h0_mat<- result_m1_sh2_presentation[[2]]

g_last <- g_mat[, 30]
h0_last <- h0_mat[, 30]
y_obs <- result_m1_sh2_presentation$y_obs
y_obs30 <- y_obs[,30]
y_true <- balldropg(t, c(9.8, 46.45))

y_pred <- matrix(NA, length(t), length(g_last), byrow = FALSE)

for (i in 1:length(g_last)) {
  theta <- c(g_last[i], h0_last[i])
  y_pred[,i] <- balldropg(t, theta)
}

y_pred1  = t(y_pred)
par(mfrow = c(1,1))
boxplot(y_pred1, col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 1)
legend(x=30, y=45, legend=c("Simulated data",
                            "True code", "Predictions"), lwd=rep(2,2), col=c("gold","blue", "orange2"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))

boxplot(y_pred1[, 1:5], col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 3)
legend(x=30, y=45, legend=c("Simulated data", "Predictions",
                            "True code"), lwd=rep(2,2), col=c("gold", "orange2","blue"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))

y_pred2 <- rowMeans(y_pred) 
y_obs2 <- rowMeans(y_obs)

df_cmp <- data.frame(
  time   = t,
  y_obs2  = y_obs2,
  y_pred2  = y_pred2,
  y_true = y_true
)
library(ggplot2)
p <- ggplot(df_cmp, aes(x = time)) +
  geom_point(aes(y = y_obs2, color = "Observed"),      size = 1) +
  geom_line( aes(y = y_pred2,  color = "Estimate"),     linetype = "dashed", size = 0.5) +
  geom_line( aes(y = y_true, color = "True"),         linetype = "solid",  size = 0.5) +
  scale_color_manual(
    name   = NULL,
    values = c("Observed" = "#800020",
               "Estimate" = "green",
               "True"     = "blue")
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x     = "Time",
    y     = ""
  ) +
  theme(
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.6)),
    legend.key           = element_rect(fill = NA)
  )

print(p)

boxplot(colMeans(g_sh1), colMeans(g_sh2), colMeans(g_sh3), colMeans(g_sh4), colMeans(g_sh5),
        names  = sheet_names,
        ylab = "g",
        col    = rainbow(5))

boxplot(colMeans(h0_sh1), colMeans(h0_sh2), colMeans(h0_sh3), colMeans(h0_sh4), colMeans(h0_sh5),
        names  = sheet_names,
        ylab = "h0",
        col    = rainbow(5))

boxplot(colMeans(sigma_sq_err_sh1), colMeans(sigma_sq_err_sh2), colMeans(sigma_sq_err_sh3), 
        colMeans(sigma_sq_err_sh4), colMeans(sigma_sq_err_sh5),
        names  = sheet_names,
        ylab = "sigma_sq_err",
        col    = rainbow(5))

boxplot(colMeans(alpha_sh1), colMeans(alpha_sh2), colMeans(alpha_sh3), colMeans(alpha_sh4), 
        colMeans(alpha_sh5),
        names  = sheet_names,
        ylab = "alpha",
        ylim = c(0.4,1),
        col    = rainbow(5))

boxplot(colMeans(psi_delta_sh1), colMeans(psi_delta_sh2), colMeans(psi_delta_sh3), 
        colMeans(psi_delta_sh4), colMeans(psi_delta_sh5),
        names  = sheet_names,
        ylab = "psi_delta",
        col    = rainbow(5))

boxplot(colMeans(k_sh1), colMeans(k_sh2), colMeans(k_sh3), colMeans(k_sh4), colMeans(k_sh5),
        names  = sheet_names,
        ylab = "k",
        col    = rainbow(5))

load("E:/Phd_Paris Saclay/resultat/Results_M1_Simulation/result_m1_sh2_presentation.RData")
load("E:Phd_Paris Saclay/resultat/Results_M1_Simulation/y_obs_m1_sh2_presentation.RData")

alpha <- result_m1_sh2_presentation[[4]]
hist(alpha)
plot(density(alpha))

library(ggplot2)

# 1) flatten the matrix to a single vector
alpha_vec <- as.vector(alpha)

# 2) put into a data.frame
df_alpha <- data.frame(alpha = alpha_vec)

# 3) ggplot density
ggplot(df_alpha, aes(x = alpha)) +
  geom_density(
    fill  = "#800020",   # burgundy fill
    color = "#550010",   # darker border
    alpha = 0.6,         # semi‐transparent
    size  = 1
  ) +
  labs(
    x     = expression(alpha),
    y     = "Density",
    title = "Posterior Density of " ~ alpha
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    plot.margin  = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

#################################################################################
#################################################################################
###################### Results for model 2 ######################################
#################################################################################
#################################################################################
set.seed(12345)

k <- 0.2
sigma_sq_err   <- 0.01
sigma_sq_delta <- sigma_sq_err / k

n_samples <- 3
n_iter    <- 20000
burn_in   <- 5000

# (g, h0), sigma, psi, k, alpha
mcmc_parameters <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
Sigma_theta     <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)

psi_vec <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

results_list <- vector("list", length(psi_vec))
names(results_list) <- paste0("psi", 1:length(psi_vec))

for (i in seq_along(psi_vec)) {
  sim_psi_delta <- psi_vec[i]
  cat("Running psi_delta =", sim_psi_delta, "\n")
  
  sigma_props <- c(NA, NA, NA, NA, 0.4, NA)
  init        <- c(9.8, 46.46, 0.01, 0.5, 0.2, 0.2)
  # if (i == 1) {
  #   # Psi1
  #   sigma_props <- c(NA, NA, NA, NA, 0.02, NA)
  # } else if (i == 10) {
  #   # Psi10
  #   sigma_props <- c(NA, NA, NA, NA, 0.3, NA)
  # } else {
  #   # Psi2 ... Psi9
  #   sigma_props <- c(NA, NA, NA, NA, 0.4, NA)
  # }
  
  g_chain     <- matrix(NA, n_iter, n_samples)
  h0_chain    <- matrix(NA, n_iter, n_samples)
  sigma_chain <- matrix(NA, n_iter, n_samples)
  alpha_chain <- matrix(NA, n_iter, n_samples)
  psi_chain   <- matrix(NA, n_iter, n_samples)
  k_chain     <- matrix(NA, n_iter, n_samples)
  loglik_mat  <- matrix(NA, n_iter, n_samples)
  delta_list  <- vector("list", n_samples)
  zeta_list   <- vector("list", n_samples)
  accept_rate <- numeric(n_samples)
  
  y_obs <- matrix(NA, n, n_samples)
  
  for (v in 1:n_samples) {
    Sigma_delta <- GP_covariance(t, sigma_sq_delta, sim_psi_delta)
    delta <- as.vector(rmvnorm(1, rep(0, n), Sigma_delta))
    y_1 <- balldropg(t, c(9.8, 46.46)) + rnorm(n, 0, sqrt(sigma_sq_err)) + delta
    y_obs[, v] <- y_1
    
    res <- mcmc_step6(
      y_1, t, n_iter, init,
      sigma_props, mcmc_parameters,
      Sigma_theta, n_burnin = burn_in
    )
    
    g_chain[, v]     <- res$theta[, 1]
    h0_chain[, v]    <- res$theta[, 2]
    sigma_chain[, v] <- res$theta[, 3]
    alpha_chain[, v] <- res$theta[, 4]
    psi_chain[, v]   <- res$theta[, 5]
    k_chain[, v]     <- res$theta[, 6]
    
    delta_list[[v]] <- res$delta
    zeta_list[[v]]  <- res$zeta
    loglik_mat[, v] <- res$loglik
    accept_rate[v]  <- res$accept_rate_psi
  }
  
  results_list[[i]] <- list(
    g_chain     = g_chain,
    h0_chain    = h0_chain,
    sigma_chain = sigma_chain,
    alpha_chain = alpha_chain,
    psi_chain   = psi_chain,
    k_chain     = k_chain,
    delta_list  = delta_list,
    zeta_list   = zeta_list,
    loglik_mat  = loglik_mat,
    accept_rate = accept_rate,
    y_obs       = y_obs,
    psi_delta   = sim_psi_delta
  )
  
  # save(results_list[[i]], file = paste0("/Users/negarsoleimani/Documents/phd/paper1/", "result_m2_sh2_psi", i, "_simple1.RData"))
}

# Define a function to plot boxplots for each parameter
plot_boxplot <- function(results_all, param_index, ylab, abline_value, ylim_value = NULL) {
  results_list <- lapply(1:10, function(i) results_all[[i]][[param_index]])
  
  par(mfrow = c(1,1), mai = c(1,1,0.5,0.5))
  boxplot(lapply(results_list, colMeans), 
          names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
          xlab = expression(psi[delta]),
          ylab = ylab,
          ylim = ylim_value,
          cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
  abline(h = abline_value)
}

# Plot for g
plot_boxplot(results_all, 1, expression(g), 9.8)

# Plot for h0
plot_boxplot(results_all, 2, expression(h), 46.45, ylim_value = c(45, 50))

# Plot for sigma_sq_err
plot_boxplot(results_all, 3, expression(sigma[err]^2), 0.1)

# Plot for alpha
plot_boxplot(results_all, 4, expression(alpha), 0.4, ylim_value = c(0, 1))

# Plot for psi_delta
plot_boxplot(results_all, 5, expression(psi[delta]), 0.4)

# Plot for k
plot_boxplot(results_all, 6, expression(k), 0.2)

