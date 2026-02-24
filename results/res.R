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
################# Results for model 0 - sheet2 ##################################
#################################################################################
#################################################################################
g_sh2 <- result_m0_sh2_classic_classic[[1]]
h0_sh2 <- result_m0_sh2_classic_classic[[2]]
sigma_sq_err_sh2 <- result_m0_sh2_classic_classic[[3]]
alpha_sh2 <- result_m0_sh2_classic_classic[[4]]
psi_delta_sh2 <- result_m0_sh2_classic_classic[[5]]
k_sh2 <- result_m0_sh2_classic_classic[[6]]
#zeta <- result_m0_sh2_classic_classic[[9]][[30]]

mean(colMeans(sigma_sq_err_sh2))
par(mfrow = c(1, 3),
    mar   = c(3, 4, 1, 1) 
)
boxplot(
  colMeans(g_sh2),
  ylab = "g",
  #xlab = "Bowling Ball",
  col  = "lightseagreen",
  main = ""
)
abline(h = 9.8, lty = 2)

boxplot(
  colMeans(h0_sh2),
  ylab = "h0",
  #xlab = "",
  col  = "lightseagreen",
  main = ""
)
abline(h = 46.45045, lty = 2)

boxplot(
  sigma_sq_err_sh2[, 50],
  ylab = expression(lambda^2),
  #ylim = c(0.01, 0.0452),
  #xlab = "",
  col  = "lightseagreen",
  main = ""
)
abline(h = 0.01, lty = 2)

boxplot(
  colMeans(psi_delta_sh2),
  ylab = expression(lambda^2),
  #ylim = c(0.01, 0.0452),
  #xlab = "",
  col  = "lightseagreen",
  main = ""
)
abline(h = 0.4, lty = 2)

alpha <- result_m0_sh2_classic_classic[[4]]

library(ggplot2)

# 1) flatten the matrix to a single vector
alpha_vec <- as.vector(alpha)

df_alpha <- data.frame(alpha = alpha_vec)

ggplot(df_alpha, aes(x = alpha)) +
  geom_density(
    fill  = "#4CCDC9",  
    color = "lightseagreen", 
    alpha = 0.6,        
    size  = 1
  ) +
  labs(
    x     = expression(alpha),
    y     = "Density",
    title = "Posterior Density of alpha"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    plot.margin  = unit(c(0.2,0.2,0.2,0.2), "cm")
  )
plot(sigma_sq_err_sh2[ , 50], type = "l", col = "#4CCDC9")
abline(h = 0.01)


g_mat <- result_m0_sh2_classic_classic[[1]]
h0_mat <- result_m0_sh2_classic_classic[[2]]

g_last <- g_mat[, 50]
h0_last <- h0_mat[, 50]
y_obs <- result_m0_sh2_classic_classic$y_obs
y_obs30 <- y_obs[,50]
y_true <- balldropg(t, c(9.8, 46.46))

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
legend(x=34, y=48, legend=c("Simulated data",
                            "True code", "Predictions"), lwd=rep(2,2), col=c("gold","blue", "orange2"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))

boxplot(y_pred1[, 1:5], col = "orange2",
        xlab = "Time", ylab = "Height")
lines(y_true, lwd = 2, col = "blue")
points(y_obs30, col = "gold", pch=20, cex = 3)
legend(x=4.2, y=47, legend=c("Simulated data", "Predictions",
                             "True code"), lwd=rep(2,2), col=c("gold", "orange2","blue"), 
       cex=1, pch=c(19,NA,NA), lty=c(0,1,1))


##-------------------------------------------------
## Empirical coverage of 50% and 95% credible CIs
##-------------------------------------------------

g_true       <- 9.8
h0_true      <- 46.45      
lambda2_true <- 0.01

# Number of simulated datasets
n_sims <- ncol(g_sh2)

# Coverage indicators (0/1) for each dataset and each CI level
cover_g_50       <- logical(n_sims)
cover_h0_50      <- logical(n_sims)
cover_lambda2_50 <- logical(n_sims)

cover_g_95       <- logical(n_sims)
cover_h0_95      <- logical(n_sims)
cover_lambda2_95 <- logical(n_sims)

ci_g_50_mat       <- matrix(NA, nrow = n_sims, ncol = 2)
ci_h0_50_mat      <- matrix(NA, nrow = n_sims, ncol = 2)
ci_lambda2_50_mat <- matrix(NA, nrow = n_sims, ncol = 2)

ci_g_95_mat       <- matrix(NA, nrow = n_sims, ncol = 2)
ci_h0_95_mat      <- matrix(NA, nrow = n_sims, ncol = 2)
ci_lambda2_95_mat <- matrix(NA, nrow = n_sims, ncol = 2)

for (v in 1:n_sims) {
  # --- 50% central credible interval: 25% and 75% quantiles ---
  ci_g_50       <- quantile(g_sh2[, v],            probs = c(0.25, 0.75))
  ci_h0_50      <- quantile(h0_sh2[, v],           probs = c(0.25, 0.75))
  ci_lambda2_50 <- quantile(sigma_sq_err_sh2[, v], probs = c(0.25, 0.75))
  
  ci_g_50_mat[v, ]       <- ci_g_50
  ci_h0_50_mat[v, ]      <- ci_h0_50
  ci_lambda2_50_mat[v, ] <- ci_lambda2_50
  
  cover_g_50[v]       <- (g_true       >= ci_g_50[1]       && g_true       <= ci_g_50[2])
  cover_h0_50[v]      <- (h0_true      >= ci_h0_50[1]      && h0_true      <= ci_h0_50[2])
  cover_lambda2_50[v] <- (lambda2_true >= ci_lambda2_50[1] && lambda2_true <= ci_lambda2_50[2])
  
  # --- 95% central credible interval: 2.5% and 97.5% quantiles ---
  ci_g_95       <- quantile(g_sh2[, v],            probs = c(0.025, 0.975))
  ci_h0_95      <- quantile(h0_sh2[, v],           probs = c(0.025, 0.975))
  ci_lambda2_95 <- quantile(sigma_sq_err_sh2[, v], probs = c(0.025, 0.975))
  
  ci_g_95_mat[v, ]       <- ci_g_95
  ci_h0_95_mat[v, ]      <- ci_h0_95
  ci_lambda2_95_mat[v, ] <- ci_lambda2_95
  
  cover_g_95[v]       <- (g_true       >= ci_g_95[1]       && g_true       <= ci_g_95[2])
  cover_h0_95[v]      <- (h0_true      >= ci_h0_95[1]      && h0_true      <= ci_h0_95[2])
  cover_lambda2_95[v] <- (lambda2_true >= ci_lambda2_95[1] && lambda2_true <= ci_lambda2_95[2])
}

# Empirical coverage across the 50 datasets
cat("50% CI coverage:\n")
cat("  g       :", mean(cover_g_50),       "\n")
cat("  h0      :", mean(cover_h0_50),      "\n")
cat("  lambda^2:", mean(cover_lambda2_50), "\n\n")

cat("95% CI coverage:\n")
cat("  g       :", mean(cover_g_95),       "\n")
cat("  h0      :", mean(cover_h0_95),      "\n")
cat("  lambda^2:", mean(cover_lambda2_95), "\n")

coverage_table <- data.frame(
  dataset        = 1:n_sims,
  cover_g_50     = cover_g_50,
  cover_h0_50    = cover_h0_50,
  cover_lambda2_50 = cover_lambda2_50,
  cover_g_95     = cover_g_95,
  cover_h0_95    = cover_h0_95,
  cover_lambda2_95 = cover_lambda2_95
)

print(coverage_table)

######################################################

get_ci <- function(samples, level = 0.95) {
  alpha <- 1 - level
  q <- quantile(samples, probs = c(alpha/2, 1 - alpha/2))
  c(lower = q[1], upper = q[2])
}

get_central_ci <- function(samples, probs = c(0.25, 0.75)) {
  q <- quantile(samples, probs = probs)
  c(lower = q[1], upper = q[2])
}

# آیا مقدار واقعی داخل بازه هست یا نه؟
is_covered <- function(true_value, ci_vec) {
  (true_value >= ci_vec[1]) && (true_value <= ci_vec[2])
}

# lower / upper
ci_g_50       <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_g_95       <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_h0_50      <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_h0_95      <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_lambda2_50 <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))
ci_lambda2_95 <- matrix(NA, nrow = n_sims, ncol = 2,
                        dimnames = list(NULL, c("lower","upper")))

cover_g_50       <- logical(n_sims)
cover_g_95       <- logical(n_sims)
cover_h0_50      <- logical(n_sims)
cover_h0_95      <- logical(n_sims)
cover_lambda2_50 <- logical(n_sims)
cover_lambda2_95 <- logical(n_sims)

for (s in 1:n_sims) {
  g_s       <- g_sh2[, s]
  h0_s      <- h0_sh2[, s]
  lambda2_s <- sigma_sq_err_sh2[, s]
  
  ## ---- CI 50
  ci_g_50[s, ]       <- get_central_ci(g_s,       probs = c(0.25, 0.75))
  ci_h0_50[s, ]      <- get_central_ci(h0_s,      probs = c(0.25, 0.75))
  ci_lambda2_50[s, ] <- get_central_ci(lambda2_s, probs = c(0.25, 0.75))
  
  ## ---- CI 95٪ 
  ci_g_95[s, ]       <- get_ci(g_s,       level = 0.95)
  ci_h0_95[s, ]      <- get_ci(h0_s,      level = 0.95)
  ci_lambda2_95[s, ] <- get_ci(lambda2_s, level = 0.95)
  
  ## ---- coverage (hit / miss) ----
  cover_g_50[s]       <- is_covered(g_true,       ci_g_50[s, ])
  cover_g_95[s]       <- is_covered(g_true,       ci_g_95[s, ])
  cover_h0_50[s]      <- is_covered(h0_true,      ci_h0_50[s, ])
  cover_h0_95[s]      <- is_covered(h0_true,      ci_h0_95[s, ])
  cover_lambda2_50[s] <- is_covered(lambda2_true, ci_lambda2_50[s, ])
  cover_lambda2_95[s] <- is_covered(lambda2_true, ci_lambda2_95[s, ])
}

len_g_50       <- ci_g_50[, "upper"]       - ci_g_50[, "lower"]
len_g_95       <- ci_g_95[, "upper"]       - ci_g_95[, "lower"]
len_h0_50      <- ci_h0_50[, "upper"]      - ci_h0_50[, "lower"]
len_h0_95      <- ci_h0_95[, "upper"]      - ci_h0_95[, "lower"]
len_lambda2_50 <- ci_lambda2_50[, "upper"] - ci_lambda2_50[, "lower"]
len_lambda2_95 <- ci_lambda2_95[, "upper"] - ci_lambda2_95[, "lower"]

g_mean       <- colMeans(g_sh2)
h0_mean      <- colMeans(h0_sh2)
lambda2_mean <- colMeans(sigma_sq_err_sh2)

bias_g       <- g_mean       - g_true
bias_h0      <- h0_mean      - h0_true
bias_lambda2 <- lambda2_mean - lambda2_true

summary_table <- data.frame(
  parameter = rep(c("g", "h0", "lambda^2"), each = 2),
  CI_level  = rep(c("50%", "95%"), times = 3),
  coverage  = c(
    mean(cover_g_50),
    mean(cover_g_95),
    mean(cover_h0_50),
    mean(cover_h0_95),
    mean(cover_lambda2_50),
    mean(cover_lambda2_95)
  ) * 100,  
  mean_CI_length = c(
    mean(len_g_50),
    mean(len_g_95),
    mean(len_h0_50),
    mean(len_h0_95),
    mean(len_lambda2_50),
    mean(len_lambda2_95)
  ),
  mean_abs_bias = c(
    mean(abs(bias_g)),
    mean(abs(bias_g)),
    mean(abs(bias_h0)),
    mean(abs(bias_h0)),
    mean(abs(bias_lambda2)),
    mean(abs(bias_lambda2))
  )
)

summary_table

library(ggplot2)

## -------- g: 50% credible intervals --------
df_ci_g <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_g_95[, "lower"],
  upper   = ci_g_95[, "upper"],
  covered = factor(cover_g_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)


# مرکز هر بازه (فقط برای رسم نقطه)
df_ci_g$center <- (df_ci_g$lower + df_ci_g$upper) / 2

par(mfrow = c(1,3))
p_g_95 <- ggplot(df_ci_g,
                 aes(y = dataset,
                     x = center,
                     xmin = lower,
                     xmax = upper,
                     color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = g_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = "g",
    y = "dataset index",
    color = "coverage",
    title = "95% credible intervals for g across 50 datasets"
  ) +
  theme_minimal()

p_g_95

## -------- h0: 50% credible intervals --------
df_ci_h0 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_h0_50[, "lower"],
  upper   = ci_h0_50[, "upper"],
  covered = factor(cover_h0_50,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_h0$center <- (df_ci_h0$lower + df_ci_h0$upper) / 2

p_h0_50 <- ggplot(df_ci_h0,
                  aes(y = dataset,
                      x = center,
                      xmin = lower,
                      xmax = upper,
                      color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = h0_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(h[0]),
    y = "dataset index",
    color = "coverage",
    title = "50% credible intervals for h0 across 50 datasets"
  ) +
  theme_minimal()

p_h0_50

## -------- h0: 95% credible intervals --------
df_ci_h0 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_h0_95[, "lower"],
  upper   = ci_h0_95[, "upper"],
  covered = factor(cover_h0_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_h0$center <- (df_ci_h0$lower + df_ci_h0$upper) / 2

p_h0_95 <- ggplot(df_ci_h0,
                  aes(y = dataset,
                      x = center,
                      xmin = lower,
                      xmax = upper,
                      color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = h0_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(h[0]),
    y = "dataset index",
    color = "coverage",
    title = "95% credible intervals for h0 across 50 datasets"
  ) +
  theme_minimal()

p_h0_95

## -------- lambda^2: 50% credible intervals --------
df_ci_lambda2 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_lambda2_50[, "lower"],
  upper   = ci_lambda2_50[, "upper"],
  covered = factor(cover_lambda2_50,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_lambda2$center <- (df_ci_lambda2$lower + df_ci_lambda2$upper) / 2

p_lambda2_50 <- ggplot(df_ci_lambda2,
                       aes(y = dataset,
                           x = center,
                           xmin = lower,
                           xmax = upper,
                           color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = lambda2_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(lambda^2),
    y = "dataset index",
    color = "coverage",
    title = expression("50% credible intervals for " * lambda^2 * " across 50 datasets")
  ) +
  theme_minimal()

p_lambda2_50

## -------- lambda^2: 95% credible intervals --------
df_ci_lambda2 <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_lambda2_95[, "lower"],
  upper   = ci_lambda2_95[, "upper"],
  covered = factor(cover_lambda2_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

df_ci_lambda2$center <- (df_ci_lambda2$lower + df_ci_lambda2$upper) / 2

p_lambda2_95 <- ggplot(df_ci_lambda2,
                       aes(y = dataset,
                           x = center,
                           xmin = lower,
                           xmax = upper,
                           color = covered)) +
  geom_pointrange() +
  geom_vline(xintercept = lambda2_true, linetype = "dashed") +
  scale_color_manual(values = c("miss" = "red", "hit" = "darkgreen")) +
  labs(
    x = expression(lambda^2),
    y = "dataset index",
    color = "coverage",
    title = expression("95% credible intervals for " * lambda^2 * " across 50 datasets")
  ) +
  theme_minimal()

p_lambda2_95
#########################################################
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

