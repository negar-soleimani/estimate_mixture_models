## ------------------------------------ Model0 ------------------------------------ ##

# =========================================================
# Figure 1, (main = template2.tex), page 14
# Simulation under M_0 (50 datasets, n = 45, Blue Basketball)
# Priors: delta ~ GP(0, Sigma), gamma_delta ~ U(0.1,1), k ~ U(0,1)
# Shared params (theta, lambda^2): Jeffreys prior
# Boxplots of the posterior means of delta_i
# =========================================================
rm(list = ls())

#load("/Users/negar/Documents/phd/Result/Model1/Classic/result_m0_sh2_classic_classic.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/result_m0_sh2_classic_classic_100.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/result_m0_sh2_classic_classic_200.RData")

#result_m0_sh2_classic_classic <- result_m0_sh2_classic_classic_100

g_sh2 <- result_m0_sh2_classic_classic[[1]]
h0_sh2 <- result_m0_sh2_classic_classic[[2]]
sigma_sq_err_sh2 <- result_m0_sh2_classic_classic[[3]]
alpha_sh2 <- result_m0_sh2_classic_classic[[4]]
psi_delta_sh2 <- result_m0_sh2_classic_classic[[5]]
k_sh2 <- result_m0_sh2_classic_classic[[6]]

# ----- (Top Left) Boxplots of the 50 posterior means of theta and lambda2 ------ #

par(
  mfrow = c(1, 3),
  mar   = c(3.7, 4.7, 1, 1),  
  cex.axis = 1.5,          
  cex.lab  = 1.8,          
  lwd      = 1.3          
)

boxplot(
  colMeans(g_sh2),
  ylab = "g",
  col  = "lightseagreen",
  main = ""
)
abline(h = 9.8, lty = 2)

boxplot(
  colMeans(h0_sh2),
  ylab = "h0",
  col  = "lightseagreen",
  main = ""
)
abline(h = 46.45045, lty = 2)

boxplot(
  colMeans(sigma_sq_err_sh2),
  ylab = expression(lambda^2),
  col  = "lightseagreen",
  main = ""
)
abline(h = 0.01, lty = 2)

# ----- (Top Right) Pooled posterior densities of alpha ------ #

alpha <- result_m0_sh2_classic_classic[[4]]

library(ggplot2)

# flatten the matrix to a single vector
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

# --------------------- (Bottom) delta --------------------- #

delta_list <- result_m0_sh2_classic_classic[[8]]
par(mfrow = c(1,1))
p <- length(delta_list[[1]][1, ])  
n_samples = 50 
delta_means <- matrix(NA, nrow = n_samples, ncol = p)

for (v in 1:n_samples) {
  delta_means[v, ] <- colMeans(delta_list[[v]])
}

df_delta <- as.data.frame(delta_means)
colnames(df_delta) <- paste0("delta_", 1:p)

par(mar = c(6, 6, 2, 2), cex.axis = 1.2, cex.lab = 1.5)

boxplot(
  df_delta,
  outline = FALSE,
  col = "lightseagreen",
  xaxt = "n",   
  ylab = expression("Posterior mean of " * delta)
)

axis(
  1,
  at = 1:p,
  labels = parse(text = paste0("delta[", 1:p, "]")),
  las = 2,   
  cex.axis = 0.7
)

# =========================================================
# Figure S.1, (supplementarymaterial.tex), page 10
# Posterior predictive vs simulated data (classical GP)
# Left: full data comparison
# Right: zoom on first 5 observations
# =========================================================

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

# =========================================================
# Table 1, (main = template2.tex), page 16
# Empirical coverage, CI length, and absolute bias
# Parameters: g, h0, lambda^2
# Based on 50 simulated datasets
# =========================================================

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

#print(coverage_table)

get_ci <- function(samples, level = 0.95) {
  alpha <- 1 - level
  q <- quantile(samples, probs = c(alpha/2, 1 - alpha/2))
  c(lower = q[1], upper = q[2])
}

get_central_ci <- function(samples, probs = c(0.25, 0.75)) {
  q <- quantile(samples, probs = probs)
  c(lower = q[1], upper = q[2])
}


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


# =========================================================
# Figure S.2, (supplementarymaterial.tex), page 20
# 95% credible intervals across 50 datasets (classical GP, m0)
# Parameters: g, h0, lambda^2
# Dots: posterior medians | Segments: CI
# Dashed line: true value | Color: coverage
# =========================================================

library(ggplot2)

## -------- g: 95% credible intervals --------
df_ci_g <- data.frame(
  dataset = factor(1:n_sims),
  lower   = ci_g_95[, "lower"],
  upper   = ci_g_95[, "upper"],
  covered = factor(cover_g_95,
                   levels = c(FALSE, TRUE),
                   labels = c("miss", "hit"))
)

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

# -------------------------------- posterior summaries for alpha -------------------------------- #

alpha <- result_m0_sh2_classic_classic[[4]]
n_samples <- 50
alpha_post_mean <- colMeans(alpha)

# Direct posterior probabilities for each dataset v
prob_alpha_lt_01 <- colMeans(alpha < 0.2)   # P(alpha < 0.1 | Y_v)
prob_alpha_gt_09 <- colMeans(alpha > 0.8)   # P(alpha > 0.9 | Y_v)

alpha_summary <- data.frame(
  dataset_id        = 1:n_samples,
  alpha_post_mean   = alpha_post_mean,
  prob_alpha_lt_01  = prob_alpha_lt_01,
  prob_alpha_gt_09  = prob_alpha_gt_09
)

print(alpha_summary)

prob_alpha_gt_09

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

boxplot(
  prob_alpha_gt_09,
  ylab = expression(hat(P)(alpha > 0.9 ~ "|" ~ Y[v])),
  col  = "lightseagreen",
  main = expression("Support for " * M[0])
)

boxplot(
  prob_alpha_lt_01,
  ylab = expression(hat(P)(alpha < 0.1 ~ "|" ~ Y[v])),
  col  = "lightseagreen",
  main = expression("Support for " * M[1])
)

## ------------------------------ Model1 (sans seuil) ------------------------------------ ##

# =========================================================
# Figure 3, (main = template2.tex), page 17
# Simulation under M1 (classical GP discrepancy)
# gamma_delta* in {0.01, 0.1, ..., 0.9}, 50 datasets each
# Boxplots of posterior means (g, h0, lambda^2, alpha, gamma_delta, k)
# Horizontal lines = true values
# =========================================================

load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi1_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi2_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi3_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi4_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi5_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi6_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi7_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi8_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi9_simple.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/result_m2_sh2_psi10_simple.RData")

# boxplot(result_m2_sh2_psi8_simple[[7]][[30]])

par(
  mfrow = c(2, 3),
  mar   = c(3.7, 4.7, 1, 1),  
  cex.axis = 1.5,          
  cex.lab  = 1.8,          
  lwd      = 1.3          
)
# ----- g ------ #
g1 <- result_m2_sh2_psi1_simple[[1]]
g2 <- result_m2_sh2_psi2_simple[[1]]
g3 <- result_m2_sh2_psi3_simple[[1]]
g4 <- result_m2_sh2_psi4_simple[[1]]
g5 <- result_m2_sh2_psi5_simple[[1]]
g6 <- result_m2_sh2_psi6_simple[[1]]
g7 <- result_m2_sh2_psi7_simple[[1]]
g8 <- result_m2_sh2_psi8_simple[[1]]
g9 <- result_m2_sh2_psi9_simple[[1]]
g10 <- result_m2_sh2_psi10_simple[[1]]

boxplot(colMeans(g1), colMeans(g2), colMeans(g3), colMeans(g4), colMeans(g5), 
        colMeans(g6), colMeans(g7), colMeans(g8), colMeans(g9), colMeans(g10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = "g",
        col = "lightseagreen")
abline(h=9.8, col = "orange", lwd = 1.5)

# ----- h0 ----- #

h01 <- result_m2_sh2_psi1_simple[[2]]
h02 <- result_m2_sh2_psi2_simple[[2]]
h03 <- result_m2_sh2_psi3_simple[[2]]
h04 <- result_m2_sh2_psi4_simple[[2]]
h05 <- result_m2_sh2_psi5_simple[[2]]
h06 <- result_m2_sh2_psi6_simple[[2]]
h07 <- result_m2_sh2_psi7_simple[[2]]
h08 <- result_m2_sh2_psi8_simple[[2]]
h09 <- result_m2_sh2_psi9_simple[[2]]
h010 <- result_m2_sh2_psi10_simple[[2]]

boxplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05),
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = "h0",
        col = "lightseagreen")
abline(h = 46.45, col = "orange", lwd = 2)

# ----- lambda2 ------ #

sigma1 <- result_m2_sh2_psi1_simple[[3]]
sigma2 <- result_m2_sh2_psi2_simple[[3]]
sigma3 <- result_m2_sh2_psi3_simple[[3]]
sigma4 <- result_m2_sh2_psi4_simple[[3]]
sigma5 <- result_m2_sh2_psi5_simple[[3]]
sigma6 <- result_m2_sh2_psi6_simple[[3]]
sigma7 <- result_m2_sh2_psi7_simple[[3]]
sigma8 <- result_m2_sh2_psi8_simple[[3]]
sigma9 <- result_m2_sh2_psi9_simple[[3]]
sigma10 <- result_m2_sh2_psi10_simple[[3]]

boxplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(lambda^2),
        col = "lightseagreen")
abline(h = 0.01, col = "orange", lwd = 1)

# ----- alpha ------ #

alpha1 <- result_m2_sh2_psi1_simple[[4]]
alpha2 <- result_m2_sh2_psi2_simple[[4]]
alpha3 <- result_m2_sh2_psi3_simple[[4]]
alpha4 <- result_m2_sh2_psi4_simple[[4]]
alpha5 <- result_m2_sh2_psi5_simple[[4]]
alpha6 <- result_m2_sh2_psi6_simple[[4]]
alpha7 <- result_m2_sh2_psi7_simple[[4]]
alpha8 <- result_m2_sh2_psi8_simple[[4]]
alpha9 <- result_m2_sh2_psi9_simple[[4]]
alpha10 <- result_m2_sh2_psi10_simple[[4]]

boxplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8", "0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(alpha),
        ylim = c(0,1),
        col = "lightseagreen")

# ----- gamma_delta ------ #

psi1 <- result_m2_sh2_psi1_simple[[5]]
psi2 <- result_m2_sh2_psi2_simple[[5]]
psi3 <- result_m2_sh2_psi3_simple[[5]]
psi4 <- result_m2_sh2_psi4_simple[[5]]
psi5 <- result_m2_sh2_psi5_simple[[5]]
psi6 <- result_m2_sh2_psi6_simple[[5]]
psi7 <- result_m2_sh2_psi7_simple[[5]]
psi8 <- result_m2_sh2_psi8_simple[[5]]
psi9 <- result_m2_sh2_psi9_simple[[5]]
psi10 <- result_m2_sh2_psi10_simple[[5]]

boxplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(gamma[delta]),
        col = "lightseagreen")

# ----- k ------ #

k1 <- result_m2_sh2_psi1_simple[[6]]
k2 <- result_m2_sh2_psi2_simple[[6]]
k3 <- result_m2_sh2_psi3_simple[[6]]
k4 <- result_m2_sh2_psi4_simple[[6]]
k5 <- result_m2_sh2_psi5_simple[[6]]
k6 <- result_m2_sh2_psi6_simple[[6]]
k7 <- result_m2_sh2_psi7_simple[[6]]
k8 <- result_m2_sh2_psi8_simple[[6]]
k9 <- result_m2_sh2_psi9_simple[[6]]
k10 <- result_m2_sh2_psi10_simple[[6]]

boxplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), colMeans(k10),
        names = c("0.01", "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(k),
        col = "lightseagreen")
abline(h = 0.1, col = "orange", lwd = 2)

# =========================================================
# Figure S.5, (supplementarymaterial.tex), page 25
# Simulation under M1 (orthogonal GP discrepancy)
# gamma_delta* in {0.01, 0.1, ..., 0.9}, 50 datasets each
# Boxplots of posterior means (g, h0, lambda^2, alpha, gamma_delta, k)
# Horizontal lines = true values
# =========================================================

load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi1_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi2_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi3_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi4_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi5_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi6_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi7_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi8_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi9_ortho.RData")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/orthogonalgp/result_m2_sh2_psi10_ortho.RData")

par(
  mfrow = c(2, 3),
  mar   = c(3.7, 4.7, 1, 1),  
  cex.axis = 1.5,          
  cex.lab  = 1.8,          
  lwd      = 1.3          
)
# ----- g ------ #
g1 <- result_m2_sh2_psi1_ortho[[1]]
g2 <- result_m2_sh2_psi2_ortho[[1]]
g3 <- result_m2_sh2_psi3_ortho[[1]]
g4 <- result_m2_sh2_psi4_ortho[[1]]
g5 <- result_m2_sh2_psi5_ortho[[1]]
g6 <- result_m2_sh2_psi6_ortho[[1]]
g7 <- result_m2_sh2_psi7_ortho[[1]]
g8 <- result_m2_sh2_psi8_ortho[[1]]
g9 <- result_m2_sh2_psi9_ortho[[1]]
g10 <- result_m2_sh2_psi10_ortho[[1]]

boxplot(colMeans(g1), colMeans(g2), colMeans(g3), colMeans(g4), colMeans(g5), 
        colMeans(g6), colMeans(g7), colMeans(g8), colMeans(g9), colMeans(g10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = "g",
        col = "lightseagreen")
abline(h=9.8, col = "orange", lwd = 1.5)

# ----- h0 ----- #

h01 <- result_m2_sh2_psi1_ortho[[2]]
h02 <- result_m2_sh2_psi2_ortho[[2]]
h03 <- result_m2_sh2_psi3_ortho[[2]]
h04 <- result_m2_sh2_psi4_ortho[[2]]
h05 <- result_m2_sh2_psi5_ortho[[2]]
h06 <- result_m2_sh2_psi6_ortho[[2]]
h07 <- result_m2_sh2_psi7_ortho[[2]]
h08 <- result_m2_sh2_psi8_ortho[[2]]
h09 <- result_m2_sh2_psi9_ortho[[2]]
h010 <- result_m2_sh2_psi10_ortho[[2]]

boxplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05),
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = "h0",
        col = "lightseagreen")
abline(h = 46.45, col = "orange", lwd = 2)

# ----- lambda2 ------ #

sigma1 <- result_m2_sh2_psi1_ortho[[3]]
sigma2 <- result_m2_sh2_psi2_ortho[[3]]
sigma3 <- result_m2_sh2_psi3_ortho[[3]]
sigma4 <- result_m2_sh2_psi4_ortho[[3]]
sigma5 <- result_m2_sh2_psi5_ortho[[3]]
sigma6 <- result_m2_sh2_psi6_ortho[[3]]
sigma7 <- result_m2_sh2_psi7_ortho[[3]]
sigma8 <- result_m2_sh2_psi8_ortho[[3]]
sigma9 <- result_m2_sh2_psi9_ortho[[3]]
sigma10 <- result_m2_sh2_psi10_ortho[[3]]

boxplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(lambda^2),
        col = "lightseagreen")
abline(h = 0.01, col = "orange", lwd = 1)

# ----- alpha ------ #

alpha1 <- result_m2_sh2_psi1_ortho[[4]]
alpha2 <- result_m2_sh2_psi2_ortho[[4]]
alpha3 <- result_m2_sh2_psi3_ortho[[4]]
alpha4 <- result_m2_sh2_psi4_ortho[[4]]
alpha5 <- result_m2_sh2_psi5_ortho[[4]]
alpha6 <- result_m2_sh2_psi6_ortho[[4]]
alpha7 <- result_m2_sh2_psi7_ortho[[4]]
alpha8 <- result_m2_sh2_psi8_ortho[[4]]
alpha9 <- result_m2_sh2_psi9_ortho[[4]]
alpha10 <- result_m2_sh2_psi10_ortho[[4]]

boxplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8", "0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(alpha),
        ylim = c(0,1),
        col = "lightseagreen")

# ----- gamma_delta ------ #

psi1 <- result_m2_sh2_psi1_ortho[[5]]
psi2 <- result_m2_sh2_psi2_ortho[[5]]
psi3 <- result_m2_sh2_psi3_ortho[[5]]
psi4 <- result_m2_sh2_psi4_ortho[[5]]
psi5 <- result_m2_sh2_psi5_ortho[[5]]
psi6 <- result_m2_sh2_psi6_ortho[[5]]
psi7 <- result_m2_sh2_psi7_ortho[[5]]
psi8 <- result_m2_sh2_psi8_ortho[[5]]
psi9 <- result_m2_sh2_psi9_ortho[[5]]
psi10 <- result_m2_sh2_psi10_ortho[[5]]

boxplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(gamma[delta]),
        col = "lightseagreen")

# ----- k ------ #

k1 <- result_m2_sh2_psi1_ortho[[6]]
k2 <- result_m2_sh2_psi2_ortho[[6]]
k3 <- result_m2_sh2_psi3_ortho[[6]]
k4 <- result_m2_sh2_psi4_ortho[[6]]
k5 <- result_m2_sh2_psi5_ortho[[6]]
k6 <- result_m2_sh2_psi6_ortho[[6]]
k7 <- result_m2_sh2_psi7_ortho[[6]]
k8 <- result_m2_sh2_psi8_ortho[[6]]
k9 <- result_m2_sh2_psi9_ortho[[6]]
k10 <- result_m2_sh2_psi10_ortho[[6]]

boxplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), colMeans(k10),
        names = c("0.01", "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(k),
        col = "lightseagreen")
abline(h = 0.1, col = "orange", lwd = 2)



# =========================================================
# Figure 4, (template2.tex), page 20
# Simulation under M_alpha (classical GP discrepancy)
# scenario_II
# Bayesian inference without thresholding and g fixed
# gamma_delta* 0.5, 50 datasets
# Pooled posterior density of alpha
# Posterior means of delta(x_i)
# pointwise posterior inclusion probabilities $\hat p_i$
# =========================================================

# ---------------------- left (scenario_II) ---------------------- #
source("data/prepare_data.R")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_II.RData")

library(ggplot2)
library(gridExtra)
library(patchwork)

res_obj <- result_scenario_II_psi0.2

## pooled posterior density of alpha across all datasets

alpha_vec <- as.vector(res_obj$alpha_chain)

df_alpha_pool <- data.frame(alpha = alpha_vec)

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(fill = "#E8C1CF", color = "#8E244D",
               alpha = 0.6, linewidth = 1) +
  labs(x = expression(alpha), y = "Density",
       title = "Pooled posterior density of alpha (Scenario II)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_alpha_pool)

## zeta == 2 means model with discrepancy

v <- 30
zeta_mat <- res_obj$zeta_list[[v]]
delta_mat <- res_obj$delta_list[[v]]
n <- ncol(delta_mat)

p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(i = 1:n, x = t, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(x = "Index i",
       y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
       title = "Pointwise inclusion probabilities (Scenario II)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta)


## ------- delta(x_i) 
ndis <- 20
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
       title = expression("Posterior summary of " * delta(x[i]) ~ "(Scenario II)")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_delta)


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
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta0)


right_col2 <- p_alpha_pool / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

print(right_col2)

right_col2 <- p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1))


# -------------------- right (scenario_IV) --------------------- #
source("data/prepare_data.R")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_IV.RData")

library(ggplot2)
library(gridExtra)
library(patchwork)

res_obj <- result_scenario_IV

## pooled posterior density of alpha across all datasets

alpha_vec <- as.vector(res_obj$alpha_chain)

df_alpha_pool <- data.frame(alpha = alpha_vec)

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(fill = "#E8C1CF", color = "#8E244D",
               alpha = 0.6, linewidth = 1) +
  labs(x = expression(alpha), y = "Density",
       title = "Pooled posterior density of alpha (Scenario IV)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_alpha_pool)

## zeta == 2 means model with discrepancy

v <- 30
zeta_mat <- res_obj$zeta_list[[v]]
delta_mat <- res_obj$delta_list[[v]]
n <- ncol(delta_mat)

p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(i = 1:n, x = t, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(x = "Index i",
       y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
       title = "Pointwise inclusion probabilities (Scenario IV)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta)


## ------- delta(x_i) 
ndis <- 20
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
  theme(plot.title = element_text(hjust = 0.5)) #print(p_delta)


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
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta0)


right_col2 <- p_alpha_pool / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

print(right_col2)

right_col2 <- p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1))

# =========================================================
# Figure S.3, (supplementarymaterial.tex), page 22
# Simulation under M_alpha (classical GP discrepancy)
# scenario_I
# Bayesian inference without thresholding and g estimated
# gamma_delta* 0.5, 50 datasets
# Pooled posterior density of alpha
# Posterior means of delta(x_i)
# pointwise posterior inclusion probabilities $\hat p_i$
# =========================================================
source("data/prepare_data.R")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_I.RData")

library(ggplot2)
library(gridExtra)
library(patchwork)

res_obj <- result_scenario_I

## pooled posterior density of alpha across all datasets

alpha_vec <- as.vector(res_obj$alpha_chain)

df_alpha_pool <- data.frame(alpha = alpha_vec)

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(fill = "#E8C1CF", color = "#8E244D",
               alpha = 0.6, linewidth = 1) +
  labs(x = expression(alpha), y = "Density",
       title = "Pooled posterior density of alpha (Scenario I)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_alpha_pool)

## zeta == 2 means model with discrepancy

v <- 30
zeta_mat <- res_obj$zeta_list[[v]]
delta_mat <- res_obj$delta_list[[v]]
n <- ncol(delta_mat)

p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(i = 1:n, x = t, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(x = "Index i",
       y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
       title = "Pointwise inclusion probabilities (Scenario I)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta)


## ------- delta(x_i) 
ndis <- 20
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
       title = expression("Posterior summary of " * delta(x[i]) ~ "(Scenario I)")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_delta)


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
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta0)


right_col2 <- p_alpha_pool / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

print(right_col2)

right_col2 <- p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1))

# =========================================================
# Figure S.4, (supplementarymaterial.tex), page 23
# Simulation under M_alpha (classical GP discrepancy)
# scenario_III
# Bayesian inference with thresholding and g estimated
# gamma_delta* 0.5, 50 datasets
# Pooled posterior density of alpha
# Posterior means of delta(x_i)
# pointwise posterior inclusion probabilities $\hat p_i$
# =========================================================

source("data/prepare_data.R")
load("/Users/negarsoleimani/Documents/phd/paper1/github/model1/threshold/result_scenario_III.RData")

library(ggplot2)
library(gridExtra)
library(patchwork)

res_obj <- result_scenario_III

## pooled posterior density of alpha across all datasets

alpha_vec <- as.vector(res_obj$alpha_chain)

df_alpha_pool <- data.frame(alpha = alpha_vec)

p_alpha_pool <- ggplot(df_alpha_pool, aes(x = alpha)) +
  geom_density(fill = "#E8C1CF", color = "#8E244D",
               alpha = 0.6, linewidth = 1) +
  labs(x = expression(alpha), y = "Density",
       title = "Pooled posterior density of alpha (Scenario III)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_alpha_pool)

## zeta == 2 means model with discrepancy

v <- 30
zeta_mat <- res_obj$zeta_list[[v]]
delta_mat <- res_obj$delta_list[[v]]
n <- ncol(delta_mat)

p_hat <- colMeans(zeta_mat == 2)

df_p <- data.frame(i = 1:n, x = t, p_hat = p_hat)

p_zeta <- ggplot(df_p, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  ylim(0, 1) +
  labs(x = "Index i",
       y = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y)),
       title = "Pointwise inclusion probabilities (Scenario III)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta)


## ------- delta(x_i) 
ndis <- 20
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
       title = expression("Posterior summary of " * delta(x[i]) ~ "(Scenario III)")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) #print(p_delta)


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
  theme(plot.title = element_text(hjust = 0.5)) #print(p_zeta0)


right_col2 <- p_alpha_pool / p_delta / p_zeta +
  patchwork::plot_layout(heights = c(1, 1.25, 1))

print(right_col2)

right_col2 <- p_delta / p_zeta0 +
  patchwork::plot_layout(heights = c(1, 1))


#################################################################################
#################################################################################
###################### Results for model 2 ######################################
#################################################################################
#################################################################################


boxplot(result_m2_sh2_psi6_simple1[[7]][[10]])

boxplot(result_m2_sh2_psi8_simple[[7]][[30]])



par(mfrow = c(2,3))


library(vioplot)

vioplot(colMeans(g1), colMeans(g2), colMeans(g3), colMeans(g4), colMeans(g5),
        colMeans(g6), colMeans(g7), colMeans(g8), colMeans(g9), colMeans(g10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "",  
        ylab = "",
        col = "#E18E96")
title(xlab = expression(gamma[delta]), line = 1.5)
title(ylab = expression(g), line = 1.5)
abline(h = 9.8, col = "#800020")




h01 <- result_m2_sh2_psi1_simple[[2]]
h02 <- result_m2_sh2_psi2_simple[[2]]
h03 <- result_m2_sh2_psi3_simple[[2]]
h04 <- result_m2_sh2_psi4_simple[[2]]
h05 <- result_m2_sh2_psi5_simple[[2]]
h06 <- result_m2_sh2_psi6_simple[[2]]
h07 <- result_m2_sh2_psi7_simple[[2]]
h08 <- result_m2_sh2_psi8_simple[[2]]
h09 <- result_m2_sh2_psi9_simple[[2]]
h010 <- result_m2_sh2_psi10_simple[[2]]

boxplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05),
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = "h0",
        col = "lightseagreen")
abline(h = 46.45, col = "orange", lwd = 2)

vioplot(colMeans(h01), colMeans(h02), colMeans(h03), colMeans(h04), colMeans(h05), 
        colMeans(h06), colMeans(h07), colMeans(h08), colMeans(h09), colMeans(h010),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(gamma[delta]), line = 1.5)
title(ylab = expression(h), line = 1.5)
abline(h = 46.45, col = "#800020")

sigma1 <- result_m2_sh2_psi1_simple[[3]]
sigma2 <- result_m2_sh2_psi2_simple[[3]]
sigma3 <- result_m2_sh2_psi3_simple[[3]]
sigma4 <- result_m2_sh2_psi4_simple[[3]]
sigma5 <- result_m2_sh2_psi5_simple[[3]]
sigma6 <- result_m2_sh2_psi6_simple[[3]]
sigma7 <- result_m2_sh2_psi7_simple[[3]]
sigma8 <- result_m2_sh2_psi8_simple[[3]]
sigma9 <- result_m2_sh2_psi9_simple[[3]]
sigma10 <- result_m2_sh2_psi10_simple[[3]]

boxplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(lambda^2),
        col = "lightseagreen")
abline(h = 0.01, col = "orange", lwd = 1)

vioplot(colMeans(sigma1), colMeans(sigma2), colMeans(sigma3), colMeans(sigma4), colMeans(sigma5), 
        colMeans(sigma6), colMeans(sigma7), colMeans(sigma8), colMeans(sigma9), colMeans(sigma10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(gamma[delta]), line = 1.5)
title(ylab = expression(lambda^2), line = 1.5)
abline(h = 0.01, col = "#800020")

alpha1 <- result_m2_sh2_psi1_simple[[4]]
alpha2 <- result_m2_sh2_psi2_simple[[4]]
alpha3 <- result_m2_sh2_psi3_simple[[4]]
alpha4 <- result_m2_sh2_psi4_simple[[4]]
alpha5 <- result_m2_sh2_psi5_simple[[4]]
alpha6 <- result_m2_sh2_psi6_simple[[4]]
alpha7 <- result_m2_sh2_psi7_simple[[4]]
alpha8 <- result_m2_sh2_psi8_simple[[4]]
alpha9 <- result_m2_sh2_psi9_simple[[4]]
alpha10 <- result_m2_sh2_psi10_simple[[4]]

boxplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8", "0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(alpha),
        ylim = c(0,1),
        col = "lightseagreen")

vioplot(colMeans(alpha1), colMeans(alpha2), colMeans(alpha3), colMeans(alpha4), colMeans(alpha5), 
        colMeans(alpha6), colMeans(alpha7), colMeans(alpha8), colMeans(alpha9), colMeans(alpha10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8", "0.9"),
        xlab = "",
        ylab = "",
        ylim = c(0,1),
        col = "#E18E96")
title(xlab = expression(gamma[delta]), line = 1.5)
title(ylab = expression(alpha), line = 1.5)

psi1 <- result_m2_sh2_psi1_simple[[5]]
psi2 <- result_m2_sh2_psi2_simple[[5]]
psi3 <- result_m2_sh2_psi3_simple[[5]]
psi4 <- result_m2_sh2_psi4_simple[[5]]
psi5 <- result_m2_sh2_psi5_simple[[5]]
psi6 <- result_m2_sh2_psi6_simple[[5]]
psi7 <- result_m2_sh2_psi7_simple[[5]]
psi8 <- result_m2_sh2_psi8_simple[[5]]
psi9 <- result_m2_sh2_psi9_simple[[5]]
psi10 <- result_m2_sh2_psi10_simple[[5]]

boxplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01", "0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(gamma[delta]),
        col = "lightseagreen")

vioplot(colMeans(psi1), colMeans(psi2), colMeans(psi3), colMeans(psi4), colMeans(psi5), 
        colMeans(psi6), colMeans(psi7), colMeans(psi8), colMeans(psi9), colMeans(psi10),
        names = c("0.01","0.1","0.2","0.3","0.4", "0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(gamma[delta]), line = 1.5)
title(ylab = expression(gamma[delta]), line = 1.5)

k1 <- result_m2_sh2_psi1_simple[[6]]
k2 <- result_m2_sh2_psi2_simple[[6]]
k3 <- result_m2_sh2_psi3_simple[[6]]
k4 <- result_m2_sh2_psi4_simple[[6]]
k5 <- result_m2_sh2_psi5_simple[[6]]
k6 <- result_m2_sh2_psi6_simple[[6]]
k7 <- result_m2_sh2_psi7_simple[[6]]
k8 <- result_m2_sh2_psi8_simple[[6]]
k9 <- result_m2_sh2_psi9_simple[[6]]
k10 <- result_m2_sh2_psi10_simple[[6]]

boxplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), colMeans(k10),
        names = c("0.01", "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = expression(gamma[delta]),
        ylab = expression(k),
        col = "lightseagreen")
abline(h = 0.1, col = "orange", lwd = 2)

vioplot(colMeans(k1), colMeans(k2), colMeans(k3), colMeans(k4), colMeans(k5), 
        colMeans(k6), colMeans(k7), colMeans(k8), colMeans(k9), colMeans(k10),
        names = c("0.01","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        xlab = "",
        ylab = "",
        col = "#E18E96")
title(xlab = expression(gamma[delta]), line = 1.5)
title(ylab = expression(k), line = 1.5)
abline(h = 0.2, col = "#800020")

res_psi8 <- result_m2_sh2_psi8_simple
g_psi8 <- result_m2_sh2_psi8_simple[[1]][,10]
h0_psi8 <- result_m2_sh2_psi8_simple[[2]][,10]
sigma_psi8 <- result_m2_sh2_psi8_simple[[3]][,10]
alpha_psi8 <- result_m2_sh2_psi8_simple[[4]][,10]
psi_psi8 <- result_m2_sh2_psi8_simple[[5]][,10]
k_psi8 <- result_m2_sh2_psi8_simple[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi8, type = "l", ylab = expression(g ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi8, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 46.46, col = "red")

plot(sigma_psi8, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.01, col = "red")

plot(alpha_psi8, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.7), xlab = "iteration")

plot(psi_psi8, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.7, col = "red")

plot(k_psi8, type = "l", ylab = expression(k ~ " with " ~ psi == 0.7), xlab = "iteration")
abline(h = 0.2, col = "red")


res_psi9 <- result_m2_sh2_psi9_simple
g_psi9 <- result_m2_sh2_psi9_simple[[1]][,10]
h0_psi9 <- result_m2_sh2_psi9_simple[[2]][,10]
sigma_psi9 <- result_m2_sh2_psi9_simple[[3]][,10]
alpha_psi9 <- result_m2_sh2_psi9_simple[[4]][,10]
psi_psi9 <- result_m2_sh2_psi9_simple[[5]][,10]
k_psi9 <- result_m2_sh2_psi9_simple[[6]][,10]

par(mfrow = c(2,3))
plot(g_psi9, type = "l", ylab = expression(g ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 9.8, col = "red")

plot(h0_psi9, type = "l", ylab = expression(h0 ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 46.46, col = "red")

plot(sigma_psi9, type = "l", ylab = expression(sigma ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.1, col = "red")

plot(alpha_psi9, type = "l", ylab = expression(alpha ~ " with " ~ psi == 0.8), xlab = "iteration")

plot(psi_psi9, type = "l", ylab = expression(psi_delta ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.8, col = "red")

plot(k_psi9, type = "l", ylab = expression(k ~ " with " ~ psi == 0.8), xlab = "iteration")
abline(h = 0.2, col = "red")




#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi1_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi2_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi3_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi4_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi5_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi6_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi7_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi8_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi9_simple.RData")
#load("/Users/negar/Documents/phd/Result/Model1/Classic/Model1/result_m2_sh2_psi10_simple.RData")

alpha_posterior_summary <- function(alpha_mat, low_cut = 0.1, high_cut = 0.9) {
  data.frame(
    dataset_id       = 1:ncol(alpha_mat),
    alpha_post_mean  = colMeans(alpha_mat),
    prob_alpha_lt_01 = colMeans(alpha_mat < low_cut),   # P(alpha < 0.1 | Y_v)
    prob_alpha_gt_09 = colMeans(alpha_mat > high_cut)   # P(alpha > 0.9 | Y_v)
  )
}

alpha1 <- result_m2_sh2_psi1_simple[[4]]
alpha2 <- result_m2_sh2_psi2_simple[[4]]
alpha3 <- result_m2_sh2_psi3_simple[[4]]
alpha4 <- result_m2_sh2_psi4_simple[[4]]
alpha5 <- result_m2_sh2_psi5_simple[[4]]
alpha6 <- result_m2_sh2_psi6_simple[[4]]
alpha7 <- result_m2_sh2_psi7_simple[[4]]
alpha8 <- result_m2_sh2_psi8_simple[[4]]
alpha9 <- result_m2_sh2_psi9_simple[[4]]
alpha10 <- result_m2_sh2_psi10_simple[[4]]

alpha_sum1  <- alpha_posterior_summary(alpha1)
alpha_sum2  <- alpha_posterior_summary(alpha2)
alpha_sum3  <- alpha_posterior_summary(alpha3)
alpha_sum4  <- alpha_posterior_summary(alpha4)
alpha_sum5  <- alpha_posterior_summary(alpha5)
alpha_sum6  <- alpha_posterior_summary(alpha6)
alpha_sum7  <- alpha_posterior_summary(alpha7)
alpha_sum8  <- alpha_posterior_summary(alpha8)
alpha_sum9  <- alpha_posterior_summary(alpha9)
alpha_sum10 <- alpha_posterior_summary(alpha10)

prob_alpha_lt_01_psi1  <- alpha_sum1$prob_alpha_lt_01
prob_alpha_lt_01_psi2  <- alpha_sum2$prob_alpha_lt_01
prob_alpha_lt_01_psi3  <- alpha_sum3$prob_alpha_lt_01
prob_alpha_lt_01_psi4  <- alpha_sum4$prob_alpha_lt_01
prob_alpha_lt_01_psi5  <- alpha_sum5$prob_alpha_lt_01
prob_alpha_lt_01_psi6  <- alpha_sum6$prob_alpha_lt_01
prob_alpha_lt_01_psi7  <- alpha_sum7$prob_alpha_lt_01
prob_alpha_lt_01_psi8  <- alpha_sum8$prob_alpha_lt_01
prob_alpha_lt_01_psi9  <- alpha_sum9$prob_alpha_lt_01
prob_alpha_lt_01_psi10 <- alpha_sum10$prob_alpha_lt_01

prob_alpha_gt_09_psi1  <- alpha_sum1$prob_alpha_gt_09
prob_alpha_gt_09_psi2  <- alpha_sum2$prob_alpha_gt_09
prob_alpha_gt_09_psi3  <- alpha_sum3$prob_alpha_gt_09
prob_alpha_gt_09_psi4  <- alpha_sum4$prob_alpha_gt_09
prob_alpha_gt_09_psi5  <- alpha_sum5$prob_alpha_gt_09
prob_alpha_gt_09_psi6  <- alpha_sum6$prob_alpha_gt_09
prob_alpha_gt_09_psi7  <- alpha_sum7$prob_alpha_gt_09
prob_alpha_gt_09_psi8  <- alpha_sum8$prob_alpha_gt_09
prob_alpha_gt_09_psi9  <- alpha_sum9$prob_alpha_gt_09
prob_alpha_gt_09_psi10 <- alpha_sum10$prob_alpha_gt_09

boxplot(
  prob_alpha_lt_01_psi1, prob_alpha_lt_01_psi2, prob_alpha_lt_01_psi3,
  prob_alpha_lt_01_psi4, prob_alpha_lt_01_psi5, prob_alpha_lt_01_psi6,
  prob_alpha_lt_01_psi7, prob_alpha_lt_01_psi8, prob_alpha_lt_01_psi9,
  prob_alpha_lt_01_psi10,
  names = c("0.01", "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
  xlab  = expression(gamma[delta]),
  ylab  = expression(hat(P)(alpha < 0.1 ~ "|" ~ Y[v])),
  ylim  = c(0,1),
  col   = "lightseagreen"
)


boxplot(
  prob_alpha_gt_09_psi1, prob_alpha_gt_09_psi2, prob_alpha_gt_09_psi3,
  prob_alpha_gt_09_psi4, prob_alpha_gt_09_psi5, prob_alpha_gt_09_psi6,
  prob_alpha_gt_09_psi7, prob_alpha_gt_09_psi8, prob_alpha_gt_09_psi9,
  prob_alpha_gt_09_psi10,
  names = c("0.01", "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
  xlab  = expression(gamma[delta]),
  ylab  = expression(hat(P)(alpha > 0.9 ~ "|" ~ Y[v])),
  ylim  = c(0,1),
  col   = "lightseagreen"
)

alpha_summary_all <- list(
  psi_0.01 = alpha_sum1,
  psi_0.1  = alpha_sum2,
  psi_0.2  = alpha_sum3,
  psi_0.3  = alpha_sum4,
  psi_0.4  = alpha_sum5,
  psi_0.5  = alpha_sum6,
  psi_0.6  = alpha_sum7,
  psi_0.7  = alpha_sum8,
  psi_0.8  = alpha_sum9,
  psi_0.9  = alpha_sum10
)


sapply(alpha_summary_all, function(x) mean(x$prob_alpha_lt_01))
sapply(alpha_summary_all, function(x) median(x$prob_alpha_lt_01))

################################################################################
################################################################################
################################################################################
###########################  seuil  ############################################
################################################################################
################################################################################
# posterior density of alpha (aggregated draws over all datasets)

alpha_vec <- as.vector(alpha_chain)
df_alpha <- data.frame(alpha = alpha_vec)

p_alpha <- ggplot(df_alpha, aes(x = alpha)) +
  geom_density(fill = "#4CCDC9", color = "lightseagreen", alpha = 0.6, linewidth = 1) +
  labs(
    x = expression(alpha),
    y = "Density",
    title = "Posterior density of α (Scenario II: thresholded allocations)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))


v <- 1
delta_chain_v <- delta_list[[v]]     # (iter x n)
zeta_chain_v  <- zeta_list[[v]]      # (iter x n)

# pointwise inclusion probabilities
p_hat <- colMeans(zeta_chain_v == 2)

# posterior summary for delta
delta_mean <- colMeans(delta_chain_v)
delta_ci   <- t(apply(delta_chain_v, 2, quantile, probs = c(0.025, 0.975)))

df_loc <- data.frame(
  i = 1:n,
  p_hat = p_hat,
  d_mean = delta_mean,
  d_low  = delta_ci[, 1],
  d_high = delta_ci[, 2]
)

library(ggplot2)
library(gridExtra)

p_delta <- ggplot(df_loc, aes(x = i, y = d_mean)) +
  geom_ribbon(aes(ymin = d_low, ymax = d_high), fill = "lightseagreen", alpha = 0.25) +
  geom_line(color = "lightseagreen", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed") +
  labs(
    x = "Index i",
    y = expression(delta(x[i])),
    title = expression("Posterior summary of " * delta(x[i]))
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

p_zeta <- ggplot(df_loc, aes(x = i, y = p_hat)) +
  geom_col(fill = "lightseagreen") +
  geom_vline(xintercept = ndis + 0.5, linetype = "dashed") +
  labs(
    x = "Index i",
    y = expression(hat(p)[i]),
    title = expression(hat(p)[i] == P(zeta[i] == 1 ~ "|" ~ Y))
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# Bottom block = delta + pointwise inclusion probabilities
bottom_block <- arrangeGrob(
  p_delta, p_zeta,
  ncol = 1,
  heights = c(1, 1)
)

grid.arrange(
  p_alpha,
  bottom_block,
  ncol = 1,
  heights = c(0.9, 2.1)
)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

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

