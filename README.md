# Bayesian Mixture Model for Computer Code Calibration
---

## Overview

This repository provides an R implementation of a **Bayesian mixture model** for calibrating computer codes in the presence of model discrepancy. The core idea is to let the data decide — via a mixture weight $\alpha$ — whether a discrepancy term $\delta(x)$ is needed to explain the gap between the physical code $f(x, \boldsymbol{\theta})$ and the field observations $y$.

The model is:
$$y_i = g(x_i)\boldsymbol{\theta} + \zeta_i\delta(x_i) + \varepsilon_i, \qquad \varepsilon_i \sim \mathcal{N}(0, \lambda^2)$$

where $\zeta_i \in \{0, 1\}$ is a latent allocation variable.

- When $\alpha \to 1$: the model collapses to $\mathcal{M}_0$ (no discrepancy).
- When $\alpha \to 0$: the model collapses to $\mathcal{M}_1$ (discrepancy).
- The pointwise inclusion probabilities $\hat{p}_i = P(\zeta_i = 1 \mid Y)$ provide a **local diagnostic**.

The discrepancy $\delta(X)$ is modelled with a **Gaussian process (GP)** prior.
The inference is carried out with a **Metropolis-within-Gibbs** MCMC sampler. All full conditionals are derived in closed form wherever possible (see the paper's supplementary material).

---

## Application

The model is applied to a **ball-drop experiment**: a ball is dropped from a height $h_0$ and its trajectory is recorded. The physical code is:

$$f(x_i, \boldsymbol{\theta}) = h_0 - \frac{1}{2}g_et_i^2$$

where $\boldsymbol{\theta} = (h_0, g_e)^\top$ are the calibration parameters (initial height and gravitational acceleration).

The dataset contains two drops for each of **11 different balls**:

| Ball | Sheet |
|---|---|
| Baseball | 1 |
| Blue Basketball | 2 |
| Bowling Ball | 3 |
| Golf Ball | 4 |
| Green Basketball | 5 |
| Orange Whiffle Ball | 6 |
| Tennis Ball | 7 |
| Volleyball | 8 |
| Whiffle Ball 1 | 9 |
| Whiffle Ball 2 | 10 |
| Yellow Whiffle Ball | 11 |

---
## Key Parameters

| Parameter | Symbol | Role |
|---|---|---|
| `g` | $g_e$ | Gravitational acceleration (calibration) |
| `h0` | $h_0$ | Initial height (calibration) |
| `sigma_sq_err` | $\lambda^2$ | Measurement noise variance |
| `alpha` | $\alpha$ | Mixture weight — probability of $\mathcal{M}_0$ |
| `psi_delta` | $\gamma_\delta$ | GP length-scale for $\delta(X)$ |
| `k` | $k$ | GP signal-to-noise ratio ($\sigma^2_\delta = \lambda^2 / k$) |

## Getting Started

### 1 — Prerequisites

Install the following R packages:

```r
install.packages(c(
  "mvtnorm",    
  "invgamma",   
  "truncnorm",  
  "readxl",     
  "ggplot2",    
  "patchwork",  
  "coda"      
))
```

### 2 — Running a simulation study

```r
# Simulation under M_0 with classical GP
source("examples/simulation_model0_CGP.R")

# Simulation under M_1 with classical GP
source("examples/simulation_model1_CGP.R")

# Simulation under M1 with orthogonal GP
source("examples/simulation_model1_OGP.R")

# Simulation under M_1 with classical GP (four scenarios)
# discrepancy is ~0 for first ndis points, then active afterwards
# Scenario (I):   seuil == FALSE, g == FALSE
# Scenario (II):  g == TRUE, seuil == FALSE 
# Scenario (III): g == FALSE, seuil == TRUE 
# Scenario (IV):  seuil == TRUE, g == TRUE
source("examples/simulation_model1_seuil_CGP.R")
```

### 3 — Running the real-data analysis

Edit `data/prepare_data.R` to select the ball (`sheet`) and the drop number, then:

```r
source("examples/real_data.R")
```

Or use the full pipeline with Gelman-Rubin convergence monitoring (see below).

---

### Initialisation vector

```r
# init = c(g, h0, sigma_sq_err, alpha, psi_delta, k)
init <- c(9.8, 46.45, 0.01, 0.5, 0.5, 0.1)
```

### Fixing parameters during sampling

Each parameter can be fixed at its initial value by setting the corresponding flag to `TRUE`:

```r
res <- mcmc_step6(
  y = y, t = t, n_iter = 20000, init = init,
  sigma_proposals = c(NA, NA, NA, NA, 0.5, NA),
  g_init      = TRUE,   # g fixed at 9.8
  h0_init     = FALSE,  # h0 estimated
  sig2er_init = FALSE,
  alpha_init  = FALSE,
  psi_init    = FALSE,
  k_init      = FALSE,
  Sigma_theta = matrix(c(0.5, 0, 0, 0.5), 2),
  n_burnin    = 2000,
  seuil       = TRUE,   # threshold mechanism active
  s           = 0.3     # threshold value
)
```

---

## Thresholding Mechanism

### Threshold-based allocation

To avoid activating the discrepancy model in regions where the estimated discrepancy is negligible, we introduce the threshold indicator

$$
\tau_i(s)=\mathbf{1}\{|\delta(x_i)|>s\},
$$

and define the allocation probabilities as

$$
P(\zeta_i=1\mid \delta,\alpha,s)=(1-\alpha)\tau_i(s),
$$

$$
P(\zeta_i=0\mid \delta,\alpha,s)=1-\tau_i(s)+\alpha\\tau_i(s).
$$

Hence,

- if $|\delta(x_i)|\le s$, then $\zeta_i=0$ deterministically (the no-discrepancy model is selected);
- if $|\delta(x_i)|>s$, the standard Bayesian mixture allocation is recovered.
---

## Diagnostics

### Gelman-Rubin convergence

Run multiple chains and check $\hat{R} < 1.1$:

```r
library(coda)
mcmc_list <- mcmc.list(
  mcmc(chain1[, c("g","h0","sigma_sq_err","alpha","psi_delta","k")]),
  mcmc(chain2[, c("g","h0","sigma_sq_err","alpha","psi_delta","k")])
)
gelman.diag(mcmc_list)
```
---

