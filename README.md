# Bayesian Mixture Model for Computer Code Calibration
---

## Overview

This repository provides an R implementation of a **Bayesian mixture model** for calibrating computer codes in the presence of model discrepancy. The core idea is to let the data decide — via a mixture weight $\alpha$ — whether a discrepancy term $\delta(x)$ is needed to explain the gap between the physical code $f(x, \boldsymbol{\theta})$ and the field observations $y$.

The model is:
$$y_i = g(x_i)\boldsymbol{\theta} + \zeta_i\,\delta(x_i) + \varepsilon_i, \qquad \varepsilon_i \sim \mathcal{N}(0, \lambda^2)$$

where $\zeta_i \in \{0, 1\}$ is a latent allocation variable.

- When $\alpha \to 1$: the model collapses to $\mathcal{M}_0$ (no discrepancy needed).
- When $\alpha \to 0$: the model collapses to $\mathcal{M}_1$ (discrepancy is present everywhere).
- The posterior of $\alpha$ quantifies the global evidence for discrepancy.
- The pointwise inclusion probabilities $\hat{p}_i = P(\zeta_i = 1 \mid Y)$ provide a **local diagnostic**.

The discrepancy $\delta(\cdot)$ is modelled with a **Gaussian process (GP)** prior. Two variants are supported:

| Prior | Description |
|---|---|
| **CGP** — Classical GP | Standard exponential kernel |
| **OGP** — Orthogonal GP | Kernel projected to be orthogonal to $g(x)\boldsymbol{\theta}$, reducing confounding |

The inference is carried out with a **Metropolis-within-Gibbs** MCMC sampler. All full conditionals are derived in closed form wherever possible (see the paper's supplementary material).

---

## Application

The model is applied to a **ball-drop experiment**: a ball is dropped from a height $h_0$ and its trajectory is recorded. The physical code is:

$$f(x_i, \boldsymbol{\theta}) = h_0 - \frac{1}{2}\,g_e\,t_i^2$$

where $\boldsymbol{\theta} = (h_0, g_e)^\top$ are the calibration parameters (initial height and gravitational acceleration). The discrepancy $\delta(x_i)$ captures systematic departures from this idealized free-fall model (e.g., air drag).

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

## Repository Structure

```
estimate_mixture_models/
│
├── data/
│   ├── Ball_drops_data.xlsx      # Real experimental data (11 balls, 2 drops each)
│   └── prepare_data.R            # Load and normalise data for one ball/drop
│
├── scripts/
│   ├── physics_model.R           # Physical model f(x, θ) = h0 - 0.5*g*t²
│   ├── helper_function_CGP.R     # GP covariance kernel (classical)
│   ├── helper_function_OGP.R     # GP covariance kernel (orthogonal)
│   ├── main_function_seuil_CGP.R # MCMC sampler — CGP with optional thresholding
│   └── main_function_OGP.R       # MCMC sampler — OGP
│
├── examples/
│   ├── simulation_model0_CGP.R   # Simulation study under M₀ (CGP)
│   ├── simulation_model0_OGP.R   # Simulation study under M₀ (OGP)
│   ├── simulation_model1_CGP.R   # Simulation study under M₁ (CGP)
│   ├── simulation_model1_CGP_seuil.R  # Simulation study — CGP with thresholding
│   ├── simulation_model1_OGP.R   # Simulation study under M₁ (OGP)
│   └── real_data.R               # Real-data analysis (single ball, single drop)
│
└── results/
    ├── real_data_seuil_gfix/     # Results — real data, g fixed, with thresholding
    │   └── figures/              # Posterior plots for each ball and drop
    └── real_data_joint_seuil_gfix/  # Results — joint analysis (two drops)
        └── figures/
```

---

## Getting Started

### 1 — Prerequisites

Install the following R packages:

```r
install.packages(c(
  "mvtnorm",    # multivariate normal sampling
  "invgamma",   # inverse-gamma sampling
  "truncnorm",  # truncated-normal proposals
  "readxl",     # Excel data import
  "ggplot2",    # plots
  "patchwork",  # composite figures
  "coda"        # MCMC diagnostics
))
```

### 2 — Running a simulation study

```r
# Simulation under M_0 with classical GP
source("examples/simulation_model0_CGP.R")

# Simulation under M_1 with classical GP
source("examples/simulation_model1_CGP_seuil.R")

# Simulation under M1 with orthogonal GP
source("examples/simulation_model1_OGP.R")
```

### 3 — Running the real-data analysis

Edit `data/prepare_data.R` to select the ball (`sheet`) and the drop number, then:

```r
source("examples/real_data.R")
```

Or use the full pipeline with Gelman-Rubin convergence monitoring (see below).

---

## Key Parameters

| Parameter | Symbol | Role |
|---|---|---|
| `g` | $g_e$ | Gravitational acceleration (calibration) |
| `h0` | $h_0$ | Initial height (calibration) |
| `sigma_sq_err` | $\lambda^2$ | Measurement noise variance |
| `alpha` | $\alpha$ | Mixture weight — probability of $\mathcal{M}_0$ |
| `psi_delta` | $\gamma_\delta$ | GP length-scale for $\delta(\cdot)$ |
| `k` | $k$ | GP signal-to-noise ratio ($\sigma^2_\delta = \lambda^2 / k$) |

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

An optional **threshold-based allocation** replaces the standard Bernoulli allocation by:

$$P(\zeta_i = 1 \mid \delta, \alpha, s) = (1 - \alpha)\,\mathbf{1}\{|\delta(x_i)| > s\}$$

When $|\delta(x_i)| \leq s$, observation $i$ is deterministically assigned to $\mathcal{M}_0$. This provides a **local diagnostic**: the pointwise inclusion probabilities $\hat{p}_i$ highlight where the discrepancy is practically significant.

Activate with `seuil = TRUE` and choose `s` on the order of magnitude of the measurement noise.

---

## Outputs

Each MCMC run returns a list with:

| Element | Dimensions | Description |
|---|---|---|
| `theta` | `n_iter × 6` | Posterior draws of $(g, h_0, \lambda^2, \alpha, \gamma_\delta, k)$ |
| `delta` | `n_iter × n` | Posterior draws of $\delta(x_1), \ldots, \delta(x_n)$ |
| `zeta` | `n_iter × n` | Posterior draws of $\zeta_1, \ldots, \zeta_n \in \{0, 1\}$ |
| `loglik` | `n_iter` | Log-likelihood at each iteration |
| `accept_rate_psi` | scalar | Metropolis acceptance rate for $\gamma_\delta$ |

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

### Key plots to inspect

- **Trace plots** — check mixing for all parameters
- **Posterior of $\alpha$** — concentrated near 1: $\mathcal{M}_0$ preferred; near 0: $\mathcal{M}_1$ preferred
- **$\delta(x_i)$ posterior mean** — identifies where discrepancy is present
- **Pointwise inclusion probabilities $\hat{p}_i$** — local allocation diagnostic

---

## Notes on Identifiability

When $g_e$ is estimated jointly with $\delta(X)$, a confounding issue arises: a shift in $g_e$ can be partially absorbed by the discrepancy function, and vice versa. Two strategies are implemented to mitigate this:

1. **Fix $g$ at its nominal value** (`g_init = TRUE`) — this is the default in the main real-data analysis.
2. **Use the OGP prior** (`helper_function_OGP.R`, `main_function_OGP.R`) — the orthogonality constraint $\langle \delta, g(x)\boldsymbol{\theta} \rangle = 0$ prevents the discrepancy from absorbing the linear component of the calibration parameters.

Note: experiments on real data show that both strategies yield similar posterior conclusions (see supplementary material, Section S.7).

---

