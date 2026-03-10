# helper function for when we want to use the "orthogonal exponential kernel function" in the covariance matrix related to the model discrepancy:
# Page 5 in the article

J_x <- function(x, psi_delta) {
  psi <- psi_delta
  psi * ( 2 - exp(-x/psi) - exp(-(1 - x)/psi) )
}

I_x <- function(x, psi_delta) {
  psi <- psi_delta
  A1 <- function(b) { psi * (1 - exp(-b/psi)) }
  A2 <- function(b) { psi^2 * (1 - exp(-b/psi)) - psi * b * exp(-b/psi) }
  A3 <- function(b) { 2*psi^3 - exp(-b/psi) * (b^2*psi + 2*b*psi^2 + 2*psi^3) }
  
  term1 <- x^2 * ( A1(x) + A1(1 - x) )
  term2 <- 2*x * ( -A2(x) + A2(1 - x) )
  term3 <- A3(x) + A3(1 - x)
  
  term1 + term2 + term3
}

A_scalar <- function(psi_delta) {
  psi <- psi_delta
  2*psi - ((2*psi^2) * (1 - exp(-1/psi)))
}

B_scalar <- function(psi_delta) {
  psi <- psi_delta
  E <- exp(-1/psi)
  ((2/3)*psi) - (psi^2) + (2*psi^3) - (4*psi^4) + E*(psi^2 + 2*psi^3 + 4*psi^4)
}

D_scalar <- function(psi_delta) {
  psi <- psi_delta
  E <- exp(-1/psi)
  ((2/5)*psi) - (psi^2) + ((4/3)*psi^3) - (8*(psi^6)) + E*(4*(psi^4) + 8*(psi^5) + 8*(psi^6))
}

h_vec_x <- function(x, psi_delta) {
  psi <- psi_delta
  a   <- t_range 
  J <- J_x(x, psi)
  I <- I_x(x, psi)
  rbind(J, -0.5 * (a^2) * I)
}

H_matrix <- function(psi_delta) {
  A <- A_scalar(psi_delta); B <- B_scalar(psi_delta); D <- D_scalar(psi_delta)
  H <- matrix(c(A, -((a^2)*B)/2, -((a^2)*B)/2, ((a^4)*D)/4), nrow = 2, byrow = TRUE)
  H <- 0.5 * (H + t(H))
  H + diag(1e-10, 2)
}


GP_correlation <- function(t, psi_delta) {
  R_se <- outer(t, t, function(ti, tj) exp(-abs(ti - tj) / psi_delta))
  
  h_mat <- t(sapply(t, function(x) h_vec_x(x, psi_delta)))
  H <- H_matrix(psi_delta)
  H_inv <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(H_inv)) return(NULL)
  
  K_star <- R_se - h_mat %*% H_inv %*% t(h_mat)
  K_star <- 0.5 * (K_star + t(K_star))  # symmetrize
  K_star <- K_star + diag(1e-10, nrow(K_star))
  K_star
}

GP_covariance_star_complete <- function(t, sigma_sq_err, k, psi_delta) {
  K_star <- GP_correlation(t, psi_delta)
  if (is.null(K_star)) {
    return(NULL)
  }
  (sigma_sq_err / k) * K_star
}

symm <- function(M) 0.5 * (M + t(M))

# Adaptive-jitter Cholesky:
# تلاش می‌کند chol(M + jitter I) را پاس کند. jitter را ضربی زیاد می‌کند.
# اگر باز هم fail شد، eigen-fallback می‌زند و SPD می‌کند.
chol_adapt <- function(M,
                       jitter0 = 1e-10,
                       mult    = 10,
                       max_tries = 12,
                       jitter_max = 1e-2) {
  M <- symm(M)
  n <- nrow(M)
  jitter <- jitter0
  
  for (i in 1:max_tries) {
    R <- tryCatch(chol(M + diag(jitter, n)),
                  error = function(e) NULL)
    if (!is.null(R)) {
      return(list(chol = R,
                  Mpd  = M + diag(jitter, n),
                  jitter = jitter,
                  tries = i,
                  fallback = FALSE))
    }
    jitter <- jitter * mult
    if (jitter > jitter_max) break
  }
  
  # Fallback: eigenvalue floor
  eig <- eigen(M, symmetric = TRUE)
  vals <- pmax(eig$values, jitter0)
  Mpd <- eig$vectors %*% diag(vals, length(vals)) %*% t(eig$vectors)
  Mpd <- symm(Mpd)
  R <- chol(Mpd)
  
  list(chol = R,
       Mpd  = Mpd,
       jitter = NA_real_,
       tries = max_tries,
       fallback = TRUE)
}

inv_from_chol <- function(cholR) chol2inv(cholR)