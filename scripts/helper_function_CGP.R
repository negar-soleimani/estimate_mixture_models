<<<<<<< HEAD
# helper function for when we want to use the simple exponential kernel function in the covariance matrix related to the model discrepancy
=======
# helper function for when we want to use the simple exponential kernel 
#function in the covariance matrix related to the model discrepancy
>>>>>>> 8450abc4c47461f1601519d1d05674d6497dd7ed

GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  n <- length(t)
  Sigma <- outer(t, t, function(ti, tj) sigma_sq_delta * exp(-abs(ti - tj) / psi_delta))
  return(Sigma)
<<<<<<< HEAD
}
=======
}
>>>>>>> 8450abc4c47461f1601519d1d05674d6497dd7ed
