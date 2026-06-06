balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range + t_min)^2) #cbind(1, -0.5 * (t * t_range)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}
