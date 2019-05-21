#' Basic Bootstrap CI
#'
#' Calculate a boostrap confidence interval using the
#' basic method
#' @export
basic_CI <- function(theta_hat, theta_hat_star, level=0.95){
  xbar <- mean(theta_hat_star)
  q <- quantile(theta_hat_star, c( (1-level)/2, 1 - (1-level)/2 ) )
  output <- c( theta_hat - (q[2] - xbar), theta_hat + (xbar - q[1]) )
  names(output) <- rev( names(output) )
  return(output)
}
