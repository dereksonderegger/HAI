#' Zero Truncated Poisson Distribution
#'
#' A poisson distribution with the zero removed.
#' @examples
#' dtpois(2, lambda=4)
#' dtpois(0, lambda=4)
#' dtpois(1:100, lambda=4) %>% sum()
#' @export
dtpois <- function(x, lambda){
  truncdist::dtrunc(x, 'pois', a=0.5, b=Inf, lambda )
}
#' @rdname dtpois
#' @export
ptpois <- function(q, lambda){
  truncdist::qtrunc(q, 'pois', a=0.5, b=Inf, lambda )
}
#' @rdname dtpois
#' @export
qtpois <- function(p, lambda){
  truncdist::ptrunc(p, 'pois', a=0.5, b=Inf, lambda )
}
#' @rdname dtpois
#' @export
rtpois <- function(n, lambda){
  truncdist::rtrunc(n, 'pois', a=0.5, b=Inf, lambda )
}





#' Probability function for Hypergeometric Truncated Poisson distribution
#'
#' @param nhat The observed number.
#' @param lambda The rate parameter for the poisson distribution.
#' @param alpha The sampling percentage.
#' @param N The population size.
#' @param log Should the log probability be returned.
#' @examples
#' dHyperTPoisson(nhat=2, lambda=2, alpha=1, N=10)
#' dHyperTPoisson(1, lambda=c(2,3), alpha=1, N=10)
#' dHyperTPoisson(1:3, lambda=2, alpha=1, N=10)
#' dHyperTPoisson(1:3, lambda=2, alpha=1, N=10, log=TRUE)
#'
#' @export
dHyperTPoisson <- Vectorize(function( nhat, lambda, alpha, N, log=FALSE){
  n <- nhat:N
  out <- sum( truncdist::dtrunc(nhat, 'hyper', a=0, params=list(m=n, n=N-n, k = round(alpha*N))) * dtpois(n, lambda) )
  if(log){ out = log(out) }
  return(out)
})

#' rdname dHyperTPoisson
#' @export
rHyperTPoisson <- function(n, lambda, alpha, N){
  nn <- rtpois( n, lambda )
  out <-
    data.frame(nn=nn) %>%
    mutate( m = truncdist::rtrunc(1, 'hyper', a=0, b=Inf, params=list(m=.$nn, n=N-.$nn, k=round(N*alpha) ) ) ) %>%
    pull(m)

  return( out )
}


