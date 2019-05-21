#' Logit functions
#' @export
logit <- Vectorize(function(p){
  if(p<0 | p > 1){
    stop(paste('logit: p must be in [0,1] - p is now: ', p))
  }
  log( p / (1-p) )
})

#' Inverse Logit functions
#' @export
ilogit <- Vectorize(function(x){
  exp(x) / (1 + exp(x))
})

