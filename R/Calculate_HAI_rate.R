#' Calculate HAI Rate
#' Given a data frame with column for cluster ID.
#' @param clusterID A vector of clusterIDs, one for each patient.
#' @param method A character string denoting which method to use.
#' @param alpha The sampling fraction.
#' @param N The population size from which the sample was taken.
#' @param bias Should we do a bias correction?
#' @param bias_reps The number of subsamples we should create for
#'   estimating the bias.
#' @param CI Should we calculate a bootstrap confidence interval
#'   in addition to the CI from the bias replicates
#' @param CI_reps The number of boostrap samples to create for
#'   calculating the bootstrap CI.
#' @param level The confidence interval level.
#'
#' @return A vector with the estimated Hospital Acquired Infection Rate
#'   and associatied confidence intervals.
#' @examples
#' pop <- generate_population(1000, lambda = .6, mix.frac = .01)
#'
#' pop %>% pull(clusterID) %>%
#'   calc_HAI( method='naive', bias=FALSE, CI=FALSE )
#'
#' clusterID <- pop %>% sample_frac(.5) %>% pull(clusterID)
#' calc_HAI(clusterID, method='naive', alpha=.5, N=nrow(pop), bias=TRUE, CI=TRUE, CI_reps=100)
#'
#'
#'
#' calc_HAI(clusterID, method='hypergeometric', alpha=.5, N=nrow(pop), bias=FALSE)
#' calc_HAI(clusterID, method='hypergeometric', alpha=.5, N=nrow(pop), bias=TRUE)
#'
#' pop %>% sample_frac(.5) %>% pull(clusterID) %>%
#'   calc_HAI( method='hypergeometric', alpha=.5, N=nrow(pop), bias=TRUE)
#'
#' pop %>% group_by(Aquired) %>% count() %>% group_by() %>%        # doing it straight
#'     mutate(percent = n / sum(n))                                # using metadata
#' pop %>% group_by(clusterID) %>% count() %>% pull(n) %>%         # via calc_HAI
#'   calc_HAI( . )
#'
#' pop %>% sample_frac( .5 ) %>% group_by(clusterID) %>% count() %>% pull(n) %>%
#'   calc_HAI( ., alpha=.5, N=nrow(pop), method = 'naive' )
#' pop %>% sample_frac( .5 ) %>% group_by(clusterID) %>% count() %>% pull(n) %>%
#'   calc_HAI( ., alpha=.5, N=nrow(pop), method = 'HypergeometricPoisson')
#' pop %>% sample_frac( .5 ) %>% group_by(clusterID) %>% count() %>% pull(n) %>%
#'   calc_HAI( ., alpha=.5, N=nrow(pop), method = 'chao')
#'
#'
#' out <- NULL
#' for( lambda in c( 0.5, 1, 1.5, 2) ){
#'   for( alpha in c( .25, .5, .75, 1) ){
#'     for( i in 1:4 ){
#'       N <- 1000
#'       pop <- generate_population(N, lambda = lambda, mix.frac = 0)
#'       HAI <- pop %>% group_by(clusterID) %>% count() %>% pull(n) %>% calc_HAI( . )
#'       for(method in c('naive','hypergeometric','chao','jack1', 'jack2', 'boot') ){
#'         Est.HAI <-pop %>% sample_frac( alpha ) %>% group_by(clusterID) %>%
#'                    count() %>% pull(n) %>%
#'                    calc_HAI( alpha=alpha, N=nrow(pop), method = method)
#'
#'         out <- rbind(out, data.frame(
#'              alpha = alpha, lambda = lambda, method = method, rep = i,
#'              HAI = HAI, Est.HAI=Est.HAI))
#'       }
#'     }
#'   }
#' }
#' ggplot( out, aes( x=HAI, y=Est.HAI, color=method )) +
#'   geom_point(alpha=.5) +
#'   facet_grid(method ~ alpha) +
#'   geom_abline(intercept = 0, slope=1, color='red')
#'
#' out %>% group_by(alpha) %>% summarize(mean=mean(HAI))
#'
#' @export
calc_HAI <- function(clusterID=NULL, method='naive', alpha=1, N=length(clusterID),
                     bias=TRUE, bias_reps=10, CI=FALSE, CI_reps=1000, level=0.95 ){

  HAI <- NULL
  method = str_to_lower(method)
  methods <- c('naive','plugin','hypergeometric','binomial1','binomial2','chao','jack1','jack2','boot')

  m <- data.frame(clusterID =clusterID ) %>% group_by(clusterID) %>% count() %>% pull(n)
  n <- length(clusterID)

  if( pmatch( method, methods ) == 1 ){  # naive
    nhat = m #/alpha
    HAI = 1 - ( length(nhat) / sum(nhat) )
  }else if( pmatch( method, methods ) == 2 ){  # plugin
    nhat = m / alpha
    HAI = 1 - ( length(nhat) / sum(nhat) )
  }else if( pmatch( method, methods ) == 3){ # Hypergeometric
    if(is.null(alpha) | is.null(N)){
      stop('alpha and N must be defined')
    }
    # function to optimize to find the MOM estimator nhat_i
    f <- Vectorize(function(ni, mi, alpha, N){
      E_mi <- alpha * ni * ( 1 - dhyper(0, ni, N-ni, round(alpha*N)) )^(-1)
      return(abs(mi - E_mi) )
    }, c('ni','mi'))

    ni <- 1:(max(m)/alpha)
    nhats <- NULL
    for( mi in unique(m) ){
      nhats <- data.frame(mi=mi, nhat_i=ni, diff = f(ni, mi, alpha, N)) %>% rbind(nhats,.)
    }
    nhats <- nhats %>% group_by(mi) %>% arrange(diff) %>% slice(1)


    HAI <-
      data.frame( clusterID = 1:length(m)) %>% mutate( mi = m ) %>%
      left_join(nhats, by='mi')%>%
      summarize( HAI = 1 - (length(nhat_i) / sum(nhat_i)) ) %>%
      pull(HAI)
  }else if( pmatch( method, methods ) == 4){ # binomial1
    nhat = m * ( 1-(1-1/N)^(alpha*N) ) / alpha
    nhat = pmax(m, nhat)
    HAI = 1 - ( length(nhat) / sum(nhat) )
  }else if( pmatch( method, methods ) == 5){ # binomial2
    nhat = m * ( 1-(1-(m/alpha)/N)^(alpha*N) ) / alpha
    nhat = pmax(m, nhat)
    HAI = 1 - ( length(nhat) / sum(nhat) )

  }else if(pmatch( method, methods )  %in% 6:9 ){
    if(is.null(alpha) | is.null(N)){
      stop('alpha and N must be defined')
    }
    comm <- data.frame(clusterID=1:length(m), m=m) %>%
      group_by(clusterID) %>% do({ merge(., data.frame(rep=1:.$m) ) }) %>%
      arrange(clusterID, rep) %>%
      mutate(Value = 1) %>%
      spread( clusterID, Value) %>% select(-m, -rep)
    comm[ is.na(comm) ] <- 0
    estimate <- vegan::specpool(comm)
    est_num_clusters <- estimate[[method]] #chao, jack1, jack2 or boot
    est_num_clusters <- max( length(m), est_num_clusters ) %>%
                        min( .,  N - sum(m -1) )
    HAI <- 1 - est_num_clusters/N
  }else{
    stop('method must be in: naive, wag, hypergeometric, chao, jack1, jack2, or boot')
  }

  # Do the Bias Correction Next. We use the logit scale to force the resulting
  # estimate to always be in [0,1]. However, we have to be careful because if we
  # get a HAI rate of 0, that results in logit(0) = -Inf and we should avoid that.
  CI_values1 <- NULL
  if( alpha == 1 & bias == TRUE ){
    # warning('No bias correction should be performed if alpha == 1')
    CI_values1 <- c(HAI, HAI)
    bias <- FALSE
  }else if( bias == TRUE ){
    ## Two ways to create a population to do bias correction
    ## 1. Bootstrap an estimated population from the sample
    # Pop <- expand.grid(clusterID=clusterID, rep=1:ceiling(1/alpha) ) %>%
    #   mutate(clusterID = str_c(clusterID,rep, sep='_')) %>% select(-rep)
    ## 2. Subsample the sample
    Pop <- tibble( clusterID = clusterID ) %>% sample_frac(alpha)

    Pop_HAI <- calc_HAI(Pop$clusterID, method='naive', N=nrow(Pop), alpha=1, bias=FALSE, CI=FALSE )
    Pop_HAI <- Pop_HAI[1]
    HAI_reps <- NULL
    for( i in 1:bias_reps ){
      temp <-
        Pop %>% sample_frac(alpha) %>% pull(clusterID) %>%
        calc_HAI(., method=method, alpha=alpha, N=length(clusterID), bias=FALSE, CI=FALSE)
      HAI_reps[i] <- temp[1]
    }
    # if we have any reps with HAI rate > 0, use those.  If not, don't do any bias adjustment.
    HAI_reps <- data.frame(x = HAI_reps) %>% filter( x > 0) %>% pull(x)
    if(length(HAI_reps) > 0 ){
      avg.diff <- mean( logit(HAI_reps) - logit(Pop_HAI[1]) )
      se <- sd(logit(HAI_reps))
      HAI <- ilogit( logit(HAI) - avg.diff )
      CI_values1 <- ilogit(  logit(HAI) + qnorm( c( (1-level)/2, 1-(1-level)/2 ) ) * se )
    }else{
      HAI <- ilogit( logit(HAI)*(2-alpha) )
      CI_values1 <- sort( ilogit( c( logit(HAI)/(2-alpha), logit(HAI)*(2-alpha) ) ) )
    }
  }

  # Calculate the CI if the user asked for it.
  if( CI == TRUE ){
    m = data.frame(clusterID=clusterID) %>% group_by(clusterID) %>% count() %>% pull(n)
    temp <- NULL
    for( i in 1:CI_reps ){
      ## bootstrap the cluster sizes
      # temp <- rbind( temp,
      #   data.frame(m=sample(m, replace=TRUE)) %>%
      #   mutate(clusterID=1:n()) %>%
      #   group_by(clusterID) %>%
      #   do({data.frame(clusterID=.$clusterID, rep=1:.$m) }) %>%
      #   pull( clusterID ) %>%
      #   calc_HAI(., method = method, alpha = alpha, N = N, bias=bias) )

      ## bootstrap the patients
      temp <- data.frame(clusterID=clusterID) %>% sample_frac(replace=TRUE) %>%
        pull( clusterID ) %>%
        calc_HAI(., method = method, alpha = alpha, N = N, bias=bias) %>%
        rbind(temp, .)
    }
    temp <- temp[,1] %>% logit() %>% discard(is.infinite)
    CI_values2 <-
      basic_CI( theta_hat=logit(HAI[1]), theta_hat_star=temp, level=level) %>% ilogit()
    # CI_values2 <- quantile(temp[,1], probs=c( (1-level)/2, 1-(1-level)/2 ) )
    names(CI_values2) <- NULL
  }


  # Return Our Estimate!
  output <- c(est=HAI, lwr1=NA, upr1=NA, lwr2=NA, upr2=NA)
  if( !is.null(CI_values1) ){
    output[2] = CI_values1[1]
    output[3] = CI_values1[2]
  }
  if( CI == TRUE ){
    output[4]=CI_values2[1]
    output[5]=CI_values2[2]
  }
  return(output)

}






