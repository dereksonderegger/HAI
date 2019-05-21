#' Create a population of HAI cases
#' @param n - Total number of patients that we could have sampled
#' @param lambda - Poisson paramter
#' @param mix.frac Mixing fraction of LogNoramal to be included
#' @return A data frame of approximately n patients, each with a ClusterID
#'         and if their infection was Community or Hospital aquired
#' @export
generate_population <- function(n, lambda, mix.frac=0.01){
    data.frame(clusterID=1:n) %>%
    mutate( type = ifelse( runif(n) > mix.frac, 'Poisson', 'logNormal') ) %>% group_by(clusterID) %>%
    mutate( clusterSize = ifelse( type == 'Poisson',
                                  truncdist2::rtrunc(100, 'pois', a=0, lambda=lambda),
                                  round(rlnorm(1, meanlog=log(35), sdlog=.25))) ) %>%
    group_by() %>% mutate( TotalPatients = cumsum(.$clusterSize) ) %>%
    filter(TotalPatients <= n) %>%
    select( -type ) %>%
    group_by( clusterID ) %>%
    do({ expand.grid( Rep = 1:.$clusterSize[1]) }) %>%
    mutate(Aquired = ifelse( Rep == 1, 'Community', 'Hospital')) %>%  # First case in the cluster is
    select(-Rep) %>% group_by() %>%                                   # assumed to be community aquired
    mutate( patientID = 1:n() ) %>%
    select(clusterID, patientID, Aquired)
}
