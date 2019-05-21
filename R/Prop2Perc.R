#' Convert proportions to percent labels
#'
#' @param input A vector of numeric values representing proportions [0,1].
#' @return A vector of factors values as percentages, arranged in assending order.
#' @export
prop2perc <- function(input){
  output <- input %>% as.character %>% as.numeric()
  output <- output * 100
  output <- output %>% factor() %>%
    factor( labels=paste(levels(.),'%', sep='') )
  return(output)
}


# prop2perc(
#   input = c(0, 0.1, .3, 1, .8, .2, .1)
#   )




