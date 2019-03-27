#' @title Generates parameters for the beta distribution # I don't think we need to show this as a separate function, 
#'                                                         but put together with gen.eff.stg1 or be called by gen.eff.stg1
#' 
#' @description Function \code{beta.ab()} returns parameters alpha and beta for generating beta r.v. (per dose) 
#' 
#' @return Vector of alpha and beta values for generating beta random variable for a dose.
#' 
#' @param m  Vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  Vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
#' 
#' @examples  
#' 
#' @export
#' 
#' @keywords internal

beta.ab <- function(m, v) {
  
  a <- seq(0.5, 20, 0.01)                            # a is a seq of alpha in beta distr.
  b <- a * (1 - m) / m
  
  vfit  <- a * b / ((a + b + 1) * (a + b)^2)
  diff  <- abs(vfit - v)
  index <- (1:length(diff))[diff == min(diff)]       # return the index of the min var.
  
  return(list(a = a[index],
              b = b[index]))                         # return alpha and beta for the min.var.
}