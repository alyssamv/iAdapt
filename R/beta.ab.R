#' @title Get Beta parameters
#' @description Function beta.ab() returns parameters alpha and beta for generating beta r.v. (per dose) 
#' @param m - mean efficacy of a dose (single value). Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v - efficacy variances of a dose (single value). Values range from 0 - 1. (e.g., 0.01)

beta.ab <- function(m, v) {
  
  a <- seq(0.5, 20, 0.01)                            # a is a seq of alpha in beta distr.
  b <- a * (1 - m) / m
  
  vfit  <- a * b / ((a + b + 1) * (a + b)^2)
  diff  <- abs(vfit - v)
  index <- (1:length(diff))[diff == min(diff)]       # return the index of the min var.
  
  return(list(a = a[index],
              b = b[index]))                         # return alpha and beta for the min.var.
}