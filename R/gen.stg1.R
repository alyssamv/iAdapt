#' @title Generate stage 1 - safety 
#' @description Function gen.y.a() uses a beta-binomial distribution to generate outcomes (Ys) corresponding 
#'  to all dose assign. in stage 1
#'  Returns: Y.alloc - vector of Ys corresponding to all dose assign. from stage 1 
#'                      (acceptable and unacceptable)
#'            d.alloc - vector of all dose assign. from stage 1
#'             Y.safe - subset of Y.alloc corresponding to only acceptable doses; to be used for 
#'                      randomization in stage 2
#'             d.safe - subset of d.alloc corresponding to only acceptable dose assign.; to be used for 
#'                      randomization in stage 2
#'          tox. safe - vector of toxicities corresponding to acceptable doses
#'                 n1 - total number of pts. allocated in stage 1
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p0  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p1  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p0 > p1
#' @param K  threshold for LR. Takes integer values: 1,2,... (recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose
#' @param m  mean efficacy of a dose (single value). Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  efficacy variances of a dose (single value). Values range from 0 - 1. (e.g., 0.01)
#' @param nbb  binomial parameter (default = 100 cells per patient)


gen.stg1 <- function(dose, dose.tox, p0, p1, K, coh.size, m, v, nbb = 100) {
  
  res <- safe.dose(dose, dose.tox, p0, p1, K, coh.size)
  d.alloc <- res$alloc.total
  val.safe <- res$alloc.safe
  
  Y.safe <- d.safe <- tox.safe <- Y.alloc <- NULL      
  n1 <- res$n1  
  
  for (i in 1:length(d.alloc)) {                            # Generate Ys for all allocations in stage 1
    ab <- beta.ab(m[d.alloc[i]]/100, v[d.alloc[i]])
    p <- rbeta(1, ab$a, ab$b)
    Y.alloc[i] <- 100*rbinom(1, nbb, p) / nbb                  
  }                          
  
  if (length(val.safe) > 2) {                               # if 2 or more acceptable doses
    d.safe <- sort(rep(val.safe[, 1], coh.size))
    tox.safe <- res$alloc.safe[, 2]
    Y.safe <- Y.alloc[1:length(d.safe)]
    
  } else if (length(val.safe) == 2) {                   # if only dose 1 acceptable
    d.safe <- sort(rep(val.safe[1], coh.size))
    tox.safe <- res$alloc.safe[2]
    Y.safe <- Y.alloc[1:length(d.safe)]
    
  } else {
    Y.safe <- d.safe  <- NULL                     # if dose 1 is unacceptable, leave it null
    tox.safe <- res$alloc.safe[,2]
  }         
  return(list(Y.safe = Y.safe,
              d.safe = d.safe,
              tox.safe = tox.safe,
              Y.alloc = Y.alloc,
              d.alloc = d.alloc,
              n1 = n1))  
}