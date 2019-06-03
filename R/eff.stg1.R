#' @title Generates efficacy outcomes for stage 1
#' 
#' @description Function \code{eff.stg1()} uses a beta-binomial distribution to generate outcomes (Ys) corresponding 
#' to acceptable dose assignments from stage 1. 
#' 
#' @return List of efficacy outcomes for subject enrolled during stage 1 (dose-escalation)
#' \itemize{
#' \item Y.safe - vector of efficacy outcomes for each subject enrolled on an acceptably toxic dose
#' \item d.safe - vector of dose allocation for each subject enrolled on an acceptably toxic dose
#' \item tox.safe - number of dose-limiting toxicities for each safe dose level
#' \item Y.alloc - vector of efficacy outcomes for all subjects enrolled on all doses (safe and unsafe)
#' \item d.alloc - vector of dose allocation for all subjects enrolled on all doses (safe and unsafe)
#' }
#' 
#'         
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
#' @param nbb  binomial parameter (default = 100 cells per patient)
#' 
#' 
#' @examples 
#' # Number of pre-specified dose levels
#' dose <- 5
#' # Vector of true toxicities associated with each dose
#' dose.tox <- c(0.05, 0.10, 0.20, 0.35, 0.45)       
#' # Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety
#' p_no <- 0.40                                     
#' p_yes <- 0.15    
#' 
#' # Likelihood-ratio (LR) threshold
#' K <- 2                                          
#' 
#' # Cohort size used in stage 1
#' coh.size <- 3 
#' 
#' # Vector of true mean efficacies per dose (here mean percent persistence per dose)
#' m <- c(5, 15, 40, 65, 80)   # MUST BE THE SAME LENGTH AS dose.tox                  
#' 
#' # Efficacy(equal) variance per dose
#' v <- rep(0.01, 5) 
#' 
#' # Total sample size (stages 1&2)                            
#' N <- 25                                        
#' 
#' # Stopping rule: if dose 1 is the only safe dose, allocate up to 9 pts.
#' stop.rule <- 9 
#' 
#' eff.stg1(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size, m, v, nbb = 100)
#' 
#' @export

 
eff.stg1 <- function(dose, dose.tox, p1, p2, K, coh.size, m, v, nbb = 100) {
  
  res <- safe.dose(dose, dose.tox, p1, p2, K, coh.size)
  d.alloc <- res$alloc.total
  val.safe <- res$alloc.safe
  
  Y.safe <- d.safe <- tox.safe <- Y.alloc <- NULL      
  n1 <- res$n1  
  
  for (i in 1:length(d.alloc)) {                            # Generate Ys for all allocations in stage 1
    ab <- beta.ab(m[d.alloc[i]]/100, v[d.alloc[i]])
    p <- stats::rbeta(1, ab$a, ab$b)
    Y.alloc[i] <- 100*stats::rbinom(1, nbb, p) / nbb                  
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
    Y.safe <- d.safe  <- NULL                           # if dose 1 is unacceptable, leave it null
    tox.safe <- res$alloc.safe[,2]
  }         
  return(list(Y.safe = Y.safe,
              d.safe = d.safe,
              tox.safe = tox.safe,
              n1 = n1,
              Y.alloc = Y.alloc,
              d.alloc = d.alloc))  
}