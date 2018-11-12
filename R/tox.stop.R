#' @title Generate DLTs and calculate the LR for each dose
#' 
#' @description Returns a 4-column matrix containing dose assign., DLTs at each dose, cohort number, and LR
#' 
#' @return 4-column matrix containing dose assign., DLTs at each dose, cohort number, and LR
#' 
#' @param dose - number of doses to be tested (scalar)
#' @param dose.tox - vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p0  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p1  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p0 > p1
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dosetox.stop <- function(dose, dose.tox, p0, p1, K, coh.size)
#' 
#' @examples
#' dose = 5                                      # dose levels
#' dose.tox <- c(0.05, 0.10, 0.15, 0.20, 0.30)   # True toxicity per dose
#' p0 = 0.40                                     # Unacceptable DLT rate
#' p1 = 0.15                                     # Acceptable DLT rate
#' K = 2                                         # Likelihood-ratio (LR) threshold
#' coh.size = 3                                  # (3 pts per dose in stage 1)
#' m <- c(5, 15, 40, 65, 80)                     # Mean efficacy per dose (here mean persistence per dose (per cents)
#' v <- rep(0.01,5)                              # Efficacy variance per dose
#' N = 25                                        # Total sample size (stages 1&2)
#' stop.rule = 9                                 # If dose 1 is the only safe dose remaining, allocate up to 9 pts. to collect more info
#' 
#' tox.stop(dose = dose, dose.tox = dose.tox, p0 = p0, p1 = p1, K = K, coh.size = coh.size)
#' 
#' @export
 
tox.stop <- function(dose, dose.tox, p0, p1, K, coh.size){ 
  dose   <- c(1:dose) # vector of counts up to number of doses given
  
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1 # current cohort corresponding to dose
    dltsi  <- rbinom(1, coh.size, dose.tox[i])	# number of DLTs for that dose based on tox prob
    
    l.p1   <- dltsi * log(p1) + (coh.size - dltsi) * log(1 - p1) # likelihood of alternative hypothesis 
    l.p0   <- dltsi * log(p0) + (coh.size - dltsi) * log(1 - p0) # likelihood of null hypothesis
    LR     <- round(exp(l.p1 - l.p0), 2)    
    
    x <- c(x, dose[i], dltsi, cohort, LR)                              			
    
    if (LR <= (1/K))                        # stop escalation
      stop <- 1
    
    if (LR > (1/K))                         # escalate to next dose
      i <- i + 1	 
    
  }       
  return(matrix(x, ncol = 4, byrow = T))
} 