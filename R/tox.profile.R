#' @title Generates DLTs and calculate the likelihood-ratio (LR) for each dose
#' 
#' @description Returns a 4-column matrix containing dose assignment, DLTs at each dose, cohort number, and LR
#' 
#' @return 4-column matrix containing dose assign., DLTs at each dose, cohort number, and LR # This is the same is descripton above. Keep only one?
#' 
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' 
#' @examples
#' dose = 5                                      # Dose levels
#' dose.tox <- c(0.05, 0.10, 0.15, 0.20, 0.30)   # True toxicity per dose
#' p1 = 0.40                                     # Unacceptable DLT rate
#' p2 = 0.15                                     # Acceptable DLT rate
#' K = 2                                         # Likelihood-ratio (LR) threshold
#' coh.size = 3                                  # Assign 3 pts per dose in stage 1

#' @export
 
tox.profile <- function(dose, dose.tox, p1, p2, K, coh.size){ 
  
  dose   <- c(1:dose)                                           # vector of counts up to number of doses given
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1                                        # current cohort corresponding to dose
    dltsi  <- stats::rbinom(1, coh.size, dose.tox[i])	# number of DLTs for that dose based on tox prob
    
    l.p2   <- dltsi * log(p2) + (coh.size - dltsi) * log(1 - p2) # likelihood of acceptable/alternative hypothesis 
    l.p1   <- dltsi * log(p1) + (coh.size - dltsi) * log(1 - p1) # likelihood of unacceptable/null hypothesis
    LR     <- round(exp(l.p2 - l.p1), 2)    
    
    x <- c(x, dose[i], dltsi, cohort, LR)                              			
    
    if (LR <= (1/K))                        # stop escalation
      stop <- 1
    
    if (LR > (1/K))                         # escalate to next dose
      i <- i + 1	 
    
  }       
  return(matrix(x, ncol = 4, byrow = T))
} 