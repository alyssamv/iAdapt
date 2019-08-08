#' @title Generates DLTs and calculate the likelihood-ratio (LR) for each dose
#' 
#' @description Gives toxicity profile (number of dose-limiting toxicities) and likelihood ratio of safety for each dose.
#' 
#' @return 4-column matrix containing dose assignment, dose-limiting toxicities at each dose, 
#' cohort number, and likelihood ratio.
#' 
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' 
#' @examples
#' # Number of pre-specified dose levels
#' dose <- 5
#' # Vector of true toxicities associated with each dose
#' dose.tox <- c(0.05, 0.10, 0.20, 0.35, 0.45)       
#' # Acceptable (p2) and unacceptable (p1) DLT rates used for establishing safety
#' p1 <- 0.40                                     
#' p2 <- 0.15    
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
#' tox.profile(dose = dose, dose.tox = dose.tox, p1 = p1, p2 = p2, K = K, coh.size = coh.size)
#' 
#' 
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