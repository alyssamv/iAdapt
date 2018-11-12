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