#' @title Calculates likelihood of safety for single dose
#' 
#' @description Function \code{LRtox()} calculates the likelihood of safety for a single dose 
#' and designates whether to escalate to the next dose (safe) or stop dose escalation and move onto stage 2 (unsafe).
#' 
#' @return List object that gives the likelihood ratio of safety and indicates whether to escalate to the 
#' next highest dose level, or stop dose escalation and move onto stage 2.
#'         
#' @param x  number of observed DLTs
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' 
#' 
#' @examples 
#' LRtox(coh.size=3,x=2,p1=0.40,p2=0.15,K=2)
#' LRtox(coh.size=3,x=1,p1=0.40,p2=0.15,K=2)
#' 
#' @export


LRtox <- function(coh.size, x, p1, p2, K = 2) {
  
  l.p2   <- x * log(p2) + (coh.size-x) * log(1 - p2) # likelihood of acceptable/alternative hypothesis 
  l.p1   <- x * log(p1) + (coh.size-x) * log(1 - p1) # likelihood of unacceptable/null hypothesis
  LR     <- round(exp(l.p2 - l.p1), 2) 
  
  if (LR > 1/K)  print ("Safe/Escalate") else print ("Unsafe/Stop")
  
  return(list(LR=LR))
}
