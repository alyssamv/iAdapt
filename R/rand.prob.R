#' @title Calculates randomization probabilities and dose allocation for next patient
#' 
#' @description Function \code{rand.prob()} calculates the updated randomization probabilities based on observed efficacies up to that point.
#' It also gives the dose allocation for the next enrolled patient based on these probabilities.
#' 
#'         
#' @param y.eff vector of all efficacy outcomes for each dose allocation
#' @param d.safe vector of dose assignment
#' 
#' @examples 
#' y.eff <- c(9,1,0,34,10,27,38,42,60,75,48,62)
#' d.safe <- c(1,1,1,2,2,2,3,3,3,4,4,4)
#' rand.prob(y.eff, d.safe)
#' 
#' @export

rand.prob <- function(y.eff, d.safe){
  
  reg <- lm(log(y.eff + 1) ~ factor(d.safe))         # Linear model with log(Y) for accept. doses 
  fit <- as.vector(reg$fitted.values)                # Fitted values for Y
  fitp <- exp(fit)
  
  dose.unique <- duplicated(d.safe)
  fitp <- fitp[dose.unique == F]
  
  rp <- fitp/sum(fitp)                                # Calculate randomization prob. for each dose                            
  rp <- ifelse(rp < 0.02, 0.02, rp) 
  
  rec.dose <- which(rp==max(rp))                      # Next dose with max rand. prob
  
  return(list(Rand.Prob=rp, Next.Dose=rec.dose))
  
}

