#' @title Visualize simulation results
#' 
#' @description Results from simulated trials (using \code{sim.trials()} function) displayed in tabular and/or graphical format
#' 
#' @return Printed tables and a list of the following objects:
#' \itemize{
#' \item pct.treated - IQR (25th percentile, median, 75th percentile) of percent of subjects treated at each dose level
#' \item efficacy - IQR of efficacy observed at each dose level
#' }
#'          
#' @param sims output from sim.trials
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
#' simulations = sim.trials(numsims = 100, dose, dose.tox, p1 = p_no, p2 = p_yes, K, 
#' coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100)
#' 
#' summary = sim.summary(simulations)
#'  
#' @export
#' 


sim.summary <- function(sims){
  sim.doses = sims$sim.d
  n.doses = max(sim.doses, na.rm = T)
  
  sim.eff = sims$sim.Y
  
  ## columns are IQR for percent of patients treated on each dose
  dose.mat.a <- matrix(NA, nrow(sim.doses), n.doses)
  for (i in 1:nrow(sim.doses)) {     
    dose.no.na <- na.omit(sim.doses[i, ])   
    dose.mat.a[i, ] <- table(factor(dose.no.na, levels = 1:n.doses))/length(dose.no.na)       
  }
  est.dose1 <- matrix(NA, n.doses, 4) 
  for (j in 1:n.doses) {
    est.dose1[j, ] <- c(j/100, quantile(dose.mat.a[, j], prob = c(0.25, 0.5, 0.75), na.rm = T))	
  }
  dose.IQR = round(est.dose1*100, 1)
  
  print(knitr::kable(dose.IQR, 
                     caption = "Percent of subjects treated on each dose level",
                     col.names = c("Dose", "25th percentile", "Median", "75th percentile")))
  
  
  ## IQR of persistence (columns) for each dose (rows)
  pers.hat.a <- matrix(NA, nrow(sim.eff), n.doses + 1)                                     # Median 
  for (i in 1:nrow(sim.eff)) {
    for (j in 1:(n.doses + 1)) {
      pers.hat.a[i,j] <- (median(sim.eff[i, sim.doses[i, ] == j - 1]))
    }
  }
  est.pers1 <- matrix(NA, (n.doses + 1), 4)  
  for (j in 1:(n.doses + 1)) {
    est.pers1[j, ] <- c((j - 1), quantile(pers.hat.a[,j], prob = c(0.25, 0.5, 0.75), na.rm = T))
  }
  
  Y = est.pers1[-1, ]
  
  print(knitr::kable(Y, 
                     caption = "Efficacy at each dose level",
                     col.names = c("Dose", "25th percentile", "Median", "75th percentile")))
  
  
  return(list(pct.treated = dose.IQR, efficacy = Y))
}
