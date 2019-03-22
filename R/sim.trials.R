#' @title Simulate full trial (both stages) x times
#' 
#' @description Results are displayed in a matrix format, where each row represents one trial simulation  
#' 
#' @return List of the following objects:
#'           sim.Y - estimated efficacy per each dose assignment 
#'           sim.d - dose assignment for each patient in the trial 
#'         
#'          
#' @param numsims  number of simulated trials
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
#' @param nbb  binomial parameter (default = 100 cells per patient)
#' @param N  total sample size for stages 1&2
#' @param stop.rule  if only dose 1 safe, allocate up to 9 (default) patients at dose 1 to collect more info
#' @param cohort ############
#' @param samedose ##########
#' 
#' @examples
#' 
#' @export


sim.trials <- function(numsims, dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule = 9, cohort = 1, samedose = T, nbb = 100){
  
  sim.yk <- sim.dk <- matrix(NA, nrow = numsims, ncol = N)    
  
  for (i in 1:numsims) {  
    
    fstudy.out <- rand.stg2(dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule, cohort, samedose, nbb)                  
    
    if (length(fstudy.out$Y.final) < N) {                           # if max sample size was not reached, fill-in with NAs        
      
      sim.yk[i,] <- c(fstudy.out$Y.final, rep(NA, N - length(fstudy.out$Y.final)))
      sim.dk[i,] <- c(fstudy.out$d.final, rep(NA, N - length(fstudy.out$d.final)))
      
    } else {
      
      sim.yk[i,] <- fstudy.out$Y.final
      sim.dk[i,] <- fstudy.out$d.final
    }
    cat(i,"\n")
  }
  return(list(sim.Y = sim.yk, sim.d = sim.dk))
}