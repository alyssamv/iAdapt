#' @title Simulate full trial (both stages) x times
#' @description Results are displayed in a matrix format, where each row represents one trial   
#' Returns: dose assignment for each patient in the trial - sim.dk   
#'          estimated persistence per each dose assignment - sim.yk
#' @param nsims - number of simulated trials
#' @param dose - number of doses to be tested (scalar)
#' @param dose.tox - vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p0 - toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p1 - toxicity under alternative (safe DLT rate). Values range from 0 - 1; p0 > p1
#' @param K - threshold for LR. Takes integer values: 1,2,... (recommended K=2)
#' @param coh.size - cohort size (number of patients) per dose
#' @param m - mean efficacy of a dose (single value). Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v - efficacy variances of a dose (single value). Values range from 0 - 1. (e.g., 0.01)
#' @param N - max sample size for stages 1&2
#' @param stop.rule - if only dose 1 safe, allocate up to 9 (default) patients at dose 1 to collect more info
#' @param nbb - binomial parameter(default=100 cells per patient)


sim.trials <-function(numsims, dose, dose.tox, p0, p1, K, coh.size, m, v, N, stop.rule = 9, cohort = 1, samedose = T, nbb = 100){
  
  sim.yk <- sim.dk <- matrix(NA, nrow=numsims, ncol=N)    
  
  for (i in 1:numsims) {  
    
    fstudy.out <- f.study.a(dose, dose.tox, p0, p1, K, coh.size, m, v, N, stop.rule=9, cohort=1, samedose=T, nbb=100)                    
    
    if (length(fstudy.out$yk.final) < N){                           # if max sample size was not reached, fill-in with NAs        
      
      sim.yk[i,] <- c(fstudy.out$yk.final, rep(NA, N-length(fstudy.out$yk.final)))
      sim.dk[i,] <- c(fstudy.out$dk.final, rep(NA, N-length(fstudy.out$dk.final)))
      
    } else {
      
      sim.yk[i,] <- fstudy.out$yk.final
      sim.dk[i,] <- fstudy.out$dk.final
    }
    cat(i,"\n")
  }
  return(list(sim.yk=sim.yk, sim.dk=sim.dk))
}