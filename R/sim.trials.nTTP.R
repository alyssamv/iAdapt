#' @title Simulate full trial (both stages) x times when using nTTP to measure toxicity
#' 
#' @description Results are displayed in a matrix format, where each row represents one trial simulation  
#' 
#' @return List of the following objects:
#' \itemize{
#' \item sim.Y - estimated efficacy per each dose assignment 
#' \item sim.d - dose assignment for each patient in the trial 
#' }
#'         
#'          
#' @param numsims  number of simulated trials
#' @param dose  number of doses to be tested (scalar)
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
#' @param nbb  binomial parameter (default = 100 cells per patient)
#' @param N  total sample size for stages 1&2
#' @param stop.rule  if only dose 1 safe, allocate up to 9 (default) patients at dose 1 to collect more info
#' @param cohort cohort size (number of patients) per dose (Stage 2). Default is 1.
#' @param samedose designates whether the next patient is allocated to the same dose as the previous patient. Default is TRUE. Function adjusts accordingly.
#' @param ntox  number (integer) of different toxicity types
#' @param W  matrix defining burden weight of each grade level for all toxicity types. The dimensions are ntox rows by 4 columns (for grades 1-4. A grade 0 event has a burden of 0). See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. Each matrix represents one toxicity type, where probabilities of each toxicity grade are specified across each dose. Each matrix has the same dimensions: n rows, representing number of doses, and 5 columns (for grades 0-4, since the probability of a grade 0 event may not be 0). Probabilities across each dose (rows) must sum to 1. See Ezzalfani et al. (2013) for details.
#' 
#' 
#' @examples
#' # Number of pre-specified dose levels
#' dose <- 6      
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
#' # Toxicity burden weight matrix
#' W = matrix(c(0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 1
#'              0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 2
#'              0.00, 0.00, 0.5, 1), ## Burden weight for grades 1-4 for toxicity 3
#'            nrow = ntox, byrow = T)
#'            
#' # Number of toxicity types
#' ntox = 3
#' 
#' # Array of toxicity event probabilities
#' TOX = array(NA, c(d, 5, ntox)) 
#' 
#' TOX[, , 1] = matrix(c(0.823, 0.152, 0.022, 0.002, 0.001,
#'                       0.791, 0.172, 0.032, 0.004, 0.001,
#'                       0.758, 0.180, 0.043, 0.010, 0.009,
#'                       0.685, 0.190, 0.068, 0.044, 0.013,
#'                       0.662, 0.200, 0.078, 0.046, 0.014,
#'                       0.605, 0.223, 0.082, 0.070, 0.020),
#'                     nrow = 6, byrow = TRUE)
#' TOX[, , 2] = matrix(c(0.970, 0.027, 0.002, 0.001, 0.000,
#'                       0.968, 0.029, 0.002, 0.001, 0.000,
#'                       0.813, 0.172, 0.006, 0.009, 0.000,
#'                       0.762, 0.183, 0.041, 0.010, 0.004,
#'                       0.671, 0.205, 0.108, 0.011, 0.005,
#'                       0.397, 0.258, 0.277, 0.060, 0.008),
#'                     nrow = 6, byrow = TRUE)
#' TOX[, , 3] = matrix(c(0.930, 0.060, 0.005, 0.001, 0.004,
#'                       0.917, 0.070, 0.007, 0.001, 0.005,
#'                       0.652, 0.280, 0.010, 0.021, 0.037,
#'                       0.536, 0.209, 0.031, 0.090, 0.134,
#'                       0.015, 0.134, 0.240, 0.335, 0.276,
#'                       0.005, 0.052, 0.224, 0.372, 0.347),
#'                     nrow = 6, byrow = TRUE)
#' 
#' sim.trials.nTTP(numsims = 10, dose = dose, p1 = p_no, p2 = p_yes, K = K, 
#' coh.size = coh.size, m = m, v = v, N = N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100,
#' W = W, TOX = TOX, ntox = ntox)
#' 
#' @export


sim.trials.nTTP <- function (numsims, dose, p1, p2, K, coh.size, m, v, 
                             N, stop.rule = 9, cohort = 1, samedose = TRUE, nbb = 100, W, TOX, ntox) 
{
  sim.yk <- sim.dk <- matrix(NA, nrow = numsims, ncol = N)
  sim.doses <- matrix(NA, nrow = numsims, ncol = dose)
  for (i in 1:numsims) {
    fstudy.out <- rand.stg2.nTTP(dose, p1, p2, K, coh.size, 
                                 m, v, N, stop.rule, cohort, samedose, nbb, W, TOX, ntox)
    n.safe <- max(fstudy.out$d.final[(fstudy.out$n1 + 1):length(fstudy.out$d.final)], 
                  na.rm = TRUE)
    sim.doses[i, ] <- c(rep(1, n.safe), rep(0, dose - n.safe))
    if (length(fstudy.out$Y.final) < N) {
      sim.yk[i, ] <- c(fstudy.out$Y.final, rep(NA, N - 
                                                 length(fstudy.out$Y.final)))
      sim.dk[i, ] <- c(fstudy.out$d.final, rep(NA, N - 
                                                 length(fstudy.out$d.final)))
    } else {
      sim.yk[i, ] <- fstudy.out$Y.final
      sim.dk[i, ] <- fstudy.out$d.final
    }
    cat(i, "\n")
  }
  return(list(sim.Y = sim.yk, sim.d = sim.dk, safe.d = sim.doses))
}
