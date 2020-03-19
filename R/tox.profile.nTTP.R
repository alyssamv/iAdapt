#' @title (nTTP toxicity endpoint) Generates nTTPs and calculate the likelihood-ratio (LR) for each dose
#' 
#' @description Gives toxicity profile (number of dose-limiting toxicities) and likelihood ratio of safety for each dose.
#' 
#' @return 4-column matrix containing dose assignment, mean nTTP at each dose, 
#' cohort number, and likelihood ratio.
#' 
#' @param dose  number of doses to be tested (scalar)
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param ntox  number (integer) of different toxicity types
#' @param W  matrix defining burden weight of each grade level for all toxicity types. The dimensions are ntox rows by 4 columns (for grades 1-4. A grade 0 event has a burden of 0). See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. Each matrix represents one toxicity type, where probabilities of each toxicity grade are specified across each dose. Each matrix has the same dimensions: n rows, representing number of doses, and 5 columns (for grades 0-4, since the probability of a grade 0 event may not be 0). Probabilities across each dose (rows) must sum to 1. See Ezzalfani et al. (2013) for details.
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
#' tox.profile(dose = dose, p1 = p1, p2 = p2, K = K, coh.size = coh.size, ntox = ntox, W = W, TOX = TOX)
#' 
#' 
#' @export

tox.profile.nTTP <- function(dose, p1, p2, K, coh.size, ntox, W, TOX){ 
  
  dose   <- c(1:dose)                                           # vector of counts up to number of doses given
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  
  # bounds for nTTP (truncated normal distribution)
  a = 0
  b = 1
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1                                        # current cohort corresponding to dose
    dltsi  <- replicate(coh.size, nTTP.sim(W = W, 
                                           TOX = TOX, 
                                           ntox = ntox, 
                                           dose = dose[i]))	# nTTPs for that dose based on tox prob
    
    l.p2   <- prod(sapply(dltsi, FUN = function(i){ dnorm((i - p2)/sigma) })) / (sigma*(pnorm((b - p2)/sigma) - pnorm((a - p2)/sigma)))^coh.size # likelihood of acceptable/alternative hypothesis 
    l.p1   <- prod(sapply(dltsi, FUN = function(i){ dnorm((i - p1)/sigma) })) / (sigma*(pnorm((b - p1)/sigma) - pnorm((a - p1)/sigma)))^coh.size # likelihood of unacceptable/null hypothesis
    LR     <- round(l.p2/l.p1, 2)    
    
    x <- c(x, dose[i], mean(dltsi), cohort, LR)                              			
    
    if (LR <= (1/K)) {                       # stop escalation
      stop <- 1
    } else if (LR > (1/K)) {                        # escalate to next dose
      i <- i + 1	 
    }
    
  }       
  return(matrix(x, ncol = 4, byrow = TRUE))
} 