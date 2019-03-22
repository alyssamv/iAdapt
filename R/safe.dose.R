#' @title Identify safe/acceptable doses  
#' 
#' @description Function safe.dose() distinguishes acceptable from unacceptable doses
#' 
#' @return List of the following objects:
#'          alloc.safe - matrix of assignments only for acceptable doses (to be used in stage 2) and their 
#'                       corresponding toxicities
#'          alloc.total - vector of all dose assignments from stage 1 
#'          n1 - total number of subjects allocated in stage 1
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


safe.dose <- function(dose, dose.tox, p1, p2, K, coh.size) {
  
  res         <- tox.profile(dose, dose.tox, p1, p2, K, coh.size)  # save output from tox.profile()
  alloc.total <- sort(rep(res[,1], coh.size))                      # sort according to dose/cohort size           
  n1          <- nrow(res)*coh.size                                # total number of patients assigned to doses              
  unsafe.dose <- which(res[, 4] <= (1/K))                          # identify which doses are unacceptably toxic
  
  if (length(unsafe.dose) == 0) {  # if-else to return only those rows for safe doses
    alloc.safe <- res[, 1:2]
  } else {
    alloc.safe <- res[res[,1] != unsafe.dose, 1:2]
  }
  return(list(alloc.safe = alloc.safe, 
              alloc.total = alloc.total, 
              n1 = n1))  # return named list
}