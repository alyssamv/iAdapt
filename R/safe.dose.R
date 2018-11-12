#' @title Identify safe doses  
#' 
#' @description Function safe.dose() distinguishes acceptable from unacceptable doses
#' 
#' @return List of the following objects:
#'          alloc.safe - matrix of assign. for only acceptable doses (to be used in stage 2) and the 
#'                       corresponding toxicities
#'          alloc.total - vector of all assign. from stage 1 
#'          n1 - total number of pts. allocated in stage 1
#'          
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p0  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p1  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p0 > p1
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose
#' 
#' @examples
#' 
#' @export


safe.dose <- function(dose, dose.tox, p0, p1, K, coh.size) {
  
  res         <- tox.stop(dose, dose.tox, p0, p1, K, coh.size)     # output from tox.stop()
  alloc.total <- sort(rep(res[,1], coh.size))                      # sort according to dose/cohort size           
  n1          <- nrow(res)*coh.size                                # total number of patients assigned to doses              
  unsafe.dose <- which(res[, 4] <= (1/K))                          # which doses are unacceptably toxic
  
  if (length(unsafe.dose) == 0) {  # if-else to return only those rows for safe doses
    alloc.safe <- res[, 1:2]
  } else {
    alloc.safe <- res[res[,1] != unsafe.dose, 1:2]
  }
  return(list(alloc.safe = alloc.safe, 
              alloc.total = alloc.total, 
              n1 = n1))  # return named list
}