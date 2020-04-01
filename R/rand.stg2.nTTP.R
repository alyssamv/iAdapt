#' @title Stage 2 Adaptive Randomization with nTTP to monitor toxicity
#' 
#' @description Function \code{rand.stg2.nTTP()} fits a linear regression for the continuous efficacy outcomes,
#' computes the randomization probabilities/dose and allocates the next patient to a dose that
#' is considered acceptably safe and has the highest efficacy. Dose safety (with nTTP) is still monitored using LR and doses
#' that become unacceptable are discarded.
#' 
#' @return List of the following objects:
#' \itemize{
#' \item Y.final - vector of all efficacy outcomes (Ys) corresponding to dose assignments (Stages 1&2)
#' \item d.final - vector of all dose assignments(Stage 1&2)
#' }
#' If no dose allocation, put NAs in d.final and y.final.
#' 
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
#' rand.stg2.nTTP(dose = dose, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size, 
#' m = m, v = v, N = N, stop.rule = stop.rule, 
#' cohort = 1, samedose = TRUE, nbb = 100, W = W, TOX = TOX, ntox = ntox) 
#' 
#' @export


rand.stg2.nTTP <- function(dose, p1, p2, K, coh.size, m, v, N, stop.rule = 9, 
                           cohort = 1, samedose = TRUE, nbb = 100, W, TOX, ntox) {
  
  res <- eff.stg1.nTTP(dose, p1, p2, K, coh.size, m, v, nbb, W, TOX, ntox)
  dose   <- c(1:dose) 
  yk.safe <- res$Y.safe
  yk.final <- res$Y.alloc
  dk.safe <- res$d.safe
  dk.final <- dk1 <- dk2 <- res$d.alloc
  toxk <- res$tox.safe
  n1 <- res$n1
  nmore <- N - n1
  nd <- length(unique(dk.safe))
  rp <- NULL
  stop <- 0
  
  if (nd == 0) {
    yk.final <- yk.final
    dk.final <- dk.final
    stop <- 1
  }
  
  if (nd == 1) {
    extra <- stop.rule - length(dk.safe)
    ab <- beta.ab(m[1]/100, v[1])
    y.extra <- 100 * stats::rbinom(extra, nbb, stats::rbeta(1, 
                                                            ab$a, ab$b))/nbb
    yk.final <- c(yk.final, y.extra)
    dk.final <- c(dk.final, rep(1, extra))
    stop <- 1
  }
  
  if (nd > 1) {
    coh.toxk <- cbind(matrix(dk.safe, ncol = coh.size, byrow = TRUE)[, 
                                                                     1], toxk)
    for (k in 1:nmore) {
      if (stop == 0) {
        reg <- stats::lm(log(yk.safe + 1) ~ factor(dk.safe))
        fit <- as.vector(reg$fitted.values)
        dose.unique <- duplicated(dk.safe)
        fitp <- exp(fit)
        fitp <- fitp[dose.unique == FALSE]
        rp <- fitp/sum(fitp)
        rp <- ifelse(rp < 0.02, 0.02, rp)
        dj <- stats::rmultinom(1, 1, prob = rp)
        if (samedose) {
          dj <- rep((1:length(dj))[dj == 1], cohort)
        } else {
          dosemat <- as.vector(dj * matrix(1:nd, ncol = cohort, 
                                           nrow = nd))
          dj <- dosemat[dosemat > 0]
        }
        ab <- beta.ab(m[dj]/100, v[dj])
        p <- stats::rbeta(1, ab$a, ab$b)
        yj <- 100 * stats::rbinom(1, nbb, p)/nbb
        toxj <- nTTP.sim(W = W, 
                         TOX = TOX, 
                         ntox = ntox, 
                         dose = dose[dj]) # stats::rbinom(1, size = 1, dose.tox[dj])
        coh.toxj <- c(dj, toxj)
        yk.safe <- c(yk.safe, yj)
        yk.final <- c(yk.final, yj)
        dk.safe <- c(dk.safe, dj)
        dk.final <- c(dk.final, dj)
        coh.toxk <- rbind(coh.toxk, coh.toxj)
        toxk <- c(toxk, toxj)
        n.obsk <- table(dk.safe)
        if (toxj == 0) {
          dk.safe <- dk.safe
          yk.safe <- yk.safe
        } else {
          LR.table.temp <- table(coh.toxk[, 1], coh.toxk[, 2])
          if (ncol(LR.table.temp) == 2) {
            LR.table <- cbind(LR.table.temp[, 2], n.obsk)
          } else {
            LR.table <- cbind(LR.table.temp[, 1], n.obsk)
          }
          loglik.p2 <- NULL
          loglik.p1 <- NULL
          lik.diff <- NULL
          accept.dose <- NULL
          for (j in 1:nrow(LR.table)) {
            loglik.p2[j] <- LR.table[j, 1] * log(p2) + 
              (LR.table[j, 2] - LR.table[j, 1]) * log(1 - 
                                                        p2)
            loglik.p1[j] <- LR.table[j, 1] * log(p1) + 
              (LR.table[j, 2] - LR.table[j, 1]) * log(1 - 
                                                        p1)
            lik.diff[j] <- exp(loglik.p2[j] - loglik.p1[j])
            accept.dose[j] <- ifelse(lik.diff[j] > (1/K), 
                                     1, 0)
          }
          if (sum(accept.dose) < max(dose)) {
            dk.safe[dk.safe >= which(accept.dose == 0)[[1]]] <- NA
          }
          new.model <- cbind(dk.safe, yk.safe)
          new.model <- stats::na.omit(new.model)
          dk.safe <- new.model[, 1]
          yk.safe <- new.model[, 2]
          yk.final <- yk.final
          dk.final <- dk.final
          if (sum(accept.dose) < max(dose)) {
            coh.toxk <- coh.toxk[!apply(coh.toxk, 1, function(x) {
              any(x >= which(accept.dose == 0)[1])
            }), ]
          } 
        }
        if (length(unique(dk.safe)) > 1) {
          dk.safe <- dk.safe
          yk.safe <- yk.safe
          dk.final <- dk.final
          yk.final <- yk.final
        } else if (length(unique(dk.safe)) == 1) {
          new.size <- nmore + length(dk2)
          length.dk1 <- length(dk.final)
          if ((length(dk.safe) < stop.rule) && (length.dk1 < 
                                                new.size)) {
            extra.one <- min(new.size - length.dk1, 
                             stop.rule - length(dk.safe))
            ab <- beta.ab(m[1]/100, v[1])
            yj.one <- 100 * stats::rbinom(extra.one, 
                                          nbb, stats::rbeta(1, ab$a, ab$b))/nbb
            yk.final <- c(yk.final, yj.one)
            dk.final <- c(dk.final, rep(1, extra.one))
            stop <- 1
          } else {
            dk.final <- dk.final
            yk.final <- yk.final
            stop <- 1
          }
        }
        if (length(unique(dk.safe)) < 1) {
          dk.final <- dk.final
          yk.final <- yk.final
          stop <- 1
        }
      } else {
        dk.final <- dk.final
        yk.final <- yk.final
      }
    }
  }
  return(list(Y.final = yk.final, d.final = dk.final, n1 = n1))
}