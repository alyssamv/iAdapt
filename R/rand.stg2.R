#' @title Stage 2 Adaptive Randomization
#' 
#' @description Function rand.stg2() fits a linear regression for the continuous efficacy outcomes,
#' computes the randomization probabilities/dose and allocates the next patient to a dose that
#' is considered acceptably safe and has the highest efficacy. Dose safety is still monitored using LR and doses
#' that become unacceptable are discarded.
#' 
#' @return List of the following objects:
#'        Y.final - vector of all efficacy outcomes (Ys) corresponding to dose assignments (Stages 1&2)
#'        d.final - vector of all dose assignments(Stage 1&2)
#' If no dose allocation, put NAs in d.final and y.final
#' 
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
#' @param cohort ##########
#' @param samedose ###########
#' 
#' @examples
#' 
#' @export


rand.stg2 <- function(dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule=9, cohort=1, samedose=T, nbb=100) {
  
  res <- eff.stg1(dose, dose.tox, p1, p2, K, coh.size, m, v, nbb)
  yk.safe <- res$Y.safe                                    
  yk.final <- res$Y.alloc                   
  dk.safe <- res$d.safe                                          # Safe doses from stage 1 used for randomization  
  dk.final <- dk1 <- dk2 <- res$d.alloc    
  toxk <- res$tox.safe                             
  n1 <- res$n1
  nmore <- N - n1                                                # nmore = max sample size - pts. used in stage 1                                                                                      
  nd <- length(unique(dk.safe))                  
  rp <- NULL
  stop <- 0                                           
  
  if (nd == 0) {                                               # If no accept. doses after stage 1, print allocation, no stage 2                   
    yk.final <- yk.final
    dk.final <- dk.final
    stop <- 1
  }
  
  if (nd == 1) {                                               # If only dose 1 safe, allocate up to 9 pts., no stage 2                
    extra <- stop.rule - length(dk.safe)
    ab <- beta.ab(m[1]/100, v[1])
    y.extra <- 100*stats::rbinom(extra, nbb, stats::rbeta(1, ab$a, ab$b) ) / nbb
    yk.final <- c(yk.final, y.extra)                          
    dk.final <- c(dk.final, rep(1,extra))     
    stop <- 1    
  } 
  
  if (nd > 1) {                                               # If 2 or more doses are accept. after stage 1, enter stage 2
    
    coh.toxk <- cbind(matrix(dk.safe, ncol = coh.size, byrow = T)[,1], toxk) # Matrix of safe dose assign. and tox. to be used for LR
    
    for (k in 1:nmore) {
      
      if (stop == 0) {                                        # As long as there are 2 or more doses in randomization
        
        reg <- stats::lm(log(yk.safe + 1) ~ factor(dk.safe))         # Linear model with log(Y) for accept. doses 
        fit <- as.vector(reg$fitted.values)                   # Fitted values for Y
        dose.unique <- duplicated(dk.safe)
        fitp <- exp(fit) 
        fitp <- fitp[dose.unique == F]
        #fitp <- ifelse(fitp > 100, 100, fitp)                 # Restrict values - %persistence can only be b/w 0 and 1
        #fitp <- ifelse(fitp < 0, 0, fitp)                              
        rp <- fitp/sum(fitp)                                  # Calculate randomization prob. for each dose
        rp <- ifelse(rp < 0.02, 0.02, rp)                   
        dj <- stats::rmultinom(1, 1, prob = rp)                      # New (next) dose assign.
        
        if (samedose) {
          dj <- rep((1:length(dj))[dj == 1], cohort)
        } else {
          dosemat <- as.vector(dj*matrix(1:nd, ncol = cohort, nrow = nd))
          dj <- dosemat[dosemat > 0]
        } 
        ab <- beta.ab(m[dj]/100, v[dj])
        p <- stats::rbeta(1, ab$a, ab$b)
        yj <- 100*stats::rbinom(1, nbb, p)/nbb                      # New Y value
        toxj <- stats::rbinom(1, size = 1, dose.tox[dj])            # New toxicity for the next patient
        
        coh.toxj <- c(dj, toxj)                              # New dose and new tox.  
        yk.safe  <- c(yk.safe, yj)
        yk.final <- c(yk.final,yj)                                    
        dk.safe  <- c(dk.safe, dj)
        dk.final <- c(dk.final,dj)  
        
        coh.toxk <- rbind(coh.toxk, coh.toxj)       
        toxk <- c(toxk,toxj)
        n.obsk <- table(dk.safe)
        
        # If no toxicities observed, keep going, else calculate the LR and establish safety
        
        if (toxj == 0) {
          dk.safe <- dk.safe
          yk.safe <- yk.safe         
          
        } else { 
          
          # Create a table with observed toxicities and total n for computing the LR
          
          LR.table.temp <- table(coh.toxk[,1], coh.toxk[,2])
          
          if (ncol(LR.table.temp) == 2) { 
            LR.table <- cbind(LR.table.temp[,2], n.obsk)                         
          }else {                        
            LR.table <- cbind(LR.table.temp[,1], n.obsk)
          }
          loglik.p2 <- NULL
          loglik.p1 <- NULL
          lik.diff <- NULL
          accept.dose <- NULL
          
          for (j in 1:nrow(LR.table)) {
            
            loglik.p2[j] <- LR.table[j, 1]*log(p2) + (LR.table[j, 2] - LR.table[j, 1])*log(1 - p2)          
            loglik.p1[j] <- LR.table[j, 1]*log(p1) + (LR.table[j, 2] - LR.table[j, 1])*log(1 - p1)          
            lik.diff[j] <- exp(loglik.p2[j] - loglik.p1[j])              
            accept.dose[j] <- ifelse(lik.diff[j] > (1/K), 1, 0)
          }              
          dk.safe[dk.safe >= which(accept.dose == 0)] <- NA           # Discard the non-safe doses and all above it by putting NAs
          
          new.model <- cbind(dk.safe,yk.safe)
          new.model <- stats::na.omit(new.model)
          dk.safe <- new.model[, 1]                              
          yk.safe <- new.model[, 2]                               
          yk.final <- yk.final
          dk.final <- dk.final
          
          coh.toxk <- coh.toxk[!apply(coh.toxk, 1, function(x){any(x >= which(accept.dose == 0))}), ]  # New cohort and tox. vector
          
        }#else LR
        
        if (length(unique(dk.safe)) > 1) {                              # continue rand. if more than 2 doses 
          
          dk.safe <- dk.safe
          yk.safe <- yk.safe
          dk.final <- dk.final
          yk.final <- yk.final
        }
        
        if (length(unique(dk.safe)) == 1) {                              # if only dose 1 left    
          new.size <- nmore + length(dk2)                                # dk2 - dose assign. from stage 1
          length.dk1 <- length(dk.final)                                 # length(dk1) - dose assign. from stage 1&2
          
          if ((length(dk.safe) < stop.rule) && (length.dk1 < new.size)) {     # if the max. sample size was not reached and less than 9 subj. at dose 1                                  
            extra.one <- min(new.size - length.dk1, stop.rule - length(dk.safe))                     
            ab <- beta.ab(m[1]/100, v[1])
            yj.one <- 100*stats::rbinom(extra.one, nbb, stats::rbeta(1, ab$a, ab$b) ) / nbb
            yk.final <- c(yk.final, yj.one)
            dk.final <- c(dk.final, rep(1, extra.one))         
            stop <- 1   
            
          } else {  
            
            dk.final <- dk.final
            yk.final <- yk.final
            stop <- 1
          }
        }
        if (length(unique(dk.safe)) < 1) {                                  # stop if no dose left
          dk.final <- dk.final
          yk.final <- yk.final
          stop <- 1
        }                         
      } else {
        dk.final <- dk.final
        yk.final <- yk.final
      }
    }# end for
  }# end(if nd>1)
  
  return(list(Y.final = yk.final, 
              d.final = dk.final, 
              n1 = n1))
}