#################################################################################################
#	              R code for implementing the adaptive two-stage early- phase design	            #
#	                 for continuous immunologic outcomes with safety constraints	                #
#################################################################################################


#############################################################################################
# Stage 1 - establish the safety profiles of the predefined doses	                          #
#	                                                                                          #
# Doses are classified as acceptable (safe) or unacceptable (unsafe).	                      #
# Based on observed toxicities/dose, calculate the likelihood-ratio(LR) and compare it to   #
# a threshold K.	                                                                          #
# If LR > 1/K => dose is acceptable, escalate to the next dose                            	#
# If LR <= 1/K => dose is unacceptable, stop allocation	                                    #
#############################################################################################


#############################################################################################
# Function tox.profile() generates dlts and calculates the LR/dose
# Returns: a 4-column matrix containing dose assignment, DLTs at each dose, cohort number, and LR

# dose number of doses to be tested (scalar)
# dose.tox vector of true toxicities for each dose. Values range from 0 - 1.
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 


tox.profile <- function(dose, dose.tox, p1, p2, K, coh.size){ 
  
  dose   <- c(1:dose)                                           # vector of counts up to number of doses given
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1                                        # current cohort corresponding to dose
    dltsi  <- rbinom(1, coh.size, dose.tox[i])	# number of DLTs for that dose based on tox prob
    
    l.p2   <- dltsi * log(p2) + (coh.size - dltsi) * log(1 - p2) # likelihood of acceptable/alternative hypothesis 
    l.p1   <- dltsi * log(p1) + (coh.size - dltsi) * log(1 - p1) # likelihood of unacceptable/null hypothesis
    LR     <- round(exp(l.p2 - l.p1), 2)    
    
    x <- c(x, dose[i], dltsi, cohort, LR)                              			
    
    if (LR <= (1/K))                        # stop escalation
      stop <- 1
    
    if (LR > (1/K))                         # escalate to next dose
      i <- i + 1	 
    
  }       
  return(matrix(x, ncol = 4, byrow = T))
} 


#####################################################################################################
#  Function safe.dose() distinguishes acceptable from unacceptable doses
#  Returns: alloc.safe - matrix of assignments for only acceptable doses (to be used in stage 2) and the
#	                       corresponding toxicities
#          alloc.total - vector of all dose assignments from stage 1 
#                   n1 - total number of subjects allocated in stage 1

# dose  number of doses to be tested (scalar)
# dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 


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


################################################################################################
#  Function beta.ab() returns parameters alpha and beta for generating beta r.v. (per dose)

# m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
# v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)


beta.ab <- function(m, v) {
  
  a <- seq(0.5, 20, 0.01)                            # a is a seq of alpha in beta distr.
  b <- a * (1 - m) / m
  
  vfit  <- a * b / ((a + b + 1) * (a + b)^2)
  diff  <- abs(vfit - v)
  index <- (1:length(diff))[diff == min(diff)]       # return the index of the min var.
  
  return(list(a = a[index],
              b = b[index]))                         # return alpha and beta for the min.var.
}



#################################################################################################
#  Function gen.eff.stg1() uses a beta-binomial distribution to generate outcomes (Ys) corresponding
#  to all dose assignments in stage 1

#  Returns: Y.safe - vector of efficacy outcomes corresponding to only acceptable doses; to be used for 
#'                   randomization in stage 2
#'          d.safe - vector of acceptable dose assignments; to be used for 
#'                   randomization in stage 2
#'       tox. safe - vector of toxicities corresponding to acceptable doses
#'              n1 - total number of subjects allocated in stage 1
#'         Y.alloc - vector of efficacy outcomes corresponding to all dose assignments from stage 1 
#'                   (acceptable and unacceptable)
#'         d.alloc - vector of all dose assignments from stage 1         
#'                 
# dose  number of doses to be tested (scalar)
# dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 
# m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
# v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
# nbb  binomial parameter (default = 100 cells per patient)

gen.eff.stg1 <- function(dose, dose.tox, p1, p2, K, coh.size, m, v, nbb = 100) {
  
  res <- safe.dose(dose, dose.tox, p1, p2, K, coh.size)
  d.alloc <- res$alloc.total
  val.safe <- res$alloc.safe
  
  Y.safe <- d.safe <- tox.safe <- Y.alloc <- NULL      
  n1 <- res$n1  
  
  for (i in 1:length(d.alloc)) {                            # Generate Ys for all allocations in stage 1
    ab <- beta.ab(m[d.alloc[i]]/100, v[d.alloc[i]])
    p <- rbeta(1, ab$a, ab$b)
    Y.alloc[i] <- 100*rbinom(1, nbb, p) / nbb                  
  }                          
  
  if (length(val.safe) > 2) {                               # if 2 or more acceptable doses
    d.safe <- sort(rep(val.safe[, 1], coh.size))
    tox.safe <- res$alloc.safe[, 2]
    Y.safe <- Y.alloc[1:length(d.safe)]
    
  } else if (length(val.safe) == 2) {                   # if only dose 1 acceptable
    d.safe <- sort(rep(val.safe[1], coh.size))
    tox.safe <- res$alloc.safe[2]
    Y.safe <- Y.alloc[1:length(d.safe)]
    
  } else {
    Y.safe <- d.safe  <- NULL                           # if dose 1 is unacceptable, leave it null
    tox.safe <- res$alloc.safe[,2]
  }         
  return(list(Y.safe = Y.safe,
              d.safe = d.safe,
              tox.safe = tox.safe,
              n1 = n1,
              Y.alloc = Y.alloc,
              d.alloc = d.alloc))  
}



#############################################################################################
#	                          Stage 2 - adaptive randomization	                              #
#	                                                                                          #
# If 2 or more doses are considered acceptable after stage 1, the remaining patients are    #
# randomized to these open doses until the total sample size N is reached.	                #
# If only dose 1 is acceptable after stage 1, allocate up to 9 patients (stop.rule=9)	      #
#	                                                                                          #
# If a DLT is observed at a certain dose, calculate the LR, compare it to K, and:	          #
# If LR > 1/K => dose is acceptable, keep in the randomization pool	                        #
# If LR <= 1/K => dose is unacceptable, discard from randomization pool	                    #
#############################################################################################

#############################################################################################
# Function rand.stg2() fits a linear regression for the continuous efficacy outcomes,
# computes the randomization probabilities/dose and allocates the next patient to a dose that
# is considered acceptably safe and has the highest randomization probability. 
# Dose safety is still monitored using LR and doses that become unacceptable are discarded.

# Returns:  Y.final - vector of all efficacy outcomes (Ys) corresponding to dose assignments (Stages 1&2)
#           d.final - vector of all dose assignments(Stage 1&2)
#  If no dose allocation, put NAs in d.final and y.final

# dose  number of doses to be tested (scalar)
# dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 
# m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
# v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
# nbb  binomial parameter (default = 100 cells per patient)
# N  total sample size for stages 1&2
# stop.rule  if only dose 1 safe, allocate up to 9 (default) patients at dose 1 to collect more info


rand.stg2 <- function(dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule=9, cohort=1, samedose=T, nbb=100) {
  
  res <- gen.eff.stg1(dose, dose.tox, p1, p2, K, coh.size, m, v, nbb)
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
    y.extra <- 100*rbinom(extra, nbb, rbeta(1, ab$a, ab$b) ) / nbb
    yk.final <- c(yk.final, y.extra)                          
    dk.final <- c(dk.final, rep(1,extra))     
    stop <- 1    
  } 
  
  if (nd > 1) {                                               # If 2 or more doses are accept. after stage 1, enter stage 2
    
    coh.toxk <- cbind(matrix(dk.safe, ncol = coh.size, byrow = T)[,1], toxk) # Matrix of safe dose assign. and tox. to be used for LR
    
    for (k in 1:nmore) {
      
      if (stop == 0) {                                        # As long as there are 2 or more doses in randomization
        
        reg <- lm(log(yk.safe + 1) ~ factor(dk.safe))         # Linear model with log(Y) for accept. doses 
        fit <- as.vector(reg$fitted.values)                   # Fitted values for Y
        dose.unique <- duplicated(dk.safe)
        fitp <- exp(fit) 
        fitp <- fitp[dose.unique == F]
        #fitp <- ifelse(fitp > 100, 100, fitp)                 # Restrict values - %persistence can only be b/w 0 and 1
        #fitp <- ifelse(fitp < 0, 0, fitp)                              
        rp <- fitp/sum(fitp)                                  # Calculate randomization prob. for each dose
        rp <- ifelse(rp < 0.02, 0.02, rp)                   
        dj <- rmultinom(1, 1, prob = rp)                      # New (next) dose assign.
        
        if (samedose) {
          dj <- rep((1:length(dj))[dj == 1], cohort)
        } else {
          dosemat <- as.vector(dj*matrix(1:nd, ncol = cohort, nrow = nd))
          dj <- dosemat[dosemat > 0]
        } 
        ab <- beta.ab(m[dj]/100, v[dj])
        p <- rbeta(1, ab$a, ab$b)
        yj <- 100*rbinom(1, nbb, p)/nbb                      # New Y value
        toxj <- rbinom(1, size = 1, dose.tox[dj])            # New toxicity for the next patient
        
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
          new.model <- na.omit(new.model)
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
            yj.one <- 100*rbinom(extra.one, nbb, rbeta(1, ab$a, ab$b) ) / nbb
            yk.final <- c(yk.final, yj.one)
            dk.final <- c(dk.final, rep(1, extra.one))         
            stop <- 1   
            
          } else {  
            
            dk.final <- dk.final
            yk.final <- yk.final
            stop <- 1
          }
        }
        if (length(unique(dk.safe)) < 1) {                               # stop if no dose left
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
  
  return(list(Y.final = yk.final, 
              d.final = dk.final, 
              n1 = n1))
}


#############################################################################################
#	                              Repeated Simulations                                        #
#############################################################################################

# Function sim.trials() simulate a full trial (both stages) x times
# Results are displayed in a matrix format, where each row represents one trial
# Returns:   sim.Y - estimated efficacy per each dose assignment 
#            sim.d - dose assignment for each patient in the trial 


# nsims  number of simulated trials
# dose  number of doses to be tested (scalar)
# dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 
# m  vector of mean efficacies per dose. Values range from 0 - 100. (e.g, T cell persistence - values b/w 5 and 80 per cent) 
# v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
# nbb  binomial parameter (default = 100 cells per patient)
# N  total sample size for stages 1&2
# stop.rule  if only dose 1 safe, allocate up to 9 (default) patients at dose 1 to collect more info


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