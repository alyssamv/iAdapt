#Function to obtain average toxicity scores at each dose level, which is the weighted sum of toxicity weights and the probability 
# of toxicity at each dose level. This computes equation #1.23 in the manuscript. 
get.thresh <- function(ntox, W, TOX){
  
  thetamax = sum(W[, 4]^2)
  
  tox_type <- seq(from = 0, to = 4, by = 1)
  tox.type <- matrix(rep(tox_type, each = ntox), ncol = ntox, byrow = TRUE)
  tox.type <- as.data.frame(tox.type)
  possible_outcomes <- expand.grid(tox.type[, 1:ntox]) #Permutes all possible AE grades into a data set
  
  # weights corresponding to AE grades
  mapped.weight <- NA
  vec <- NULL
  for (i in 1:ntox) {
    for (k in 1:nrow(possible_outcomes)) {
      mapped.weight[k] <- ifelse(possible_outcomes[k, i] > 0, W[i, possible_outcomes[k, i]], 0)
    }
    vec <- cbind(vec, mapped.weight)
  }
  
  
  mapped_data <- matrix(vec, ncol = 3, byrow = FALSE)
  scores <- apply(mapped_data^2, 1, sum) # Calculate toxicity scores for each profile in data set
  normalized_scores <- sqrt(scores/thetamax) # Calculate nTTP
  
  prob.data <- array(NA, c(nrow(possible_outcomes), ntox, d))
  
  # probabilities corresponding to AE grades
  for (j in 1:d) {
    for (i in 1:nrow(possible_outcomes)) {
      for (k in 1:ntox) {
        tox.prob <- TOX[j, possible_outcomes[i, k] + 1, k] #ifelse(possible_outcomes[i, k] > 0, TOX[j, possible_outcomes[i, k] + 1, k], TOX[j, 1, k])
        prob.data[i, k, j] <- tox.prob
      }
    }
  }
  
  prob.tox <- apply(prob.data[,,], c(1, 3), prod)
  prob.scores <- normalized_scores*prob.tox
  thresh <- apply(prob.scores, 2, sum)
  return(thresh)
}


#Function to caculate the probabiity of a DLT at each dose level using 
#the individual toxicity probabilities across grades and doses
# grade.thresh is a vector of the grade at which each toxicity type qualifies as a DLT
dlt.prob <- function(TOX, ntox, grade.thresh){ 
  
  cp <- function(p) 
  {
    ev <- do.call(expand.grid, replicate(length(p), 0:1, simplify = FALSE)) # all permutations (rows) of DLT for each tox type
    pe <- apply(ev, 1, function(x) prod(p*(x == 1) + (1 - p)*(x == 0))) # Law of total prob ; DLT occurs with probability p 
    tapply(pe,rowSums(ev), sum) # sums probabilities for each group (total # DLTs across tox types)
  }
  
  ## probability of DLT for each dose at each tox type
  dlt.matrix <- matrix(0, nrow = nrow(TOX), ncol = d)
  for (i in 1:ntox) {
    xxx <- TOX[, -(1:grade.thresh[i]), i] # for the dose and tox type, give only the cols corresponding to the grades that qualify as DLT
    xxx <- cbind(xxx, 0) # ensures that xxx is a matrix (instead of vector)
    dlt.matrix[i, ] <- apply(xxx, 1, sum) # why sum?
  }
  ## probability of DLT for each dose across tox types
  ptox <- NULL
  for (i in 1:d) {
    ptox[i] <- 1 - cp(dlt.matrix[1:ntox, i])[1] #Ptox is probability of DLT at each dose level
  } 
  
  return(ptox)
}


# Function to write the function that calculates all possible nTTP from weight matrix W
# grade.thresh is a vector of the grades at which each toxicity type qualifies as a DLT

nTTP.calc.fn <- function(ntox, W, grade.thresh){
  hyp.plot = "function(){
  W = cbind(rep(0, ntox), W)
  thetamax = sum(W[, ncol(W)]^2)
  tt = list(h2 = c(), h1 = c()) \n
  "
  forloop = ""
  cond = "("
  summ = "sum("
  vars = letters[(1:ntox) + 8]
  for (xx in 1:ntox) {
    
    forloop = paste0(forloop, "for (", vars[xx], " in 1:5) {\n")
    
    if (xx %in% c(1, ntox - 1)) {
      cond = paste0(cond, vars[xx], " >= ", grade.thresh[xx] + 1, " || ") # add 1 to grade thresh bc first grade is 0, but we have for (1:5)
      summ = paste0(summ, "W[", xx, ", ", vars[xx], "]^2, ")
    } else {
      cond = paste0(cond, vars[xx], " >= ", grade.thresh[xx] + 1, ")")
      summ = paste0(summ, "W[", xx, ", ", vars[xx], "]^2)")
    }
    
  }
  
  closure = paste(rep("}", ntox), collapse = "")
  code = paste0(forloop,
                "if ", cond, "{\n",
                'tt[["h1"]] = append(tt[["h1"]], 
                sqrt(', summ, "/thetamax))\n",
                "} else {",
                'tt[["h2"]] = append(tt[["h2"]], 
                sqrt(', summ, "/thetamax))\n",
                "}", closure)
  
  fn.text = paste0(hyp.plot, code, "\nreturn(tt)}")
  
  fn = eval(parse(text = fn.text))
  
  return(fn)
}

# # Function that calls nTTP.calc.fn, plots distributions, and gives suggested hypothesis values
# # grade.thresh is a vector of the grade at which each toxicity type qualifies as a DLT
# nTTP.hypothesis <- function(ntox, W, grade.thresh){
#   
#   ### Calculate all possible nTTP values ###
#   ### DLT defined as grade >= 3 in any toxicity category ###
#   nTTP.calc = eval(parse(text = "nTTP.calc.fn(ntox, W, grade.thresh)"))
#   tt = nTTP.calc()
#   
#   ## Plot distributions
#   par(mfrow = c(2, 1))
#   hist(tt[["h2"]], xlim = c(0, 1), main = "No DLTs")
#   abline(v = mean(tt[["h2"]]), lty = 3, lwd = 2)
#   hist(tt[["h1"]], xlim = c(0, 1), main = ">0 DLTs")
#   abline(v = mean(tt[["h1"]]), lty = 3, lwd = 2)
#   
#   
#   ## Suggested hypotheses
#   sug = lapply(tt, mean) # hypothesis values
#   
#   return(list(suggested.H = sug, all.nTTP = tt))
# }

#################################################################################
################################## STAGE 1 ######################################
#################################################################################

nTTP.indiv.sim <- function(W, TOX, ntox, dose){

  # Simulate grade observed for each toxicity type, based on TOX probabilities
  Tox <- NA
  for (k in 1:ntox) {
    # Tox vector of observed toxicity's grades for each toxicity given the current dose
    Tox[k] = sample(0:4, size = 1, replace = TRUE, prob = TOX[dose, , k])
  }

  # Corresponding weights
  toxscores <- NA
  for (k in 1:ntox) {
    # ifelse statement bc W does not have zero weights for grade 0
    toxscores[k] <- ifelse(Tox[k] == 0, 
                           0, # weight for grade 0
                           W[k, Tox[k]]) # weight of AE given grade, tox type
  }
  
  # nTTP calculation
  thetamax = sum(W[, 4]^2)
  nTTP <- sqrt(sum(toxscores^2) / thetamax) # The observed nTTP for the patient 
  
  return(nTTP)
}


# Escalation scheme for stage 1, based on observed nTTP values at each dose (calculate LR)
tox.profile.nTTP <- function(dose, p1, p2, K, coh.size, ntox, W, TOX, sigma = 0.15){ 
  
  dose   <- c(1:dose) # vector of counts up to number of doses given
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  nttp  <- NULL
  
  # bounds for nTTP (truncated normal distribution)
  a = 0
  b = 1
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1                                        # current cohort corresponding to dose
    
    # nTTPs for all patients (size coh.size) on dose i
    coh.nttp  <- replicate(coh.size, nTTP.indiv.sim(W = W, 
                                                    TOX = TOX, 
                                                    ntox = ntox, 
                                                    dose = dose[i]))	# nTTPs for that dose based on tox prob
    
    # Calculate LR
    l.p2   <- prod(sapply(coh.nttp, FUN = function(i){ dnorm((i - p2)/sigma) })) / 
      (sigma*(pnorm((b - p2)/sigma) - pnorm((a - p2)/sigma)))^coh.size # likelihood of acceptable/alternative hypothesis 
    l.p1   <- prod(sapply(coh.nttp, FUN = function(i){ dnorm((i - p1)/sigma) })) / 
      (sigma*(pnorm((b - p1)/sigma) - pnorm((a - p1)/sigma)))^coh.size # likelihood of unacceptable/null hypothesis
    LR     <- round(l.p2/l.p1, 2)    
    
    x <- c(x, dose[i], mean(coh.nttp), cohort, LR)                              			
    
    # list of observed nTTP
    nttp = append(nttp, coh.nttp) 
    
    if (LR <= (1/K)) {       # stop escalation
      stop <- 1
    } else if (LR > (1/K)) { # escalate to next dose i + 1
      i <- i + 1	 
    }
    
  }       
  return(list(mnTTP = matrix(x, ncol = 4, byrow = TRUE), 
              all_nTTP = nttp))
} 



safe.dose.nTTP <- function(dose, p1, p2, K, coh.size, W, TOX, ntox, sigma = 0.15) {
  
  # simulate 
  tox <- tox.profile.nTTP(dose = dose, p1 = p1, p2 = p2, K = K, 
                          coh.size = coh.size, ntox = ntox, W = W, TOX = TOX, sigma = sigma)
  res <- tox$mnTTP
  all_nttp <- tox$all_nTTP # all observed nTTP in stage 1
  
  alloc.total <- sort(rep(res[, 1], coh.size)) # dose allocation for all patients (stage 1)
  n1 <- nrow(res) * coh.size # sample size in stage 1
  
  unsafe.dose <- which(res[, 4] <= (1/K))
  if (length(unsafe.dose) == 0) {
    alloc.safe <- res[, 1:2]
  } else {
    alloc.safe <- res[res[, 1] != unsafe.dose, 1:2]
  }
  return(list(alloc.safe = alloc.safe, alloc.total = alloc.total, 
              n1 = n1, all_nttp = all_nttp))
}


eff.stg1.nTTP <- function(dose, p1, p2, K, coh.size, m, v, nbb = 100, W, TOX, ntox, sigma = 0.15) {
  
  res <- safe.dose.nTTP(dose, p1, p2, K, coh.size, W, TOX, ntox, sigma = sigma)
  d.alloc <- res$alloc.total
  val.safe <- res$alloc.safe
  Y.safe <- d.safe <- tox.safe <- Y.alloc <- NULL
  n1 <- res$n1
  all_nttp <- res$all_nttp
  
  # Simulate efficacy outcomes for all doses that enrolled patients
  for (i in 1:length(d.alloc)) {
    ab <- beta.ab(m[d.alloc[i]]/100, v[d.alloc[i]])
    p <- stats::rbeta(1, ab$a, ab$b)
    Y.alloc[i] <- 100 * stats::rbinom(1, nbb, p)/nbb
  }
  
  if (length(val.safe) > 2) { # if more than one dose found safe
    d.safe <- sort(rep(val.safe[, 1], coh.size))
    tox.safe <- res$alloc.safe[, 2]
    Y.safe <- Y.alloc[1:length(d.safe)]
  } else if (length(val.safe) == 2) { # if exactly one dose (the first) found safe
    d.safe <- sort(rep(val.safe[1], coh.size))
    tox.safe <- res$alloc.safe[2]
    Y.safe <- Y.alloc[1:length(d.safe)]
  } else { # if no doses found safe (first dose is too toxic)
    Y.safe <- d.safe <- NULL
    tox.safe <- res$alloc.safe[, 2]
  }
  
  return(list(Y.safe = Y.safe, d.safe = d.safe, tox.safe = tox.safe, 
              n1 = n1, Y.alloc = Y.alloc, d.alloc = d.alloc, all_nttp = all_nttp))
}

#################################################################################
################################# STAGE 1+2 #####################################
#################################################################################

rand.stg2.nTTP <- function(dose, p1, p2, K, coh.size, m, v, N, stop.rule = 9, 
                           cohort = 1, samedose = TRUE, nbb = 100, W, TOX, ntox, sigma = 0.15) {
  
  res <- eff.stg1.nTTP(dose, p1, p2, K, coh.size, m, v, nbb = 100, 
                       W = W, TOX = TOX, ntox = ntox, sigma = sigma) # stage 1
  dose   <- c(1:dose) 
  yk.safe <- res$Y.safe # efficacy of all safe doses
  yk.final <- res$Y.alloc # efficacy of all patients enrolled in stage 1
  dk.safe <- res$d.safe # dose allocation for all safe doses
  dk.final <- dk1 <- dk2 <- res$d.alloc # dose allocation of all patients enrolled in stage 1
  toxk <- res$tox.safe # mnTTP at each dose level in stage 1
  n1 <- res$n1 # number of patients enrolled in stage 1
  nmore <- N - n1 # number of subjects left to randomize
  all_nttp <- res$all_nttp
  nd <- length(unique(dk.safe)) # number of safe doses
  rp <- NULL
  stop <- 0
  
  if (nd == 0) { # if no doses found safe
    yk.final <- yk.final
    dk.final <- dk.final
    stop <- 1
  }
  if (nd == 1) { # if exactly one dose found safe
    extra <- stop.rule - length(dk.safe)
    ab <- beta.ab(m[1]/100, v[1])
    y.extra <- 100 * stats::rbinom(extra, nbb, stats::rbeta(1, ab$a, ab$b))/nbb
    yk.final <- c(yk.final, y.extra)
    dk.final <- c(dk.final, rep(1, extra))
    stop <- 1
  }
  if (nd > 1) { # if more than one dose found safe
    coh.toxk <- cbind(matrix(dk.safe, ncol = coh.size, byrow = TRUE)[, 1], toxk)
    for (k in 1:nmore) {
      if (stop == 0) {
        reg <- stats::lm(log(yk.safe + 1) ~ factor(dk.safe))
        fit <- as.vector(reg$fitted.values)
        dose.unique <- duplicated(dk.safe)
        fitp <- exp(fit)
        fitp <- fitp[dose.unique == FALSE]
        rp <- fitp/sum(fitp)
        rp <- ifelse(rp < 0.02, 0.02, rp) #################################### why?
        dj <- stats::rmultinom(1, 1, prob = rp)
        if (samedose == TRUE) {
          dj <- rep((1:length(dj))[dj == 1], cohort)
        } else {
          dosemat <- as.vector(dj * matrix(1:nd, ncol = cohort, 
                                           nrow = nd))
          dj <- dosemat[dosemat > 0]
        }
        ab <- beta.ab(m[dj]/100, v[dj])
        p <- stats::rbeta(1, ab$a, ab$b)
        yj <- 100 * stats::rbinom(1, nbb, p)/nbb
        
        # simulate nTTP for next individual
        toxj <- nTTP.indiv.sim(W = W, 
                               TOX = TOX, 
                               ntox = ntox, 
                               dose = dose[dj]) 
        
        coh.toxj <- c(dj, toxj)
        yk.safe <- c(yk.safe, yj)
        yk.final <- c(yk.final, yj)
        dk.safe <- c(dk.safe, dj)
        dk.final <- c(dk.final, dj)
        all_nttp <- append(all_nttp, toxj)
        coh.toxk <- data.frame(dose = dk.final, 
                               nttp = all_nttp) 
        
        # calculate new LR for all patients on the given dose
        l.p2 <- NULL
        l.p1 <- NULL
        LR <- NULL
        accept.dose <- NULL
        a = 0
        b = 1
        for (j in 1:max(coh.toxk[, 1])) {
          
          rows = which(coh.toxk[,1] == j)
          nttp = coh.toxk[rows, 2]
          nj = length(rows)
          
          # likelihood of acceptable/alternative hypothesis 
          l.p2[j]   <- prod(sapply(nttp, FUN = function(i){ dnorm((i - p2)/sigma) })) / 
            (sigma*(pnorm((b - p2)/sigma) - pnorm((a - p2)/sigma)))^nj 
          # likelihood of unacceptable/null hypothesis
          l.p1[j]   <- prod(sapply(nttp, FUN = function(i){ dnorm((i - p1)/sigma) })) / 
            (sigma*(pnorm((b - p1)/sigma) - pnorm((a - p1)/sigma)))^nj 
          # likelihood ratio
          LR[j]     <- round(l.p2[j]/l.p1[j], 2)
          
          
          accept.dose[j] <- ifelse(LR[j] > (1/K), 1, 0)
        }
        
        unsafe_dose = which(accept.dose == 0)
        if (length(unsafe_dose) > 0) {
          smallest_unsafe_dose = unsafe_dose[1]
          
          dk.safe[dk.safe >= smallest_unsafe_dose] <- NA
          coh.toxk <- coh.toxk[!apply(coh.toxk, 1, function(x) {
            any(x >= smallest_unsafe_dose)
          }), ]
        }
        new.model <- cbind(dk.safe, yk.safe)
        new.model <- stats::na.omit(new.model)
        dk.safe <- new.model[, 1]
        yk.safe <- new.model[, 2]
        yk.final <- yk.final
        dk.final <- dk.final
        
        
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
        } else if (length(unique(dk.safe)) < 1) {
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



#################################################################################
################################## WRAPPER ######################################
#################################################################################


sim.trials.nTTP <- function (numsims, dose, p1, p2, K, coh.size, m, v, 
                             N, stop.rule = 9, cohort = 1, samedose = TRUE, nbb = 100, 
                             W, TOX, ntox, sigma = 0.15) 
{
  sim.yk <- sim.dk <- matrix(NA, nrow = numsims, ncol = N)
  sim.doses <- matrix(NA, nrow = numsims, ncol = dose)
  for (i in 1:numsims) {
    fstudy.out <- rand.stg2.nTTP(dose = dose, p1 = p1, p2 = p2, K = K, coh.size = coh.size, 
                                 m = m, v = v, N = N, stop.rule = stop.rule, 
                                 cohort = cohort, samedose = samedose, nbb = nbb, W = W, 
                                 TOX = TOX, ntox = ntox, sigma = sigma)
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




sim.summary <- function(sims, print = TRUE){
  sim.doses = sims$sim.d
  n.doses = max(sim.doses, na.rm = TRUE)
  sim.eff = sims$sim.Y
  dose.mat.a <- matrix(NA, nrow(sim.doses), n.doses)
  for (i in 1:nrow(sim.doses)) {
    dose.no.na <- na.omit(sim.doses[i, ])
    dose.mat.a[i, ] <- table(factor(dose.no.na, levels = 1:n.doses))/length(dose.no.na)
  }
  est.dose1 <- matrix(NA, n.doses, 5)
  
  for (j in 1:n.doses) {
    est.dose1[j, ] <- c(j/100, 
                        quantile(dose.mat.a[, j], 
                                 prob = c(0.25, 0.5, 0.75), na.rm = TRUE),
                        round(mean(dose.mat.a[, j], na.rm = TRUE), 2))
  }
  dose.IQR = round(est.dose1 * 100, 1)
  
  pers.hat.a <- matrix(NA, nrow(sim.eff), n.doses + 1)
  for (i in 1:nrow(sim.eff)) {
    for (j in 1:(n.doses + 1)) {
      pers.hat.a[i, j] <- (median(sim.eff[i, sim.doses[i, 
                                                       ] == j - 1]))
    }
  }
  est.pers1 <- matrix(NA, (n.doses + 1), 5)
  for (j in 1:(n.doses + 1)) {
    est.pers1[j, ] <- c((j - 1), 
                        quantile(pers.hat.a[, j], 
                                 prob = c(0.25, 0.5, 0.75), na.rm = TRUE),
                        round(mean(pers.hat.a[, j], na.rm = TRUE), 2))
  }
  Y = est.pers1[-1, ]
  
  if (print == TRUE) {
    print(knitr::kable(dose.IQR, caption = "Percent allocation per dose level", 
                       col.names = c("Dose", "25th percentile", "Median", 
                                     "75th percentile", "Mean")))
    print(knitr::kable(Y, caption = "Estimated efficacy per dose level", 
                       col.names = c("Dose", "25th percentile", "Median", 
                                     "75th percentile", "Mean")))
  }
  
  return(list(pct.treated = dose.IQR, efficacy = Y))
}
