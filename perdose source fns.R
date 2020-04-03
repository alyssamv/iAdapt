######################################################
############### Functions ############################
######################################################


#Function to obtain average toxicity scores at each dose level, which is the weighted sum of toxicity weights and the probability 
# of toxicity at each dose level. This computes equation #1.23 in the manuscript. 
get.thresh <- function(ntox, W, TOX){
  
  W = W / sum(W[,4])
  
  tox_type <- seq(from = 0, to = 4,by = 1)
  tox.type <- matrix(rep(tox_type, each = ntox), ncol = ntox, byrow = TRUE)
  tox.type <- as.data.frame(tox.type)
  possible_outcomes <- expand.grid(tox.type[, 1:ntox]) #Permutes all possible toxicity profiles into a data set
  #map toxicity profiles to weights
  mapped.weight <- NA
  v <- NULL
  
  for (i in 1:ntox) {
    for (k in 1:nrow(possible_outcomes)) {
      mapped.weight[k] <- ifelse(possible_outcomes[k, i] > 0, W[i, possible_outcomes[k, i]], 0)
    }
    v <- cbind(v, mapped.weight)
  }
  
  mapped_data <- matrix(v, ncol = 3, byrow = FALSE)
  scores <- apply(mapped_data, 1, sum) #Calulate toxicity scores for each profile in data set
  prob.data <- array(NA, c(nrow(possible_outcomes), ntox, d))
  
  for (j in 1:d) {
    for (i in 1:nrow(possible_outcomes)) {
      for (k in 1:ntox) {
        tox.prob <- ifelse(possible_outcomes[i, k] > 0, TOX[j, possible_outcomes[i, k] + 1, k], TOX[j, 1, k])
        prob.data[i, k, j] <- tox.prob
      }
    }
  }
  
  prob.tox <- apply(prob.data[,,], c(1, 3), prod)
  prob.scores <- scores*prob.tox
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


# Functions to get hypotheses
# Function to write th efunction that calculates all possible nTTP from weight matrix W
# grade.thresh is a vector of the grade at which each toxicity type qualifies as a DLT
nTTP.calc.fn <- function(ntox, W, grade.thresh){
  hyp.plot = "function(){
  W = cbind(rep(0, ntox), W)
  thetamax = sum(W[,ncol(W)])
  tt = list(h2 = c(), h1 = c()) \n
  "
  forloop = ""
  cond = "("
  summ = "sum("
  vars = letters[(1:ntox) + 8]
  for (xx in 1:ntox) {
    
    forloop = paste0(forloop, "for (", vars[xx], " in 1:5) {\n")
    
    if (xx %in% c(1, ntox - 1)) {
      cond = paste0(cond, vars[xx], " >= ", grade.thresh[xx], " || ")
      summ = paste0(summ, "W[", xx, ", ", vars[xx], "]^2, ")
    } else {
      cond = paste0(cond, vars[xx], " >= ", grade.thresh[xx], ")")
      summ = paste0(summ, "W[", xx, ", ", vars[xx], "]^2)")
    }
    
  }
  
  closure = paste(rep("}", ntox), collapse = "")
  code = paste0(forloop,
                "if ", cond, "{\n",
                'tt[["h1"]] = append(tt[["h1"]], 
                sqrt(', summ, ")/thetamax)\n",
                "} else {",
                'tt[["h2"]] = append(tt[["h2"]], 
                sqrt(', summ, ")/thetamax)\n",
                "}", closure)
  
  fn.text = paste0(hyp.plot, code, "\nreturn(tt)}")
  
  fn = eval(parse(text = fn.text))
  
  return(fn)
}

# Function that calls nTTP.calc.fn, plots distributions, and gives suggested hypothesis values
# grade.thresh is a vector of the grade at which each toxicity type qualifies as a DLT
nTTP.hypothesis <- function(ntox, W, grade.thresh){
  
  ### Calculate all possible nTTP values ###
  ### DLT defined as grade >= 3 in any toxicity category ###
  nTTP.calc = eval(parse(text = "nTTP.calc.fn(ntox, W, grade.thresh)"))
  tt = nTTP.calc()
  
  ## Plot distributions
  par(mfrow = c(2, 1))
  hist(tt[["h2"]], xlim = c(0, 0.8), main = "No DLTs")
  abline(v = mean(tt[["h2"]]), lty = 3, lwd = 2)
  hist(tt[["h1"]], xlim = c(0, 0.8), main = ">0 DLTs")
  abline(v = mean(tt[["h1"]]), lty = 3, lwd = 2)
  
  
  ## Suggested hypotheses
  sug = lapply(tt, mean) # hypothesis values
  
  return(list(suggested.H = sug, all.nTTP = tt))
}
