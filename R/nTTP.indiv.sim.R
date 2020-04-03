# W weight matrix of toxicity burden for each type of toxicity
# TOX list of matrices giving probability of observing given toxicity grade for each toxicity type
# ntox number of toxicity types
# dose dose to simulate nTTP for

nTTP.indiv.sim <- function(W, TOX, ntox, dose){
  Tox <- NA
  random.tox <- runif(ntox) #draw random number from uniform(0,1)
  comb.curr <- dose #Define current dose level
  
  W.1 <- W / sum(W[,4]) # normalize toxicity burden matrix (equivalent of tox.weights from Ezzalfani)
  
  for (k in 1:ntox) {
    Tox[k] <- ifelse(random.tox[k] < TOX[comb.curr, 1, k] || TOX[comb.curr,1,k] == 0, 0, 
                     ifelse(random.tox[k] < TOX[comb.curr, 2, k] + TOX[comb.curr, 1, k], 1, 
                            ifelse(random.tox[k] < TOX[comb.curr, 2, k] + TOX[comb.curr, 1, k] + TOX[comb.curr, 3, k], 2,
                                   ifelse(random.tox[k] < TOX[comb.curr, 1, k] + TOX[comb.curr, 2, k] + TOX[comb.curr, 3, k] + TOX[comb.curr, 4, k], 3, 4))))
  }
  #Tox #vector of observed toxicity's grades for each toxicity given the current dose
  toxscores <- NA
  for (k in 1:ntox) {
    toxscores[k] <- ifelse(max(Tox[k]) == 0, 0, W.1[k, max(Tox[k])]) 
  }
  
  nTTP <- sum(toxscores) #The observed toxicity score the patient 
  return(nTTP)
}
