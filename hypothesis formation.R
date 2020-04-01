
source(file.path(getwd(), "source fns.R"))


nTTP.hypothesis <- function(ntox, W) {
  thetamax = sum(W[,ncol(W)])
  # for (i in 1:ntox) {
  #   W <- tox.weights(W) 
  # }
  
  
  ### Calculate all possible nTTP values ###
  ### DLT defined as grade >= 3 in any toxicity category ###
  tt = list(h2 = c(),
            h1 = c())
  for (i in 1:5) {
    for (j in 1:5) {
      for (k in 1:5) {
        if (i >= 4 || j >= 4 || k >= 4 ){
          tt[["h1"]] = append(tt[["h1"]], 
                               sqrt(sum(W[1, i]^2, W[2, j]^2, W[3, k]^2))/thetamax)
        } else {
          tt[["h2"]] = append(tt[["h2"]], 
                                sqrt(sum(W[1, i]^2, W[2, j]^2, W[3, k]^2))/thetamax)
        }
      }
    }
  }
  
  ## Plot distributions
  par(mfrow = c(2, 1))
  hist(tt[[1]], xlim = c(0, 0.8), main = "No DLTs")
  abline(v = mean(tt[[1]]), lty = 3, lwd = 2)
  hist(tt[[2]], xlim = c(0, 0.8), main = ">0 DLTs")
  abline(v = mean(tt[[2]]), lty = 3, lwd = 2)
  
  
  ## Suggested hypotheses
  sug = lapply(tt, mean) # hypothesis values
  
  return(sug)
}

ntox <- 3 #Number of unique toxicities
#### Define the weight Matrix ####
W <- matrix(c(0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 1
              0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 2
              0, 0.00, 0.00, 0.5, 1), ## Burden weight for grades 1-4 for toxicity 3
            nrow = ntox, byrow = T)


nTTP.hypothesis(ntox, W)
