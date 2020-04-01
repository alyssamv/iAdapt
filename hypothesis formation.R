
source(file.path(getwd(), "perdose source fns.R"))

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

ntox <- 3 #Number of unique toxicities
grade.thresh = c(4, 4, 4)
#### Define the weight Matrix ####
W <- matrix(c(0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 1
              0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 2
              0.00, 0.00, 0.5, 1), ## Burden weight for grades 1-4 for toxicity 3
            nrow = ntox, byrow = T)

ok = nTTP.hypothesis(ntox, W, grade.thresh)

W <- matrix(c(0.1, 0.35, 0.7, 1.00, #Burden weight for grades 1-4 for toxicity 1
              0.08, 0.23,  0.6, 0.80, #Burden weight for grades 1-4 for toxicity 2
              0.00, 0.15, 0.45, 0.80), ##Burden weight for grades 1-4 for toxicity 3
            nrow = ntox, byrow = T)

nTTP.hypothesis(ntox, W, grade.thresh)
