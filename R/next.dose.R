
# Functions that can be used in ongoing trials to inform about the next dose allocation
# In Stage 1 - these can be used for escalation based on toxicity (DLT) info
# In Stage 2 - these can be used to compute the randomization probabilities and also to check safety at each dose



# Calculate the likelihood ratio per dose for observed toxicities 

# coh.size - cohort size (number of patients) per dose
# x - number of observed DLTs per dose
# p1 - toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2 - toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K - threshold for LR. Takes integer values: 1,2,...(recommended K=2)

LRtox <- function(coh.size, x, p1, p2, K) {
  
  l.p2   <- x * log(p2) + (coh.size-x) * log(1 - p2) # likelihood of acceptable/alternative hypothesis 
  l.p1   <- x * log(p1) + (coh.size-x) * log(1 - p1) # likelihood of unacceptable/null hypothesis
  LR     <- round(exp(l.p2 - l.p1), 2) 
  
  if (LR > 1/K)  print ("Safe/Escalate") else print ("Unsafe/Stop")
  
  return(list(LR=LR))
}

# Test
# LRtox(coh.size=3,x=2,p1=0.40,p2=0.15,K=2)
# LRtox(coh.size=3,x=1,p1=0.40,p2=0.15,K=2)


# Fit a linear model, calculate the randomization probabilities for safe doses 
# and output the next dose allocation

# y.eff <- vector of all efficacy outcomes for each dose allocation
# d.safe <- vector of dose assignment

rand.prob <- function(y.eff, d.safe){
  
  reg <- lm(log(y.eff + 1) ~ factor(d.safe))         # Linear model with log(Y) for accept. doses 
  fit <- as.vector(reg$fitted.values)                # Fitted values for Y
  fitp <- exp(fit)
  
  dose.unique <- duplicated(d.safe)
  fitp <- fitp[dose.unique == F]
  
  rp <- fitp/sum(fitp)                                # Calculate randomization prob. for each dose                            
  rp <- ifelse(rp < 0.02, 0.02, rp) 
  
  rec.dose <- which(rp==max(rp))                      # Next dose with max rand. prob
  
  return(list(Rand.Prob=rp, Next.Dose=rec.dose))
  
}

#Test
# y.eff <- c(9,1,0,34,10,27,38,42,60,75,48,62)
# d.safe <- c(1,1,1,2,2,2,3,3,3,4,4,4)
# rand.prob(y.eff, d.safe)


# y.eff <- c(9,1,0,34,10,27,38,42,60,75,48,62,90,89,98,100)
# d.safe <- c(1,1,1,2,2,2,3,3,3,4,4,4,2,2,3,3)
# rand.prob(y.eff, d.safe)
