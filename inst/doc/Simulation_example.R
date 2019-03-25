## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(iAdapt)

## ------------------------------------------------------------------------
# Number of pre-specified dose levels
dose <- 5 

# Vector of true toxicities associated with each dose
dose.tox <- c(0.05, 0.10, 0.15, 0.20, 0.30)       

# Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety
p_no <- 0.40                                     
p_yes <- 0.15    

# Likelihood-ratio (LR) threshold
K <- 2                                          

# Cohort size used in stage 1
coh.size <- 3 

# Vector of true mean efficacies per dose (here mean percent persistence per dose)
m <- c(5, 15, 40, 65, 80)                     

# Efficacy(equal) variance per dose
v <- rep(0.01,5) 

# Total sample size (stages 1&2)                            
N <- 25                                        

# Stopping rule: if dose 1 is the only safe dose, allocate up to 9 pts.
stop.rule <- 9   

## ------------------------------------------------------------------------
tox.profile(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size)

## ------------------------------------------------------------------------
safe.dose(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size) 

## ------------------------------------------------------------------------
eff.stg1(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size, m, v, nbb = 100)

## ------------------------------------------------------------------------
rand.stg2(dose, dose.tox, p_no, p_yes, K, coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = T, nbb = 100) 

