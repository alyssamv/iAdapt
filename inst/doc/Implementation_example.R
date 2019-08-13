## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(iAdapt)

# Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety
p_no <- 0.40                                     
p_yes <- 0.15    

# Likelihood-ratio (LR) threshold
K <- 2                                          

# Cohort size used in stage 1
coh.size <- 3 

# number of observed DLTs
x <- 1


## ------------------------------------------------------------------------
LRtox(coh.size, x, p_no, p_yes, K)

## ------------------------------------------------------------------------
LRtox(coh.size, x = 2, p_no, p_yes, K)

## ------------------------------------------------------------------------
y.eff <- c(9, 1, 0, 34, 10, 27, 38, 42, 60, 75, 48, 62)
d.safe <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
rand.prob(y.eff, d.safe)

