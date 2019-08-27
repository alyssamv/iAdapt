## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ------------------------------------------------------------------------
library(iAdapt)

# Number of pre-specified dose levels
dose <- 5 

# Vector of true toxicities associated with each dose
dose.tox <- c(0.05, 0.10, 0.20, 0.35, 0.45)       

# Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety
p_no <- 0.40                                     
p_yes <- 0.15    

# Likelihood-ratio (LR) threshold
K <- 2                                          

# Cohort size used in stage 1
coh.size <- 3 

# Vector of true mean efficacies per dose (here mean T-cell persistence per dose (%))
m <- c(5, 15, 40, 65, 80)   # MUST BE THE SAME LENGTH AS dose.tox                  

# Efficacy (equal) variance per dose
v <- rep(0.01, 5) 

# Total sample size (stages 1&2)                            
N <- 25                                        

# Stopping rule: if dose 1 is the only safe dose, allocate up to 9 pts before ending the trial to collect more information
stop.rule <- 9   

## ---- echo=F, fig.width=6, fig.height=3----------------------------------
par(mfrow = c(1,2))
plot(x = 1:5, y = dose.tox, ylim = c(0, p_no + 0.1),
     type = 'b',
     xlab = "Dose level", ylab = "True toxicity rate",
     main = "(a) Dose-toxicity")
lines(x = 0:6, y = rep(p_no, 7), lty = 3) # threshold toxicity value
points(x = 4, y = dose.tox[4], pch = 15, col = "green") # highlight best dose

plot(x = 1:5, y = m, ylim = c(0, 100),
     type = 'b',
     xlab = "Dose level", ylab = "True efficacy (T-cell % persistence)",
     main = "(b) Dose-efficacy")
points(x = 4, y = m[4], pch = 15, col = "green") # highlight best dose

dev.off()

## ------------------------------------------------------------------------
set.seed(3)
tox.profile(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size)

## ------------------------------------------------------------------------
set.seed(3)
safe.dose(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size) 

## ------------------------------------------------------------------------
set.seed(3)
eff.stg1(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size, m, v, nbb = 100)

## ------------------------------------------------------------------------
set.seed(3)
rand.stg2(dose, dose.tox, p_no, p_yes, K, coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100) 

## ---- results='hide'-----------------------------------------------------
numsims = 100

set.seed(3)
simulations = sim.trials(numsims = numsims, dose, dose.tox, p1 = p_no, p2 = p_yes, K, coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100)

## ------------------------------------------------------------------------
head(simulations$safe.d)
head(simulations$sim.Y)
head(simulations$sim.d)

## ------------------------------------------------------------------------
colSums(simulations$safe.d) / numsims

## ------------------------------------------------------------------------
sim.tables = sim.summary(simulations)

## ------------------------------------------------------------------------
# Vector of true mean efficacies per dose (here mean percent persistence per dose)
m <- c(15, 35, 80, 60, 40)   # MUST BE THE SAME LENGTH AS dose.tox                  


## ---- echo=F, fig.width=6, fig.height=3----------------------------------
par(mfrow = c(1,2))
plot(x = 1:5, y = dose.tox, ylim = c(0, p_no + 0.1),
     type = 'b',
     xlab = "Dose level", ylab = "True toxicity rate",
     main = "(a) Dose-toxicity")
lines(x = 0:6, y = rep(p_no, 7), lty = 3) # threshold toxicity value
points(x = 3, y = dose.tox[3], pch = 15, col = "green") # highlight most effecive dose
points(x = 4, y = dose.tox[4], pch = 15, col = "red") # highlight MTD

plot(x = 1:5, y = m, ylim = c(0, 100),
     type = 'b',
     xlab = "Dose level", ylab = "True efficacy (T-cell % persistence)",
     main = "(b) Dose-efficacy")
points(x = 3, y = m[3], pch = 15, col = "green") # highlight best dose
points(x = 4, y = m[4], pch = 15, col = "red") # highlight MTD

dev.off()

## ------------------------------------------------------------------------
set.seed(1)
tox.profile(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size)

## ------------------------------------------------------------------------
set.seed(1)
safe.dose(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size) 

## ------------------------------------------------------------------------
set.seed(1)
eff.stg1(dose = dose, dose.tox = dose.tox, p1 = p_no, p2 = p_yes, K = K, coh.size = coh.size, m, v, nbb = 100)

## ------------------------------------------------------------------------
set.seed(1)
rand.stg2(dose, dose.tox, p_no, p_yes, K, coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100) 

## ---- results='hide'-----------------------------------------------------
numsims = 100

set.seed(1)
simulations = sim.trials(numsims = numsims, dose, dose.tox, p1 = p_no, p2 = p_yes, K, coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100)

## ------------------------------------------------------------------------
head(simulations$safe.d)
head(simulations$sim.Y)
head(simulations$sim.d)

## ------------------------------------------------------------------------
colSums(simulations$safe.d) / numsims

## ------------------------------------------------------------------------
sim.tables = sim.summary(simulations)

