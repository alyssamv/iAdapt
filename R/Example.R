
# One trial example on how to generate toxicity, efficacy data, and dose allocations for each stage of the adaptive design
# Function sim.trials can be used for repeated trials

#######################################################################
#                           Input parameters                          #
#######################################################################

# Number of pre-specified dose levels
dose <- 5 

# Vector of true toxicities associated with each dose
dose.tox <- c(0.05, 0.10, 0.15, 0.20, 0.30)       

# Acceptable (p2) and unacceptable (p1) DLT rates used for establishing safety
p1 <- 0.40                                     
p2 <- 0.15    

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

#########################################################################
#                              Stage 1                                  #
#########################################################################

# Function to generate and tabulate toxicities per dose level

tox.profile(dose, dose.tox, p1, p2, K, coh.size)
# Output: dose assignment(col1), toxicities per dose(col2), coh.number(col3), LR value(col4)
#     [,1] [,2] [,3] [,4]
#[1,]    1    0    1 2.84
#[2,]    2    0    2 2.84
#[3,]    3    0    3 2.84
#[4,]    4    0    4 2.84
#[5,]    5    2    5 0.20                          

# Dose 5 is unacceptably toxic (LR ≤ 1/2), thus it won't be used in stage 2

# Function to select only the acceptable toxic doses 
 safe.dose(dose, dose.tox, p1, p2, K, coh.size) 
# Output: safe dose assignment(col1), toxicities per dose(col2)

#$alloc.safe                                       
#       [,1] [,2]
#[1,]    1    0
#[2,]    2    0
#[3,]    3    0
#[4,]    4    0

# Other output
#$alloc.total
#[1] 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5           # complete dose allocation vector used in stage 1

#$n1
#[1] 15                                       # no. subjects used in stage 1

# Function to generate efficacy outcomes (here percent persistence) for each dose
gen.eff.stg1(dose, dose.tox, p1, p2, K, coh.size, m, v, nbb=100)

# Selected output
# Efficacy outcomes only for safe doses to be used in stage 2
#$Y.safe
#[1]  9  1  0 34 10 27 38 42 60 75 48 62          

# Safe doses to be used in stage 2
#$d.safe
#[1] 1 1 1 2 2 2 3 3 3 4 4 4                      

# Number of toxicities for safe doses
#$tox.safe
#[1] 0 0 0 0   

# Other output not shown here: vector of all outcomes and dose assignments for unacceptable + safe doses

#############################################################################
#                                 Stage 2                                   #
#############################################################################

# Function to fit a linear regression for the continuous efficacy outcomes,
# compute the randomization probabilities per dose and allocate the next subject 
# to an acceptable safe dose that has the highest randomization probability 
# (cohort of size 1) 
rand.stg2(dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule=9, cohort=1, samedose=T, nbb=100) 

# Output: vector of all efficacy outcomes (N=25)
#$Y.final
#[1]  9  1  0 34 10 27 38 42 60 75 48 62 90 89 73 50 47 55 68 64 48 80 45 28 57        

# 9 1 0 34 10 27 38 42 60 75 48 62 90 89 73 – efficacy outcomes from stage 1 
# Note: only efficacy from acceptable doses are used to compute randomization probabilites
# Here: 9  1  0 34 10 27 38 42 60 75 48 62 (see output from safe.dose above)

# 50 47 55 68 64 48 80 45 28 57 – efficacy outcomes from stage 2

# Output: vector of all dose allocations (N=25)
#$d.final
#[1] 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 3 4 4 4 4 3 4 4 2 3                                 

# 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 - dose allocation from stage 1 (cohorts of 3 subj. per dose)
# 3 4 4 4 4 3 4 4 2 3 - dose allocation from stage 2 (cohorts of 1) - no dose 5, b/c it was considered unacceptably toxic in stage 1

#The complete vectors of dose allocations and efficacy outcomes can be used to compute the operating characteristics of the design in repeated simulations.


