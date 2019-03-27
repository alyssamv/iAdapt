# iAdapt
R package for phase I two-stage adaptive dose-finding clinical trial design.

This package provides software based on the phase I trial design by [Chiuzan et al. (2018)](https://www.tandfonline.com/doi/abs/10.1080/19466315.2018.1462727). Stage 1 is safety-driven dose-escalation, and Stage 2 employs efficacy-driven randomization while continuing to monitor dose safety. This design is relevant in the face of a non-monotonous dose-toxicity relationship - a phenomenon that most often seen in immunologic therapies [REFERENCE]. Additionally, the probabilistic nature (as opposed to rule-based) provides an advantage in identifying the optimal dose to carry forward in development, by allowing more than one dose to be examined for efficacy. 

The design uses a likelihood paradigm to determine the acceptability (safety) of a dose. e.g. In Stage 1, when the likelihood ratio for a dose is greater than a prespecified threshold, the dose is considered acceptable and subsequent patients are enrolled on the next dose level. Conversely, if the likelihood ratio is less than or equal to the threshold, escalation is stopped.

The function of this package is two-fold:
(i) Produce trial outcomes through simulation for an inputted scenario, and
(ii) Implement the design in a real trial.

Vignettes are provided to walk the user through each function of the package. 