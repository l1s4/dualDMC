# dualDMC

Extension of the Diffusion Model for Conflict Tasks (DMC) for two sources of irrelevant information

Parameters: \
$\mu_c$: drift rate controlled process \
$b$: upper decision boundary / threshold (-b corresponds to lower decision boundary) \
$A_1$: amplitude first automatic activation \
$A_2$: amplitude second automatic activation \
$\tau_1$: scale parameter gamma function first automatic process \
$\tau_2$: scale parameter gamma function second automatic process \
$\sigma$: diffusion constant superimposed process \
Note: $a_1 = a_2 = 2$


## BayesFlow

Code used for specifying the experiment, the priors and the simulation settings is based on "Amortized Bayesian Workflow for Modeling Congruency Effects Using the Diffusion Model for Conflict Tasks." Schaefer, S.B., Radev, S.T., Göttmann, J. et al. Comput Brain Behav (2026). https://doi.org/10.1007/s42113-026-00266-y \
https://github.com/simschaefer/amortized-dmc   

BayesFlow models were trained using the computational resource bwUniCluster funded by the Ministry of Science, Research and the Arts Baden-Württemberg and the Universities of the State of Baden-Württemberg, Germany, within the framework program bwHPC.
