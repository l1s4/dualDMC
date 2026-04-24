library(lattice)
library(ggplot2)
library(dplyr)
library(microbenchmark)
library(Rcpp)
source('helpers.R')
Rcpp::sourceCpp('../src/dualDMC.cpp')


# just for testing #############################################################
source('../src/dualDMC.R')
result_cpp <- simDDMCtrial(0.5, 50, 30, 10, 30, 150, 0.01, 0.0, -1, 1)
print(result_cpp$rt)
print(result_cpp$dec)
plot(result_cpp$XTraj)

result_R <- MDMC_T(1500, list(sigma = 0.00, dt = 0.01, mu_c = 0.5, 
b = 50, tau1 = 30, tau2 = 150, A1 = 30, A2 = 10, a1 = 2, a2 = 2, 
automatic1 = "incongruent", automatic2 = "congruent"))
print(result_R$rt)
print(result_R$dec)

# varying parameters
tau1_vals <- seq(20, 100, 40)
tau2_vals <- seq(20, 100, 40)
A1_vals <- seq(10, 25, 10)
A2_vals <- seq(10, 25, 10)

param_grid <- expand.grid(tau1 = tau1_vals, tau2 = tau2_vals, 
                          A1 = A1_vals, A2 = A2_vals)

# constant parameters
param_grid$a1    <- 2
param_grid$a2    <- 2
param_grid$b     <- 60
param_grid$sigma <- 0
param_grid$dt    <- 0.01
param_grid$mu_c  <- 0.6

res <- simDDMC(param_grid, 10)
resR <- MDMC_Sim(10, 1000, param_grid)
df_all <- bind_rows(resR, .id = "param_set")
df_all$param_set <- as.integer(df_all$param_set)
df_all$congruency <- paste(df_all$first_pr, df_all$second_pr, sep = "_")



# DDMC #########################################################################
# varying parameters
tau1_vals <- seq(20, 200, 40)   # for testing
tau2_vals <- seq(20, 200, 40)   # for testing
A1_vals <- seq(15, 25, 5)
A2_vals <- seq(15, 25, 5)
#tau1_vals <- seq(20, 200, 20)
#tau2_vals <- seq(20, 200, 20)
#A1_vals <- seq(10, 25, 5)
#A2_vals <- seq(10, 25, 5)

#b_vals <- seq(50, 70, 10)
#mu_c_vals <- seq(0.4, 0.6, 0.1)

param_grid <- expand.grid(tau1 = tau1_vals, tau2 = tau2_vals, 
                          A1 = A1_vals, A2 = A2_vals)

# constant parameters
param_grid$a1    <- 2
param_grid$a2    <- 2
param_grid$b     <- 60
param_grid$sigma <- 4 
param_grid$dt    <- 0.01
param_grid$mu_c  <- 0.6

N_sim <- 500		# number of simulations per parameter set

datDDMC <- simDDMC(param_grid, N_sim)


# save data
write.csv(datDDMC, "out/data/datDDMC.csv")
datDDMC <- read.csv("out/data/datDDMC.csv")



# Plots ########################################################################
# Note: only varying A1, A2, tau1, tau2

# get parameter sets with A1 == A2
df_As_eq <- datDDMC[datDDMC$A1 == datDDMC$A2, ]
df_As_eq <- mk_congruency(df_As_eq)   # add congruency column
df_As_lst <- split(df_As_eq, df_As_eq$A1)

# get parameter sets with tau1 == tau2
df_taus_eq <- datDDMC[datDDMC$tau1 == datDDMC$tau2, ]
df_taus_eq <- mk_congruency(df_taus_eq)   # add congruency column
df_taus_lst <- split(df_taus_eq, df_taus_eq$tau1)

# Mean RTs
#TODO: add info about other parameters somewhere in the plot
lapply(df_As_lst, plt_var_taus_rt, corr_only = T)   # plot for every level of As
lapply(df_taus_lst, plt_var_As_rt, corr_only = T)   # plot for every level of taus


# Mean ERs
lapply(df_As_lst, plt_var_taus_er)
lapply(df_taus_lst, plt_var_As_er)

# CDFs
#TODO: correct trials only? 
lapply(df_As_lst, plt_cdfs)
lapply(df_taus_lst, plt_cdfs)


# Delta plot
#TODO: correct trials only? 
lapply(df_As_lst, plt_delta, 
       cond1 = "congruent_congruent", cond2 = "congruent_incongruent")
lapply(df_taus_lst, plt_delta, 
       cond1 = "congruent_congruent", cond2 = "congruent_incongruent")
