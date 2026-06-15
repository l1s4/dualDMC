library(lattice)
library(latticeExtra)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(microbenchmark)
library(Rcpp)
source('helpers.R')
Rcpp::sourceCpp('../src/dualDMC.cpp')



# varying parameters
#tau1_vals <- seq(20, 200, 40)
#tau2_vals <- seq(20, 200, 40)
#A1_vals <- seq(15, 25, 5)
#A2_vals <- seq(15, 25, 5)
tau1_vals <- seq(30, 230, 50)
tau2_vals <- seq(30, 230, 50)
A1_vals <- seq(10, 25, 5)
A2_vals <- seq(10, 25, 5)
#b_vals <- seq(50, 70, 10)
#mu_c_vals <- seq(0.4, 0.6, 0.1)

param_grid <- expand.grid(tau1 = tau1_vals, tau2 = tau2_vals, 
                          A1 = A1_vals, A2 = A2_vals)

# constant parameters
param_grid$a1     <- 2
param_grid$a2     <- 2
param_grid$b      <- 60
param_grid$sigma  <- 4 
param_grid$dt     <- 0.01
param_grid$mu_c   <- 0.6
param_grid$ndt_m  <- 300
param_grid$ndt_sd <- 30

N_sim <- 2000		# number of simulations per parameter set

datDDMC <- simDDMC(param_grid, N_sim)


save(datDDMC, file = "out/data/datDDMC.RData")      # save data
load("out/data/datDDMC.RData")
#write.csv(datDDMC, "out/data/datDDMC.csv")
#datDDMC <- read.csv("out/data/datDDMC.csv")




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
pdf("out/plots/RT_plt_vary_taus.pdf")
lapply(df_As_lst, plt_var_taus_rt, corr_only = T)   # plot for every level of As
dev.off()
pdf("out/plots/RT_plt_vary_As.pdf")
lapply(df_taus_lst, plt_var_As_rt, corr_only = T)   # plot for every level of taus
dev.off()

# Mean ERs
pdf("out/plots/ER_plt_vary_taus.pdf")
lapply(df_As_lst, plt_var_taus_er)
dev.off()
pdf("out/plots/ER_plt_vary_As.pdf")
lapply(df_taus_lst, plt_var_As_er)
dev.off()

# CDFs
pdf("out/plots/CDF_plt_vary_taus.pdf")
lapply(df_As_lst, plt_cdfs_var_taus)
dev.off()
pdf("out/plots/CDF_plt_vary_As.pdf")
lapply(df_taus_lst, plt_cdfs_var_As)
dev.off()

# CAFs
pdf("out/plots/CAF_plt_vary_taus.pdf")
lapply(df_As_lst, plt_cafs_var_As, n_bins = 4)
dev.off()
pdf("out/plots/CAF_plt_vary_As.pdf")
lapply(df_taus_lst, plt_cafs_var_taus, n_bins = 4)
dev.off()

# Delta plots
pdf("out/plots/delta_plt_vary_taus.pdf")
lapply(df_As_lst, plt_delta_vary_taus)
dev.off()
pdf("out/plots/delta_plt_vary_As.pdf")
lapply(df_taus_lst, plt_delta_vary_As)
dev.off()

# Plot densities
pdf("out/plots/dens_plt_vary_taus.pdf")
lapply(df_As_lst, rt_denss_vary_tau)
dev.off()
pdf("out/plots/dens_plt_vary_As.pdf")
lapply(df_taus_lst, rt_denss_vary_A)
dev.off()
