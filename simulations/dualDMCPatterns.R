library(lattice)
library(latticeExtra)
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

## constant parameters
#param_grid$a1     <- 2
#param_grid$a2     <- 2
#param_grid$b      <- 60
#param_grid$sigma  <- 4 
#param_grid$dt     <- 0.01
#param_grid$mu_c   <- 0.6#0.5
#param_grid$ndt_m  <- 330
#param_grid$ndt_sd <- 30


#N_sim <- 5000		# number of simulations per parameter set
#datDDMC <- simDDMC(param_grid, N_sim)


#save(datDDMC, file = "out/data/datDDMC6.RData")      # save data
load("out/data/datDDMC1.RData")
#load("out/data/datDDMC2.RData")
#load("out/data/datDDMC3.RData")
#load("out/data/datDDMC4.RData")
#load("out/data/datDDMC5.RData")
#load("out/data/datDDMC6.RData")

dat_nr <- "_1"    # set plot suffix

# Descriptive ##################################################################
#dat <- mk_congruency(datDDMC)
#dat$congruency <- as.factor(dat$congruency)
#dat$error <- ifelse(dat$dec == 1, 0, 1)
#
## RT
#mean(dat$rt) |> round(2)
#sd(dat$rt) |> round(2)
#aggregate(rt ~ error, data = dat, FUN = mean) |> round(2)
#aggregate(rt ~ error, data = dat, FUN = sd) |> round(2)
#
#mean(dat$error) * 100 # ER
#sd(dat$error) *100    # ER
#
#
## per congruency condition
#aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = mean)
#aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = mean)$error * 100
#aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = sd)
#aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = sd)$error * 100



## calculate congruency effects and pattern #####################################
subset <- dat[dat$A1 == 25 & dat$A2 == 15 & dat$tau1 == 30 & dat$tau2 == 30, ]

aggregate(rt ~ congruency, FUN=mean, data=subset)
check_ceff <- function(dat, A1, A2, tau1, tau2) {
  df <- dat[
    dat$A1 == A1 & dat$A2 == A2 & dat$tau1 == tau1 & dat$tau2 == tau2 & 
      dat$dec == 1, ]
  
  meanrts <- aggregate(rt ~ congruency, FUN = mean, df)$rt
  ceff1 <- meanrts[2] - meanrts[1]
  ceff2 <- meanrts[4] - meanrts[3]
  diff_diffs <- ceff1 - ceff2
#  print(paste(
#    "A1 =", A1, "A2 =", A2, "tau1 =", tau1, "tau2 =", tau2, ":", 
#    "Diff (ci-cc) - (ii-ic) =", diff_diffs))
  
  return(diff_diffs)
}
check_ceff(dat, 10, 25, 30, 30)

param_grid$diffs <- NA
for (i in 1:nrow(param_grid)) {
  print(paste("processing set :", i))
  param_grid[i, ]$diffs <- check_ceff(dat, param_grid[i,]$A1, param_grid[i,]$A2, 
                                 param_grid[i,]$tau1, param_grid[i,]$tau2)
  
}
min(param_grid$diffs)
max(param_grid$diffs)
param_grid$pattern <- ifelse(param_grid$diffs < -30, "OA", 
                             ifelse(param_grid$diffs > 30, "UA", "A"))

param_grid[param_grid$pattern == "OA", ]
param_grid[param_grid$pattern == "UA", ]

param_grid[
  param_grid$A1 == param_grid$A2 & param_grid$A1 == 15 & 
    param_grid$pattern == "UA", ]

param_grid[
  param_grid$tau1 == param_grid$tau2 & param_grid$tau1 == 30 & 
    param_grid$pattern == "UA", ]


## Plots ########################################################################
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
pdf(paste0("out/plots/RT_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_var_taus_rt, corr_only = T)   # plot for every level of As
dev.off()
pdf(paste0("out/plots/RT_vary_As", dat_nr, ".pdf"))
lapply(df_taus_lst, plt_var_As_rt, corr_only = T)   # plot for every level of taus
dev.off()

# Mean ERs
pdf(paste0("out/plots/ER_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_var_taus_er)
dev.off()
pdf(paste0("out/plots/ER_vary_As", dat_nr, ".pdf"))
lapply(df_taus_lst, plt_var_As_er)
dev.off()

# CAFs
pdf(paste0("out/plots/CAF_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_cafs_var_As, n_bins = 4)
dev.off()
pdf(paste0("out/plots/CAF_vary_As", dat_nr, ".pdf"))
lapply(df_taus_lst, plt_cafs_var_taus, n_bins = 4)
dev.off()

# Delta plots
pdf(paste0("out/plots/delta_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_unc_delta_vary_taus)
dev.off()
pdf(paste0("out/plots/delta_vary_As", dat_nr, ".pdf")) 
lapply(df_taus_lst, plt_unc_delta_vary_As)
dev.off()

## CDFs
#pdf(paste0("out/plots/CDF_vary_taus", dat_nr, ".pdf"))
#lapply(df_As_lst, plt_cdfs_var_taus)
#dev.off()
#pdf(paste0("out/plots/CDF_vary_As", dat_nr, ".pdf"))
#lapply(df_taus_lst, plt_cdfs_var_As)
#dev.off()

## Plot densities
#pdf(paste0("out/plots/dens_vary_taus", dat_nr, ".pdf"))
#lapply(df_As_lst, rt_denss_vary_tau)
#dev.off()
#pdf(paste0("out/plots/dens_vary_As", dat_nr, ".pdf"))
#lapply(df_taus_lst, rt_denss_vary_A)
#dev.off()