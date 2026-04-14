library(lattice)
library(dplyr)
source('../MDMC-Functions.R')
source('helpers.R')

# MDMC #########################################################################
# varying parameters
tau1_vals <- seq(30, 100, 70)
tau2_vals <- seq(30, 100, 70)
A1_vals <- seq(10, 50, 20)
A2_vals <- seq(10, 50, 20)
#b_vals <- seq(50, 70, 10)
#mu_c_vals <- seq(0.4, 0.6, 0.1)

param_grid <- expand.grid(tau1 = tau1_vals, tau2 = tau2_vals, 
                          A1 = A1_vals, A2 = A2_vals)

# constant parameters
param_grid$a1    <- 2
param_grid$a2    <- 2
param_grid$b     <- 50
param_grid$sigma <- 4
param_grid$dt    <- 0.01
param_grid$mu_c  <- 0.5

N_sim <- 500		# number of simulations per parameter set
N_time <- 500		# max time

dat_all <- MDMC_Sim(N_sim, N_time, param_grid)      # simulate
means <- get_mean_dvs_MDMC(dat_all, TRUE)                # get mean RT and mean ER

df_rts <- means$df_rts
df_ers <- means$df_ers

df_rts$rt <- dfrt$rt + 300		# add ndt

# save mean rts and ers of simulated data
write.csv(df_rts, "sim-data/df_rts.csv")
write.csv(df_ers, "sim-data/df_ers.csv")

df_rts <- read.csv("sim-data/df_rts.csv")
df_ers <- read.csv("sim-data/df_ers.csv")


df_rts$first_pr <- as.factor(df_rts$first_pr)
df_rts$second_pr <- as.factor(df_rts$second_pr)


# Plots
df_1 <- df_rts[df_rts$A1 == 10 & df_rts$A2 == 10, ]
plt_A_tau(df_1, "RT")

df_2 <- df_rts[df_rts$tau1 == 30 & df_rts$tau2 == 30, ]
plt_tau_A(df_2, "RT")

# MRDM #########################################################################
delta1_vals <- seq(0.3, 0.8, 0.2)
delta2_vals <- seq(0.3, 0.8, 0.2)
P1_vals <- seq(5, 15, 5)
P2_vals <- seq(5, 15, 5)
param_grid_l <- expand.grid(delta1 = delta1_vals, delta2 = delta2_vals, 
                            P1 = P1_vals, P2 = P2_vals)

## constant parameters
param_grid_l$b <- 50
param_grid_l$sigma <- 4
param_grid_l$dt <- 1
param_grid_l$mu_c <- 0.25


dat_all_l <- MRDM_Sim(N_sim, N_time, param_grid_l)

means_l <- get_mean_dvs_MDDM(dat_all_l, TRUE)

df_rts_l <- means_l$df_rts
df_ers_l <- means_l$df_ers

# save mean rts and ers of simulated data
write.csv(df_rts_l, "sim-data/df_rts_l.csv")
write.csv(df_ers_l, "sim-data/df_ers_l.csv")

df_rts_l <- read.csv("sim-data/df_rts_l.csv")
df_ers_l <- read.csv("sim-data/df_ers_l.csv")

df_rts_l$first_pr <- as.factor(df_rts_l$first_pr)
df_rts_l$second_pr <- as.factor(df_rts_l$second_pr)
df_ers_l$first_pr <- as.factor(df_ers_l$first_pr)
df_ers_l$second_pr <- as.factor(df_ers_l$second_pr)

df_rts_l <- df_rts_l[df_rts_l$rt != Inf, ]    # exclude potential Infs

# Plots
df_1 <- df_rts_l[df_rts_l$delta1 == 0.3 & df_rts_l$delta2 == 0.3, ]
plt_delta_P(df_1, "RT")

df_2 <- df_ers_l[df_ers_l$delta1 == 0.3 & df_ers_l$delta2 == 0.3, ]
plt_delta_P(df_2, "ER")

df_3 <- df_rts_l[df_rts_l$P1 == 5 & df_rts_l$P2 == 5, ]
plt_P_delta(df_3, "RT")

df_4 <- df_ers_l[df_ers_l$P1 == 5 & df_ers_l$P2 == 5, ]
plt_P_delta(df_4, "ER")
