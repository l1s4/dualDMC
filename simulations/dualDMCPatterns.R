library(lattice)
library(ggplot2)
library(dplyr)
library(microbenchmark)
library(Rcpp)
source('../dualDMC.R')
source('helpers.R')
Rcpp::sourceCpp('../dualDMC.cpp')


# just for testing
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
write.csv(datDDMC, "sim-data/datDDMC.csv")
datDDMC <- read.csv("sim-data/datDDMC.csv")



# mean RT, ER ##################################################################

# loop parameter sets and extract subsets with A1 == A2
df_As_eq <- datDDMC[datDDMC$A1 == datDDMC$A2, ]
df_As_lst <- split(df_As_eq, df_As_eq$A1)

#TODO: move to function: below for each dataframe 

# get means per tau values and  congruency conditions
means <- aggregate(rt ~ auto1 + auto2 + tau1 + tau2, FUN = mean, data = df_As_lst[[1]])
means$auto1 <- factor(means$auto1, levels = c(1, -1), 
                      labels = c("congruent", "incongruent"))
means$auto2 <- factor(means$auto2, levels = c(1, -1), 
                      labels = c("congruent", "incongruent"))


#TODO: auto-detect if bigger than 5x5 => split into two plots?
#TODO: add info about other parameters somewhere in the plot
xyplot(rt ~ auto1 | factor(tau1) * factor(tau2), groups = auto2, data = means, 
       type = c("p", "a"), main = "Mean RTs", 
       strip = strip.custom(strip.names = c(TRUE, TRUE), 
                            var.name = c("tau1", "tau2")), 
       auto.key = list(title = "second automatic process", space = "right", 
                       cex = .7),
)



# loop parameter sets and extract subsets with tau1 == tau2
df_taus_eq <- datDDMC[datDDMC$tau1 == datDDMC$tau2, ]
df_taus_lst <- split(df_taus_eq, df_taus_eq$tau1)




################################################################################
# remove below after double-checking above solution
# MDMC
means <- get_mean_dvs_DDMC(, TRUE)           # get mean RT and mean ER

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
df_ers$first_pr <- as.factor(df_ers$first_pr)
df_ers$second_pr <- as.factor(df_ers$second_pr)

# Plots
df_1 <- df_rts[df_rts$A1 == 10 & df_rts$A2 == 10, ]
plt_A_tau(df_1, "RT")

df_2 <- df_rts[df_rts$tau1 == 30 & df_rts$tau2 == 30, ]
plt_tau_A(df_2, "RT")

df_3 <- df_ers[df_ers$tau1 == 30 & df_ers$tau2 == 30, ]
plt_tau_A(df_3, "ER")

# MRDM 
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







# CDFs #########################################################################

# MDMC
df1 <- df_all[df_all$A1 == 10 & df_all$A2 == 10 , ]
df2 <- df_all[df_all$A1 == 30 & df_all$A2 == 30 , ]

# CDFs
plt_cdfs_mdmc(df2)

# Delta plot
df_cc <- df2[df2$congruency == "congruent_congruent", ]
df_ci <- df2[df2$congruency == "congruent_incongruent", ]
df_ic <- df2[df2$congruency == "incongruent_congruent", ]
df_ii <- df2[df2$congruency == "incongruent_incongruent", ]

plt_delta_mdmc(df_cc, df_ci)








