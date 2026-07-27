library(lattice)
library(latticeExtra)
library(dplyr)
library(purrr)
library(tidyr)
library(microbenchmark)
library(Rcpp)
source('helpers.R')
Rcpp::sourceCpp('../src/dualDMC.cpp')

## plot activation #############################################################
N <- 400    # number of timepoints
b <- 60
activ_cc <- simDDMCactivation(N, mu_c = 0.5, b, 20, 20, 30, 30, 1, 1)
activ_ci <- simDDMCactivation(N, mu_c = 0.5, b, 20, 20, 30, 30, 1, -1)
activ_ic <- simDDMCactivation(N, mu_c = 0.5, b, 20, 20, 30, 30, -1, 1)
activ_ii <- simDDMCactivation(N, mu_c = 0.5, b, 20, 20, 30, 30, -1, -1)

plot_activ <- function(activ, N, b, congruency) {
  plot(activ$cont_traj, type = "n", xlab = "RT (ms)", ylab = "E[X(t)]", 
       main = congruency, xlim = c(0, N-100), ylim = c(-b+30, b+20))
  abline(h = 0)
  lines(activ$cont_traj, col = "black")
  lines(activ$auto1_traj, col = "red")
  lines(activ$auto2_traj, col = "blue")
  lines(activ$super_traj, col = "green")
  abline(h = c(-b, b), v = which(activ$super_traj >= 60)[1], lty = 2)
  legend("topright", legend = c("controlled", "auto 1", "auto 2", "superimposed"), 
         col = c("black", "red", "blue", "green"), lty = 1, bg = "white")
}

pdf("../../MA/images/UA_Activation.pdf", 12, 8)
par(mfrow = c(2, 2))
plot_activ(activ_cc, N, b, "congruent-congruent")
plot_activ(activ_ci, N, b, "congruent-incongruent")
plot_activ(activ_ic, N, b, "incongruent-congruent")
plot_activ(activ_ii, N, b, "incongruent-incongruent")
dev.off()


## plot CAFs for A1 = A2 = 25 & tau1 = tau2 = 30 for each simulation ###########
dat1 <- get(load("out/data/datDDMC1.RData"))
dat2 <- get(load("out/data/datDDMC2.RData"))
dat3 <- get(load("out/data/datDDMC3.RData"))
dat4 <- get(load("out/data/datDDMC4.RData"))
dat5 <- get(load("out/data/datDDMC5.RData"))
dat6 <- get(load("out/data/datDDMC6.RData"))

dat1$simulation <- 1
dat2$simulation <- 2
dat3$simulation <- 3
dat4$simulation <- 4
dat5$simulation <- 5
dat6$simulation <- 6

datAll <- rbind(dat1, dat2, dat3, dat4, dat5, dat6)
subset <- datAll[
  datAll$A1 == 25 & datAll$A2 == 25 & datAll$tau1 == 30 & datAll$tau2 == 30, 
  ]
df <- mk_congruency(subset)

n_bins = 5
df$error <- ifelse(df$dec == 1, 0, 1)
df$acc <- 1 - df$error
  
caf_data <- df %>% 
  group_by(simulation, congruency) %>% 
  mutate(rt_bin = ntile(rt, n_bins)) %>%
  group_by(simulation, congruency, rt_bin) %>% 
  summarise(mean_acc = mean(acc), mean_rt = mean(rt), .groups = "drop")

pdf("../../MA/images/CAFs_As25_taus30_sims.pdf", 8, 4)
par(mar = c(1, 1, 1, 1))
xyplot(mean_acc ~ as.numeric(rt_bin) | 
         factor(simulation, levels = c("2", "4", "6", "1", "3", "5")), 
       groups = factor(congruency, 
          labels = c("congruent-congruent", "congruent-incongruent", 
          "incongruent-congruent", "incongruent-incongruent")), 
       data = caf_data, xlab = "RT Bin", ylab = "Accuracy", type = c("b", "g"), 
       ylim = c(0, 1.1), 
       auto.key = list(cex = 0.8, title = "Condition", space = "right"),
       strip = strip.custom(strip.names = T, var.name = "Simulation", sep = " ")
)
dev.off()

## Simulate data ###############################################################
# varying parameters
tau1_vals <- seq(30, 230, 50)
tau2_vals <- seq(30, 230, 50)
A1_vals <- seq(10, 20, 5)
A2_vals <- seq(10, 20, 5)

muc_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
b_vals <- c(30, 40, 50, 60, 70, 80, 90)

#param_grid <- expand.grid(tau1 = tau1_vals, tau2 = tau2_vals, 
#                          A1 = A1_vals, A2 = A2_vals)
param_grid <- expand.grid(mu_c = muc_vals, b = b_vals)

# constant parameters
param_grid$a1     <- 2
param_grid$a2     <- 2
#param_grid$b      <- 60#50#70
param_grid$sigma  <- 4 
param_grid$dt     <- 0.01
#param_grid$mu_c   <- 0.6#0.5
param_grid$ndt_m  <- 330
param_grid$ndt_sd <- 30

param_grid$tau1 <- 80
param_grid$tau2 <- 80
param_grid$A1 <- 15
param_grid$A2 <- 15


N_sim <- 5000		# number of simulations per parameter set

datDDMC <- simDDMC(param_grid, N_sim)

save(datDDMC, file = "out/data/datDDMCx.RData")      # save data


## Descriptive #################################################################
load("out/data/datDDMC1.RData")
dat <- mk_congruency(datDDMC)
dat$congruency <- as.factor(dat$congruency)
dat$error <- ifelse(dat$dec == 1, 0, 1)

# RT
mean(dat$rt) |> round(2)
sd(dat$rt) |> round(2)
aggregate(rt ~ error, data = dat, FUN = mean) |> round(2)
aggregate(rt ~ error, data = dat, FUN = sd) |> round(2)

mean(dat$error) * 100 # ER
sd(dat$error) *100    # ER

# per congruency condition
aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = mean)
aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = mean)$error * 100
aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = sd)
aggregate(cbind(rt, error) ~ congruency, data = dat, FUN = sd)$error * 100


### calculate congruency effects and pattern ###################################
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


## Plots vary one parameter ####################################################

dat1 <- get(load("out/data/datDDMC1.RData"))
datx <- get(load("out/data/datDDMCx.RData"))

datx <- mk_congruency(datx)
dat1 <- mk_congruency(dat1)

dat_b <- datx[datx$A1 == 15 & datx$A2 == 15 & 
                datx$tau1 == 80 & datx$tau2 == 80 & 
                datx$mu_c == 0.6, ]
dat_muc <- datx[datx$A1 == 15 & datx$A2 == 15 & 
                  datx$tau1 == 80 & datx$tau2 == 80 & 
                  datx$b == 60, ]
dat_A1 <- dat1[dat1$A2 == 15 & 
                 dat1$tau1 == 80 & dat1$tau2 == 80 & 
                 dat1$b == 60 & dat1$mu_c == 0.6, ]
dat_A2 <- dat1[dat1$A1 == 15 & 
                 dat1$tau1 == 80 & dat1$tau2 == 80 & 
                 dat1$b == 60 & dat1$mu_c == 0.6, ]
dat_tau1 <- dat1[dat1$A1 == 15 & dat1$A2 == 15 & 
                   dat1$tau2 == 80 & 
                   dat1$b == 60 & dat1$mu_c == 0.6, ]
dat_tau2 <- dat1[dat1$A1 == 15 & dat1$A2 == 15 & 
                   dat1$tau1 == 80 & 
                   dat1$b == 60 & dat1$mu_c == 0.6, ]

# mean RTs
pdf("../../MA/images/meanRT_var_bs.pdf", 8, 4)
plt_means_single(df = dat_b, v = "b", type = "RT", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/meanRT_var_mucs.pdf", 8, 4)
plt_means_single(df = dat_muc, v = "mu_c", type = "RT", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/meanRT_var_A1s.pdf", 8, 4)
plt_means_single(df = dat_A1, v = "A1", type = "RT", n_rows = 1, n_cols = 4)
dev.off()
pdf("../../MA/images/meanRT_var_tau1s.pdf", 8, 4)
plt_means_single(df = dat_tau1, v = "tau1", type = "RT", n_rows = 1, n_cols = 5)
dev.off()

# Mean ERs
pdf("../../MA/images/meanER_var_bs.pdf", 8, 4)
plt_means_single(df = dat_b, v = "b", type = "ER", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/meanER_var_mucs.pdf", 8, 4)
plt_means_single(df = dat_muc, v = "mu_c", type = "ER", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/meanER_var_A1s.pdf", 8, 4)
plt_means_single(df = dat_A1, v = "A1", type = "ER", n_rows = 1, n_cols = 4)
dev.off()
pdf("../../MA/images/meanER_var_tau1s.pdf", 8, 4)
plt_means_single(df = dat_tau1, v = "tau1", type = "ER", n_rows = 1, n_cols = 5)
dev.off()

# CAFs
pdf("../../MA/images/CAF_var_bs.pdf", 8, 4)
plt_cafs_single(df = dat_b, n_bins = 4, v = "b", n_rows = 1, n_cols = 7, )
dev.off()
pdf("../../MA/images/CAF_var_mucs.pdf", 8, 4)
plt_cafs_single(df = dat_muc, n_bins = 4, v = "mu_c", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/CAF_var_A1s.pdf", 8, 4)
plt_cafs_single(df = dat_A1, n_bins = 4, v = "A1", n_rows = 1, n_cols = 4)
dev.off()
pdf("../../MA/images/CAF_var_tau1s.pdf", 8, 4)
plt_cafs_single(df = dat_tau1, n_bins = 4, v = "tau1", n_rows = 1, n_cols = 5)
dev.off()

# Delta uncond.
pdf("../../MA/images/Delta_var_bs.pdf", 8, 4)
plt_delta_single(df = dat_b, v = "b", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/Delta_var_mucs.pdf", 8, 4)
plt_delta_single(df = dat_muc, v = "mu_c", n_rows = 1, n_cols = 7)
dev.off()
pdf("../../MA/images/Delta_var_A1s.pdf", 8, 4)
plt_delta_single(df = dat_A1, v = "A1", n_rows = 1, n_cols = 4)
dev.off()
pdf("../../MA/images/Delta_var_tau1s.pdf", 8, 4)
plt_delta_single(df = dat_tau1, v = "tau1", n_rows = 1, n_cols = 5)
dev.off()

# Delta cond.
pdf("../../MA/images/Delta_var_bs.pdf", 8, 4)
plt_delta_single(df = dat_b, v = "b", n_rows = 1, n_cols = 7, conditional = TRUE)
dev.off()
pdf("../../MA/images/Delta_var_mucs.pdf", 8, 4)
plt_delta_single(df = dat_muc, v = "mu_c", n_rows = 1, n_cols = 7, conditional = TRUE)
dev.off()
pdf("../../MA/images/Delta_var_A1s.pdf", 8, 4)
plt_delta_single(df = dat_A1, v = "A1", n_rows = 1, n_cols = 4, conditional = TRUE)
dev.off()
pdf("../../MA/images/Delta_var_tau1s.pdf", 8, 4)
plt_delta_single(df = dat_tau1, v = "tau1", n_rows = 1, n_cols = 5, conditional = TRUE)
dev.off()


## Plots vary A1 x A2 / tau1 x tau2 ############################################

load("out/data/datDDMC1.RData")
dat_nr <- "_6"    # set plot suffix

# get parameter sets with A1 == A2
df_As_eq <- datDDMC[datDDMC$A1 == datDDMC$A2, ]
df_As_eq <- mk_congruency(df_As_eq)   # add congruency column
df_As_lst <- split(df_As_eq, df_As_eq$A1)

# get parameter sets with tau1 == tau2
df_taus_eq <- datDDMC[datDDMC$tau1 == datDDMC$tau2, ]
df_taus_eq <- mk_congruency(df_taus_eq)   # add congruency column
df_taus_lst <- split(df_taus_eq, df_taus_eq$tau1)

## taus = 30
#tmp <- df_taus_lst[[1]]
#pdf("../../MA/images/taus30_sim1.pdf", 6, 6)
#plt_var_As_clean_rt(tmp, corr_only = "T")
#dev.off()
#
## taus = 230
#tmp <- df_taus_lst[[5]]
#pdf("../../MA/images/taus230_sim1.pdf", 6, 6)
#plt_var_As_clean_rt(tmp, corr_only = "T")
#dev.off()

# Mean RTs
pdf(paste0("out/plots/rt/RT_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_means, v1 = "tau1", v2 = "tau2", type = "RT")
dev.off()
pdf(paste0("out/plots/rt/RT_vary_As", dat_nr, ".pdf"))
lapply(df_taus_lst, plt_means, v1 = "A1", v2 = "A2", type = "RT")
dev.off()

# Mean ERs
pdf(paste0("out/plots/er/ER_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_means, v1 = "tau1", v2 = "tau2", type = "ER")
dev.off()
pdf(paste0("out/plots/er/ER_vary_As", dat_nr, ".pdf"))
lapply(df_taus_lst, plt_means, v1 = "A1", v2 = "A2", type = "ER")
dev.off()

# CAFs
pdf(paste0("out/plots/caf/CAF_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_cafs, n_bins = 4, "tau1", "tau2")
dev.off()
pdf(paste0("out/plots/caf/CAF_vary_As", dat_nr, ".pdf"))
lapply(df_taus_lst, plt_cafs, n_bins = 4, "A1", "A2")
dev.off()

# Delta plots
pdf(paste0("out/plots/delta/delta_vary_taus", dat_nr, ".pdf"))
lapply(df_As_lst, plt_delta, "tau1", "tau2")
dev.off()
pdf(paste0("out/plots/delta/delta_vary_As", dat_nr, ".pdf")) 
lapply(df_taus_lst, plt_delta, "A1", "A2")
dev.off()

pdf(paste0("out/plots/delta/delta_vary_taus_cond", dat_nr, ".pdf"))
lapply(df_As_lst, plt_delta, "tau1", "tau2", TRUE)
dev.off()
pdf(paste0("out/plots/delta/delta_vary_As_cond", dat_nr, ".pdf")) 
lapply(df_taus_lst, plt_delta, "A1", "A2", TRUE)
dev.off()

## CDFs
#pdf(paste0("out/plots/CDF_vary_taus", dat_nr, ".pdf"))
#lapply(df_As_lst, plt_cdfs, v1 = "tau1", v2 = "tau2")
#dev.off()
#pdf(paste0("out/plots/CDF_vary_As", dat_nr, ".pdf"))
#lapply(df_taus_lst, plt_cdfs, v1 = "A1", v2 = "A2")
#dev.off()

## Plot densities
#pdf(paste0("out/plots/dens_vary_taus", dat_nr, ".pdf"))
#lapply(df_As_lst, plt_denss, v1 = "tau1", v2 = "tau2")
#dev.off()
#pdf(paste0("out/plots/dens_vary_As", dat_nr, ".pdf"))
#lapply(df_taus_lst, plt_denss, v1 = "A1", v2 = "A2")
#dev.off()