#' ---
#' title: "MDMC-dRiftDM"
#' author: "Lisa Fischer"
#' output: "pdf_document"
#' ---

options(browser = "/usr/bin/firefox")
library(dRiftDM)
source("dRiftDMHelper.R")

mu_mdmc <- function(prms_model, prms_solve, t_vec, one_cond, ddm_opts) {
	muc = prms_model[["muc"]]

	A1 = prms_model[["A1"]]
	A2 = prms_model[["A2"]]
	tau1 = prms_model[["tau1"]]
	tau2 = prms_model[["tau2"]]

	mu_a1 <- A1 / tau1 * exp(1 - t_vec / tau1) * (1 - t_vec / tau1)
	mu_a2 <- A2 / tau2 * exp(1 - t_vec / tau2) * (1 - t_vec / tau2)
	return (muc + mu_a1 + mu_a2)
}

mdmc_dm <- function() {
	# define parameters
	prms_model <- c(
		muc = 4, alpha = 4, b = .6, 
		tau1 = .02, tau2 = .17, 
		non_dec = .3, sd_non_dec = .02,
		A1 = .15, A2 = .45 
	)
	
	# define conditions
	conds <- c("compcomp", "compincomp", "incompcomp", "incompincomp")
#	conds <- c("comp", "incomp")

	# get access to pre-built component functions
	comps <- component_shelf()

	# call drift_dm() function 
	model <- drift_dm(
		prms_model = prms_model, 
		conds = conds, 
		subclass = "mdmc", 
		mu_fun = mu_mdmc, 
		x_fun = comps$x_beta, 
		b_fun = comps$b_constant, 
		dt_b_fun = comps$dt_b_fun,
		nt_fun = comps$nt_truncated_normal
	)

	# 'A' in incomp is negative of 'A' in comp for both automatic processes
	instructions <- "
	A1 ~ incompcomp == -(A1 ~ compcomp)
	A2 ~ incompcomp == (A2 ~ compcomp)

	A1 ~ compincomp == (A1 ~ compcomp)
	A2 ~ compincomp == -(A2 ~ compcomp)
	
	A1 ~ incompincomp == -(A1 ~ compcomp)
	A2 ~ incompincomp == -(A2 ~ compcomp)
	"
#	instructions <- "
#	A1 ~ incomp == -(A1 ~ comp)
#	"

	model <- modify_flex_prms(model, instr = instructions)
}

mdmc <- mdmc_dm()


### Parameter recovery ########################################################
#prms_solve(mdmc)["dx"] <- .001
#prms_solve(mdmc)["dt"] <- .001

# different parameter representations: convert ms to s
prms_ms <- c(muc = 0.7, alpha = 0, b = 80, tau = 20, non_dec = 330, 
	     sd_non_dec = 30, A = 50)
convert_prms(prms_ms, sigma_old = 4, sigma_new = 1, t_from_to = "ms->s")

lower_sim_bnd <- c(
  muc = 3,
  b = 0.35,
  non_dec = 0.275,
  sd_non_dec = .01,
  tau1 = 0.025,
  tau2 = 0.125,
  A1 = 0.05,
  A2 = 0.15,
  alpha = 2.5
)
upper_sim_bnd <- c(
  muc = 6.5,
  b = 0.70,
  non_dec = 0.425,
  sd_non_dec = .04,
  tau1 = 0.200,
  tau2 = 0.300,
  A1 = 0.25,
  A2 = 0.35,
  alpha = 6
)

data_prms <- simulate_data(mdmc, n = 500, k = 100, lower = lower_sim_bnd, upper = upper_sim_bnd)
synth_data <- data_prms$synth_data

# Fit model
lower_bnd <- c(
	muc = 0.50, 
  b = .15,  
	non_dec = .15, 
  sd_non_dec = .0075,
	tau1 = .015, 
  tau2 = .015, 
	A1 = .005, 
  A2 = .005,
  alpha = 2
)
upper_bnd <- c(
	muc = 9.00, 
  b = 1.20,  
	non_dec = .60, 
  sd_non_dec = .10,
	tau1 = .25, 
  tau2 = .25, 
	A1 = .30, 
  A2 = .30,
  alpha = 8.00
)

fit_result <- estimate_dm(
	drift_dm_obj = mdmc, 
	obs_data = synth_data, 
	optimizer = "DEoptim", 
	n_cores = 5, 
	lower = lower_bnd,
	upper = upper_bnd, 
	control = list(trace = TRUE, itermax = 2000, steptol = 400), 
	approach = "agg_c"		# Aggregated data
)

save(fit_result, file="fit_result1203")

fit_stats <- calc_stats(object = fit_result, level = "group", type = c("fit_stats"))
fit_stats

# quantiles, cafs
check_fit_q <- calc_stats(object = fit_result, level = "group", resample = FALSE, type = c("quantiles"))
check_fit_c <- calc_stats(object = fit_result, level = "group", resample = FALSE, type = c("cafs"))
plot(check_fit_q)
plot(check_fit_c)

check_fit_qr <- calc_stats(object = fit_result, level = "group", resample = TRUE, type = c("quantiles"))
check_fit_cr <- calc_stats(object = fit_result, level = "group", resample = TRUE, type = c("cafs"))
plot(check_fit_qr)
plot(check_fit_cr)

# delta functions
check_fit_delta_cc <- calc_stats(object = fit_result, level = "group", resample = FALSE, type = c("delta_funs"), minuends = rep("compcomp", 3), subtrahends = c("compincomp", "incompcomp", "incompincomp"))
check_fit_delta_ci <- calc_stats(object = fit_result, level = "group", resample = FALSE, type = c("delta_funs"), minuends = rep("compincomp", 2), subtrahends = c("incompcomp", "incompincomp"))
check_fit_delta_ic <- calc_stats(object = fit_result, level = "group", resample = FALSE, type = c("delta_funs"), minuends = rep("incompcomp"), subtrahends = c("incompincomp"))

plot(check_fit_delta_cc)
plot(check_fit_delta_ci)
plot(check_fit_delta_ic)

# densities
max_rt <- fit_result$obs_data_ids$RT |> max()
check_fit_dens <- calc_stats(object = fit_result, level = "group", resample = FALSE, type = "densities", t_max = max_rt + 0.01)
plot(check_fit_dens)

par(mfrow = c(2, 2))
plot(check_fit_dens, conds = "compcomp")
plot(check_fit_dens, conds = "compincomp")
plot(check_fit_dens, conds = "incompcomp")
plot(check_fit_dens, conds = "incompincomp")

par(mfrow = c(1, 1))


## Individual data

#data_prms <- simulate_data(mdmc, n = 250, k = 5, lower = lower_sim_bnd, upper = upper_sim_bnd)
#synth_data <- data_prms$synth_data
#
#fit_result_ind <- estimate_dm(
#	drift_dm_obj = mdmc, 
#	obs_data = synth_data, 
#	optimizer = "DEoptim", 
#	n_cores = 5, 
#	lower = lower_sim_bnd,
#	upper = upper_sim_bnd, 
#	control = list(trace = TRUE), 
#	verbose = 1
#)
#
#
## Correlations
#prms <- c("muc", "b", "non_dec", "sd_non_dec", "tau1", "tau2", "A1", "A2", "alpha")
#recov_prms <- coef(fit_result)
#orig_prms <- data_prms$prms
#
#recov_prms
#orig_prms
#
#sapply(prms, function(one_prm) {
#	cor_val <- cor(recov_prms[[one_prm]], orig_prms[[one_prm]])
#	round(cor_val, 3)
#})
#
#check_fit_ind <- calc_stats(object = fit_result_ind, level = "group", resample = TRUE, type = c("quantiles", "cafs"))
#plot(check_fit_ind)
