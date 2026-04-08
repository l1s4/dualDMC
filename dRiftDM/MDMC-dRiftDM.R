library(dRiftDM)
library(ggplot2)
source("dRiftDMHelper.R")

mu_mdmc <- function(prms_model, prms_solve, t_vec, one_cond, ddm_opts) {
	muc   = prms_model[["muc"]]
	A1    = prms_model[["A1"]]
	A2    = prms_model[["A2"]]
	tau1  = prms_model[["tau1"]]
	tau2  = prms_model[["tau2"]]

	mu_a1 <- A1 / tau1 * exp(1 - t_vec / tau1) * (1 - t_vec / tau1)
	mu_a2 <- A2 / tau2 * exp(1 - t_vec / tau2) * (1 - t_vec / tau2)

	return (muc + mu_a1 + mu_a2)
}

mdmc_dm <- function(prms_model) {
	
	# define conditions
	conds <- c("compcomp", "compincomp", "incompcomp", "incompincomp")

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

  # set A1 and A2 in relation to compcomp conditions
	# 'A' in incomp is negative of 'A' in comp for both automatic processes
	instructions <- "
	A1 ~ incompcomp == -(A1 ~ compcomp)
	A2 ~ incompcomp == (A2 ~ compcomp)

	A1 ~ compincomp == (A1 ~ compcomp)
	A2 ~ compincomp == -(A2 ~ compcomp)
	
	A1 ~ incompincomp == -(A1 ~ compcomp)
	A2 ~ incompincomp == -(A2 ~ compcomp)
	"

	model <- modify_flex_prms(model, instr = instructions)
}



# different parameter scalings 
prms_ms <- c(muc = 0.5, alpha = 2, b = 50, tau = 30, 
	     non_dec = 330, sd_non_dec = 30, A = 19)
convert_prms(prms_ms, sigma_old = 4, sigma_new = 1, t_from_to = "ms->s")

prms_s <- c(muc = 5.5, alpha = 2, b = 0.4, tau = .05, 
	    non_dec = 0.3, sd_non_dec = 0.03, A = .5)
convert_prms(prms_s, sigma_old = 1, sigma_new = 4, t_from_to = "s->ms")


# define parameters Ulrich2015 Eriksen (1) + Simon (2)
prms_model <- c(
	muc = 5.5, alpha = 2, b = .4, 
	tau1 = .05, tau2 = .05, 
	non_dec = .0, sd_non_dec = .03,
	A1 = .5, A2 = .5
)

mdmc <- mdmc_dm(prms_model)

print(mdmc)
summary(mdmc)
plot(mdmc)

# traces
exp_behavior <- simulate_traces(object = mdmc, k = 1, sigma = 0)
traces <- simulate_traces(object = mdmc, k = 5)

par(mfrow = c(1, 2))
plot(exp_behavior, col = c("green", "red", "blue", "pink"))
plot(traces, col = c("green", "red", "blue", "pink"))
