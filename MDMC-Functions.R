#library(DMCfun)
#dmc = dmcSim(amp = 0, nTrlData = 1, fullData = TRUE)
#plot(dmc, figType = "trials")
#dmcSimApp()

MDMC_AV <- function(N, parameters, automatic1, automatic2) {
  mu_c <- parameters$mu_c
  tau1 <- parameters$tau1
  tau2 <- parameters$tau2
  a1 <- parameters$a1
  a2 <- parameters$a2
  A1 <- parameters$A1
  A2 <- parameters$A2
  
  t <- 1:N
  controlled <- cumsum(rep(mu_c, N))
  
  auto_traj <- function(A, tau, a, t, congruency) {
    sign <- ifelse(congruency == "congruent", 1, -1)
    sign * A * exp(-t / tau) * ((t * exp(1)) / ((a - 1) * tau))^(a - 1)
  }
  
  trajectory_a1 <- auto_traj(A1, tau1, a1, t, automatic1)
  trajectory_a2 <- auto_traj(A2, tau2, a2, t, automatic2)
  
  superimposed <- controlled + trajectory_a1 + trajectory_a2
  
  list(
    superimposed = superimposed,
    trajectory_a1 = trajectory_a1,
    trajectory_a2 = trajectory_a2,
    controlled = controlled
  )
}

MDMC_T <- function(N, parameters) {
  #  Purpose: simulate one DMC trial 
  #  Input: 
  #      - N = max Time
  #      - Parameters: tau1, tau2, a1, a2, A1, A2, mu_c, b, sigma, dt, 
  #	   automatic1, automatic2
  #	 	- automatic1 & automatic2: "congruent" / "incongruent"
  #  Return: X(t) rt, decision, type of first and second automatic process

  mu_c 	<- parameters$mu_c
  tau1 	<- parameters$tau1
  tau2 	<- parameters$tau2
  a1 		<- parameters$a1
  a2 		<- parameters$a2
  sigma <- parameters$sigma
  A1 		<- parameters$A1
  A2 		<- parameters$A2
  dt 		<- parameters$dt
  b 		<- parameters$b
  automatic1 	<- parameters$automatic1
  automatic2 	<- parameters$automatic2
  
  t <- seq(dt, N, by = dt)
  
  auto_drift <- function(A, tau, a, t, congruency) {
    sign <- ifelse(congruency == "congruent", 1, -1)
    sign * A * exp(-t / tau) * ((t * exp(1)) / ((a - 1) * tau))^(a - 1) *
      (((a - 1) / t) - (1 / tau))
  }

  # Drifts of automatic processes
  mu_a1 <- auto_drift(A1, tau1, a1, t, automatic1)
  mu_a2 <- auto_drift(A2, tau2, a2, t, automatic2)
  # Drift of superimposed process
  mu_t <- mu_c + mu_a1 + mu_a2

  # Simulate noisy process
  dX <- mu_t * dt + sigma * sqrt(dt) * rnorm(length(t))
  X <- cumsum(dX)     # X(t)
  
  # Get reaction time and decision
  pb <- suppressWarnings(min(which(X > b), na.rm = TRUE))
  nb <- suppressWarnings(min(which(X < -b), na.rm = TRUE))
  
  if (is.infinite(pb)) pb <- Inf
  if (is.infinite(nb)) nb <- Inf
 	

  if (pb < nb) {
    rt <- pb * dt
    decision <- "pb"
  } else {
    rt <- nb * dt
    decision <- "nb"
  }

  list(
    mut = mu_t, 
    dX = dX, 
    TrajX = X,
    autom1 = automatic1,
    autom2 = automatic2,
    rt = rt,
    dec = decision
  )
}

#test <- MDMC_T(1500, list(sigma = 0.00, dt = 1, mu_c = 0.5, 
#b = 50, tau1 = 30, tau2 = 30, A1 = 20, A2 = 20, a1 = 2, a2 = 2, 
#automatic1 = "congruent", automatic2 = "congruent"))
#print(test$mut)
#print(test$dX)
#print(test$TrajX)
#print(test$rt)
#print(test$dec)

MDMC_Sim <- function(N_sim, N_time, param_grid) {
	# Purpose: simulate multiple DMC-trials per parameter set
	# Input: Number of Simulations N_sim per parameter set in param_grid
	# Return: mean RTs per condition
	dat_all <- list()

	# Loop parameter sets
	for (j in seq_len(nrow(param_grid))) {
		cat("Parameter Set:", j, "\n", unlist(param_grid[j, ]), "\n")

		first_pr <- character(N_sim)
		second_pr <- character(N_sim)
		rt <- numeric(N_sim)
		decision <- character(N_sim)

		for (i in 1:N_sim) {
	
			# randomly choose congruency
			first <- sample(c("congruent", "incongruent"), 1)
			second <- sample(c("congruent", "incongruent"), 1)

			Sim <- MDMC_T(N = N_time, 
				      parameters = as.list(c(param_grid[j, ], 
							     automatic1 = first, 
							     automatic2 = second)
				      )
			)
			first_pr[i]   <- first
			second_pr[i]  <- second
			rt[i]         <- Sim$rt
			decision[i]   <- Sim$dec
		}

		df <- data.frame(first_pr, second_pr, rt, decision)

		df$error <- ifelse(df$decision == "nb", 1, 0)

		df$A1   <- param_grid$A1[j]
		df$A2   <- param_grid$A2[j]
		df$tau1 <- param_grid$tau1[j]
		df$tau2 <- param_grid$tau2[j]
		df$b    <- param_grid$b[j]
		df$mu_c <- param_grid$mu_c[j]

		dat_all[[j]] <- df

	}

	return(dat_all)
}

#tau1_vals <- c(30, 60)
#tau2_vals <- c(160, 180)
#A1_vals <- 0
#A2_vals <- 50
#param_grid <- expand.grid(tau1 = tau1_vals, tau2 = tau2_vals, A1 = A1_vals, A2 = A2_vals)
## constant parameters
#param_grid$a1 <- 2
#param_grid$a2 <- 2
#param_grid$b <- 50
#param_grid$sigma <- 4
#param_grid$dt <- 1
#param_grid$mu_c <- 0.25
#
#
#test <- MDMC_Sim(50, 2000, param_grid)
#str(test)


#Sim <- MDMC_Sim(N_sim = 500, N_time = input$N, param_grid = data.frame(
#mu_c  = 0.4, sigma = 4, tau1  = 100, tau2  = 170, a1    = 2
#a2    = 2, A1    = 30,A2    = 30, dt    = 1, b     = 55))[[1]]
#print(Sim)


MDDM_T <- function(N, parameters) {
	dt 	  <- parameters$dt
	sigma <- parameters$sigma
	mu_c 	<- parameters$mu_c
	b	    <- parameters$b
	P1     <- parameters$P1
	P2     <- parameters$P2

  delta_1 <- parameters$delta_1
  delta_2 <- parameters$delta_2

  automatic1 <- parameters$automatic1
  automatic2 <- parameters$automatic2
  auto1 <- ifelse(automatic1 == "congruent", 1, 0)
  auto2 <- ifelse(automatic2 == "congruent", 1, 0)

	t <- seq(dt, N, by = dt)


	dX_c  <- mu_c * dt + sigma * sqrt(dt) * rnorm(length(t))
	dX_a1 <- (delta_1 * auto1) * dt + sigma * sqrt(dt) * rnorm(length(t))
	dX_a2 <- (delta_2 * auto2) * dt + sigma * sqrt(dt) * rnorm(length(t))
	X_c   <- cumsum(dX_c)     	# X(t) controlled
	X_a1  <- cumsum(dX_a1)     	# X(t) auto 1
	X_a2  <- cumsum(dX_a2)     	# X(t) auto 2

	# Get RT for each process
  
	# lower b by process-specific parameter P
	if(auto1  == 1) {b <- b - P1}
	if(auto2 == 1) {b <- b - P2}
	
	pb_c <- suppressWarnings(min(which(X_c > b), na.rm = TRUE))
	if (is.infinite(pb_c)) pb_c <- Inf

	pb_a1 <- suppressWarnings(min(which(X_a1 > b), na.rm = TRUE))
	if (is.infinite(pb_a1)) pb_a1 <- Inf
	
	pb_a2 <- suppressWarnings(min(which(X_a2 > b), na.rm = TRUE))
	if (is.infinite(pb_a2)) pb_a2 <- Inf

	fastest <- min(pb_c, pb_a1, pb_a2)

	# RT / decision those of fastest process
	rt <- fastest * dt

	decision <- ifelse(fastest == pb_c, "c", 
			   ifelse(fastest == pb_a1, "a1", "a2"))

	list(TrajXc = X_c, 
	     TrajXa1 = X_a1, 
	     TrajXa2 = X_a2, 
	     rt = rt, 
	     dec = decision
	)

}

A <- MDDM_T(500, list(dt = 0.1, mu_c = 0.5, delta_1 = 0.2, delta_2 = 0.8, 
                      sigma = 4, b = 50, P1 = 5, P2 = 5,
                      automatic1 = "congruent", automatic2 = "incongruent"))
print(A)

MDDM_Sim <- function(N_sim, N_time, param_grid) {
  dat_all <- list()
  
  # Loop parameter sets
	for (j in seq_len(nrow(param_grid))) {
		cat("Parameter Set:", j, "\n", unlist(param_grid[j, ]), "\n")

		first_pr <- character(N_sim)
		second_pr <- character(N_sim)
		rt <- numeric(N_sim)
		decision <- character(N_sim)

		for (i in 1:N_sim) {
	
			# randomly choose congruency
			first <- sample(c("congruent", "incongruent"), 1)
			second <- sample(c("congruent", "incongruent"), 1)

			Sim <- MDDM_T(N = N_time, 
				      parameters = as.list(c(param_grid[j, ], 
							     automatic1 = first, 
							     automatic2 = second)
				      )
			)
			first_pr[i] <- first
			second_pr[i] <- second
			rt[i] <- Sim$rt
			decision[i] <- Sim$dec
		}

		df <- data.frame(first_pr, second_pr, rt, decision)
    
    df$error <- ifelse(decision == "a1" & first_pr == "congruent", 0, 
                  ifelse(decision == "a2" & second_pr == "congruent", 0, 
                  ifelse(decision == "c", 0, 1)))

		df$mu_c <- param_grid$mu_c[j]
		df$delta_1 <- param_grid$delta_1[j]
		df$delta_2 <- param_grid$delta_2[j]
		df$b <- param_grid$b[j]

		dat_all[[j]] <- df
	}
	return(dat_all)
}

#delta_1_vals <- c(0.3, 0.6)
#delta_2_vals <- c(0.3, 0.6)
#param_grid <- expand.grid(delta_1 = delta_1_vals, delta_2 = delta_2_vals) 
## constant parameters
#param_grid$b <- 50
#param_grid$sigma <- 4
#param_grid$dt <- 1
#param_grid$mu_c <- 0.25
#test <- MDDM_Sim(50, 1000, param_grid)
#str(test)
#test[[1]]
#
#Sim <- MDMC_Sim(N_sim = 500, N_time = input$N, param_grid = data.frame(
#mu_c  = 0.4, sigma = 4, tau1  = 100, tau2  = 170, a1    = 2
#a2    = 2, A1    = 30,A2    = 30, dt    = 1, b     = 55))[[1]]
#print(Sim)
