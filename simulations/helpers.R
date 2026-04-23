# Get mean RTs and ERs for each parameter set for MDMC and MDDM
get_mean_dvs_DDMC <- function(data, corr_only){
  # get mean RT and mean ER values for every parameter combination across 
  # conditions
  # data: list of dataframes containing simulated trials for each combination
  # corr_only: bool, if TRUE only correct trials are included in mean rts
  meanrts <- list()
  meaners <- list()
  for (i in 1:length(data)) {
    d <- data[[i]]
    
    A1    <- d$A1[1]
    A2    <- d$A2[1]
    tau1  <- d$tau1[1]
    tau2  <- d$tau2[1]
    b     <- d$b[1]
    mu_c  <- d$mu_c[1]
    
    d_rt <- d     # include errors
    if (corr_only) d_rt <- d[d$decision == "pb", ]  # correct trials only for rt
    
    # correct trials only for rt
    dfrt <- aggregate(rt ~ first_pr + second_pr, FUN = mean, data = d_rt)
    dfrt$A1   <- A1
    dfrt$A2   <- A2
    dfrt$tau1 <- tau1
    dfrt$tau2 <- tau2
    dfrt$b    <- b
    dfrt$mu_c <- mu_c
    
    dfer <- aggregate(error ~ first_pr + second_pr, FUN = mean, data = d)
    dfer$A1   <- A1
    dfer$A2   <- A2
    dfer$tau1 <- tau1
    dfer$tau2 <- tau2
    dfer$b    <- b
    dfer$mu_c <- mu_c
    
    meanrts[[i]] <- dfrt
    meaners[[i]] <- dfer
  }
  
  out <- list(df_rts = bind_rows(meanrts), df_ers = bind_rows(meaners))
  return (out)
}

get_mean_dvs_MDDM <- function(data, corr_only){
  # get mean RT and mean ER values for every parameter combination across 
  # conditions
  # data: list of dataframes containing simulated trials for each combination
  # corr_only: bool, if TRUE only correct trials are included in mean rts
  meanrts <- list()
  meaners <- list()
  for (i in 1:length(data)) {
    d <- data[[i]]
    
    delta1  <- d$delta1[1]
    delta2  <- d$delta2[1]
    P1      <- d$P1[1]
    P2      <- d$P2[1]
    b       <- d$b[1]
    mu_c    <- d$mu_c[1]
    
    d_rt <- d     # include errors
    if (corr_only) d_rt <- d[d$error == 0, ]  # correct trials only for rt
    
    # correct trials only for rt
    dfrt <- aggregate(rt ~ first_pr + second_pr, FUN = mean, data = d_rt)
    dfrt$P1   <- P1
    dfrt$P2   <- P2
    dfrt$delta1 <- delta1
    dfrt$delta2 <- delta2
    dfrt$b    <- b
    dfrt$mu_c <- mu_c
    
    dfer <- aggregate(error ~ first_pr + second_pr, FUN = mean, data = d)
    dfer$P1   <- P1
    dfer$P2   <- P2
    dfer$delta1 <- delta1
    dfer$delta2 <- delta2
    dfer$b    <- b
    dfer$mu_c <- mu_c
    
    meanrts[[i]] <- dfrt
    meaners[[i]] <- dfer
  }
  
  out <- list(df_rts = bind_rows(meanrts), df_ers = bind_rows(meaners))
  return (out)
}


# MDMC plot functions
plt_A_tau <- function(df, type) {
  # df: constant A, varying combinations of taus
  ymin <- 0
  ymax <- ifelse(type == "ER", 1, max(df$rt) + 200)
  
  if(type == "ER") {
    ymax = 1
    xyplot(error ~ first_pr | as.factor(tau1) * as.factor(tau2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           main = bquote(atop("Mean ER", scriptstyle(paste("A1 = A2 = ", .(df$A1[1]))))),
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("tau1", "tau2")),
           auto.key = list(title = "second automatic process",
                           space = "right", cex = .7)
    )
  } else {
    ymax = max(df$rt) + 200
    xyplot(rt ~ first_pr | as.factor(tau1) * as.factor(tau2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           main = bquote(atop(paste("Mean ", .(type)), 
                         scriptstyle(paste("A1 = A2 = ", .(df$A1[1])))
           )),
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("tau1", "tau2")),
           auto.key = list(title = "second automatic process",
                           space = "right", cex = .7)
    )
  }
}

plt_tau_A <- function(df, type) {
  # df: constant tau, varying combinations of As
  ymin <- 0
  ymax <- ifelse(type == "ER", 1, max(df$rt) + 200)
  
  if (type == "ER") {
    ymax <- 1
    xyplot(error ~ first_pr | as.factor(A1) * as.factor(A2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           main = bquote(atop("Mean ER", scriptstyle(paste("tau1 = tau2 = ", .(df$A1[1]))))),
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("A1", "A2")),
           auto.key = list(title = "second automatic process"),
           space = "right", cex = .7)
  } else {
    xyplot(rt ~ first_pr | as.factor(A1) * as.factor(A2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           main = bquote(atop("Mean RT", scriptstyle(paste("tau1 = tau2 = ", .(df$A1[1]))))),
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("A1", "A2")),
           auto.key = list(title = "second automatic process"),
           space = "right", cex = .7)
  }
}


# MRDM plot functions
plt_delta_P <- function(df, type) {
  # df: constant delta, varying combinations of Ps
  ymin <- 0
  ymax <- ifelse(type == "ER", 1, max(df$rt) + 200)
  
  if(type == "ER") {
    ymax = 1
    xyplot(error ~ first_pr | as.factor(P1) * as.factor(P2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           main = bquote(atop("Mean ER", scriptstyle(paste("delta1 = delta2 = ", .(df$delta1[1]))))),
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("P1", "P2")),
           auto.key = list(title = "second automatic process"),
           space = "right", cex = .7)
    
  } else {
    ymax <- max(df$rt) + 200
    xyplot(rt ~ first_pr | as.factor(P1) * as.factor(P2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           main = bquote(atop("Mean RT", scriptstyle(paste("delta1 = delta2 = ", .(df$delta1[1]))))),
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("P1", "P2")),
           auto.key = list(title = "second automatic process"),
           space = "right", cex = .7)
  }
}

plt_P_delta <- function(df, type) {
  # df: constant P, varying combinations of deltas
  ymin <- 0
  
  if(type == "ER") {
    ymax <- 1
    xyplot(error ~ first_pr | as.factor(delta1) * as.factor(delta2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           main = bquote(atop(paste("Mean ", .(type)), scriptstyle(paste("P1 = P2 = ", .(df$P1[1]))))),
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("delta1", "delta2")),
           auto.key = list(title = "second automatic process"),
           space = "right", cex = .7)
  } else {
    ymax <- max(df$rt) + 200
    xyplot(rt ~ first_pr | as.factor(delta1) * as.factor(delta2), 
           groups = second_pr, data = df, type = c("p", "a"), 
           xlab = "first automatic process", ylab = type, ylim = c(ymin, ymax), 
           main = bquote(atop(paste("Mean ", .(type)), scriptstyle(paste("P1 = P2 = ", .(df$P1[1]))))),
           strip = strip.custom(strip.names = c(TRUE, TRUE), 
                                var.name = c("delta1", "delta2")),
           auto.key = list(title = "second automatic process"),
           space = "right", cex = .7)
  }
}


# MDMC plot ECDFs
plt_cdfs_mdmc <- function(df) {
  ggplot(df, aes(rt, colour = congruency)) + 
    stat_ecdf(geom = "point") + theme_classic()
}

# Delta plots MDMC
plt_delta_mdmc <- function(df1, df2) {
  probs <- seq(0.1, 0.9, by = 0.1)      # quantile probabilities
  q1 <- quantile(df1$rt, probs = probs)
  q2 <- quantile(df2$rt, probs = probs)
  delta <- q1 - q2
  mean_rt <- (q1 + q2) / 2
  
  # Plot
  plot(mean_rt, delta, type = "b", pch = 16, 
       xlab = "Mean RT (ms)", ylab = "delta", main = "Delta Plot")
}