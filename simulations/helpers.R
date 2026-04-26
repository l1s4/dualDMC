plt_var_taus_rt <- function(data, corr_only) {
  # Input: 
  #   - data: dataframe containing columns auto1, auto2 (-1 incongruent, 1 congruent), 
  #           A1, A2, tau1, tau2, rt, dec (1 correct, -1 incorrect). 
  #           data aggregated across A1, A2 (i.e. A1 == A2 is expected in data)
  #   - corr_only: bool, if TRUE only correct trials are used
  # Output: lattice plot of mean RTs for combinations of tau1 x tau2
  
  if (corr_only) {data <- data[data$dec == 1, ]}    # correct trials only
  means <- aggregate(rt ~ auto1 + auto2 + tau1 + tau2, FUN = mean, 
                     data = data)
  means$auto1 <- factor(means$auto1, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  means$auto2 <- factor(means$auto2, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  
  xyplot(rt ~ auto1 | factor(tau1) * factor(tau2), groups = auto2, data = means, 
         type = c("p", "a"), main = "Mean RTs", 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("tau1", "tau2")), 
         auto.key = list(title = "second automatic process", space = "top", 
                         cex = .7),
  ) |> print()
}

plt_var_As_rt <- function(data, corr_only) {
  # Input: 
  #   - data: dataframe containing columns auto1, auto2 (-1 incongruent, 1 congruent), 
  #           A1, A2, tau1, tau2, rt, dec (1 correct, -1 incorrect). 
  #           data aggregated across tau1, tau2 (i.e. tau1 == tau2 is expected in data)
  #   - corr_only: bool, if TRUE only correct trials are used
  # Output: lattice plot of mean RTs for combinations of A1 x A2
  
  if (corr_only) {data <- data[data$dec == 1, ]}    # correct trials only
  means <- aggregate(rt ~ auto1 + auto2 + A1 + A2, FUN = mean, 
                     data = data)
  means$auto1 <- factor(means$auto1, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  means$auto2 <- factor(means$auto2, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  
  xyplot(rt ~ auto1 | factor(A1) * factor(A2), groups = auto2, data = means, 
         type = c("p", "a"), main = "Mean RTs", 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("A1", "A2")), 
         auto.key = list(title = "second automatic process", space = "top",
                         cex = .7), 
#         auto.key = list(title = "second automatic process", corner = c(1, 1), 
#                         x = 0.95, y = 0.1, cex = 0.7)
  ) |> print()
}

plt_var_taus_er <- function(data) {
  # Input: 
  #   - data: dataframe containing columns auto1, auto2 (-1 incongruent, 1 congruent), 
  #           A1, A2, tau1, tau2, rt, dec (1 correct, -1 incorrect). 
  #           data aggregated across A1, A2 (i.e. A1 == A2 is expected in data)
  #   - corr_only: bool, if TRUE only correct trials are used
  # Output: lattice plot of mean ERs for combinations of tau1 x tau2
  data$error <- ifelse(data$dec == -1, 1, 0)
  means <- aggregate(error ~ auto1 + auto2 + tau1 + tau2, FUN = mean, 
                     data = data)
  means$auto1 <- factor(means$auto1, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  means$auto2 <- factor(means$auto2, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  
  xyplot(error ~ auto1 | factor(tau1) * factor(tau2), groups = auto2, data = means, 
         type = c("p", "a"), main = "Mean ERs", ylim = c(0, 1), 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("tau1", "tau2")), 
         auto.key = list(title = "second automatic process", space = "top", 
                         cex = .7),
  ) |> print()
}

plt_var_As_er <- function(data) {
  # Input: 
  #   - data: dataframe containing columns auto1, auto2 (-1 incongruent, 1 congruent), 
  #           A1, A2, tau1, tau2, rt, dec (1 correct, -1 incorrect). 
  #           data aggregated across tau1, tau2 (i.e. tau1 == tau2 is expected in data)
  #   - corr_only: bool, if TRUE only correct trials are used
  # Output: lattice plot of mean ERs for combinations of A1 x A2
  data$error <- ifelse(data$dec == -1, 1, 0)
  means <- aggregate(error ~ auto1 + auto2 + A1 + A2, FUN = mean, 
                     data = data)
  means$auto1 <- factor(means$auto1, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  means$auto2 <- factor(means$auto2, levels = c(1, -1), 
                        labels = c("congruent", "incongruent"))
  
  xyplot(error ~ auto1 | factor(A1) * factor(A2), groups = auto2, data = means, 
         type = c("p", "a"), main = "Mean ERs", ylim = c(0, 1), 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("A1", "A2")), 
         auto.key = list(title = "second automatic process", space = "top", 
                         cex = .7),
  ) |> print()
}

# Create congruency column from auto1 and auto2 columns
mk_congruency <- function(df) {
  # auto1, auto2 congruent if == 1, else incongruent
  df$congruency <- ifelse(df$auto1 == 1 & df$auto2 == 1, "congruent_congruent", 
                   ifelse(df$auto1 == 1 & df$auto2 == -1, "congruent_incongruent", 
                   ifelse(df$auto1 == -1 & df$auto2 == 1, "incongruent_congruent", 
                   "incongruent_incongruent")))
  df
}

# ECDF plot
plt_cdfs <- function(df) {
  # Input: dataframe containing reaction times in different congruency 
  #        conditions
  # Output: plot of ECDF per congruency condition
  ggplot(df, aes(rt, colour = congruency)) + 
    stat_ecdf(geom = "point", size = 0.7) + 
    theme_classic()
}

# Delta plots MDMC
plt_delta <- function(df, cond1, cond2) {
  # cond1, cond2: congruency conditions to be used
  probs <- seq(0.1, 0.9, by = 0.1)        # quantile probabilities
  
  df1 <- df[df$congruency == cond1, ]
  df2 <- df[df$congruency == cond2, ]
  
  q1 <- quantile(df1$rt, probs = probs)
  q2 <- quantile(df2$rt, probs = probs)
  
  delta <- q1 - q2
  mean_rt <- (q1 + q2) / 2
  
  # Plot
  plot(mean_rt, delta, type = "b", pch = 16, 
       xlab = "Mean RT (ms)", ylab = "delta", main = "Delta Plot")
}