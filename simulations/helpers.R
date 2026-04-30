# Create congruency column from auto1 and auto2 columns
mk_congruency <- function(df) {
  # auto1, auto2 congruent if == 1, else incongruent
  df$congruency <- ifelse(df$auto1 == 1 & df$auto2 == 1, "congruent_congruent", 
                          ifelse(df$auto1 == 1 & df$auto2 == -1, "congruent_incongruent", 
                                 ifelse(df$auto1 == -1 & df$auto2 == 1, "incongruent_congruent", 
                                        "incongruent_incongruent")))
  df
}

# mean RT plots
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
         type = "b", main = "Mean RTs", sub = paste("As:", data$A1[1]), 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("tau1", "tau2")), 
         auto.key = list(title = "second automatic process", space = "top", 
                         cex = .7) 
  ) # |> print()
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
         type = "b", main = "Mean RTs", xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("A1", "A2")), 
         auto.key = list(title = "second automatic process", space = "top",
                         cex = .7), 
#         auto.key = list(title = "second automatic process", corner = c(1, 1), 
#                         x = 0.95, y = 0.1, cex = 0.7)
  )# |> print()
}

# mean ER plots
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
         type = "b", main = "Mean ERs", ylim = c(0, 1), 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("tau1", "tau2")), 
         auto.key = list(title = "second automatic process", space = "top", 
                         cex = .7),
  ) #|> print()
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
         type = "b", main = "Mean ERs", ylim = c(0, 1), 
         xlab = "first automatic process", 
         strip = strip.custom(strip.names = c(TRUE, TRUE), 
                              var.name = c("A1", "A2")), 
         auto.key = list(title = "second automatic process", space = "top", 
                         cex = .7),
  ) #|> print()
}



# ECDF plots
plt_cdfs_var_taus <- function(df) {
  df_ecdf <- df %>%
    group_by(tau1, tau2, congruency) %>%
    arrange(rt, .by_group = TRUE) %>%
    mutate(ecdf_y = (row_number() - 1) / (n() - 1)) %>%
    ungroup()

  xyplot(ecdf_y ~ rt | factor(tau1) + factor(tau2), data = df_ecdf,
         groups = congruency, type = "l",
         xlab = "RT", ylab = "ECDF", xlim = c(0, 1200), 
         strip = strip.custom(strip.names = TRUE, var.name = c("tau1", "tau2")),
         auto.key = list(title = "condition", cex = 0.7, columns = 2)
  ) #|> print()
}

plt_cdfs_var_As <- function(df) {
  df_ecdf <- df %>%
    group_by(A1, A2, congruency) %>%
    arrange(rt, .by_group = TRUE) %>%
    mutate(ecdf_y = (row_number() - 1) / (n() - 1)) %>%
    ungroup()
  
  xyplot(ecdf_y ~ rt | factor(A1) + factor(A2), data = df_ecdf,
         groups = congruency, type = "l",
         xlab = "RT", ylab = "ECDF", xlim = c(0, 1200), 
         strip = strip.custom(strip.names = TRUE, var.name = c("A1", "A2")),
         auto.key = list(title = "condition", cex = 0.7, columns = 2)
  ) #|> print()
}

# CAF plots
plt_nxn_cafs_As <- function(df, n_bins) {
  
  df$error <- ifelse(df$dec == 1, 0, 1)
  df$acc <- 1 - df$error
  
  df$tau1 <- factor(df$tau1)
  df$tau2 <- factor(df$tau2)
  
  df_lst <- split(df, list(df$tau1, df$tau2), drop = TRUE)
  
  binned_lst <- lapply(df_lst, function(x) {
    x$rt_bin <- cut(
      x$rt,
      breaks = quantile(x$rt, probs = seq(0, 1, length.out = n_bins + 1)),
      include.lowest = TRUE,
      labels = FALSE
    )
    aggregate(acc ~ rt_bin + congruency, FUN = mean, data = x)
  })
  
  data <- bind_rows(binned_lst, .id = "list_name") %>%
    separate(list_name, into = c("tau1", "tau2"), sep = "\\.", convert = TRUE)
  
  xyplot(
    acc ~ as.numeric(rt_bin) | factor(tau1) + factor(tau2), data = data, 
    groups = congruency, type = "b", main = "CAFs", ylim = c(0, 1.1),
    xlab = "RT bin", ylab = "Accuracy",
    auto.key = list(title = "condition", cex = .7, columns = 2, space = "top"),
    strip = strip.custom(strip.names = TRUE, var.name = c("tau1", "tau2"))
  ) #|> print()
}

plt_nxn_cafs_taus <- function(df, n_bins) {
  
  df$error <- ifelse(df$dec == 1, 0, 1)
  df$acc <- 1 - df$error
  
  df$tau1 <- factor(df$A1)
  df$tau2 <- factor(df$A2)
  
  df_lst <- split(df, list(df$A1, df$A2), drop = TRUE)
  
  binned_lst <- lapply(df_lst, function(x) {
    
    x$rt_bin <- cut(
      x$rt,
      breaks = quantile(x$rt, probs = seq(0, 1, length.out = n_bins + 1)),
      include.lowest = TRUE,
      labels = FALSE
    )
    
    aggregate(acc ~ rt_bin + congruency, FUN = mean, data = x)
  })
  
  data <- bind_rows(binned_lst, .id = "list_name") %>%
    separate(list_name, into = c("A1", "A2"), sep = "\\.", convert = TRUE)
  
  xyplot(
    acc ~ as.numeric(rt_bin) | factor(A1) + factor(A2),
    data = data, groups = congruency, type = "b", main = "CAFs", 
    xlab = "RT bin", ylab = "Accuracy",
    auto.key = list(title = "condition", cex = .7, columns = 2, space = "top"),
    strip = strip.custom(strip.names = TRUE, var.name = c("A1", "A2")),
    #col = c("green", "blue", "purple", "red"),
    ylim = c(0, 1.1)
  )# |> print()
}


# Delta plots DDMC
plt_delta <- function(df, cond1, cond2, cond3, cond4) {
  probs <- seq(0.1, 0.9, by = 0.1)        # quantile probabilities
  
  df1 <- df[df$congruency == cond1, ]
  df2 <- df[df$congruency == cond2, ]
  df3 <- df[df$congruency == cond3, ]
  df4 <- df[df$congruency == cond4, ]
  
  q1 <- quantile(df1$rt, probs = probs)
  q2 <- quantile(df2$rt, probs = probs)
  q3 <- quantile(df3$rt, probs = probs)
  q4 <- quantile(df4$rt, probs = probs)
  
  delta1 <- q1 - q2
  delta2 <- q1 - q3
  delta3 <- q1 - q4
  mean_rt <- (q1 + q2 + q3 + q4) / 4
  
  print(delta1)
  
  
  # Plot
  plot(mean_rt, delta1, type = "n", pch = 16, 
       xlim = c(300, 500), ylim = c(-10, 100), 
       xlab = "Mean RT (ms)", ylab = "delta", main = "Delta Plot")
  points(mean_rt, delta1, type = "b", col = "green")
  points(mean_rt, delta2, type = "b", col = "blue")
  points(mean_rt, delta3, type = "b", col = "red")
}


plt_nxn_delta_taus <- function(df) {
  
  df$tau1 <- factor(df$A1)
  df$tau2 <- factor(df$A2)
  
  probs = seq(0.1, 0.9, by = 0.1)
  
  df_lst <- split(df, list(df$A1, df$A2))
  
  delta_lst <- lapply(df_lst, function(x) {
    
    qs <- tapply(x$rt, x$congruency, quantile, probs = probs)

    q_names <- names(qs)
    ref <- qs[[4]]    # reference category: incongruent_incongruent
   
    deltas <- lapply(qs, function(q) ref - q)
    mean_rt <- Reduce("+", qs) / length(qs)
   
    do.call(rbind, lapply(names(deltas), function(x) {
      data.frame(mean_rt = mean_rt, delta = deltas[[x]], cond = x)
    }))
  })
  
  data <- bind_rows(delta_lst, .id = "list_name") |>
    separate(list_name, into = c("A1", "A2"), sep = "\\.", convert = TRUE)
  
  # drop reference category
  data <- data[data$cond != "incongruent_incongruent", ]
  
  xyplot(
    delta ~ mean_rt | factor(A1) + factor(A2), groups = cond, data = data,
    type = "l", main = "Delta Plots", xlim = c(100, 700), ylim = c(-50, 200), 
    xlab = "Mean RT (ms)", ylab = "Delta (ms)",
    auto.key = list(title = "condition", cex = 0.7, columns = 2, space = "top"),
    strip = strip.custom(strip.names = TRUE, var.name = c("A1", "A2"))
  )# |> print()
}

plt_nxn_delta_As <- function(df) {
  
  df$tau1 <- factor(df$tau1)
  df$tau2 <- factor(df$tau2)
  
  probs = seq(0.1, 0.9, by = 0.1)
  
  df_lst <- split(df, list(df$tau1, df$tau2))
  
  delta_lst <- lapply(df_lst, function(x) {
    
    qs <- tapply(x$rt, x$congruency, quantile, probs = probs)
    
    q_names <- names(qs)
    ref <- qs[[4]]    # reference category: incongruent_incongruent
    
    deltas <- lapply(qs, function(q) ref - q)
    mean_rt <- Reduce("+", qs) / length(qs)
    
    do.call(rbind, lapply(names(deltas), function(x) {
      data.frame(mean_rt = mean_rt, delta = deltas[[x]], cond = x)
    }))
  })
  
  data <- bind_rows(delta_lst, .id = "list_name") |>
    separate(list_name, into = c("tau1", "tau2"), sep = "\\.", convert = TRUE)
  
  # drop reference category
  data <- data[data$cond != "incongruent_incongruent", ]
  
  xyplot(
    delta ~ mean_rt | factor(tau1) + factor(tau2), groups = cond, data = data,
    type = "l", main = "Delta Plots", xlim = c(100, 700), ylim = c(-50, 200), 
    xlab = "Mean RT (ms)", ylab = "Delta (ms)",
    auto.key = list(title = "condition", cex = 0.7, columns = 2, space = "top"),
    strip = strip.custom(strip.names = TRUE, var.name = c("tau1", "tau2"))
  )# |> print()
}