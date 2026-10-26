library(docstring)

mk_congruency <- function(df) {
  #' Add congruency column to df
  #' @param df data frame. Data frame to create congruency column for. 
  #' Must include columns named auto1 (-1 incongruent, 1 congruent) and 
  #' auto2 (-1 incongruent, 1 congruent)
  #' @return modified dataframe including congruency column 'congruency' 
  
  df$congruency <- ifelse(
    df$auto1 == 1 & df$auto2 == 1, "congruent_congruent", 
      ifelse(df$auto1 == 1 & df$auto2 == -1, "congruent_incongruent", 
        ifelse(df$auto1 == -1 & df$auto2 == 1, "incongruent_congruent", 
          "incongruent_incongruent")))
  df
}

mk_subtitle <- function(df) {
  #' Generate subtitle containing parameter values 
  #' @param df data frame. Data frame to create subtitle for
  #' @return string containing unique values for different parameters in df

  A1_values <- paste(sort(unique(df$A1)), collapse = " / ")
  A2_values <- paste(sort(unique(df$A2)), collapse = " / ")
  tau1_values <- paste(sort(unique(df$tau1)), collapse = " / ")
  tau2_values <- paste(sort(unique(df$tau2)), collapse = " / ")
  b_values <- paste(sort(unique(df$b)), collapse = " / ")
  muc_values <- paste(sort(unique(df$mu_c)), collapse = " / ")
  sigma_values <- paste(sort(unique(df$sigma)), collapse = " / ")
  dt_values <- paste(sort(unique(df$dt)), collapse = " / ")
  ndtm_values <- paste(sort(unique(df$ndt_m)), collapse = " / ")
  ndtsd_values <- paste(sort(unique(df$ndt_sd)), collapse = " / ")
  
  subtitle <- bquote(
    b == .(b_values) * "," ~
    mu[c] == .(muc_values) * "," ~ 
    A[1] == .(A1_values) * "," ~ 
    A[2] == .(A2_values) * "," ~ 
    tau[1] == .(tau1_values) * "," ~ 
    tau[2] == .(tau2_values) * "," ~ 
    sigma == .(sigma_values) * "," ~ 
    dt == .(dt_values) * "," ~ 
    ndt[M] == .(ndtm_values) * "," ~ 
    ndt[SD] == .(ndtsd_values)
  )
  
  subtitle
}

mk_label <- function(v) {
  switch (v,
          "b" = "b", 
          "mu_c" = expression(mu[c]), 
          "tau1" = expression(tau[1]), 
          "tau2" = expression(tau[2]), 
          "A1" = expression(A[1]), 
          "A2" = expression(A[2]), 
          "ndt_m" = expression(ndt[M]), 
          "ndt_sd" = expression(ndt[SD]), 
          v
  )
}

plt_means <- function(df, v1, v2 = NA, type, corr_only = TRUE, subt = TRUE, 
                      title = TRUE, ylim_custom = NA, legend_pos = "top", 
                      xlab_name = "First Automatic Process", 
                      lines_name = "Second Automatic Process", 
                      small_panel_text = FALSE, auto_key = TRUE) {
  #' Lattice interaction plots for means per condition
  #' @param df data frame. Must include columns auto1, auto2 (-1 incongruent, 
  #' 1 congruent), v1, v2, rt, dec (1 correct, -1 incorrect). 
  #' @param v1 character. First parameter to aggregate over
  #' @param v2 character. Second parameter to aggregate over
  #' @param type character. Determines dependent variable, either 'RT' or 'ER'
  #' @param corr_only boolean. If TRUE, only correct trials are used 
  #' for type = 'RT'
  #' @param subt boolean. If TRUE, a subtitle with parameter values is added
  #' @param title boolean. If TRUE, a title is added
  #' @param ylim_custom vector, c(ymin, ymax). To set custom y-axis limits
  #' @param xlab_name character. Name for the x-axis.
  #' @param lines_name character. Name for the grouping-factor
  #' @param small_panel_text boolean. If TRUE, the panel headers are .7
  #' @param legend_pos position of the legend. Must be lattice-compatible
  
  if (type == "RT") {
    if (corr_only) {df <- df[df$dec == 1, ]}    # correct trials only
    
    if (!is.na(v2)) {
      form <- as.formula(paste("rt ~ auto1 + auto2 +", v1, "+", v2))
      plot_form <- as.formula(
        paste0("rt ~ auto1 | factor(", v1, ") + factor(", v2, ")")
      )
    } else {
      form <- as.formula(paste("rt ~ auto1 + auto2 +", v1))
      plot_form <- as.formula(paste0("rt ~ auto1 | factor(", v1, ")"))
    }
    ylims <- if (any(is.na(ylim_custom))) c(360, 550) else ylim_custom
  } else {
    df$error <- ifelse(df$dec == -1, 1, 0)
    if (!is.na(v2)){
      form <- as.formula(paste("error ~ auto1 + auto2 +", v1, "+", v2))
      plot_form <- as.formula(
        paste0("error ~ auto1 | factor(", v1, ") + factor(", v2, ")")
      )
    } else {
      form <- as.formula(paste("error ~ auto1 + auto2 +", v1))
      plot_form <- as.formula(paste0("error ~ auto1 | factor(", v1, ")"))
    }
    ylims <- if (any(is.na(ylim_custom))) c(0, 1) else ylim_custom
  }
  
  means <- aggregate(form, FUN = mean, data = df)
  means$auto1 <- factor(means$auto1, levels = c(1, -1), 
                        labels = c("Cong", "Inc"))
  means$auto2 <- factor(means$auto2, levels = c(1, -1), 
                        labels = c("Cong", "Inc"))
  
  label1 <- mk_label(v1)
  label2 <- mk_label(v2)
  subtitle <- mk_subtitle(df)
  
  lattice::xyplot(
    plot_form, groups = auto2, data = means, type = c("b", "g"), 
    ylim = ylims, ylab = type, xlab = xlab_name, 
    main = if(title) {paste0("Mean ", type, "s")}, 
    sub = if(subt) {list(subtitle, cex = .7)}, 
    par.strip.text = if (small_panel_text) {par.strip.text = list(cex = .7)}, 
    strip = strip.custom(strip.names = TRUE, sep = " = ", var.name = c(label1, label2)), 
    auto.key = if (auto_key) {list(title = lines_name, space = legend_pos, cex = .7)}
  )
}

plt_cafs <- function(df, v1, v2 = NA, n_bins = 5, subt = TRUE, title = TRUE, 
                     auto_key = TRUE) {
  #' Lattice plot for CAFs per condition
  #' @param df data frame. Must include columns auto1, auto2 (-1 incongruent, 
  #' 1 congruent), v1, v2, rt, dec (1 correct, -1 incorrect). 
  #' @param n_bins integer. Number of RT bins
  #' @param v1 character. First parameter to aggregate over
  #' @param v2 character. Second parameter to aggregate over
  #' @param subt boolean. If TRUE, a subtitle with parameter values is added
  #' @param title boolean. If TRUE, a 'CAFs' is added as a title
  
  df$error <- ifelse(df$dec == 1, 0, 1)
  df$acc <- 1 - df$error
  
  if (!is.na(v2)) {
    caf_data <- df %>% 
      group_by(.data[[v1]], .data[[v2]], congruency) %>% 
      mutate(rt_bin = ntile(rt, n_bins)) %>%
      group_by(.data[[v1]], .data[[v2]], congruency, rt_bin) %>% 
      summarise(mean_acc = mean(acc), mean_rt = mean(rt), .groups = "drop")
    plot_form <- as.formula(
      paste0("mean_acc ~ rt_bin | factor(", v1, ") + factor(", v2, ")"))
  } else {
    caf_data <- df %>% 
      group_by(.data[[v1]], congruency) %>% 
      mutate(rt_bin = ntile(rt, n_bins)) %>%
      group_by(.data[[v1]], congruency, rt_bin) %>% 
      summarise(mean_acc = mean(acc), mean_rt = mean(rt), .groups = "drop")
    plot_form <- as.formula(
      paste0("mean_acc ~ rt_bin | factor(", v1, ")"))
  }
  
  label1 <- mk_label(v1) 
  label2 <- mk_label(v2)
  subtitle <- mk_subtitle(df)
  lattice::xyplot(
    plot_form, groups = congruency, data = caf_data, type = c("b", "g"), 
    ylim = c(0, 1.1), xlab = "RT bin", ylab = "Accuracy",
    sub = if(subt) {list(subtitle, cex = .7)}, main = if(title) {"CAFs"},
    strip = strip.custom(strip.names = TRUE, sep = " = ", 
                         var.name = c(label1, label2)),
    auto.key = if (auto_key) {list(title = "Condition", cex = .7, columns = 2, 
                    space = "top")}, 
    scales = list(x = list(at = seq(1, n_bins, by=1)))
  )
}



# Delta plots 
plt_delta <- function(df, v1, v2 = NA, conditional = FALSE, subt = TRUE, 
                      title = TRUE, corr_only = TRUE) {
  #' Lattice plot for delta functions of the first process
  #' Calculated as difference between (auto1 == 1) - (auto1 == -1)
  #' @param df data frame. Must include columns auto1, auto2 (-1 incongruent, 
  #' 1 congruent), v1, v2, rt, dec (1 correct, -1 incorrect). 
  #' @param v1 character. First parameter to aggregate over
  #' @param v2 character. Second parameter to aggregate over
  #' @param conditional boolean. If TRUE, two delta functions will be 
  #' calculated. One is for auto2 == -1, one for auto2 == 1
  #' @param subt boolean. If TRUE, a subtitle with parameter values is added
  #' @param title boolean. If TRUE, a 'CAFs' is added as a title
  #' @param corr_only boolean. If TRUE, only correct trials are used 
  
  if (corr_only) {df <- df[df$dec == 1, ]}    # correct trials only
  probs <- seq(0.1, 0.9, by = 0.1)
  
  if(conditional) {
    df_lst <- split(df, list(df[[v1]], df[[v2]], df$auto2), sep = ":")
  } else {
    df_lst <- split(df, list(df[[v1]], df[[v2]]), sep = ":")
  }
  
  delta_dat <- do.call(rbind, lapply(names(df_lst), function(nm) {
    df <- df_lst[[nm]]
    q_pos <- quantile(df$rt[df$auto1 == 1], probs = probs)
    q_neg <- quantile(df$rt[df$auto1 == -1], probs = probs)
    parts <- strsplit(nm, ":")[[1]]
    v1 <- as.numeric(parts[1])
    v2 <- as.numeric(parts[2])
    if (conditional) {auto2 <- as.numeric(parts[3])}
    data.frame(
      v1 = v1, v2 = v2, mean_rt = (q_pos + q_neg) / 2,
      delta = q_neg - q_pos, bin = seq_along(probs), 
      auto2 = if (conditional) {auto2} else {NA}
    )
  }))
  subtitle <- mk_subtitle(df)
  
  if (all(df[[v1]] == df[[v2]])) {
    label <- switch (
      v1, 
      "tau1" = expression(tau[1] == tau[2]), 
      "A1" = expression(A[1] == A[2]), 
      v1
    )
    n_cols <- length(unique(delta_dat$v1))
    
    lattice::xyplot(
      delta ~ mean_rt | factor(v1), data = delta_dat, 
      groups = if (conditional) {auto2}, main = if(title) {"Delta Plots"}, 
      sub = if(subt) {list(subtitle, cex = .7)}, type = c("b", "g"), 
      xlim = c(250, 700), ylim = c(-40, 120), layout = c(n_cols, 1), 
      xlab = "Mean RT (ms)", ylab = "Delta (ms)",
      strip = strip.custom(strip.names = TRUE, sep = " = ", var.name = label)
    )
  } else {
    label1 <- mk_label(v1)
    label2 <- mk_label(v2)
    lattice::xyplot(
      delta ~ mean_rt | factor(v1) + factor(v2), data = delta_dat, 
      groups = if (conditional) {auto2}, main = if(title) {"Delta Plots"}, 
      sub = if(subt) {list(subtitle, cex = .7)}, type = c("b", "g"), 
      xlim = c(250, 700), ylim = c(-40, 130), 
      xlab = "Mean RT (ms)", ylab = "Delta (ms)",
      strip = strip.custom(strip.names = TRUE, sep = " = ", var.name = c(label1, label2))
    )
  }
}

plt_cdfs <- function(df, v1, v2 = NA) {
  #' Lattice plot for CDFs per condition
  #' @param df data frame. Must include columns auto1, auto2 (-1 incongruent, 
  #' 1 congruent), v1, v2, rt, dec (1 correct, -1 incorrect). 
  #' @param v1 character. First parameter to aggregate over
  #' @param v2 character. Second parameter to aggregate over
  
  if (!is.na(v2)) {
    plot_form <- as.formula(
      paste0("~ rt | factor(", v1, ") + factor(", v2, ")"))
  } else {
    plot_form <- as.formula(paste0("~ rt | factor(", v1, ")"))
  }
  label1 <- mk_label(v1)
  label2 <- mk_label(v2)
  subtitle <- mk_subtitle(df)

  latticeExtra::ecdfplot(
    plot_form, groups = congruency, data = df, type = c("l", "g"), 
    xlim = c(0, 1000), xlab = "RT", ylab = "ECDF",  main = "CDFs", 
    sub = list(subtitle, cex = .7), 
    strip = strip.custom(strip.names = TRUE, sep = " = ", 
                         var.name = c(label1, label2)), 
    auto.key = list(title = "condition", cex = 0.7, columns = 2)
  )
}

plt_denss <- function(df, v1, v2 = NA) {
  #' Lattice plot for densities per condition
  #' @param df data frame. Must include columns auto1, auto2 (-1 incongruent, 
  #' 1 congruent), v1, v2, rt, dec (1 correct, -1 incorrect). 
  #' @param v1 character. First parameter to aggregate over
  #' @param v2 character. Second parameter to aggregate over

  if (!is.na(v2)) {
    plot_form <- as.formula(
      paste0("~ rt | factor(", v1, ") + factor(", v2, ")"))
  } else {
    plot_form <- as.formula(paste0("~ rt | factor(", v1, ")"))
  }
  
  label1 <- mk_label(v1)
  label2 <- mk_label(v2)
  subtitle <- mk_subtitle(df)
  
  densityplot(plot_form, groups = factor(congruency), data = df, 
              plot.points = FALSE, xlab = "RT", main = "Densities", 
              sub = list(subtitle, cex = .7),
              auto.key = list(title = "condition", cex = 0.7, columns = 2, 
                              space = "top"), 
              strip = strip.custom(strip.names = TRUE, sep = " = ", 
                                   var.name = c(label1, label2)), 
  )
}