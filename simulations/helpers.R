library(docstring)

check_ceff <- function(dat, A1, A2, tau1, tau2, mu_c, b) {
  #' Calculate D score between condition: (CI - CC) - (II - IC)
  #' @param dat data frame. Must include columns A1, A2, tau1, tau2, rt, congruency
  #' @param A1 int, A1 value
  #' @param A2 int, A2 value
  #' @param tau1 int, tau1 value
  #' @param tau2 int, tau2 value
  #' @return diff_diffs, int D score
  df <- dat[dat$b == b & dat$mu_c == mu_c & dat$dec == 1 & 
              dat$A1 == A1 & dat$A2 == A2 & 
              dat$tau1 == tau1 & dat$tau2 == tau2, ]
  
  meanrts <- aggregate(rt ~ congruency, FUN = mean, df)$rt
  ceff1 <- meanrts[2] - meanrts[1]
  ceff2 <- meanrts[4] - meanrts[3]
  diff_diffs <- ceff1 - ceff2
  
  diff_diffs
}

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
  #' @param corr_only boolean. If TRUE, only correct trials are used. 
  #' Only applies if type is 'RT'
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
    ylims <- if (any(is.na(ylim_custom))) c(-0.1, 1) else ylim_custom
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
    par.strip.text = if (small_panel_text) {par.strip.text = list(cex = .8)}, 
    strip = strip.custom(strip.names = TRUE, sep = " = ", 
                         var.name = c(label1, label2)), 
    auto.key = if (auto_key) {list(title = lines_name, space = legend_pos, 
                                   cex = .7)}
  )
}

plt_cafs <- function(df, v1, v2 = NA, n_bins = 5, subt = TRUE, title = TRUE, 
                     auto_key = TRUE, small_panel_text = FALSE) {
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
      dplyr::group_by(.data[[v1]], .data[[v2]], congruency) %>% 
      dplyr::mutate(rt_bin = ntile(rt, n_bins)) %>%
      dplyr::group_by(.data[[v1]], .data[[v2]], congruency, rt_bin) %>% 
      dplyr::summarise(mean_acc = mean(acc), .groups = "drop")
    plot_form <- as.formula(
      paste0("mean_acc ~ rt_bin | factor(", v1, ") + factor(", v2, ")")
    )
  } else {
    caf_data <- df %>% 
      dplyr::group_by(.data[[v1]], congruency) %>% 
      dplyr::mutate(rt_bin = ntile(rt, n_bins)) %>%
      dplyr::group_by(.data[[v1]], congruency, rt_bin) %>% 
      dplyr::summarise(mean_acc = mean(acc), .groups = "drop")
    plot_form <- as.formula(paste0("mean_acc ~ rt_bin | factor(", v1, ")"))
  }

  label1 <- mk_label(v1) 
  label2 <- mk_label(v2)
  subtitle <- mk_subtitle(df)
  lattice::xyplot(
    plot_form, data = caf_data, type = c("b", "g"), 
    groups = factor(
      congruency, labels = c("congruent-congruent", "congruent-incongruent", 
                             "incongruent-congruent", "incongruent-incongruent")
    ), ylim = c(0, 1.1), xlab = "RT bin", ylab = "Accuracy",
    sub = if(subt) {list(subtitle, cex = .7)}, main = if(title) {"CAFs"},
    scales = list(x = list(at = seq(1, n_bins, by=1))), 
    par.strip.text = if (small_panel_text) {par.strip.text = list(cex = .8)}, 
    strip = strip.custom(strip.names = TRUE, sep = " = ", 
                         var.name = c(label1, label2)),
    auto.key = if (auto_key) {
      list(title = "Condition", cex = .7, columns = 2, space = "top")
    } 
  )
}

plt_delta <- function(df, v1, v2, conditional = FALSE, subt = TRUE, 
                      title = TRUE, corr_only = TRUE, auto_key = FALSE){
  #' Lattice plot for delta functions of the first task
  #' Calculated as difference between (auto1 == -1) - (auto1 == 1)
  #' @param df data frame. Must include columns auto1, auto2 (-1 incongruent, 
  #' 1 congruent), v1, v2, rt, dec (1 correct, -1 incorrect). 
  #' @param v1 character. First parameter to aggregate over
  #' @param v2 character. Second parameter to aggregate over
  #' @param conditional boolean. If TRUE, two delta functions will be 
  #' calculated. One is for auto2 == -1, one for auto2 == 1
  #' @param subt boolean. If TRUE, subtitle with parameter values in df is added
  #' @param title boolean. If TRUE, 'CAFs' is added as a title
  #' @param corr_only boolean. If TRUE, only correct trials are used 
  
  if (corr_only) {df <- df[df$dec == 1, ]}    # correct trials only
  probs = seq(0.1, 0.9, 0.1)
  
  df$auto1 <- factor(df$auto1, levels =c(-1, 1) ,labels = c("inc", "con"))
  delta_dat <- df %>% 
    {if (conditional) group_by(., .data[[v1]], .data[[v2]], auto1, auto2) else 
      group_by(., .data[[v1]], .data[[v2]], auto1)} %>%
    # calculate rt quantiles per parameter combination and auto 1/2 congruencies
    summarise(q = list(quantile(rt, probs = probs)), .groups = "drop") %>%
    tidyr::unnest_longer(q, indices_to = "bin") %>%
    {if (conditional) tidyr::pivot_wider(
      ., id_cols = c(all_of(v1), all_of(v2), auto2, bin), names_from = auto1, 
      values_from = q
    ) else tidyr::pivot_wider(
        ., id_cols = c(all_of(v1), all_of(v2), bin), names_from = auto1, 
        values_from = q
    )} %>% mutate(mean_rt = (inc + con) / 2, delta = inc - con) %>% 
    as.data.frame()
  
  names(delta_dat)[(which(colnames(delta_dat) == v1))] <- "v1"
  names(delta_dat)[(which(colnames(delta_dat) == v2))] <- "v2"
  
  subtitle <- mk_subtitle(df)
  
  if (all(df[[v1]] == df[[v2]])) {
    label <- switch (v1, "tau1" = expression(tau[1] == tau[2]), 
                     "A1" = expression(A[1] == A[2]), v1)
    n_cols <- length(unique(delta_dat$v1))
    plt_form <- as.formula(paste0("delta ~ mean_rt | factor(v1)"))
  } else {
    label1 <- mk_label(v1)
    label2 <- mk_label(v2)
    plt_form <- as.formula(paste0("delta ~ mean_rt | factor(v1) + factor(v2)"))
    n_cols <- NA
  }
    
  lattice::xyplot(
    plt_form, data = delta_dat, xlim = c(250, 700), ylim = c(-40, 130), 
    groups = if (conditional) {auto2}, main = if(title) {"Delta Plots"}, 
    sub = if(subt) {list(subtitle, cex = .7)}, type = c("b", "g"), 
    xlab = "Mean RT (ms)", ylab = "Delta (ms)", layout = if(!is.na(n_cols)) {c(n_cols, 1)}, 
    strip = strip.custom(
      strip.names = TRUE, sep = " = ", 
      var.name = if (!all(df[[v1]] == df[[v2]])) c(label1, label2) else label
    ), 
    auto.key = if (auto_key) {list(title = "auto2", space = "top", cex = .7)}
  )
}