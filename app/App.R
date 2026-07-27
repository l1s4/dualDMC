library(shiny)
library(Rcpp)
Rcpp::sourceCpp("../src/dualDMC.cpp")

ui <- fluidPage(
  sidebarPanel(width = 2, style = "margin: 10px; overflow-y:scroll; max-height: 10%; font-size:14px",
    selectInput("automProcess1", "Type of first automatic process",  c("congruent", "incongruent")),
    selectInput("automProcess2", "Type of second automatic process", c("congruent", "incongruent")),
    sliderInput("mu_c",          "mu_c [drift rate controlled]",     0, 1, 0.5),
    sliderInput("b",             "b [decision boundary]",            0, 100, 60),
    sliderInput("tau1",          "tau1 [scale parameter auto1]",     0, 250, 80), 
    sliderInput("tau2",          "tau2 [scale parameter auto2]",     0, 250, 80),
    sliderInput("A1",            "A1 [amplitude auto1]",             0, 50, 15),
    sliderInput("A2",            "A2 [amplitude auto2]",             0, 50, 15),
    sliderInput("ndt_m",         "mean ndt [non-decision time]",     0, 700, 330),
    sliderInput("ndt_sd",        "sd of ndt [non-decision time]",    0, 100, 30),
    sliderInput("N",             "N [number of timepoints]",         5, 2000, 500),
    sliderInput("dt",            "dt [step size]",                   0.1, 1, 1),
    sliderInput("sigma",         "sigma [SD Wiener process]",        0, 10, 4),
    sliderInput("nSim",          "nSim [number of simulations]",     300, 10000, 5000),
    checkboxInput("YlimFixed",   "keep plot y-axis constant at [-50, +50]", value = FALSE),
    actionButton("SimTrial",     "Simulate trial") 
  ),
  mainPanel(width = 5, wellPanel(plotOutput("AVplot"))),
  mainPanel(width = 5, wellPanel(plotOutput("MeanRTPlot"))),
  mainPanel(width = 5, wellPanel(plotOutput("DeltaPlot"))),
  mainPanel(width = 5, wellPanel(plotOutput("DeltaPlotCond"))),
  mainPanel(width = 5, wellPanel(plotOutput("MeanERPlot"))),
  mainPanel(width = 5, wellPanel(plotOutput("DensityPlot"))),
  mainPanel(width = 5, wellPanel(plotOutput("XPlot")))
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(trials = list(), ap1 = list(), ap2 = list())     # to store simulated trials
  
  # change of parameter values
  observeEvent(input$change, {
    updateSliderInput(session, "N")
    updateSliderInput(session, "mu_c")
    updateSliderInput(session, "tau1")
    updateSliderInput(session, "tau2")
    updateSliderInput(session, "A1")
    updateSliderInput(session, "A2")
    updateSliderInput(session, "dt")
    updateSliderInput(session, "b")
    updateSliderInput(session, "nSim")
    updateSliderInput(session, "ndt_m")
    updateSliderInput(session, "ndt_sd")
  })
  
  # simulate trial button 
  observeEvent(input$SimTrial, {
    auto1 <- ifelse(input$automProcess1 == "congruent", 1, -1)
    auto2 <- ifelse(input$automProcess2 == "congruent", 1, -1)
    M <- simDDMCtrial(
      mu_c  = input$mu_c,
      b     = input$b, 
      A1    = input$A1,
      A2    = input$A2, 
      tau1  = input$tau1,
      tau2  = input$tau2,
      dt    = input$dt, 
      sigma = input$sigma,
      ndt_m = input$ndt_m, 
      ndt_sd= input$ndt_sd,
      auto1 = auto1,
      auto2 = auto2
      
    )
    
    trial_id <- paste0("trial_", length(rv$trials) + 1)
    rv$trials[[trial_id]] <- M$XTraj
    rv$ap1[[trial_id]] <- M$auto1
    rv$ap2[[trial_id]] <- M$auto2
  })
  
  # plot activation functions
  output$AVplot <- renderPlot({
    M <- simDDMCactivation(
      N     = input$N, 
      mu_c  = input$mu_c, 
      b     = input$b, 
      A1    = input$A1, 
      A2    = input$A2, 
      tau1  = input$tau1, 
      tau2  = input$tau2,
      auto1 <- ifelse(input$automProcess1 == "congruent", 1, -1),
      auto2 <- ifelse(input$automProcess2 == "congruent", 1, -1)
    )
    
    ymin <- if (input$YlimFixed) -50 else -input$b - 10 
    ymax <- if (input$YlimFixed) +50 else input$b + 10 

    cross_pb <- suppressWarnings(min(which(M$superimposed > input$b)))
    cross_mb <- suppressWarnings(min(which(M$superimposed < -input$b)))

    plot(M$cont_traj[1:input$N],
      type = "n", ylim = c(ymin, ymax),
      ylab = "Mean Activation", xlab = "t [ms]"
    )
    abline(h = 0, lty = 3)
    abline(h = c(-input$b, input$b), lty = 2)
    lines(M$cont_traj[1:input$N], col = "black")
    lines(M$auto1_traj[1:input$N], col = "green")
    lines(M$auto2_traj[1:input$N], col = "blue")
    lines(M$super_traj[1:input$N], col = "red")
    legend("bottomright",
      col = c("black", "green", "blue", "red"), lty = 1,
      legend = c("controlled", "automatic 1", "automatic 2", "superimposed")
    )
  })
  
  # plot simulated trials
  output$XPlot <- renderPlot({
    ymin <- if (input$YlimFixed) -50 else -input$b - 20
    ymax <- if (input$YlimFixed) +50 else input$b + 20

    plot(isolate(rv[['trial_1']]), type = "n", ylab = "X(t)", xlab = "t [ms]", 
         ylim = c(ymin, ymax), xlim = c(0, input$N)) 
    abline(h = c(-input$b, input$b), col = "black", lty = 2)


    for (trial_id in names(rv$trials)) {
      trial <- rv$trials[[trial_id]]

      # check which boundary was crossed first and set line color
      cross_pb <- which(trial > input$b)
      cross_mb <- which(trial < -input$b)
      pb <- suppressWarnings(min(cross_pb))
      pm <- suppressWarnings(min(cross_mb))
      if (pb < pm) line_col <- "green" else line_col <- "red"

      lines(trial, col = line_col)
    }
  })
  
  # simulate data for nSim trials, recompute everytime input changes
  # used in delta plots, meanRT, meanER, density plot
  SimData <- reactive({
    simDDMC(
      df = data.frame(
        mu_c   = input$mu_c,
        b      = input$b, 
        A1     = input$A1,
        A2     = input$A2, 
        tau1   = input$tau1,
        tau2   = input$tau2,
        dt     = input$dt, 
        sigma  = input$sigma,
        ndt_m  = input$ndt_m, 
        ndt_sd = input$ndt_sd), 
      input$nSim)
  })
  # plot mean RTs
  output$MeanRTPlot <- renderPlot({
    Sim <- SimData()
    
    s_dfrt <- aggregate(rt ~ auto1 + auto2, FUN = mean, data = Sim[Sim$dec == 1, ])
    s_dfrt$rt <- s_dfrt$rt
    
    # reorder factor levels and add labels
    s_dfrt$auto1 <- factor(s_dfrt$auto1, levels = c(1, -1), 
                          labels = c("congruent", "incongruent"))
    s_dfrt$auto2 <- factor(s_dfrt$auto2, levels = c(1, -1), 
                          labels = c("congruent", "incongruent"))

    ymin <- if (input$YlimFixed) 0 else min(s_dfrt$rt) - 200
    ymax <- if (input$YlimFixed) 1000 else max(s_dfrt$rt) + 200
      
    interaction.plot(
      x.factor = s_dfrt$auto1, trace.factor = s_dfrt$auto2, 
      response = s_dfrt$rt, ylim = c(ymin, ymax), col = c("red", "blue"), 
      xlab = "first dimension", trace.label = "second dimension", 
      ylab = "mean RT [correct trials]", legend = T, 
      main = "mean RT [correct trials] per condition"
    )
  })
  
  # plot mean ERs
  output$MeanERPlot <- renderPlot({
    Sim <- SimData()

    Sim$error <- ifelse(Sim$dec == -1, 1, 0)
    s_dfer <- aggregate(error ~ auto1 + auto2, FUN = mean, data = Sim)
    
    # reorder factor levels and add labels
    s_dfer$auto1 <- factor(s_dfer$auto1, levels = c(1, -1), 
                          labels = c("congruent", "incongruent"))
    s_dfer$auto2 <- factor(s_dfer$auto2, levels = c(1, -1), 
                          labels = c("congruent", "incongruent"))

    interaction.plot(
      x.factor = s_dfer$auto1, trace.factor = s_dfer$auto2, 
      response = s_dfer$error, ylim = c(0, 1), col = c("red", "blue"), 
      xlab = "first dimension", trace.label = "second dimension", 
      legend = T, ylab = "mean ER", main = "mean ER per condition"
    )
  })
  
  output$DeltaPlot <- renderPlot({
    Sim <- SimData()
    
    d_uc1 <- Sim[Sim$auto1 == 1, ]
    d_ui1 <- Sim[Sim$auto1 == -1, ]
    
    probs <- seq(0.1, 0.9, by = 0.1)
    
    q_uc1 <- quantile(d_uc1$rt, probs)
    q_ui1 <- quantile(d_ui1$rt, probs)
    
    delta_u1 <- q_ui1 - q_uc1
    mean_rts <- (q_uc1+q_ui1)/2
    
    ymin <- min(delta_u1) - 20
    ymax <- max(delta_u1) + 20
    xmin <- min(mean_rts) - 25
    xmax <- max(mean_rts) + 25 
    
    plot(mean_rts, delta_u1, type = "b", pch = 16, 
         xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = "blue", 
         xlab = "Mean RT (ms)", ylab = "delta", 
         main = "Delta Plot Task 1 [across Task 2 congruency]")
  })
  
  output$DeltaPlotCond <- renderPlot({
    Sim <- SimData()
    
    d_1c2c <- Sim[Sim$auto1 == 1 & Sim$auto2 == 1, ]
    d_1c2i <- Sim[Sim$auto1 == 1 & Sim$auto2 == -1, ]
    d_1i2c <- Sim[Sim$auto1 == -1 & Sim$auto2 == 1, ]
    d_1i2i <- Sim[Sim$auto1 == -1 & Sim$auto2 == -1, ]
    
    probs <- seq(0.1, 0.9, by = 0.1)
    
    q_1c2c <- quantile(d_1c2c$rt, probs)
    q_1c2i <- quantile(d_1c2i$rt, probs)
    q_1i2c <- quantile(d_1i2c$rt, probs)
    q_1i2i <- quantile(d_1i2i$rt, probs)
    
    
    delta_2c <- q_1i2c - q_1c2c
    delta_2i <- q_1i2i - q_1c2i
    
    mean_rts_2c <- (q_1i2c + q_1c2c) / 2
    mean_rts_2i <- (q_1i2i + q_1c2i) / 2
    
    ymin <- min(delta_2c, delta_2i) - 20
    ymax <- max(delta_2c, delta_2i) + 20
    xmin <- min(mean_rts_2c, mean_rts_2i) - 25
    xmax <- max(mean_rts_2c, mean_rts_2i) + 25 
    
    plot(mean_rts_2c, delta_2c, type = "n", pch = 16, 
         xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
         xlab = "Mean RT (ms)", ylab = "delta", main = "Delta Plot Task 1")
    points(mean_rts_2c, delta_2c, type = "b", pch = 16, col = "blue")
    points(mean_rts_2i, delta_2i, type = "b", pch = 16, col = "red")
    legend("topright", title = "Task two", legend = c("congruent", "incongruent"), 
           col = c("blue", "red"), lty = 1, pch = 16)
  })
  
  output$DensityPlot <- renderPlot({
    Sim <- SimData()
    
    d_cc <- Sim[Sim$auto1 == 1 & Sim$auto2 == 1, ]
    d_ci <- Sim[Sim$auto1 == 1 & Sim$auto2 == -1, ]
    d_ic <- Sim[Sim$auto1 == -1 & Sim$auto2 == 1, ]
    d_ii <- Sim[Sim$auto1 == -1 & Sim$auto2 == -1, ]
    
    plot(density(d_cc$rt), type = "n", main = "Density", xlab = "RT")
    lines(density(d_cc$rt), col = "green")
    lines(density(d_ci$rt), col = "blue")
    lines(density(d_ic$rt), col = "orange")
    lines(density(d_ii$rt), col = "red")
    legend("topright", title = "Condition", legend = c("cc", "ci", "ic", "ii"), 
           col = c("green", "blue", "orange", "red"), lty = 1)
  })
}

shinyApp(ui, server)
