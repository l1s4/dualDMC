library(shiny)
library(Rcpp)
Rcpp::sourceCpp("../dualDMC.cpp")

ui <- fluidPage(
  sidebarPanel(width = 2, style = "margin: 10px; overflow-y:scroll; max-height: 10%; font-size:14px",
    selectInput("automProcess1", "Type of first automatic process", c("congruent", "incongruent")),
    selectInput("automProcess2", "Type of second automatic process", c("congruent", "incongruent")),
    sliderInput("mu_c", "mu_c [drift rate controlled]", 0, 1, 0.5),
    sliderInput("b", "b [decision boundary]",           0, 100, 50),
    sliderInput("sigma", "sigma [SD Wiener process]",   0, 10, 4),
    sliderInput("tau1", "tau1 [scale parameter auto1]", 0, 200, 20), 
    sliderInput("tau2", "tau2 [scale parameter auto2]", 0, 200, 30),
    sliderInput("A1", "A1 [amplitude auto1]",           0, 100, 20),
    sliderInput("A2", "A2 [amplitude auto2]",           0, 100, 20),
    sliderInput("N", "N [number of timepoints]",        5, 2000, 500),
    sliderInput("dt", "dt [step size]",                 0.1, 1, 1),
    sliderInput("ndt", "ndt [non-decision time [ms]]",  0, 700, 300),
    sliderInput("nSim", "nSim [number of simulations]", 500, 20000, 1000),
    checkboxInput("YlimFixed", "keep plot y-axis constant at [-50, +50]", value = FALSE),
    actionButton("SimTrial", "Simulate trial") 
  ),
  mainPanel(width = 5, wellPanel(plotOutput("AVplot"))),
  mainPanel(width = 5, wellPanel(plotOutput("XPlot"))), 
  mainPanel(width = 5, wellPanel(plotOutput("MeanRTPlot"))),
  mainPanel(width = 5, wellPanel(plotOutput("MeanERPlot")))
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(trials = list(), ap1 = list(), ap2 = list())     # to store simulated trials
  
  # change of parameter values
  observeEvent(input$change, {
    updateSliderInput(session, "N")
    updateSliderInput(session, "mu_c", max = N)
    updateSliderInput(session, "tau1", max = N)
    updateSliderInput(session, "tau2", max = N)
    updateSliderInput(session, "A1", max = N)
    updateSliderInput(session, "A2", max = N)
    updateSliderInput(session, "dt", max = N)
    updateSliderInput(session, "b", max = N)
    updateSliderInput(session, "nSim")
    updateSliderInput(session, "ndt")
  })
  
  # simulate Trial button 
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
      auto1 = auto1,
      auto2 = auto2
    )
    
    trial_id <- paste0("trial_", length(rv$trials) + 1)
    rv$trials[[trial_id]] <- M$XTraj
    rv$ap1[[trial_id]] <- M$auto1
    rv$ap2[[trial_id]] <- M$auto2
  })

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

    tr_len <- input$N - 1 # trim vector to avoid ugly plot
    plot(M$cont_traj[1:tr_len],
      type = "n", ylim = c(ymin, ymax),
      ylab = "Mean Activation", xlab = "t"
    )
    abline(h = 0, lty = 3)
    abline(h = c(-input$b, input$b), lty = 2)
    abline(v = ifelse(cross_pb < cross_mb, cross_pb, cross_mb))
    lines(M$cont_traj[1:tr_len], col = "black")
    lines(M$auto1_traj[1:tr_len], col = "green")
    lines(M$auto2_traj[1:tr_len], col = "blue")
    lines(M$super_traj[1:tr_len], col = "red")
    legend("bottomright",
      col = c("black", "green", "blue", "red"), lty = 1,
      legend = c("controlled", "automatic 1", "automatic 2", "superimposed")
    )
  })

  output$XPlot <- renderPlot({
    ymin <- if (input$YlimFixed) -50 else -input$b - 20
    ymax <- if (input$YlimFixed) +50 else input$b + 20

    plot(isolate(rv[['trial_1']]), type = "n", ylab = "X(t)", xlab = "t", 
         ylim = c(ymin, ymax), xlim = c(0, input$N)) 
    abline(h = input$b, col = "black", lty = 2)
    abline(h = -input$b, col = "black", lty = 2)


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
  
  output$MeanRTPlot <- renderPlot({
    Sim <- simDDMC(
      df = data.frame(
        mu_c  = input$mu_c,
        b     = input$b, 
        A1    = input$A1,
        A2    = input$A2, 
        tau1  = input$tau1,
        tau2  = input$tau2,
        dt    = input$dt, 
        sigma = input$sigma),
      input$nSim)
    
    s_dfrt <- aggregate(rt ~ auto1 + auto2, FUN = mean, data = Sim[Sim$dec == 1, ])
    s_dfrt$rt <- s_dfrt$rt + input$ndt    # add non decision time
    
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

  output$MeanERPlot <- renderPlot({
    Sim <- simDDMC(
      df = data.frame(
        mu_c  = input$mu_c,
        b     = input$b, 
        A1    = input$A1,
        A2    = input$A2, 
        tau1  = input$tau1,
        tau2  = input$tau2,
        dt    = input$dt, 
        sigma = input$sigma),
      input$nSim)

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
}

shinyApp(ui, server)
