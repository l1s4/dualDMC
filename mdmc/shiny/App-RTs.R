library(shiny)
source("../MDMC-Functions.R")

ui <- fluidPage(
  sidebarPanel(
    sliderInput("mu_c", "mu_c [Constant drift rate of controlled process]", 0, 1, 0.5),
    sliderInput("sigma", "sigma [SD of Wiener process]", 0, 10, 4),
    sliderInput("tau1", "tau1 [Scale parameter first gamma activation]", 0, 200, 20), 
    sliderInput("tau2", "tau2 [Scale parameter second gamma activation]", 0, 200, 30),
    sliderInput("A1", "A1 [Amplitude of first automatic activation]", 0, 100, 20),
    sliderInput("A2", "A2 [Amplitude of second automatic activation]", 0, 100, 20),
    sliderInput("b", "b [Decision boundary]", 0, 150, 50),
    sliderInput("a1", "a1 [Shape parameter first gamma activation]", 2, 10, 2),
    sliderInput("a2", "a2 [Shape parameter second gamma activation]", 2, 10, 2),
    sliderInput("N", "N [Number of Timepoints]", 5, 2000, 500),
    sliderInput("dt", "dt [Step size]", 0.001, 1, 0.2),
    sliderInput("ndt", "ndt [non-decision time [ms]]", 0, 700, 300),
    sliderInput("nSim", "nSim", 500, 20000, 1000),
    checkboxInput("YlimFixed", "Fixed axis for mean RT plot [0, 1000]", value = FALSE)
  ),
  mainPanel(wellPanel(plotOutput("AVplot"))),
  mainPanel(wellPanel(plotOutput("MeanRTPlot"))), 
  mainPanel(wellPanel(plotOutput("MeanERPlot")))
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(trials = list(), ap1 = list(), ap2 = list(), dfrt = list())     # to store simulated trials
  
  observeEvent(input$change, {
    N <- input$N
    updateSliderInput(session, "mu_c", max = N)
    updateSliderInput(session, "tau1", max = N)
    updateSliderInput(session, "tau2", max = N)
    updateSliderInput(session, "a1", max = N)
    updateSliderInput(session, "a2", max = N)
    updateSliderInput(session, "A1", max = N)
    updateSliderInput(session, "A2", max = N)
    updateSliderInput(session, "dt", max = N)
    updateSliderInput(session, "b", max = N)
    updateSliderInput(session, "nSim")
  })


  output$AVplot <- renderPlot({
    var_cond <- expand.grid(type_auto_1 = c("congruent", "incongruent"), 
              type_auto_2 = c("congruent", "incongruent")
    )
    AVSims <- list()
    for(cond in seq_len(nrow(var_cond))) {
      AVSim <- MDMC_AV(N = input$N, 
      parameters = list(
        mu_c = input$mu_c,
        tau1 = input$tau1,
        tau2 = input$tau2,
        a1   = input$a1,
        a2   = input$a2,
        A1   = input$A1,
        A2   = input$A2
      ),
      automatic1 = var_cond[cond, ]$type_auto_1, 
      automatic2 = var_cond[cond, ]$type_auto_2)
      AVSims[[cond]] <- AVSim
    }

    ymin <- -input$b - 20 
    ymax <- input$b + 20

    par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
    for(i in 1:4) {
      AVs <- AVSims[[i]]

      plot(AVs$trajectory_a1, type = "n", ylim = c(ymin, ymax), 
		     main = paste(unlist(var_cond[i, ]), collapse = "-"), 
		     xlab = "t[ms]"
      )
      lines(AVs$trajectory_a1, col = "red")
      lines(AVs$trajectory_a2, col = "blue")
      lines(AVs$controlled, col = "black")
      lines(AVs$superimposed, col = "green")

      abline(h = c(-input$b, input$b), lty = 2)

      legend(x = 400, y = -10, lty = 1, 
		       c("auto1", "auto2", "mu_c", "superimposed"), 
		       col = c("red", "blue", "black", "green"))
	
      # RT: check which boundary is crossed first by superimposed 
      p <- suppressWarnings(min(which(AVs$superimposed > input$b)))
      n <- suppressWarnings(min(which(AVs$superimposed < -input$b)))
      abline(v = ifelse(p < n, p, n))
    }
  })

  output$MeanRTPlot <- renderPlot({
    Sim <- MDMC_Sim(
      N_sim = input$nSim, 
      N_time = input$N, 
      param_grid = data.frame(
        mu_c  = input$mu_c,
        sigma = input$sigma,
        tau1  = input$tau1,
        tau2  = input$tau2,
        a1    = input$a1,
        a2    = input$a2,
        A1    = input$A1,
        A2    = input$A2, 
        dt    = input$dt, 
        b     = input$b 
      )
    )[[1]]

    s_dfrt <- aggregate(rt ~ first_pr + second_pr, FUN = mean, data = Sim[Sim$decision == "pb", ])
    s_dfrt$rt <- s_dfrt$rt + input$ndt    # add non decision time

    ymin <- if (input$YlimFixed) 0 else min(s_dfrt$rt) - 200
    ymax <- if (input$YlimFixed) 1000 else max(s_dfrt$rt) + 200
      
    interaction.plot(
      x.factor = s_dfrt$first_pr, trace.factor = s_dfrt$second_pr, 
      response = s_dfrt$rt, ylim = c(ymin, ymax), 
      xlab = "first", trace.label = "second", col = c("red", "blue"), 
      ylab = "mean RT [correct trials]", legend = T, 
      main = "mean RT [correct trials] per condition"
    )
  })

  output$MeanERPlot <- renderPlot({
    Sim <- MDMC_Sim(
      N_sim = 1000, 
      N_time = input$N, 
      param_grid = data.frame(
        mu_c  = input$mu_c,
        sigma = input$sigma,
        tau1  = input$tau1,
        tau2  = input$tau2,
        a1    = input$a1,
        a2    = input$a2,
        A1    = input$A1,
        A2    = input$A2, 
        dt    = input$dt, 
        b     = input$b 
      )
    )[[1]]

    s_dfer <- aggregate(error ~ first_pr + second_pr, FUN = mean, data = Sim)

    interaction.plot(
      x.factor = s_dfer$first_pr, trace.factor = s_dfer$second_pr, 
      response = s_dfer$error, ylim = c(0, 1), 
      xlab = "first", ylab = "mean ER", col = c("red", "blue"), 
      trace.label = "second", legend = T, 
      main = "mean ER per condition"
    )
  })

}

shinyApp(ui, server)