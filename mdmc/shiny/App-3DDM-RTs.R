library(shiny)
source("../MDMC-Functions.R")

ui <- fluidPage(
  sidebarPanel(
    sliderInput("mu_c", "mu_c [Constant drift rate of controlled process]", 0, 1, 0.5),
    sliderInput("sigma", "sigma [SD of Wiener process]", 0, 10, 4),
    sliderInput("delta1", "delta1 [weight first automatic process]", 0, 1, 0.5), 
    sliderInput("delta2", "delta2 [weight first automatic process]", 0, 1, 0.5), 
    sliderInput("b", "b [Decision boundary]", 0, 150, 50),
    sliderInput("N", "N [Number of Timepoints]", 5, 2000, 500),
    sliderInput("dt", "dt [Step size]", 0.001, 1, 0.2),
    sliderInput("ndt", "ndt [non-decision time [ms]]", 0, 700, 300),
    sliderInput("nSim", "nSim", 500, 20000, 1000),
    checkboxInput("YlimFixed", "Fixed axis for mean RT plot [0, 1000]", value = FALSE)
  ),
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

  output$MeanRTPlot <- renderPlot({
    Sim <- MDDM_Sim(
      N_sim = input$nSim, 
      N_time = input$N, 
      param_grid = data.frame(
        mu_c  = input$mu_c,
        sigma = input$sigma,
        delta_1  = input$delta1,
        delta_2  = input$delta2,
        dt    = input$dt, 
        b     = input$b 
      )
    )[[1]]

    s_dfrt <- aggregate(rt ~ first_pr + second_pr, FUN = mean, data = Sim[Sim$error == 0, ])
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
    Sim <- MDDM_Sim(
      N_sim = 1000, 
      N_time = input$N, 
      param_grid = data.frame(
        mu_c  = input$mu_c,
        sigma = input$sigma,
        delta_1  = input$delta1,
        delta_2  = input$delta2,
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