library(shiny)
source("../MDMC-Functions.R")

ui <- fluidPage(
  sidebarPanel(
    sliderInput("mu_c", "mu_c [Constant drift rate of controlled process]", 0, 1, 0.5),
    sliderInput("sigma", "sigma [SD of Wiener process]", 0, 10, 4),
    sliderInput("tau1", "tau1 [Scale parameter first Gamma activation function]", 0, 200, 20), 
    sliderInput("tau2", "tau2 [Scale parameter second Gamma activation function]", 0, 200, 30),
    sliderInput("A1", "A1 [Amplitude of first automatic activation]", 0, 100, 20),
    sliderInput("A2", "A2 [Amplitude of second automatic activation]", 0, 100, 20),
    sliderInput("b", "b [Decision boundary]", 0, 100, 50),
    sliderInput("a1", "a1 [Shape parameter of first Gamma function]", 2, 10, 2),
    sliderInput("a2", "a2 [Shape parameter of second Gamma function]", 2, 10, 2),
    sliderInput("N", "N [Number of Timepoints]", 5, 2000, 500),
    sliderInput("dt", "dt [Step size]", 0.1, 1, 1),
    selectInput(
      "automProcess1", "Type of first automatic process",
      c("congruent", "incongruent")
    ),
    selectInput(
      "automProcess2", "Type of second automatic process",
      c("congruent", "incongruent")
    ),
    checkboxInput("YlimFixed", "keep plot y-axis constant at [-50, +50]", value = FALSE),
    actionButton("SimTrial", "Simulate trial")
  ),
  mainPanel(wellPanel(plotOutput("AVplot"))),
  mainPanel(wellPanel(plotOutput("XPlot"))), 
  mainPanel(wellPanel(plotOutput("MeanRTPlot")))
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(trials = list(), ap1 = list(), ap2 = list())     # to store simulated trials
  
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
  })


  output$AVplot <- renderPlot({
    M <- MDMC_AV(
      N = input$N,
      parameters = list(
        mu_c = input$mu_c,
        tau1 = input$tau1,
        tau2 = input$tau2,
        a1   = input$a1,
        a2   = input$a2,
        A1   = input$A1,
        A2   = input$A2
      ),
      automatic1 = input$automProcess1,
      automatic2 = input$automProcess2
    )

    ymin <- if (input$YlimFixed) -50 else -input$b - 10 #min(M$trajectory_a1, M$trajectory_a2) - 30
    ymax <- if (input$YlimFixed) +50 else input$b + 10 #max(M$trajectory_a1, M$trajectory_a2) + 30

    cross_pb <- suppressWarnings(min(which(M$superimposed > input$b)))
    cross_mb <- suppressWarnings(min(which(M$superimposed < -input$b)))

    tr_len <- input$N - 1 # trim vector to avoid ugly plot
    plot(M$controlled[1:tr_len],
      type = "n", ylim = c(ymin, ymax),
      ylab = "Mean Activation", xlab = "t"
    )
    abline(h = 0, lty = 3)
    abline(h = c(-input$b, input$b), lty = 2)
    abline(v = ifelse(cross_pb < cross_mb, cross_pb, cross_mb))
    lines(M$controlled[1:tr_len], col = "black")
    lines(M$trajectory_a1[1:tr_len], col = "green")
    lines(M$trajectory_a2[1:tr_len], col = "blue")
    lines(M$superimposed[1:tr_len], col = "red")
    legend("bottomright",
      col = c("black", "green", "blue", "red"), lty = 1,
      legend = c("controlled", "automatic 1", "automatic 2", "superimposed")
    )
  })

  observeEvent(input$SimTrial, {

    M <- MDMC_T(
      N = input$N,
      parameters = list(
        mu_c  = input$mu_c,
        sigma = input$sigma,
        tau1  = input$tau1,
        tau2  = input$tau2,
        a1    = input$a1,
        a2    = input$a2,
        A1    = input$A1,
        A2    = input$A2, 
        dt    = input$dt, 
        b     = input$b, 
        automatic1 = input$automProcess1,
        automatic2 = input$automProcess2
      )
    )

    trial_id <- paste0("trial_", length(rv$trials) + 1)
    rv$trials[[trial_id]] <- M$TrajX
    rv$ap1[[trial_id]] <- M$autom1
    rv$ap2[[trial_id]] <- M$autom2
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

      # Check which boundary was crossed first and set line color
      cross_pb <- which(trial > input$b)
      cross_mb <- which(trial < -input$b)
      pb <- suppressWarnings(min(cross_pb))
      pm <- suppressWarnings(min(cross_mb))
      if (pb < pm) line_col <- "green" else line_col <- "red"

      lines(trial, col = line_col)
    }
  })
}

shinyApp(ui, server)
