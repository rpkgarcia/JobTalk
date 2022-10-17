#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(MASS)  # for mvnorm



generate_data = function(rho_y=.7, e_sd=.9, big_T=200){
  
  # Correlation matrix 
  R = matrix(0, nrow = 3, ncol = 3)
  diag(R) = rho_y
  Sigma = matrix(c(e_sd, e_sd^2, e_sd^3, 
                   e_sd^2, e_sd, e_sd^2, 
                   e_sd^3, e_sd^2, e_sd), nrow = 3, ncol =3)
  
  # Generate the data.  Initial value is 0
  sim_data = matrix(0, nrow = 3, ncol = big_T)
  
  # The rest of the values
  for(t in 2:big_T){
    sim_data[,t] = R%*%sim_data[,c(t-1)] + mvrnorm(1, rep(0, 3), 
                                                   Sigma = Sigma)
  }
  
  return(sim_data)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("rho",
                        "Choose correlation strength (rho):",
                        min = 0,
                        max = 1,
                        value = .7, 
                        step = .01), 
            numericInput("big_T", 
                         "Choose sample size:", 
                         min = 50, max = 10000, 
                         step = 50, value = 200),
            actionButton("action1", label = "(re)Generate")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel( 
            tabPanel("Plot", plotlyOutput("the_plot")),
            tabPanel(
                     p("Model", style = "color:black"),   
                     p("Example with a 3x1 vector autoregressive model with order 1 (VAR(1)). ", style = "color:black"),
                     withMathJax(), 
                     helpText('The Model: $$ \\begin{bmatrix}
                      y_{1,t}\\\\ 
                      y_{2,t}\\\\
                      y_{3,t}\\\\
                      \\end{bmatrix}=\\begin{bmatrix}
                      \\rho & 0 & 0 \\\\ 
                      0 & \\rho& 0 \\\\
                      0  & 0& \\rho \\\\
                      \\end{bmatrix}\\begin{bmatrix}
                      y_{1,t-1}\\\\ 
                      y_{2,t-1}\\\\
                       y_{3,t-1}\\\\
                      \\end{bmatrix} + \\begin{bmatrix}
                      e_{1,t}\\\\ 
                      e_{2,t}\\\\
                      e_{3,t}\\\\
                      \\end{bmatrix}$$', style = "color:black"),
                     helpText("where", style = "color:black"),
                     helpText("$$e_{t} \\overset{iid}{\\sim} N_3(0, \\Sigma_e)$$", style = "color:black"),
                     helpText("$$\\Sigma_e = \\begin{bmatrix}
                      .9 & .9^2 & .9^3 \\\\ 
                      .9^2  & .9& .9^2 \\\\
                      .9^3  & .9^2& .9^3 \\\\
                      \\end{bmatrix}$$", style = "color:black"))
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white')
    # })
  
  observeEvent(input$action1, {
    
    output$the_plot <- renderPlotly({
      
      # Multivariate 3x1 VAR(1) process
      # generate 3xbig_T matrix of data
      y = generate_data(rho_y = input$rho, 
                        big_T = input$big_T)
      
      plot1 <- plot_ly(
          x = 1:input$big_T,
          y = y[1, ], 
          type = 'scatter',
          mode = 'lines')%>%
          layout(showlegend = F)
      
      plot2 <- plot_ly(
        x = 1:input$big_T,
        y = y[2,], 
        type = 'scatter',
        mode = 'lines')%>%
        layout(showlegend = F)
      
      plot3 <- plot_ly(
        x = 1:input$big_T,
        y = y[3,], 
        type = 'scatter',
        mode = 'lines')%>%
        layout(showlegend = F)
      
      fig <- subplot(plot1, plot2, plot3, 
                     nrows =3, titleY = TRUE) %>% layout(
                     plot_bgcolor='white', 
                     hovermode="x unified")
      })
  }) 
}

# Run the application 
shinyApp(ui = ui, server = server)
