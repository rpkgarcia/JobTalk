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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("n",
                        "Sample Size (n):",
                        min = 1,
                        max = 10000,
                        value = 30, 
                        step = 1), 
            numericInput("seed",
                         "Set Seed:",
                         min = 1,
                         max = 10000,
                         value = 62, 
                         step = 1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlotly({
        # generate bins based on input$bins from ui.R
        
        sims <- 5000
        
        set.seed(input$seed)
        get_z <- function(estimate=TRUE){
            x    <- rexp(input$n, rate = 1)
            sigma <- sqrt(1/(input$n))
            
            if(estimate == TRUE){
                sigma <- sqrt(var(x)/(input$n-1))
            }
            z <- (mean(x)- 1)/sigma 
            return(z)
        }
        
        sim_means_Est <- replicate(sims, get_z(TRUE))
        sim_means <-  replicate(sims, get_z(FALSE))
        # draw the histogram with the specified number of bins
        the_values = seq(-3, 3, length.out = sims)
        trace_z <- dnorm(the_values)
        trace_t <- dt(the_values, df=(input$n-1))
        
        the_data <- data.frame(the_values, trace_z, trace_t)
        
        fig <- plot_ly(
                xlab = "x", ylab = "p(x)", 

                alpha = 0.6) 
        fig <- fig %>% add_histogram(x = ~sim_means,histnorm = "probability", 
                                     name = "True Sigma")
        fig <- fig %>%  add_histogram(x = ~sim_means_Est,histnorm = "probability", 
                                      name = "Estimated Sigma")
        fig <- fig %>% layout(barmode = "overlay")
        
        # fig <- fig %>% add_trace(the_data, y =~trace_t,, x = ~the_values, 
        #                          name = "t-dist", mode = "lines", 
        #                          type = "scatter")
        # fig <- fig %>% add_trace(the_data, y =~trace_z, x = ~the_values,
        #                          name = "Z-dist", mode = "lines", 
        #                          type = "scatter")
        fig
        
        # hist(z, col = 'darkgray', border = 'white')
        # curve(dnorm(x), from =-3, to =3, ylab = "p(x)")
        # curve(dt(x, df =(input$n -1)), from =-3, to =3, add = T, col = "blue")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
