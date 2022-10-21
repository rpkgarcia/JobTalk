#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Sidebar with a slider input for number of bins 
    
     fluidRow(column(2), column(4,
            numericInput("n",
                         "Sample Size (n):",
                         min = 1,
                         max = 10000,
                         value = 5, 
                         step = 1))), 

        # Show a plot of the generated distribution
        fluidRow(
           plotOutput("distPlot")
        )
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R

        curve(dnorm(x), from =-4, to =4, ylab = "", lwd = 2,  axes=FALSE, xlab = "")
        curve(dt(x, df =input$n), from =-4, to =4, add = T, col = "blue", lwd = 2, lty =2)
        legend("topright", border = "white", legend = c("Z-distr", "t-distr"), 
               lty = c(1, 2), col = c("black", "blue"), lwd = 2, bty = "n")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
