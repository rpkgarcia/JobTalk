#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
bartlett = function(x){
    y = 1- abs(x) 
    return(y)
}


qs = function(x){
    p1 = sin(6*pi*x/5)/(6*pi*x/5) - cos(6*pi*x/5)
    p2 = 25/(12*pi^2*x^2)
    y = p1*p2
}

parzen = function(x){
    if(abs(x) <0.5){
        y = 1 - 6*x^2 + 6*abs(x)^3
    } else{
        y = 2*(1-abs(x))^3
    }
    return(y)
}

kernel_weight = function(method, x_range){
    if(method == "Bartlett"){
        y = sapply(x_range, bartlett)
    } else if(method == "Parzen"){
        y = sapply(x_range, parzen)
    } else if(method == "Quadratic Spectral"){
        y = sapply(x_range, qs)
    } else if(method == "Truncated"){
        y = rep(1, length(x_range))
    }
    return(y)
}
library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput(inputId = "kernels", 
                               label = "Select Kernel",
                               choices = c("Bartlett", 
                                           "Parzen", 
                                           "Quadratic Spectral", 
                                           "Truncated"), 
                               selected = "Bartlett"), 
            selectInput("x_range", 
                        "X-Axis Range", 
                        c("Full View", "Half View"))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        if(input$x_range == "Full View"){
            x = seq(-1, 1, length.out = 200)
        } else{
            x = seq(0, 1, length.out = 200)
        }
        
        # Plot extra curves as needed
        the_colors = c("black", "red", "green", "blue")
        possible_kernels = c("Bartlett","Parzen", "Quadratic Spectral", "Truncated") 
        keep = which(possible_kernels %in% input$kernels)
        the_colors = the_colors[keep]
        the_lty = c(1:4) 
        the_lty = the_lty[keep]
        
        plot(x, y = rep(0, length(x)), ylim = c(0,1), col = "white", 
             ylab = "Weight", xlab = "Lags")
        
        if(length(input$kernels) != 0 ){
            for(i in 1:length(input$kernels)){
                y = kernel_weight(input$kernels[i], x)
                lines(x, y, 
                      col = the_colors[i], 
                      lwd = 2, 
                      lty = the_lty[i])
            }
            legend("topright", 
                   legend = c(input$kernels), 
                   lty =  the_lty, 
                   lwd = 2, 
                   col = the_colors)   
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
