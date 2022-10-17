# This app is to illustrate the difference between the different lugsails



# Support Functions -------------------------------------------------------

bartlett = function(x){
  y = 1- abs(x) 
  return(y)
}

lugsail_bartlett = function(x, r, c){
  
  # Actual lugsail 
  y1 = bartlett(x)/(1-c) 
  y2 = 0 
  
  if(abs(x) < 1/r){
    y2 = bartlett(x*r)*c/(1-c) 
  }
  y = y1- y2
  

  return(y)
}



# Recommended b using Andrews rule
lugsail_parameters = function(big_T = 200, b = .071, q = 1, method = "Zero"){
  
  if(method == "Over"){
    r = 3
    c = 2/(1+r^q)
    
  } else if(method == "Adaptive*"){
    r = 2
    M  = big_T * b
    c_num = (log(big_T) - log(M) + 1)
    c_den = r^q*(log(big_T) - log(M)) + 1
    c = c_num/c_den 
    
  } else {
    # Zero or Manual lugsail
    r = 2
    c = r^(-q)
    
  }
  parameters = list(r = r, c = round(c, 2))
  return(parameters)
}


# Shiny App ---------------------------------------------------------------





  
  ui <- fluidPage(
    
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        radioButtons("lugsailOptions", "Select Lugsail Type",
                     c("Zero", "Adaptive*", "Over")),
        numericInput("r_value", "Select r:", value = 2,
                     min = 1, max = 3, step = .01), 
        numericInput("c_value", "Select c:", value = .5,
                     min = .01, max = 3, step = .01), 
        p("*Sample Size = 200, and b = 0.071")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("distPlot")
      )
    )

  )
  
  server <- function(input, output, session) {
    # observe({
    #   x <- input$lugsailOptions
    #   
    #   # Can use character(0) to remove all choices
    #   if (is.null(x))
    #     x <- character(0)
    #   
    #   # Can also set the label and select items
    #   updateSelectInput(session, "inSelect",
    #                     label = paste("Select input label", length(x)),
    #                     choices = x,
    #                     selected = tail(x, 1)
    #   )
    # })
    
    observe({
      x <- input$lugsailOptions
      

      r = lugsail_parameters(method = x)$r
      c = lugsail_parameters(method = x)$c
  
      #   x <- character(0)
      
      # Can also set the label and select items
      updateNumericInput(session, "r_value",
                         "Select r:", value = r,
                         min = .01, max = 3, step = .01
      )
      
      # Can also set the label and select items
      updateNumericInput(session, "c_value",
                         "Select c:", value = c,
                         min = .01, max = 3, step = .01
      )
    })
    
    
    output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x = seq(-1, 1, length.out = 200)
      y = sapply(x, bartlett)
      y_lugsail = sapply(x,lugsail_bartlett, 
                         r = input$r_value, 
                         c = input$c_value)
      
      plot(x, y, xlab = "x", ylab = "weight", type = "l", 
           ylim = range(c(y, y_lugsail)), lwd = 2)
      lines(x, y_lugsail, col = "red", lwd = 2)
      legend("topright", 
             c("Bartlett", "Lugsail"), 
             lty = c(1, 1), 
             lwd = c(2, 2), 
             col = c("black", "red"))
    })
    
  }
  
  shinyApp(ui, server)
