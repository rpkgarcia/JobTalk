# This app is to illustrate the difference between the different lugsails

library(stringr)

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
    
  } else if(method == "Adaptive"){
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
        sliderInput("big_T_value", "Select Sample Size:", 
                     value = 200,
                     min=100, 
                     max = 100000, 
                     step = 100),
        sliderInput("b_value", "Select Bandwidth:", 
                     value = .071,
                     min = .02,
                     max = 1, 
                     step = .02),
        checkboxGroupInput("kernelOptions", 
                           "Select Kernels",
                     c("Zero", "Adaptive", "Over"))
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
  
    
    output$distPlot <- renderPlot({

      x = seq(-1, 1, length.out = 200)
      y = sapply(x, bartlett)
      
      # Get Maximum values for y-axis
      y_max = c(1)
      if(any(str_detect(input$kernelOptions, "Adaptive"))){
        r = lugsail_parameters(big_T = input$big_T_value, 
                               b = input$b_value, 
                               method = "Adaptive")$r
        c = lugsail_parameters(big_T = input$big_T_value, 
                               b = input$b_value,
                               method ="Adaptive")$c
        y_lugsail = sapply(x,lugsail_bartlett,
                           r = r,
                           c = c)
        y_max = c(y_max, max(y_lugsail))
      } 
      
      if(any(str_detect(input$kernelOptions, "Over"))){
        r = lugsail_parameters(method = "Over")$r
        c = lugsail_parameters(method ="Over")$c
        y_lugsail = sapply(x,lugsail_bartlett,
                           r = r,
                           c = c)
        y_max = c(y_max, max(y_lugsail))
      } 
      
      plot(x, y, xlab = "x",
           ylab = "weight", 
           type = "l", ylim = c(0, max(y_max)),
           lwd = 2)
      
      # Plot extra curves as needed
      the_colors = rainbow(3)
      
      keep = which(c("Zero","Adaptive", "Over") %in%input$kernelOptions)
      the_colors = the_colors[keep]
      the_lwd = 2:4
      the_lwd = the_lwd[keep]
    
      if(length(input$kernelOptions) != 0 ){
        for(i in 1:length(input$kernelOptions)){

          option = input$kernelOptions[i]
          option = gsub(" Lugsail", "", option)
          r = lugsail_parameters(big_T = input$big_T_value, 
                                 b = input$b_value, 
                                 method = option)$r
          c = lugsail_parameters(big_T = input$big_T_value, 
                                 b = input$b_value,
                                 method = option)$c
          y_lugsail = sapply(x,lugsail_bartlett,
                             r = r,
                             c = c)
        
          lines(x, y_lugsail, col = the_colors[i], lwd = 2, lty = the_lwd[i])
        }

      }


      legend("topright", 
             c("Bartlett", input$kernelOptions), 
             lty = c(1, the_lwd), 
             lwd = 2, 
             col = c("black", the_colors))
    })
    
  }
  
  shinyApp(ui, server)
