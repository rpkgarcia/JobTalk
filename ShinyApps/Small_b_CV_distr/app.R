#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# Support Functions -------------------------------------------------------


library(shiny)
library(stringr)

# Get Critical Values 
fitted_model <- function(cv_matrix, alpha_level=0.05, m=1){
    chisq_cv <- qchisq(1-alpha_level, df = m)/m 
    try_b <- cv_matrix[,1]
    try_b <- as.numeric(gsub("b=", "", try_b))
    specific_cvs = cv_matrix[,m+1] - chisq_cv
    fit = lm(specific_cvs ~ 0+  poly(try_b, 3, raw = T))     
    return(fit)
}

# Get specific fitted/predicted value 
fitted_value <- function(fit, b, alpha = 0.05, m = 1){
    chisq_cv <- qchisq(1-alpha, df = m)/m
    fit_b <- sum(fit$coefficients*poly(b, 3, raw = T))
    fit_b <- chisq_cv + fit_b
    return(fit_b)
}


# Shiny App ---------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Sidebar with a slider input for number of bins 
    fluidRow(column(2), 
             column(4, selectInput("b",
                        "Select b:",
                        list("0.2" = 0.2,
                             "0.07" = 0.07,
                             "value" = 30))
                    ), 
            column(4, selectInput("alpha",
                        "Select Significance Level:",
                        list("0.01" = "01",
                             "0.025" = "025",
                             "0.05" = "05", 
                             "0.10" = "10"))
                   )
    ), 

    fluidRow(
        column(6, plotOutput("distPlot")), 
        column(6, plotOutput("cvPlot"))
    )

    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$cvPlot <-renderPlot({
        
        # Alot of this came from "fitted_fixed_b.R"
        alpha <- as.numeric(paste(".", input$alpha, sep = ""))
        core_link <- "https://raw.githubusercontent.com/rpkgarcia/JobTalk/main/ShinyApps/Small_b_CV_distr/Bartlet_CV/"
        bartlett_file <- paste(core_link, 
                               "Bartlett_", 
                               input$alpha, ".csv", sep = "")
        bartlett_table <- read.csv(bartlett_file)
        bartlett_lug_file <- paste(core_link, 
                                   "Bartlett_Lugsail_", 
                                   input$alpha, ".csv", sep = "")
        bartlett_lug_table <- read.csv(bartlett_lug_file)
        
        
        m <- 1
        chisq_cv <- qchisq(1-alpha, df = m)/m 
        try_b <- bartlett_table[,1]
        try_b <- as.numeric(gsub("b=", "", try_b))
        fit_bartlett <- fitted_model(bartlett_table, alpha_level = alpha)
        fit_bartlett_lug <- fitted_model(bartlett_lug_table, alpha_level = alpha)
        
        bartlett_fit <- fitted_value(fit_bartlett, 
                                     as.numeric(input$b), alpha)
        bartlett_lug_fit <- fitted_value(fit_bartlett_lug, 
                                         as.numeric(input$b), alpha) 
        y_max <- max(10, bartlett_fit, bartlett_lug_fit)
        
        plot(c(0, try_b), 
             c(chisq_cv, fit_bartlett$fitted.values + chisq_cv), 
             type = "l", 
             xlim = c(0, 0.5), ylim = c(0, y_max),
             col = "red", lwd = 4, xlab = "Bandwidth (b)", 
             ylab = "Critical Value", 
             main = c("CV vs b"))
        mtext( cex=1, side = 3,
              "Given: kernel, m, significance level (alpha)")
        lines(c(0, try_b), 
              c(chisq_cv, fit_bartlett_lug$fitted.values + chisq_cv), 
              col = "blue", lwd = 4)
        
        points(x = c(input$b, input$b), 
               y = c(bartlett_fit,
                     bartlett_lug_fit))
        # text(x = c(input$b, input$b), 
        #      y = c(bartlett_fit,
        #            bartlett_lug_fit), 
        #      labels = round(c(bartlett_fit,
        #                       bartlett_lug_fit), 3))
        
        # My observed values 
        abline(h = chisq_cv, lwd = 4, lty = 2)
        abline(v = input$b, lwd = 4, lty = 3, col = "grey")
        
        legend("bottomright", 
               legend = c("Small-b CV", 
                          "Fixed-b: Bartlett", 
                          "Fixed-b: Lugsail", 
                          "CVs at b = 0.07"), 
               lty = c(2, 1, 1, 3), 
               col = c("black", "red", "blue", "grey"), 
               lwd = 4)
    })

    output$distPlot <- renderPlot({
        CV_link <- "https://raw.githubusercontent.com/rpkgarcia/JobTalk/main/ShinyApps/Small_b_CV_distr/Bartlet_CV/"
        distr_link <-"https://raw.githubusercontent.com/rpkgarcia/JobTalk/main/ShinyApps/Small_b_CV_distr/distr_est/"
        
        file_name <- paste(distr_link, 
                           "distr_est_b", 
                           input$b, ".csv", sep = "")
        alpha <- as.numeric(paste(".", input$alpha, sep = ""))
        dist_keep <- read.csv(file_name)
        
        bartlett_density <- density(dist_keep$Bartlett, from = 0)
        bartlett_lug_density <- density(dist_keep$Bartlett_lug, from = 0)
        
        # Bartlett Fits (need this for plot bounds)
        file_name <- paste(CV_link, 
                           "Bartlett_", 
                           input$alpha, ".csv", sep = "")
        bart_table <- read.csv(file_name)
        bart_fit <- fitted_model(bart_table, alpha_level = alpha)
        bart_cv <- fitted_value(bart_fit, as.numeric(input$b), alpha)
        
        # Bartlett Lugsail Fits (need this for plot bounds)
        file_name <- paste(CV_link, 
                           "Bartlett_Lugsail_", 
                           input$alpha, ".csv", sep = "")
        bart_lug_table <- read.csv(file_name)
        bart_lug_fit <- fitted_model(bart_lug_table, alpha_level = (alpha))
        bart_lug_cv <- fitted_value(bart_lug_fit, as.numeric(input$b), alpha)
        
        
        # Plot density curves
        x_max <-max(10, bart_cv, bart_lug_cv)
        plot(bartlett_density, col = "red", xlim = c(0.01, x_max*1.03), 
             main = "Distribution of Test Statistics", 
             xlab = "X", lwd = 4, ylim = c(0, 1))
        mtext( cex=1, side = 3,
               "Given: kernel, m, bandwidth (b)")
        curve(dchisq(x, 1), add = T, lwd = 4, lty = 2)
        lines(bartlett_lug_density, col = "blue", lwd = 4)
        
        
        # Standard criticval value 
        cv <- qchisq(1-alpha, 1)
        lines(c(cv, cv), c(0, .20), col = "black", lty = 3, lwd = 4)
        text(cv , .22, round(cv, 2))
        
        # Bartlett CV
        lines(c(bart_cv, bart_cv), c(0, .24), col = "red", lty = 3, lwd = 4)
        text(bart_cv , .26, round(bart_cv, 2), col = "red")

        # Bartlett_lugsail CV
        lines(c(bart_lug_cv, bart_lug_cv), 
              c(0, .20), col = "blue", lty = 3, lwd = 4)
        text(bart_lug_cv , .22, round(bart_lug_cv, 2), col = "blue")
        
        # Legend 
        legend("topright", 
               legend = c("Small-b", 
                          "Fixed-b: Bartlett", 
                          "Fixed-b: Lugsail"), 
               lty = c(2, 1, 1), 
               col = c("black", "red", "blue"), 
               lwd = 4)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
