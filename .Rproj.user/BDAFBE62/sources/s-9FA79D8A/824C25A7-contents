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



# needed before app -------------------------------------------------------

# Simulations Settings
big_T = 1000
rho = .7
num_sim = 10

# Function to generate simulation 
generate_data = function(rho_y=.7, e_sd=1, big_T=200){
    # Generate the data.  Initial value is 0
    sim_data = rep(0, big_T)
    
    # The rest of the values
    for(t in 2:big_T){
        sim_data[t] = rho_y*sim_data[c(t-1)] + rnorm(1, 0, e_sd)
    }
    
    return(sim_data)
}

# True variance matrix 
TRUE_LRV = (1/(1-rho^2))*(1+rho)/(1-rho)


# Estimate autocovariance for a given lag
R = function(h, one_sim_data){
    big_T = length(one_sim_data)
    index = 1:(big_T -h)
    est = (one_sim_data[index])*(one_sim_data[(index+h)])
    
    autocov_s = sum(est)/big_T
    
    if(h!=0){
        autocov_s = 2*autocov_s
    }
    
    return(autocov_s)
}

# All estimated autocovariances 
all_R = function(one_sim_data){
    big_T = length(one_sim_data)
    all_auto_cov <- sapply(0:(big_T-1), R, one_sim_data= one_sim_data)
    return(all_auto_cov)
}



# bartlett function
bartlett = function(x){
    if(abs(x)<1){
        k_x = 1-abs(x)
    } else{
        k_x = 0 
    }
    return(k_x)
}


# Lugsail Bartlett 
lugsail_bartlett = function(x, r = 2, c = 0.5){
    y1 = bartlett(x)/(1-c) 
    y2 = 0 
    
    if(abs(x) < 1/r){
        y2 = bartlett(x*r)*c/(1-c) 
    }
    y = y1- y2
    return(y)
}


# Get estimator for covariance
get_estimator = function(method = "bartlett", sim_data, all_autocovariances, try_b){
    W = matrix(0, nrow = length(try_b), ncol = nrow(all_autocovariances))
    
    for(i in 1:length(try_b)){
        b = try_b[i]
        M = b*nrow(sim_data)
        if(method =="bartlett"){
            new_weights = sapply(0:(M)/(M+1), bartlett) 
        } else{
            new_weights = sapply(0:(M)/(M+1), lugsail_bartlett)
        }
        W[i, 1:length(new_weights)] = new_weights
    }
    rownames(W) = paste("b=", round(try_b, 3), sep = "")
    colnames(W) = paste("Lag=", 0:(nrow(sim_data)-1), sep = "")
    
    
    omega= W %*% all_autocovariances
    colnames(omega) = paste("Sim", 1:num_sim)
    rownames(omega) = paste("b=", round(try_b, 3))
    
    return(omega)
}


# Get Critical values 
lugsail_cv = function(b){
    cv = 2.73 + 47.575*b - 57.656*b^2 + 88.866*b^3
}

bartlett_cv = function(b){
    cv = 2.8398 + 11.5731*b + 13.5167*b^2 - 4.6680*b^3
}




# Cov_list should be 
Make_CIs = function(Cov_list, b, the_means, big_T){
    
    keep = which(Cov_list$lugsail_est>=0)[1:5]
    keep_means = the_means[keep]
    bartlett_est = sqrt(Cov_list$bartlett_est[keep]/big_T)
    lugsail_est = sqrt(Cov_list$lugsail_est[keep]/big_T)
    
    LB_truth = keep_means - sqrt(Cov_list$truth/big_T)*sqrt(qchisq(0.95, 1))
    UB_truth = keep_means + sqrt(Cov_list$truth/big_T)*sqrt(qchisq(0.95, 1))
    truth_ci = data.frame(LB= LB_truth, UB= UB_truth)
    
    LB_bartlett = keep_means - bartlett_est*sqrt(bartlett_cv(b))
    UB_bartlett = keep_means + bartlett_est*sqrt(bartlett_cv(b))
    bartlett_ci = data.frame(LB= LB_bartlett, UB= UB_bartlett)
    
    LB_lugsail = keep_means - lugsail_est*sqrt(lugsail_cv(b))
    UB_lugsail = keep_means+ lugsail_est*sqrt(lugsail_cv(b))
    lugsail_ci = data.frame(LB= LB_lugsail, UB= UB_lugsail)
    
    ci = list(truth = truth_ci, bartlett = bartlett_ci, lugsail = lugsail_ci)
    return(ci)
}

plot_ci = function(ci){
    x_max = max(abs(range(unlist(ci))*1.1))
    x_range = c(-x_max, x_max)
    plot(x= 0 , y =0, ylim = c(0, 6), xlim = x_range, yaxt= "n", 
         ylab = "Simulation", col = "white")
    abline(v = 0, col = "black", lty = 2)
    
    for(i in 1:5){
        lines(ci$truth[i,], c(c(i,i)-.95), lwd = 2)
    }
    
    for(i in 1:5){
        lines(ci$bartlett[i,], c(c(i,i)-.85), col = "red", lwd = 2)
    }
    
    for(i in 1:5){
        lines(ci$lugsail[i,], c(c(i,i)-.75), col = "blue", lwd = 2)
    }
    
    legend("topleft", 
           c("True Mean", 
             "True Confidence Interval", 
             "Bartlett (fixed b)", 
             "Lugsail (fixed b)"), 
           lty = c(2, 1, 1, 1), 
           col = c("black", "black", "red", "blue"), 
           lwd = c(1, 2, 2, 2))
}


put_it_together = function(ci){
    plot_ci(ci)
}



# this is app -------------------------------------------------------------




# Define UI for application that draws a histogram
ui <- fluidPage(
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("b_value",
                        "Choose number of lags (M/T):",
                        min = 0,
                        max = .5,
                        value = .032, 
                        step = .005), 
            numericInput("big_T", 
                         "Choose sample size (T):", 
                         min = 50, max = 10000, 
                         step = 50, value = 1000),
            actionButton("action1", label = "(re)Generate")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel( 
                tabPanel(
                    p("Model", style = "color:black"),   
                    p("Example with autoregressive model with order 1 (AR(1)). ", style = "color:black"),
                    withMathJax(), 
                    helpText('The Model: 
                        $$ y_{t}=\\rho y_{t-1}+e_{t}$$', 
                             style = "color:black"),
                    helpText("where", style = "color:black"),
                    helpText("$$e_{t} \\overset{iid}{\\sim} N(0, 1)$$", style = "color:black"), 
                    helpText("$$\\rho = 0.7$$", style = "color:black"), 
                    helpText("Recommended to use $$\\frac{M}{T} \\approx T^{-1/2}$$", style = "color:black")),
                tabPanel("Confidence Intervals", plotOutput("ci_plot")),
                tabPanel("Data Visual", plotlyOutput("data_plot"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    observeEvent(input$action1, {
        b = input$b_value
        
        sim_data = replicate(num_sim, generate_data(big_T = input$big_T, 
                                                    rho_y = rho))
        orig_sim_data = sim_data
        the_means = colMeans(sim_data)
        sim_data = apply(sim_data, 1, function(row) row - the_means)
        sim_data = t(sim_data)
        
        big_T = nrow(sim_data)
        
        all_autocovariances <- apply(sim_data, 2, all_R) 
        rownames(all_autocovariances) = paste("Lag=", 0:(nrow(sim_data)-1), sep = "")
        colnames(all_autocovariances) = paste("Sim ", 1:(ncol(sim_data)), sep = "")
        
        bartlett_est = get_estimator(sim_data = sim_data, 
                                     all_autocovariances = all_autocovariances, 
                                     try_b = b)
        
        lugsail_est = get_estimator("Lugsail",
                                    sim_data = sim_data, 
                                    all_autocovariances = all_autocovariances, 
                                    try_b = b)
        
        Cov_list = list(truth = TRUE_LRV, 
                        bartlett_est = bartlett_est, 
                        lugsail_est = lugsail_est)
        
        keep = which(Cov_list$lugsail_est>=0)[1:5]
        
        ci = Make_CIs(Cov_list, b, the_means, big_T)
        
        
        output$ci_plot <- renderPlot({
            
            output$ci_plot <- renderPlot({
                put_it_together(ci)
            })
        })
        
        
        output$data_plot <- renderPlotly({
            plot1 <- plot_ly(
                x = 1:nrow(sim_data),
                y = sim_data[, keep[1]], 
                type = 'scatter',
                mode = 'lines')%>%
                layout(showlegend = F)
            
            plot2 <- plot_ly(
                x = 1:nrow(sim_data),
                y = sim_data[,keep[2]], 
                type = 'scatter',
                mode = 'lines')%>%
                layout(showlegend = F)
            
            plot3 <- plot_ly(
                x = 1:nrow(sim_data),
                y = sim_data[,keep[3]], 
                type = 'scatter',
                mode = 'lines')%>%
                layout(showlegend = F)
            plot4 <- plot_ly(
                x = 1:nrow(sim_data),
                y = sim_data[,keep[4]], 
                type = 'scatter',
                mode = 'lines')%>%
                layout(showlegend = F)
            
            plot5 <- plot_ly(
                x = 1:nrow(sim_data),
                y = sim_data[,keep[5]], 
                type = 'scatter',
                mode = 'lines')%>%
                layout(showlegend = F)
            
            fig <- subplot(plot1, plot2, plot3, plot4, plot5,
                           nrows =5, titleY = TRUE) %>% layout(
                               plot_bgcolor='white', 
                               hovermode="x unified")
        })
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
