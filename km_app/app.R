library(shiny)
library(ggplot2)
library(MASS)
source("clustering_functions.R")

theta.init <- list(0.5, 0.5, c(0,0), c(1,1), diag(2), diag(2)) # initializations for em
Sigma <- matrix(c(10, 7, 7, 10), nrow = 2)
# Define UI ----
ui <- fluidPage(
  
  # Title ----
  titlePanel("Comparing K-Means and EM Clustering"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "iters",
                  label = "Number of iterations:",
                  min = 0,
                  max = 10,
                  value = 0),
      # Choice for mean of cluster 1    
      sliderInput(inputId = "mu11",
                  label = "Mean for class 1 (horizontal):",
                  min = -10,
                  max = 10,
                  value = 0.1
                   ),
      sliderInput(inputId = "mu12",
                  label = "Mean for class 1 (vertical):",
                  min = -10,
                  max = 10,
                  value = 0.1
                  ),
      # Choice for mean of cluster 2
      sliderInput(inputId = "mu21",
                  label = "Mean for class 2 (horizontal):",
                  min = -10,
                  max = 10,
                  value = 11.0
      ),
      sliderInput(inputId = "mu22",
                  label = "Mean for class 2 (vertical):",
                  min = -10,
                  max = 10,
                  value = 0.1
      ),
      
      # Choice for Sigma of cluster 1  
      sliderInput(inputId = "sigma11",
                  label = "Covariance for class 1 (main diagonal):",
                  min = 1,
                  max = 10,
                  value = 10
      ),
      sliderInput(inputId = "sigma12",
                  label = "Covariance for class 1 (off diagonal):",
                  min = 1,
                  max = 10,
                  value = 7
      ),
      # Choice for Sigma of cluster 2 
      sliderInput(inputId = "sigma21",
                  label = "Covariance for class 2 (main diagonal):",
                  min = 1,
                  max = 10,
                  value = 10
      ),
      sliderInput(inputId = "sigma22",
                  label = "Covariance for class 2 (off diagonal):",
                  min = 1,
                  max = 10,
                  value = 7
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("",
                 fluidRow(
                   column(8, plotOutput("km")),
                   column(8, plotOutput("em"))
                 ))
      )
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  output$km <- renderPlot({
    set.seed(101)  
    # Run KM for widget input
    Sigma_1 <- matrix(c(input$sigma11, input$sigma12, input$sigma12, input$sigma11), nrow = 2)
    Sigma_2 <- matrix(c(input$sigma21, input$sigma22, input$sigma22, input$sigma21), nrow = 2)
    sample <- mv_data(mu_1 = c(input$mu11, input$mu12),
                      mu_2 = c(input$mu21, input$mu22),
                      Sigma_1 = Sigma_1, Sigma_2 = Sigma_2)
    X <- t(sample)
    R <- matrix(0, nrow=ncol(X), ncol=2)
    centroid <- matrix(c(5,5,10,-10), ncol=2,byrow=TRUE)
    
    r <- km_(X, centroid, R, max.iter = input$iters) #input from widget
    
    data <- data.frame("X1"=X[1,], "X2"=X[2,], "Class"=factor(r$cluster))
    ggplot(data, aes(x=X1, y=X2, col=Class)) + geom_point() +
      geom_point(aes(x=r$centroid[1,1],y=r$centroid[1,2]),
                 size=5, shape=3, color="red", stroke=2) +
      geom_point(aes(x=r$centroid[2,1],y=r$centroid[2,2]),
                 size=5, shape=3, color="red", stroke=2) +
      theme_bw() + scale_color_brewer(palette="Dark2") +
      ggtitle("K-Means")
    
  })
  
  output$em <- renderPlot({
    
    # Run EM for widget input
    set.seed(101)
    Sigma_1 <- matrix(c(input$sigma11, input$sigma12, input$sigma12, input$sigma11), nrow = 2)
    Sigma_2 <- matrix(c(input$sigma21, input$sigma22, input$sigma22, input$sigma21), nrow = 2)
    sample <- mv_data(mu_1 = c(input$mu11, input$mu12),
                      mu_2 = c(input$mu21, input$mu22),
                      Sigma_1 = Sigma_1, Sigma_2 = Sigma_2)
    
    if (input$iters == 0) {
      pred <- numeric(400)
      data <- data.frame(X1 = sample[,1], X2 = sample[,2], "Class" = as.factor(pred))
      ggplot(data, aes(X1, X2, color = Class)) +
        geom_point() +
        theme_bw() + scale_color_brewer(palette="Dark2") +
        ggtitle("Expectation Maximization")
    } else {
      r <- EM(sample, theta.init, max.iter = input$iters)
      pred <- apply(r$responsibilities, 1, function(row) which.max(row))
      data <- data.frame(X1 = sample[,1], X2 = sample[,2], "Class" =  as.factor(pred))
      ggplot(data, aes(X1, X2, color = Class)) +
        geom_point() + stat_ellipse(type = "norm") +
        theme_bw() + scale_color_brewer(palette="Dark2") +
        ggtitle("Expectation Maximization")
    }
    
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)