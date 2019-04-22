library(shiny)
library(ggplot2)
library(MASS)
source("clustering_functions.R")

set.seed(101)
# Parameters
mu_1 <- c(-3, 0.1); mu_2 <- c(6.0, 0.1)
Sigma1 <- matrix(c(10, 7, 7, 10), nrow = 2)
# Generate two samples
sample_1 <- mvrnorm(n=200, mu=mu_1, Sigma = Sigma)
sample_2 <- mvrnorm(n=200, mu=mu_2, Sigma = Sigma)

sample <- rbind(sample_1, sample_2) # matrix format
colnames(sample) <- c("X1", "X2")
classes <- c(rep(1, 200), rep(2, 200))

X <- sample
X <- t(X)
R <- matrix(0, nrow=ncol(X), ncol=2)
centroid <- matrix(c(-5,5,15,0),ncol=2,byrow=TRUE)

theta.init <- list(0.5, 0.5, c(0,0), c(1,1), diag(2), diag(2)) # initializations for em
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
                  max = 15,
                  value = 0,
                  animate = animationOptions(600, loop=TRUE))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data",
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