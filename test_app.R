library(shiny)
library(ggplot2)
library(tidyverse)

bcl <- read.csv("bcl-data.csv", stringsAsFactors = FALSE)


ui <- fluidPage(
  titlePanel("Looking at BCL"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("priceInput", "Price", 
                  min = 0, max = 100, 
                  value = c(25, 40), pre = "$"),
      radioButtons("typeInput", "Product type",
                   choices = c("BEER", "REFRESHMENT","SPIRITS","WINE"),
                   selected = "WINE"),
      uiOutput("countryOutput")
      
      ),
    mainPanel("outputs will go here!",
              plotOutput("coolplot"),
              tableOutput("results"))
  )
  
)
server <- function(input, output) {
  output$countryOutput <- renderUI({
    selectInput(
    "countryInput", "Country",
    sort(unique(bcl$Country)),
    selected = "CANADA"
    )
  })
  filtered <- reactive({
    if (is.null(input$countryInput)) {
      return(NULL)
    }
    
    bcl %>%
      filter(Price >= input$priceInput[1],
             Price <= input$priceInput[2],
             Type == input$typeInput,
             Country == input$countryInput
      )
  })
  
  output$coolplot <- renderPlot({
    if (is.null(filtered())) {
      return()
    }
    ggplot(filtered(), aes(Alcohol_Content)) +
      geom_histogram()
  })
  output$results <- renderTable({
    filtered()
  })
  
}


shinyApp(ui = ui, server = server)