library(tidyverse)
library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)

panel_structure <- read_csv("panel_structure.csv")


widget_maker <- function(name, label, type, values, default_value, boundaries){
  if (type == "factor") {
    choices <- unlist(str_split(values, ","))
    default_value <- if_else(is.na(default_value), "NA", default_value)
    radioButtons(inputId = name, label = label, choices = choices, selected = default_value)
  } else {
    limits <- unlist(str_split(boundaries, "\\-"))
    numericInput(inputId = name, label, min = limits[1], max = limits[2], default_value)
  }
  
}


shinyUI(fluidPage(
  titlePanel("HSCT predictions"),
  sidebarLayout(
    sidebarPanel(h2("Enter variables"),
            
                 
      numericInput("priceInput", "Price", 
                  min = 0, max = 100, value = 0),
      
      pmap(panel_structure[,2:ncol(panel_structure)], widget_maker),
      
      
    ),
    mainPanel(h2("Prediction"))
  )
))
