library(tidyverse)
library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)

panel_structure <- read_csv("panel_structure.csv")

wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px;"


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
      wellPanel(
        actionLink("showClinical", 					
                   tags$div("Clinical variables", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
        uiOutput("expandClinical"),
        
    ),
    mainPanel(h2("Prediction"))
  )
))

