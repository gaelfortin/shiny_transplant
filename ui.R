library(tidyverse)
library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)

# panel_structure <- read_csv("panel_structure.csv")

wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px;"


shinyUI(fluidPage(
  titlePanel("HSCT predictions"),
  sidebarLayout(
    sidebarPanel(h2("Enter variables"),
      #Clinical data panel
      wellPanel(
        actionLink("showClinical", 					
                   tags$div("Clinical variables", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
        uiOutput("expandClinical"),
      #Cytogenetics data panel
      wellPanel(
        actionLink("showCytogenetics", 					
                   tags$div("Cytogenetic variables", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
      uiOutput("expandCytogenetics"),
      #Genetics (core) panel
      wellPanel(
        actionLink("showGenetics_core", 					
                   tags$div("Genetics (core) variables", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
      uiOutput("expandGenetics_core"),
      #Genetics (myeloid) panel
      wellPanel(
        actionLink("showGenetics_myeloid", 					
                   tags$div("Genetics (myeloid) variables", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
      uiOutput("expandGenetics_myeloid"),
      #Genetics (rare) panel
      wellPanel(
        actionLink("showGenetics_rare", 					
                   tags$div("Genetics (rare) variables", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
      uiOutput("expandGenetics_rare"),
      #ELN2017 classification
      wellPanel(
        actionLink("showeln17", 					
                   tags$div("ELN 2017 classification", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
      uiOutput("expandeln17"),
      #NPM1 MRD
      wellPanel(
        actionLink("showmrd", 					
                   tags$div("NPM1 MRD", HTML("&#9662;")),
                   style = "color:rgb(0,0,0);")),
      uiOutput("expandmrd"),
      
    ),
    
    
    
    
    
    
    mainPanel(h2("Prediction"))
  )
))

