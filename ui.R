library(tidyverse)
library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)




shinyUI(fluidPage(
  titlePanel("HSCT prediction"),
  fluidRow(
    column(3,
           style = "border: 1px solid rgb(205,205,205);border-radius: 10px; background:rgb(245,245,245);",
           h4("Enter variables"),
        #Clinical data panel
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showClinical", 					
                     tags$div("Clinical variables", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
          uiOutput("expandClinical"),
        #Cytogenetics data panel
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showCytogenetics", 					
                     tags$div("Cytogenetic variables", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
        uiOutput("expandCytogenetics"),
        #Genetics (core) panel
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showGenetics_core", 					
                     tags$div("Genetics (core) variables", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
        uiOutput("expandGenetics_core"),
        #Genetics (myeloid) panel
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showGenetics_myeloid", 					
                     tags$div("Genetics (myeloid) variables", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
        uiOutput("expandGenetics_myeloid"),
        #Genetics (rare) panel
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showGenetics_rare", 					
                     tags$div("Genetics (rare) variables", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
        uiOutput("expandGenetics_rare"),
        #ELN2017 classification
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showeln17", 					
                     tags$div("ELN 2017 classification", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
        uiOutput("expandeln17"),
        #NPM1 MRD
        wellPanel(
          style = "background:rgb(255,255,255);",
          actionLink("showmrd", 					
                     tags$div("NPM1 MRD", HTML("&#9662;")),
                     style = "color:rgb(0,0,0);")),
        uiOutput("expandmrd"),
      
      #action compute
      actionButton("compute", tags$b("Compute prediction"), class="btn btn-primary", style = "margin-bottom:20px"),
      
    ),
    
    
    
    
    
    
    column(7,
           h4("Prediction result"),
           useShinyalert(),
            tableOutput("resultsdata"),
              textOutput("hsct_prediction")
              )
  )
))

