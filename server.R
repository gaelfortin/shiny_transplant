library(tidyverse)
library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)

#######Gerstung loading

load("multistage.RData", envir=globalenv())
cr <<- cr
set1 <- brewer.pal(8, "Set1")
pastel1 <- brewer.pal(8, "Pastel1")
s <- !crGroups %in% c("Nuisance","GeneGene")  & ! names(crGroups) %in% c("ATRA","VPA")
VARIABLES <- names(crGroups)[s]  ## nom des variables du modele
rg <- c("Fusions"=5, "CNA"=4,"Genetics"=3, "Clinical"=7, "Demographics"=8, "Treatment"=6)
o <- order(rg[crGroups[s]],((coef(coxRFXPrdTD)^2/diag(coxRFXPrdTD$var2) + coef(coxRFXNrdTD)^2/diag(coxRFXNrdTD$var2) + coef(coxRFXRelTD)^2/diag(coxRFXRelTD$var2)) * apply(data[names(crGroups)], 2, var))[VARIABLES], decreasing=TRUE)
VARIABLES <- VARIABLES[o]
NEWGRP <- c(0,diff(as.numeric(as.factor(crGroups))[s][o])) != 0
names(NEWGRP) <- VARIABLES
INTERACTIONS <- names(crGroups)[crGroups %in% "GeneGene"] 
NUISANCE <- names(crGroups)[crGroups %in% "Nuisance" | names(crGroups) %in%  c("ATRA","VPA")] 

SCALEFACTORS<- rep(1, length(VARIABLES))
names(SCALEFACTORS) <- VARIABLES
w <- crGroups[VARIABLES] %in% c("Demographics","Clinical")
r <- regexpr("(?<=_)[0-9]+$", VARIABLES[w], perl=TRUE)
SCALEFACTORS[w][r!=-1] <- as.numeric(regmatches(VARIABLES[w],r)) # redefini les scales factors, notamment pour les variables quantitatives (sont marques dans le nom des variables)

#######App loading
oldata <- read.csv("donneesAML.txt", h = T, sep = "\t", stringsAsFactor = T) ## donnees mises a disposition par Gerstung
panel_structure <- read_csv("panel_structure.csv")

widget_maker <- function(name, label, type, values, default_value, boundaries){
  if (type == "factor") {
    choices <- unlist(str_split(values, ","))
    default_value <- if_else(is.na(default_value), "NA", default_value)
    radioButtons(inputId = name, label = label, choices = choices, selected = default_value)
  } else {
    limits <- unlist(str_split(boundaries, "\\-"))
    numericInput(inputId = name, label, min = limits[1], max = limits[2], default_value)
  }}


eln_widget <- as.character(panel_structure %>% 
                filter(name == "eln17"))
mrd_widget <- as.character(panel_structure %>% 
                filter(name == "MRD"))

wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px;"

# as_tibble(t(unlist(input)))


shinyServer(function(input, output) {
  # GUI menus to input variables
    output$expandClinical <- renderUI({
      conditionalPanel(
        condition = 'input.showClinical % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Clinical Data") %>% 
                 select(-panel), widget_maker),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
      })
    output$expandCytogenetics <- renderUI({
      conditionalPanel(
        condition = 'input.showCytogenetics % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Cytogenetics") %>% 
                 select(-panel), widget_maker),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
    })
    output$expandGenetics_core<- renderUI({
      conditionalPanel(
        condition = 'input.showGenetics_core % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Genetics (core)") %>% 
                 select(-panel), widget_maker),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
    })
    output$expandGenetics_myeloid<- renderUI({
      conditionalPanel(
        condition = 'input.showGenetics_myeloid % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Genetics (myeloid panels)") %>% 
                 select(-panel), widget_maker),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
    })
    output$expandGenetics_rare<- renderUI({
      conditionalPanel(
        condition = 'input.showGenetics_rare % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Genetics (rare mutations)") %>% 
                 select(-panel), widget_maker),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
    })
    output$expandeln17<- renderUI({
      conditionalPanel(
        condition = 'input.showeln17 % 2',
        wellPanel(
          widget_maker(eln_widget[2], eln_widget[3], eln_widget[4], eln_widget[5], eln_widget[6], eln_widget[7]),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
    })
    output$expandmrd<- renderUI({
      conditionalPanel(
        condition = 'input.showmrd % 2',
        wellPanel(
          widget_maker(mrd_widget[2], mrd_widget[3], mrd_widget[4], mrd_widget[5], mrd_widget[6], mrd_widget[7]),
          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        ))
    })
    
    ##### Merge all inputs for preparation
    observeEvent(input$compute, {
      results <- 
        as_tibble(t(unlist(reactiveValuesToList(input)))) %>% 
          select(any_of(panel_structure$name)) %>% 
          write_csv("results.csv")
          
      
      
      output$resultsdata <- renderTable({
        results
      })
    })
    
    
    
    
    
    
    
})
