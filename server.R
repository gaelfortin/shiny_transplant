library(tidyverse)
library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)

oldata <- read.csv("donneesAML.txt", h = T, sep = "\t", stringsAsFactor = T) ## donnees mises a disposition par Gerstung
panel_structure <- read_csv("panel_structure.csv")



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




shinyServer(function(input, output) {
    observe({ print(input$AOD_10) })

})
