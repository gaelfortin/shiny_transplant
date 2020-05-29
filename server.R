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

wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px; margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset"



shinyServer(function(input, output) {
  # GUI menus to input variables
    output$expandClinical <- renderUI({
      conditionalPanel(
        condition = 'input.showClinical % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Clinical Data") %>% 
                 select(-panel), widget_maker),
          style = wellStyle
        ))
      })
    output$expandCytogenetics <- renderUI({
      conditionalPanel(
        condition = 'input.showCytogenetics % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Cytogenetics") %>% 
                 select(-panel), widget_maker),
          style = wellStyle
        ))
    })
    output$expandGenetics_core<- renderUI({
      conditionalPanel(
        condition = 'input.showGenetics_core % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Genetics (core)") %>% 
                 select(-panel), widget_maker),
          style = wellStyle
        ))
    })
    output$expandGenetics_myeloid<- renderUI({
      conditionalPanel(
        condition = 'input.showGenetics_myeloid % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Genetics (myeloid panels)") %>% 
                 select(-panel), widget_maker),
          style = wellStyle
        ))
    })
    output$expandGenetics_rare<- renderUI({
      conditionalPanel(
        condition = 'input.showGenetics_rare % 2',
        wellPanel(
          pmap(panel_structure %>% 
                 filter(panel == "Genetics (rare mutations)") %>% 
                 select(-panel), widget_maker),
          style = wellStyle
        ))
    })
    output$expandeln17<- renderUI({
      conditionalPanel(
        condition = 'input.showeln17 % 2',
        wellPanel(
          widget_maker(eln_widget[2], eln_widget[3], eln_widget[4], eln_widget[5], eln_widget[6], eln_widget[7]),
          style = wellStyle
        ))
    })
    output$expandmrd<- renderUI({
      conditionalPanel(
        condition = 'input.showmrd % 2',
        wellPanel(
          widget_maker(mrd_widget[2], mrd_widget[3], mrd_widget[4], mrd_widget[5], mrd_widget[6], mrd_widget[7]),
          style = wellStyle
        ))
    })
    
    ##### Merge all inputs for preparation
    observeEvent(input$compute, {
      patient_data <- 
        as_tibble(t(unlist(reactiveValuesToList(input)))) %>% 
          select(any_of(panel_structure$name)) 
          
      ##### Transform inputs
      patient_data[patient_data=="Absent"]<-"0"
      patient_data[patient_data=="Wildtype"]<-"0"
      patient_data[patient_data=="Present"]<-"1"
      patient_data[patient_data=="Mutated"]<-"1"
      patient_data[patient_data=="Male"]<-"1"
      patient_data[patient_data=="Favorable"]<-"1"
      patient_data[patient_data=="Female"]<-"2"
      patient_data[patient_data=="Intermediate"]<-"2"
      patient_data[patient_data=="Unfavorable"]<-"3"
      patient_data[patient_data=="< 4 log"]<-"bad"
      patient_data[patient_data=="> 4 log"]<-"good"
      mrd <- patient_data$MRD[1]
      
      ##### Encode additional AML_type variables
      AML <- patient_data$AML_type[1] #de novo,secondary,therapy-related,other,N/A
      patient_data$oAML[1] <- case_when(AML == "de novo" ~ 0,
                                        AML == "secondary" ~ 0,
                                        AML == "therapy-related" ~ 0,
                                        AML == "other" ~ 1,
                                        AML == NA ~ NA_real_)
      patient_data$tAML[1] <- case_when(AML == "de novo" ~ 0,
                                        AML == "secondary" ~ 0,
                                        AML == "therapy-related" ~ 1,
                                        AML == "other" ~ 0,
                                        AML == NA ~ NA_real_)
      patient_data$sAML[1] <- case_when(AML == "de novo" ~ 0,
                                        AML == "secondary" ~ 1,
                                        AML == "therapy-related" ~ 0,
                                        AML == "other" ~ 0,
                                        AML == NA ~ NA_real_)

      ###### Transform numeric values
      patient_data <- mutate_all(patient_data, as.numeric)
      patient_data$AOD_10[1] <- patient_data$AOD_10[1]/10
      patient_data$LDH_1000[1] <- patient_data$LDH_1000[1]/1000
      patient_data$wbc_100[1] <- patient_data$wbc_100[1]/100
      patient_data$platelet_100[1] <- patient_data$platelet_100[1]/100
      patient_data$PB_Blasts_100[1] <- patient_data$PB_Blasts_100[1]/100
      patient_data$BM_Blasts_100[1] <- patient_data$BM_Blasts_100[1]/100
      
      patient_data$MRD[1] <- mrd #Put MRD character value back
      
      ##### Add treatment value (=transplantation)
      patient_data$transplantCR1 <- 0
      patient_data$transplantRel <- 1
      
      patient_data <- patient_data %>% 
        select(-AML_type)
      
      write_csv(patient_data, "patient_data.csv")
      
      ##### Output input values (temporary output)
      output$resultsdata <- renderTable({
        patient_data
      })
      
      
      #####-------COMPUTE GERSTUNG SCORE-------#####
      data2 <- patient_data %>% select(-eln17, -MRD)
      names(oldata)[!(names(oldata) %in% names(data2))]
      # remettre les colonnes dans le bon ordre
      data2 = data2[,names(oldata)]
      # transformation en numerique
      data2 <- data2 %>% 
        mutate(across(everything(), as.numeric))
      
      f_calcul = function(newdata,x, models, timeSearch) #scoring function
      {
        for(i in 1:ncol(newdata))
        {
          if(SCALEFACTORS[colnames(newdata)[i]]==1)
          {
            newdata[, i] = round(newdata[, i])
          }
        }
        
        
        ######################### RAJOUT DES INTERACTIONS DANS LE JEU DE DONNEES D'ANALYSE #################
        ## rajoute les interactions et les nuisances dans le newdata frame
        l = list()
        for(n in INTERACTIONS){
          s <- strsplit(n, ":")[[1]]
          if(length(s) == 2)
          {
            l[[n]] <- newdata[s[1]] * newdata[s[2]]
          }
          else
          {
            l[[n]] <- newdata[s[1]] * newdata[s[2]] * newdata[s[3]]
          }
          names(l[[n]]) = ""
        }
        temp = data.frame(t(unlist(l)))
        names(temp) = names(l)
        newdata = cbind.data.frame(newdata,temp)
        
        l = list()
        for(n in NUISANCE)
        {
          l[[n]] <- NA
          names(l[[n]]) = ""
        }
        temp = data.frame(t(unlist(l)))
        names(temp) = names(l)
        newdata = cbind.data.frame(newdata,temp)
        
        #################################### FONCTIONS POUR LES CALCULS ####################################
        
        
        computeIncidence <- function(coxRFX, r, x) {
          #r=PredictRiskMissing(coxRFX, data, var="var2")
          if(!is.null(coxRFX$na.action))
            coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
          #r <- PredictRiskMissing(coxRFX, data, var="var2")
          
          H0 <- summary(survfit(coxRFX,se.fit = FALSE),times = x)
          hazardDist <- splinefun(H0$time, -log(H0$surv), method="monoH.FC")
          r0 <- coxRFX$means %*% coef(coxRFX)
          lambda0 <- -log(H0$surv)*exp(-r0)
          if(length(lambda0) < length(x)){
            lambda0 <- c(lambda0, hazardDist(x[(length(lambda0)+1):length(x)]))
          }
          inc <- exp(-lambda0* exp(r[,1]))
          ciup2 <- exp(-lambda0*exp( rep(r[,1] + 2*sqrt(r[,2]) * c(1), each=length(x))))
          cilo2 <- exp(-lambda0*exp( rep(r[,1] + 2*sqrt(r[,2]) * c(-1), each=length(x))))
          ciup <- exp(-lambda0*exp( rep(r[,1] + sqrt(r[,2]) * c(1), each=length(x))))
          cilo <- exp(-lambda0*exp( rep(r[,1] + sqrt(r[,2]) * c(-1), each=length(x))))
          #p <- PartialRisk(coxRFX, dataImputed)
          return(list(inc=inc, r=r, x=x, hazardDist=hazardDist, r0 = r0, ciup=ciup, cilo=cilo, ciup2=ciup2, cilo2=cilo2))
        }
        
        
        
        riskMissing <- function(newdata)
        {
          sapply(models, function(m){
            fit <- get(paste("coxRFX",m,"TD", sep=""))
            if(!is.null(fit$na.action))
              fit$Z <- fit$Z[-fit$na.action,]
            PredictRiskMissing(fit, newdata,  var="var2")}, simplify = FALSE)
        }
        
        
        crAdjust <- function(x, y, time=x$x) {
          xadj <- .crAdjust(x$inc, y$inc, time)
        }
        
        .crAdjust <- function(inc1, inc2, time) {
          cumsum(c(1,diff(inc1) * splinefun(time, inc2)(time[-1])))
        }
        
        
        
        
        
        
        survPredict <- function(surv){
          s <- survfit(surv~1)
          splinefun(s$time, s$surv, method="monoH.FC")
        }
        prsP <- survPredict(Surv(prdData$time1, prdData$time2, prdData$status))(x) # Baseline Prs (measured from relapse)
        
        coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])) 
        tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=x[-1]))) ## Hazard (function of CR length)	
        
        coxphOs <- coxph(Surv(time1,time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))) 
        tdOsBaseline <- exp(predict(coxphOs, newdata=data.frame(time0=x[-1])))	 ## Hazard (function of induction length), only for OS (could do CIR,NRM,PRS seperately)
        
        
        
        
        computeAbsoluteProbabilities <- function(i, newdata)
        {
          newdata = newdata[i,]
          
          
          ## KM incidence of NCD and CR
          kmNcd <- computeIncidence(coxRFX = coxRFXNcdTD, r = riskMissing(newdata)[["Ncd"]], x=x)
          kmCr <- computeIncidence(coxRFX = coxRFXCrTD, r = riskMissing(newdata)[["Cr"]], x=x)
          
          ## Correct KM estimate for competing risk
          ncd <- crAdjust(x= kmNcd, time=x, y=kmCr) ## Correct KM estimate for competing risk
          cr <- crAdjust(x= kmCr, time=x, y=kmNcd) ## Correct KM estimate for competing risk
          
          ## KM incidence of Relapse and NRD
          kmRel <- computeIncidence(coxRFX = coxRFXRelTD, r = riskMissing(newdata)[["Rel"]], x=x)
          kmNrd <-  computeIncidence(coxRFX = coxRFXNrdTD, r = riskMissing(newdata)[["Nrd"]], x=x)
          
          ## Correct KM estimate for competing risk
          relCr <- crAdjust(x= kmRel, time=x, y=kmNrd) ## Correct KM estimate for competing risk
          nrsCr <- crAdjust(x = kmNrd, time = x, y = kmRel)
          
          ## KM incidence of PRS
          kmPrs <-  computeIncidence(coxRFX = coxRFXPrdTD, r = riskMissing(newdata)[["Prd"]], x=x)
          
          ## Outcome after Remission
          rsCr <- computeHierarchicalSurvival(x = x, diffS0 = diff(relCr), S1Static = prsP, haz1TimeDep = tdPrmBaseline * exp(kmPrs$r[,1]-kmPrs$r0))
          osCr <- 1-(1-nrsCr)-(1-rsCr)
          
          ## Outcome from diagnosis
          osDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = osCr, haz1TimeDep = tdOsBaseline) - (1-ncd)
          nrsDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = nrsCr, haz1TimeDep = tdOsBaseline)
          rsDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = rsCr, haz1TimeDep = tdOsBaseline)
          relDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = relCr, haz1TimeDep = tdOsBaseline)
          
          
          absolutePredictions=data.frame(x=x,  osDiag=osDiag, nrsDiag=nrsDiag, rsDiag=rsDiag, relDiag=relDiag, osCr=osCr, nrsCr=nrsCr, relCr=relCr, rsCr=rsCr, ncd = ncd)
          
          return(absolutePredictions[x == timeSearch,])
          
          
          
          
        }
        
        
        
        cppFunction('NumericVector computeHierarchicalSurvival(NumericVector x, NumericVector diffS0, NumericVector S1Static, NumericVector haz1TimeDep) {
              int xLen = x.size();
              double h;
              NumericVector overallSurvival(xLen);
              for(int i = 0; i < xLen; ++i) overallSurvival[i] = 1;
              for(int j = 1; j < xLen; ++j){
              if(diffS0[j-1] != 0){
              h = haz1TimeDep[j-1];
              for(int i = j; i < xLen; ++i){
              overallSurvival[i] += diffS0[j-1] * (1-pow(S1Static[i-j], h));
              }
              }
              }
              return overallSurvival;
}')
        
        
        
        
        res = sapply(1:nrow(newdata), computeAbsoluteProbabilities, newdata = newdata)
        colnames(res) = rownames(newdata)
        
        return(res)
      }
      
      #####------5 year survival and HSCT prediction-----#####
      # imputation des donnees manquantes
      set.seed(123)
      data2i = data.frame(ImputeMissing(oldata, data2))
      
      # 5 ans allRel
      timeSearch = 1825 # temps retenu pour la prediction (en jours)
      x <- seq(0,timeSearch,1)#0:2000
      models <- c("Ncd","Cr","Rel","Nrd","Prd")
      res = f_calcul(data2i, x = x, models = models, timeSearch = timeSearch)				
      
      ###la valeur qui m'interesse est osDiag, sauf erreur c'est
      
      osDiag <- res[2,1]
      
      ##le cutoff est 
      
      gerstung <- ifelse(osDiag < 0.30, 1, 0)
      
      
      ### et apres, on n' a plus qu'a lancer l'operateur booleen
      
      
      eln17_candidates <- case_when(
        patient_data$MRD=="good" ~0,
        patient_data$MRD=="bad" ~1,
        patient_data$eln17==1 ~0,
        patient_data$eln17==2 ~1,
        patient_data$eln17==3 ~1,
        missing=NULL
      )
      
      
      
      HSCT <- if_else(gerstung==1 & eln17_candidates==1,"candidate","not candidate", missing=NULL)
      
      
      ###et on retourne vers l'app pour afficher
      output$hsct_prediction <- renderText({
        paste0("According to Fenwarth et al, the patient is ", HSCT , " for HSCT in first Complete Remission")
        })
      
      })
    
    
    
    
    
    
    
})
