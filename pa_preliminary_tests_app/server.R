#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)
library(dplyr)


load("../multistage.RData", envir=globalenv())
cr <<- cr
set1 <- brewer.pal(8, "Set1")
pastel1 <- brewer.pal(8, "Pastel1")
s <- !crGroups %in% c("Nuisance","GeneGene")  & ! names(crGroups) %in% c("ATRA","VPA")
VARIABLES <- names(crGroups)[s] 
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
SCALEFACTORS[w][r!=-1] <- as.numeric(regmatches(VARIABLES[w],r))

CATEGORIES <- sapply(VARIABLES, function(x){
    if(length(unique(data[,x])) <= 10){
        c <- min(data[,x]):max(data[,x])
        if(all(c %in% 0:1))
            names(c) <- c("absent","present")
        else if(x =="gender")
            names(c) <- c("male","female")
        return(c)
    }
    else
        NULL
})

LABELS <- sapply(VARIABLES, function(x){
    r <- round(range(data[,x]*SCALEFACTORS[x], na.rm=TRUE),1)
    i <- paste0(" [",r[1],"-",r[2],"]")
    paste0(sub(paste0("_",SCALEFACTORS[x],"$"),"",x), ifelse(is.null(CATEGORIES[[x]]),i,""))
})

LIMITS <- sapply(VARIABLES, function(x){
    r <- round(range(data[,x], na.rm=TRUE),1)})

LABELS["AOD_10"] <- sub("AOD", "Age at diagnosis (yr)", LABELS["AOD_10"])
LABELS["LDH_1000"] <- sub("LDH", "Lactic Acid Dehydrogenase (units/l)", LABELS["LDH_1000"])
LABELS["wbc_100"] <- sub("wbc", "White cell count (1e-9/l)", LABELS["wbc_100"])
LABELS["HB_10"] <- sub("HB", "Hemoglobin (g/dl)", LABELS["HB_10"])
LABELS["BM_Blasts_100"] <- sub("BM_Blasts", "Bone marrow blasts (%)", LABELS["BM_Blasts_100"])
LABELS["PB_Blasts_100"] <- sub("PB_Blasts", "Peripheral blood blasts (%)", LABELS["PB_Blasts_100"])
LABELS["platelet_100"] <- sub("platelet", "Platelet count (1e-9/l)", LABELS["platelet_100"])
LABELS["VPA"] <- "VPA (Valproic acid)"
LABELS["transplantCR1"] <- "Allograft in CR1"
LABELS["transplantRel"] <- "Allograft after Relapse"
LABELS["gender"] <- "Gender"
LABELS <- sub("t_*([a-z,0-9]+)_([a-z,0-9]+)", "t(\\1;\\2)", LABELS)
LABELS[crGroups[VARIABLES] %in% c("Fusions","CNA")] <- gsub("_","/",LABELS[crGroups[VARIABLES] %in% c("Fusions","CNA")])
LABELS <- sub("plus","+",LABELS)
LABELS <- sub("minus|^mono","-",LABELS)
LABELS <- sub("(_|/)*other"," (other)", LABELS)
LABELS <- sub("_([0-9a-zA-Z]+)"," (\\1)", LABELS)


COMPVAR <- list(`Allogeneic HSCT`=c(none="none", `in first CR`="transplantCR1", `after relapse`="transplantRel"), `AML type`=c(primary='AML', secondary='sAML',`therapy-related`='tAML',other='oAML')) ## Compound variables (factors)
COMPIDX <- numeric(length(VARIABLES))
names(COMPIDX) <- VARIABLES
COMPIDX[c("transplantRel","oAML")] <- 1 ## Index of last elements for display
VAR2COMP <- unlist(sapply(names(COMPVAR), function(n) rep(n, length(COMPVAR[[n]])))) 
names(VAR2COMP) <- unlist(COMPVAR)

message(VARIABLES[crGroups[VARIABLES]])

shinyServer(function(input, output) {
    

    
    
    output$contents <- renderTable({
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, it will be a data frame with 'name',
        # 'size', 'type', and 'datapath' columns. The 'datapath'
        # column will contain the local filenames where the data can
        # be found.
        inFile <- input$file1
        
        if (is.null(inFile))
            return(NULL)
        
        read.csv(inFile$datapath, header = input$header)
        
        data2 <- output$contents[,1:(ncol(output$contents)-2)]
        
        ###on importe le dataset qui permet d'imputer les donn?es manquantes
        oldata = read.csv("../donneesAML.txt", h = T, sep = "\t", stringsAsFactor = T) ## donnees mises a disposition par Gerstung
        
        data2 <- getData2()
        
        names(oldata)[!(names(oldata) %in% names(data2))]
        
        # remettre les colonnes dans le bon ordre
        data2 = data2[,names(oldata)]
        
        
        
        ###on importe le dataset qui permet d'imputer les donn?es manquantes
        oldata = read.csv("../donneesAML.txt", h = T, sep = "\t", stringsAsFactor = T) ## donnees mises a disposition par Gerstung
        
        names(oldata)[(!names(oldata) %in% names(data2))]
        
        message(names(oldata) %in% names(data2))
        
        # remettre les colonnes dans le bon ordre
        data2 = data2[names(oldata),]
        
        return(data2)
    })
    
    
    
    

    
    
    # transformation en numerique
    
    #data2 <- data2 %>% 
        #mutate(across(everything(), as.numeric))
    
    
        
    
    
    ################### FONCTION DE CALCUL DES SURVIES SELON GERSTUNG ##################################
    
    
    f_calcul = function(newdata,x, models, timeSearch)
    {
        
        # arrondir les valeurs imputees quand necessaire
        load("../multistage.RData", envir=globalenv())
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
            
            
            ## Confidence intervals
            # osLoDiag <- osUpDiag <- rep(NA, length(osDiag))
            
            # if("analytical" %in% input$ciType){
            # PlogP2 <- function(x) {(x * log(x))^2}
            # errOsCr <- kmNrd$r[,2] * PlogP2(kmNrd$inc) * (1-(1-kmRel$inc) * (1-kmPrs$inc))^2 + kmRel$r[,2] * PlogP2(kmRel$inc) * (1-kmPrs$inc)^2* kmNrd$inc^2 +  kmPrs$r[,2] * PlogP2(kmPrs$inc) * (1-kmRel$inc)^2* kmNrd$inc^2 
            # errOsCr <- sqrt(errOsCr / PlogP2(osCr))
            # osUpCr <- osCr ^ exp(2* errOsCr)
            # osLoCr <- osCr ^ exp(-2*errOsCr)
            #segments(z, osLo[z+1] ,z,osUp[z+1], col=1, lwd=2)
            # }
            # if("simulated" %in% input$ciType){
            # Simulate CI
            # nSim <- 200
            # osCrMc <- sapply(1:nSim, function(i){
            # r <- exp(rnorm(5,0,sqrt(c(kmRel$r[,2],kmNrd$r[,2],kmPrs$r[,2], kmNcd$r[,2], kmCr$r[,2]))))
            # nrsCr <- .crAdjust(kmNrd$inc^r[2],  kmRel$inc^r[1], time=x) ## Correct KM estimate for competing risk
            # diffCir <- diff(kmRel$inc^r[1]) * kmNrd$inc[-1]^r[2] ## Correct KM estimate for competing risk							
            # rsCr <- computeHierarchicalSurvival(x = x, diffS0 = diffCir, S1Static = prsP, haz1TimeDep = tdPrmBaseline * exp(kmPrs$r[,1]-kmPrs$r0+log(r[3])))
            # osCr <- 1-(1-nrsCr)-(1-rsCr)
            
            # cr <- .crAdjust(kmCr$inc^r[5], kmNcd$inc^r[4], time=x)
            # ncd <- .crAdjust(kmNcd$inc^r[4], kmCr$inc^r[5], time=x)
            # osDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = osCr, haz1TimeDep = tdOsBaseline)
            
            # return(cbind(osCr, osDiag - (1-ncd)))
            # }, simplify="array")
            # osCrMcQ <- apply(osCrMc,1:2,quantile, c(0.025,0.975))
            # osLoCr <- osCrMcQ[1,,1]
            # osUpCr <- osCrMcQ[2,,1]
            # osLoDiag <- osCrMcQ[1,,2]
            # osUpDiag <- osCrMcQ[2,,2]
            # }
            
            #absolutePredictions=data.frame(x=x, cr=cr, ncd=ncd, osDiag=osDiag, nrsDiag=nrsDiag, rsDiag=rsDiag, relDiag=relDiag, osCr=osCr, nrsCr=nrsCr, relCr=relCr, rsCr=rsCr, osUpCr=osUpCr, osLoCr=osLoCr, osLoDiag=osLoDiag, osUpDiag=osUpDiag)
            
            
            absolutePredictions=data.frame(x=x,  osDiag=osDiag, nrsDiag=nrsDiag, rsDiag=rsDiag, relDiag=relDiag, osCr=osCr, nrsCr=nrsCr, relCr=relCr, rsCr=rsCr, ncd = ncd)
            
            return(absolutePredictions[x == timeSearch,])
            
            #return(absolutePredictions)
            
            
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
    
    

    
    output$text_out <- renderText({
        # imputation des donnees manquantes
        
        inFile <- input$file1
        
        if (is.null(inFile))
            return(NULL)
        
        patient.file <- read.csv(inFile$datapath, header = input$header)
        
        data2 <- output$patient.file[,1:(ncol(output$patient.file)-2)]
        
        ###on importe le dataset qui permet d'imputer les donn?es manquantes
        oldata = read.csv("../donneesAML.txt", h = T, sep = "\t", stringsAsFactor = T) ## donnees mises a disposition par Gerstung
        
        data2 <- getData2()
        
        names(oldata)[!(names(oldata) %in% names(data2))]
        
        # remettre les colonnes dans le bon ordre
        data2 = data2[,names(oldata)]
        
        
        
        ###on importe le dataset qui permet d'imputer les donn?es manquantes
        oldata = read.csv("../donneesAML.txt", h = T, sep = "\t", stringsAsFactor = T) ## donnees mises a disposition par Gerstung
        
        names(oldata)[(!names(oldata) %in% names(data2))]
        
        message(names(oldata) %in% names(data2))
        
        # remettre les colonnes dans le bon ordre
        data2 = data2[names(oldata),]
        
        set.seed(123)
        

        data2 <- output$contents
        
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
            patient$MRD=="good" ~0,
            patient$MRD=="bad" ~1,
            patient$eln17==1 ~0,
            patient$eln17==2 ~1,
            patient$eln17==3 ~1,
            missing=NULL
        )
        
        
        
        HSCT <- if_else(gerstung==1 & eln17_candidates==1,"candidate","not candidate", missing=NULL)
        
        
        ###et on retourne vers l'app pour afficher
        text_out <- paste0("According to Fenwarth et al, the patient is ", HSCT , " for HSCT in first Complete Remission")
        return(text_out)
    })
    
    
    dynamicWellExpand <- function(condition, variables, style="", before=NULL){
        conditionalPanel(condition = condition,
                         wellPanel(
                             before,
                             lapply(variables, makeMenu),
                             style = style
                         )
        )
    }
    
    
    output$expandClinical <- renderUI({
        dynamicWellExpand('input.showClinical % 2', 
                          variables=VARIABLES[crGroups[VARIABLES] %in% c("Clinical","Demographics")], 
                          before=tags$em(tags$b(crGroups[VARIABLES[1]])),
                          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset"))
    })
    output$expandDrivers <- renderUI({
        dynamicWellExpand('input.showDrivers % 2', 
                          variables = VARIABLES[crGroups[VARIABLES] %in% c("Genetics","Fusions","CNA")],
                          style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        )
    })
    output$expandTreatment <- renderUI({
        dynamicWellExpand('input.showTreatment % 2', 
                          variables=VARIABLES[crGroups[VARIABLES] == "Treatment"],
                          style = paste(wellStyle,"margin-bottom:0px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
        )
    })

})
