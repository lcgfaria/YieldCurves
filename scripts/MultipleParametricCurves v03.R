
# Identification ----------------------------------------------------------
# 
# FGV - Fundacao Getulio Vargas
# EESP - Escola de Economia São Paulo
# MPEF - Mestrado Profissional em Economia e Finanças
#
# Computacao Aplicada (Prof.Ernesto Colla)
# End-course project
#
# YIELD CURVES COMPARISON
# Autor: Luis Giovanni Faria
#
# Date: 13/09/2022
#

# Initialization ----------------------------------------------------------

# Remove working set
rm(list=ls())
graphics.off()

# Working Directory
sysInfo <- Sys.info()
if (sysInfo[[4]] == "GIO-YOGA") {
  workingDirectory <- "C:/Users/giovanni/OneDrive/$FGV/Yield Curves/R/dev"
} else if (sysInfo[[4]] == "DESKTOP-G3EH8EA") {
  workingDirectory <- "C:/Users/lcgfa/OneDrive/$FGV/Yield Curves/R/dev"
} else {
  workingDirectory <- getwd()
}
setwd(workingDirectory)
getwd()

inputDirectory <- "../database/"
outputDirectory <- "../output/"

# Load R and Python libraries
library(readxl)
library(xlsx)
library(reticulate)
library(Matrix)
library(ggplot2)
library(grid)
library(gridExtra)

source("../lib/yc_lib.r")
source("../lib/yc_loaddata.r")
source_python("../lib/kr_model.py")
source_python("../lib/kr_utils.py")
source_python("../lib/yc_utils.py")

nsimDE <- 5 # number of simulations differential evolution
nsimSD <- 5 # number of simulations steepest descent
faceValue <- 1000
daysYear <- 365
datePlots <- list() # list of date reference curve plots 

# Parametric models to be run with optimization methods DE and SD
# NS = Nelson-Siegel
# BL = Bliss
# NSS= Svensson
# DP = DePooter
# BC = Bjork-Christensen
models <- c("NS","BL","NSS","DP","BC") # Models to optimize

##############################################################################
## Load data

# run data from tesouro diretoto
bondData <- ReadTesouroDireto(inputDirectory,coupon=0.1)

#bondSet$DateMaturity <- as.Date(bondSet$DateMaturity)
bondData$YieldMid <- (bondData$RateBuy+bondData$RateSell)/2
bondData <- bondData[c(1,2,3,8,9)]


# vector of reference dates
workingDates <- unique(bondData$DateReference)

# Data frame with results from each model
rmse <- data.frame(DateReference=as.Date(character()),Model=character(),RMSE=double())

workingDates = c(as.Date("2021-09-30"),as.Date("2021-10-29"),as.Date("2021-11-30"),as.Date("2021-12-30"))

# main loop to compute all reference dates
for (k in 1:length(workingDates)) {
  dateReference <- workingDates[k]
  
  cat("=> Working with reference date",format(dateReference,"%Y-%m-%d"),sprintf(" Date %d/%d\n",k,length(workingDates)),"\n")
  
  # Set to generate price vector and cashfow matrix (subset of BondPrinces) 
  bondSet <- subset(bondData, DateReference==dateReference)

  # view working set
  print(bondSet)
  
  ##############################################################################
  # Compute YTM for bonds
  
  ## Calculate yield to maturity
  bondSet$Yield <- c(0)

  for (i in 1:nrow(bondSet)) {
    bondSet$Yield[i] <- YieldToMaturity(price=bondSet$PricePU[i],x0=0.00,
                                        days=as.numeric(rev(seq(from=bondSet$DateMaturity[i],to=dateReference, by="-6 month"))-dateReference),
                                        faceValue=faceValue,
                                        coupon=bondSet$Coupon[i]/2,
                                        daysYear=daysYear,epson=1e-6,nmax=100)
  } # end-for
  
  # view working set  
  print(bondSet)
  
  ###############################################################################
  ## Build Constrained Cubic Spline
  ##
  
  # Initialize yield vectors and CCS parameters
  bondSet$CCS <- c(0)
  # sort by maturity
  bondSort <- bondSet[order(bondSet$DateMaturity,decreasing=FALSE),]  
  
  x <- as.numeric(bondSort$DateMaturity-bondSort$DateMaturity[1])/nrow(bondSet)
  y <- bondSort$Yield
  
  # Compute fitted CCS yields for maturity dates
  for (i in 1:nrow(bondSet)) {
    xstar <- (x[length(x)]-x[1])*as.numeric(bondSet$DateMaturity[i]-bondSet$DateMaturity[1])/
      as.numeric(bondSet$DateMaturity[nrow(bondSet)]-bondSet$DateMaturity[1])
    bondSet$CCS[i] <- CCSpline(x,y,xstar)
  }
  
  # Generate CCS curve
  numberOfPoints <- 100
  maxMaturity <- max(bondSet$DateMaturity)
  CCScurve <- data.frame(matrix(ncol=3, nrow=numberOfPoints))
  colnames(CCScurve) <- c("X", "DateMaturity", "CCS")
  class(CCScurve$DateMaturity) <- "Date"
  
  for (i in 1:numberOfPoints) {
    xstar <- x[1]+(x[length(x)]-x[1])/numberOfPoints*(i-1)
    CCScurve[i,1] <- xstar
    CCScurve[i,2] <- as.Date(bondSet$DateMaturity[1]+
                               as.numeric(bondSet$DateMaturity[nrow(bondSet)]-bondSet$DateMaturity[1])/numberOfPoints*(i-1))
    CCScurve[i,3] <- CCSpline(x,y,xstar)
  } # end-for
  
  # store errors
  rmse[nrow(rmse)+1,1] <- dateReference
  rmse[nrow(rmse),2] <- "CSS"
  rmse[nrow(rmse),3] <- sqrt(sum((bondSet$Yield-bondSet$CCS)^2))
  
  ###############################################################################
  ## Build parametric curves
  ##
  
  ## Define model functions as a list
  ModelFunction = list();

  ModelFunction[["NS"]] = function(params,t) {
    ## Description: Compute Nelson-Siegel (1987) function
    return(params["Beta1"]+
             params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
             params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"])))
  } # end-function
  
  ModelFunction[["DPR"]] = function(params,t) {
    ## Description: Compute Diebold-Piazzesi-Rudesbusch (2005) function
    return(params["Beta1"]+
             params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])))
  } # end-function  

  ModelFunction[["BL"]] = function(params,t) {
    ## Description: Compute Bliss (1997) function
    return(params["Beta1"]+
             params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
             params["Beta3"]*((1-exp(-t/params["Lambda2"]))/(t/params["Lambda2"])-exp(-t/params["Lambda2"])))
  } # end-function
  
  ModelFunction[["NSS"]] = function(params,t) {
    ## Description: Compute Nelson-Siegel-Svensson (1994) function
    return(params["Beta1"]+
             params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
             params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"]))+
             params["Beta4"]*((1-exp(-t/params["Lambda2"]))/(t/params["Lambda2"])-exp(-t/params["Lambda2"])))
  } # end-function
  
  ModelFunction[["DP"]] = function(params,t) {
    ## Description: Compute De Pooter (2007) function
    return(params["Beta1"]+
             params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
             params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"]))+
             params["Beta4"]*((1-exp(-t/params["Lambda2"]))/(t/params["Lambda2"])-exp(-2*t/params["Lambda2"])))
  } # end-function
  
  ModelFunction[["BC"]] = function(params,t) {
    ## Description: Compute Björk-Christensen (1999) function
    return(params["Beta1"]+
             params["Beta2"]*t/(2*params["Lambda1"])+
             params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
             params["Beta4"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"]))+
             params["Beta5"]*((1-exp(-2*t/params["Lambda1"]))/(2*t/params["Lambda1"])))
  } # end-function
  

  ###############################################################################
  ## Run Differential Evolution for  parametric curves
  ##
  
  NP <- 100 # population size
  nmaxG = 500 # number of generations
  set.seed(0)
  lower = c(Beta1=0,Beta2=-0.05,Beta3=-0.5,Beta4=-0.5,Beta5=-0.5,Lambda1=0,Lambda2=0);
  upper = c(Beta1=1,Beta2=0.05,Beta3=0.5,Beta4=0.5,Beta5=-0.5,Lambda1=1.25,Lambda2=1.25);
  
  ### Restrictionsestrições
  gx = list();
  gx[[1]] = function(params){return(-params["Beta1"]-params["Beta2"])}
  
  nsim <- nsimDE # number of simulations
  
  for (m in models) {
    cat("Reference date ",format(dateReference,"%Y-%m-%d")," Calculating Model - Differential Evolution = ",m,"\n")
    # data structure to store results
    paramsOpt = matrix(NA,nrow=7,ncol=nsim)
    errors <- c()
    rownames(paramsOpt) = c("Beta1","Beta2","Beta3","Beta4","Beta5","Lambda1","Lambda2")
    
    for (i in 1:nsim) {
      cat(sprintf("Simulation %d/%d\n",i,nsim));
      # Initial population
      paramsG0 <- t(as.matrix(data.frame(Beta1 = runif(NP,lower["Beta1"],upper["Beta1"]),
                                         Beta2 = runif(NP,lower["Beta2"],upper["Beta2"]),
                                         Beta3 = runif(NP,lower["Beta3"],upper["Beta3"]),
                                         Beta4 = runif(NP,lower["Beta4"],upper["Beta4"]),
                                         Beta5 = runif(NP,lower["Beta5"],upper["Beta5"]),
                                         Lambda1 = runif(NP,lower["Lambda1"],upper["Lambda1"]),
                                         Lambda2 = runif(NP,lower["Lambda2"],upper["Lambda2"]))))
      # get optimum parameters
      paramsOpt[,i] <- DifferentialEvolution(RootMeanSquareError,ModelFunction[[m]],paramsG0,
                                             t=as.numeric(as.Date(bondSet$DateMaturity)-dateReference)/daysYear,
                                             yield=bondSet$Yield,
                                             F = 1,
                                             CR = 0.5,
                                             lower=lower,
                                             upper=upper,
                                             gx = gx,
                                             nmaxG=500)
      # get errors
      errors[i] <- RootMeanSquareError(ModelFunction[[m]],paramsOpt[,i],
                                       t=as.numeric(as.Date(bondSet$DateMaturity)-dateReference)/daysYear,
                                       yield=bondSet$Yield)   
    } # end-simulations
    # get fitted models
    rmse[nrow(rmse)+1,1] <- dateReference
    rmse[nrow(rmse),2] <- m
    rmse[nrow(rmse),3] <- min(errors)
    
    bondSet[,ncol(bondSet)+1] <- ModelFunction[[m]](params=paramsOpt[,which.min(errors)],
                                                    t=as.numeric(as.Date(bondSet$DateMaturity)-dateReference)/daysYear)
    colnames(bondSet)[ncol(bondSet)] <- m
  } # end-for
  
  print(rmse)
  print(bondSet)
  
  ###############################################################################
  ## Run Gradient method for  parametric curves
  ##
  
  ## Number of simulations
  nsim <- nsimSD
  
  ## Set seed
  set.seed(0)
  
  for (m in models) {
    cat("Reference date ",format(dateReference,"%Y-%m-%d")," Calculating Model - Steepest Descent = ",m,"\n")
    ## Intervals
    params0 <- as.matrix(data.frame(Beta1 = runif(nsim,0,1),
                                    Beta2 = runif(nsim,-1,1),
                                    Beta3 = runif(nsim,-1,1),
                                    Beta4 = runif(nsim,-1,1),
                                    Beta5 = runif(nsim,-1,1),
                                    Lambda1 = runif(nsim,1,5),
                                    Lambda2 = runif(nsim,1,5)));
    ## data structure with parameters in each interaction
    params <- matrix(NA,nrow(params0),ncol(params0));
    colnames(params) = colnames(params0);
    ## data structure with errors
    errors <- c();
    ## data structure with fitted curves
    curves <- matrix(NA,nrow(bondSet),nsim);
    for (i in 1:nsim){
      cat(sprintf("Simulation %d/%d\n",i,nsim));
      #RmseFunc=paste("RootMeanSquareError",m,sep="");
      params[i,] <- OptimizeSteepestDescent(RootMeanSquareError,ModelFunction[[m]],
                                            x0=params0[i,],
                                            t=as.numeric(as.Date(bondSet$DateMaturity)-dateReference)/daysYear,
                                            yield=bondSet$Yield,
                                            lower = c(0,-Inf,-Inf,-Inf,1,1),
                                            upper = c(Inf, Inf, Inf, Inf, Inf, Inf))
      errors[i] <- RootMeanSquareError(ModelFunction[[m]],
                                       params=params[i,],
                                       t=as.numeric(as.Date(bondSet$DateMaturity)-dateReference)/daysYear,
                                       yield=bondSet$Yield)
      curves[,i] <- ModelFunction[[m]](params=params[i,],
                                       t=as.numeric(as.Date(bondSet$DateMaturity)-dateReference)/daysYear)
    } # end-for
    # get minimum error
    rmse[nrow(rmse)+1,1] <- dateReference
    rmse[nrow(rmse),2] <- paste(m,"_G",sep="")
    rmse[nrow(rmse),3] <- min(errors)
    
    # get fitted yiels for minimum error
    bondSet[,ncol(bondSet)+1] <- curves[,which.min(errors)]
    colnames(bondSet)[ncol(bondSet)] <- paste(m,"_G",sep="")
  } # end for
  
  ############################################################################
  ## KERNEL RIDGE MODEL
  
  # Generate price vector

  cat("Reference date ",format(dateReference,"%Y-%m-%d")," Calculating Model - Kernel Ridge\n")
    
  # Generate price vector
  priceVector <- as.matrix(bondSet$PricePU)
  
  # Generate cashflow matrix
  maxYearsToMaturity <- as.integer(trunc(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear+1))
  cashflowMatrix <- matrix(0, nrow=nrow(bondSet),ncol=maxYearsToMaturity*daysYear)
  
  ## Load bonds cashflow
  for (i in 1:nrow(bondSet)) {
    dateCoupon <- rev(seq(from=bondSet$DateMaturity[i],to=dateReference, by="-6 month"))
    timeToCoupon <- as.numeric(dateCoupon-dateReference)
    for (j in 1:length(dateCoupon)) {
      cashflowMatrix[i,timeToCoupon[j]] <- (trunc(j/length(dateCoupon))+bondSet$Coupon[i]/2)*faceValue  
    } #end-for
  } # end-for
  
  # View reduced cashflow matrix
  print(as.vector(priceVector))
  print(ReduceSparseMatrix(cashflowMatrix))
  
  # converte objetos de R para python
  priceVector <- r_to_py(priceVector)
  cashflowMatrix <- r_to_py(cashflowMatrix)
  
  # Compute fitted KR curve in python
  g_hat <- KRcurve(dateReference,priceVector,cashflowMatrix,maxYearsToMaturity)[[1]]
  bondSet$KR <- FitKR(cashflowMatrix,g_hat,dateReference)[[2]]
  
  rmse[nrow(rmse)+1,1] <- dateReference
  rmse[nrow(rmse),2] <- "KR"
  rmse[nrow(rmse),3] <- sqrt(sum((bondSet$KR-bondSet$Yield)^2))
  
  ###############################################################################
  
  # plot results
  print(rmse)
  write.xlsx(bondSet,file=paste(outputDirectory,"yc_rmse_",format(dateReference,"%Y-%m-%d"),".xlsx",sep=""),
             sheetName="BondSet",append=FALSE)
  write.xlsx(rmse,file=paste(outputDirectory,"yc_rmse_",format(dateReference,"%Y-%m-%d"),".xlsx",sep=""),
             sheetName="RMSE",append=TRUE)
  # Plot model curves curve
  p <- ggplot(data=bondSet,aes(x=DateMaturity))+
    geom_point(data=bondSet,aes(y=Yield,colour="Market")) +  geom_line(data=bondSet,aes(y=Yield,colour="Market"))+
    #geom_point(data=CCScurve,aes(y=YieldCCS,colour="Constrained Cubic Spline")) + 
    geom_line(data=CCScurve,aes(y=CCS,colour="Constrained Cubic Spline"))+
    geom_point(data=bondSet,aes(y=NS,colour="Nelson-Siegel"))+geom_line(aes(y=NS,colour="Nelson-Siegel"))+
    #geom_point(data=bondSet,aes(y=DPR,colour="Diebold-Piazzesi-Rudebush"))+geom_line(aes(y=DPR,colour="Diebold-Piazzesi-Rudebush"))+
    geom_point(data=bondSet,aes(y=BL,colour="Bliss"))+geom_line(aes(y=BL,colour="Bliss"))+
    geom_point(data=bondSet,aes(y=DP,colour="DePooter"))+geom_line(aes(y=DP,colour="DePooter"))+
    geom_point(data=bondSet,aes(y=BC,colour="Björk-Christensen"))+geom_line(aes(y=BC,colour="Björk-Christensen"))+
    geom_point(data=bondSet,aes(y=NSS,colour="Nelson-Siegel-Svensson"))+geom_line(aes(y=NSS,colour="Nelson-Siegel-Svensson"))+
    geom_point(data=bondSet,aes(y=KR,colour="Kernel-Ridge"))+geom_line(aes(y=KR,colour="Kernel-Ridge"))+
    theme(legend.position = c(0.8, 0.3),
          legend.direction = "vertical") +
    scale_colour_manual(
      breaks = c("Market","Constrained Cubic Spline","Nelson-Siegel","Diebold-Piazzesi-Rudebush",
                 "Bliss","DePooter","Björk-Christensen","Nelson-Siegel-Svensson","Kernel-Ridge"),
      values = c("darkgreen","brown","pink","orange","yellow","grey","lightgreen","blue","red"))+
    labs(title = paste("Yield Curves Comparison ",format(dateReference,"%Y-%m-%d")),
         y="Yield",
         x="Maturity [year]")
  # Save plot
  ggsave(filename = paste(outputDirectory,"yc_plot_",format(dateReference,"%Y-%m-%d"),".png",sep=""),
         units = "in",width = 8, height = 6,dpi = 100)
  print(p)
  datePlots[[k]] <- p

  # skip to next date
} # end-for

# plot dates
grid.arrange(datePlots[[1]],datePlots[[2]],datePlots[[3]],datePlots[[4]],
               nrow=2,ncol=2,top=textGrob("Multiple Dates Comparison",gp=gpar(fontsize=15,font=1)),
               bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))

print(rmse)

# print stat results
modelList <- unique(rmse$Model)
modelSummary <- data.frame(Model=character(),Mean=double(),SD=double())
for (m in modelList) {
  modelSet <- subset(rmse,Model==m)
  modelSummary[nrow(modelSummary)+1,] <- c(m,mean(modelSet$RMSE),sd(modelSet$RMSE))
} # end-for
print(modelSummary)
write.xlsx(modelSummary,file=paste(outputDirectory,"yc_summary_",format(workingDates[1],"%Y-%m-%d"),"_to_",
                                   format(workingDates[length(workingDates)],"%Y-%m-%d"),".xlsx",sep=""),
           sheetName=format(dateReference,"%Y-%m-%d"), append = FALSE)


####################################################################################################
