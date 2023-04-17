
# Identification ----------------------------------------------------------
# 
# FGV - Fundacao Getúlio Vargas
# EESP - Escola de Economia São Paulo
# MPEF - Mestrado Profissional em Economia e Finanças
#
# Projeto de dissertação
#
# YIELD CURVES COMPARISON
# Autor: Luis Giovanni Faria
#
# Date: 23/01/2023
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
cat("Working directory =>",getwd())

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
library(dplyr)
#library(bizdays)
#library(lubridate)

source("../lib/yc_lib.r")
source("../lib/yc_loaddata.r")
source_python("../lib/kr_model.py")
source_python("../lib/kr_utils.py")
source_python("../lib/yc_utils.py")
#source_python("../lib/yc_krdate.py")

# set working flags
#cc <- TRUE # use continuous capitalization for YTM
#ac <- FALSE # consider accrued coupons on first payment
#pg <- TRUE # print intermediate plots (not batch mode)
#wd <- FALSE # define if time considers working days (wd=TRUE Brazil 252 days only)
#br <- TRUE # define if working with Tesouro Direto input files
#kr <- FALSE # run KR example set

nsimDE <- 10 # number of simulations differential evolution
#nsimSD <- 10 # number of simulations steepest descent

faceValue <- 1000
daysYear <- 365
results <- data.frame() # full data
resultsF <- data.frame() # filtered data
datePlots <- list() # list of date reference curve plots 
rmsePlots <- list() # list of RMSE plots
curvas <- data.frame(DateReference=as.Date(character())) # summary of curves by reference date
curvasF <- data.frame(DateReference=as.Date(character())) # summary of curves by reference date filtered

###############################################################################
## Define parametric curves
##

# Parametric models to be run with optimization methods DE and SD
# NS = Nelson-Siegel
# BL = Bliss
# NSS= Svensson
# DP = DePooter
# BC = Bjork-Christensen
#models <- c("NS","BL","NSS","DP","BC") # Models to optimize
models <- c("NSS") # Models to optimize

## Define model functions as a list
ModelFunction = list();

ModelFunction[["NS"]] = function(params,t) {
  ## Description: Compute Nelson-Siegel (1987) function
  return(params["Beta1"]+
           params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
           params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"])))
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

ModelFunction[["BC"]] = function(params,t) {
  ## Description: Compute Björk-Christensen (1999) function
  return(params["Beta1"]+
           params["Beta2"]*t/(2*params["Lambda1"])+
           params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
           params["Beta4"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"]))+
           params["Beta5"]*((1-exp(-2*t/params["Lambda1"]))/(2*t/params["Lambda1"])))
} # end-function


##############################################################################
## Load data

# run data from tesouro direto
bondData <- ReadTesouroDireto(inputDirectory,coupon=0.1)

bondData$YieldMid <- (bondData$RateBuy+bondData$RateSell)/2
bondData <- bondData[c(1,2,3,8,9)]
bondData$Years <- as.numeric((bondData$DateMaturity-bondData$DateReference)/daysYear)
bondData <- cbind(bondData[c(1,2,3,4)],bondData[c(6)],bondData[c(5)])
bondData <- subset(bondData,DateReference!=DateMaturity)

#q=as.data.frame(table(bondData$DateReference)) # conta títulos por data

df <- bondData[c(1,5)] %>%
  group_by(DateReference) %>%
  summarize(max_value=max(Years))

# Plot model curves
bondYears <- bondData[c(1,3)]
bondYears$Years <- as.numeric((bondYears$DateMaturity-bondYears$DateReference)/daysYear)

ggplot(data=bondYears,aes(x=DateReference))+
  geom_point(data=bondYears,aes(y=Years,colour="Vencimento"),size=0.3) +
  geom_line(data=df,aes(y=max_value,colour="Vencimento Máximo"),size=0.5,linetype="dashed")+
  theme(legend.position = "bottom", #c(0.8, 0.3),
        legend.direction = "horizontal") +
  scale_colour_manual(name=NULL,
    breaks = c("Vencimento","Vencimento Máximo"),
    values = c("blue","red"))+
  scale_linetype_manual(values=c("Vencimento" = "dashed", "Vencimento Máximo" = "solid"))+
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold")) +
  labs(title = "Faixas de Vencimento",
       y="Vencimento Máximo [anos]",
       x="Data de Referência [Ano Calendário]")

# Save plot
ggsave(filename = paste(outputDirectory,"yc_maturity",".png",sep=""),
       units = "in",
       width = 8, height = 6,
       dpi = 100)


# temp working data frame of unique reference dates
dates <- as.data.frame(unique(bondData$DateReference))
colnames(dates) <- c("date")
dates$yearmonth <- substr(as.character(dates$date),1,7)

# vector of end-of-month reference dates
uniqueYearMonth <- unique(dates$yearmonth)

workingDates <- c()

for (i in 1:length(uniqueYearMonth)) {
  subdates <- subset(dates,yearmonth==uniqueYearMonth[i])
  #n <- nrow(subdates)
  n <- trunc(nrow(subdates)/2)
  #workingDates <- append(workingDates,sort(subdates$date,partial=n-1)[n-1])
  workingDates <- append(workingDates,sort(subdates$date,partial=n)[n])
}

#workingDates=c(as.Date("2021-10-15"),as.Date("2021-11-16"),as.Date("2021-12-14"),as.Date("2022-01-14"));k=1 # TESTE

nDateInitial <- 1
nDateFinal <- length(workingDates)

# main loop to compute all reference dates
for (k in nDateInitial:nDateFinal) {
  dateReference <- workingDates[k]
  
  cat("\n","=> Working with reference date",format(dateReference,"%Y-%m-%d"),sprintf(" Date %d/%d\n",k-nDateInitial+1,nDateFinal-nDateInitial+1),"\n")
  
  # Set to generate price vector and cashflow matrix (subset of BondPrices) 
  bondSet <- subset(bondData,DateReference==dateReference)

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
  df <- data.frame(matrix(ncol=4, nrow=nrow(bondSet)))
  names(df) <- c("CCSinPrice","CCSinYield","CCSoutPrice","CCSoutYield")
  bondSet <- cbind(bondSet,df)

  # sort by maturity
  bondSort <- bondSet[order(bondSet$DateMaturity,decreasing=FALSE),]  
  
  x <- as.numeric(bondSort$DateMaturity-bondSort$DateMaturity[1])/nrow(bondSet)
  y <- bondSort$Yield
  
  # Compute fitted CCS yields for maturity dates
  for (i in 1:nrow(bondSet)) {
    # in-sample
    xstar <- (x[length(x)]-x[1])*as.numeric(bondSet$DateMaturity[i]-bondSet$DateMaturity[1])/
      as.numeric(bondSet$DateMaturity[nrow(bondSet)]-bondSet$DateMaturity[1])
    bondSet$CCSinYield[i] <- CCSpline(x,y,xstar)
    bondSet$CCSinPrice[i] <- PricePU(bondSet$CCSinYield[i],
                                     as.numeric(rev(seq(from=bondSet$DateMaturity[i],to=dateReference, by="-6 month"))-dateReference),
                                     faceValue,
                                     bondSet$Coupon[i]/2)
    # out-of-sample
    xstar <- (x[length(x)]-x[1])*as.numeric(bondSet$DateMaturity[i]-bondSet$DateMaturity[1])/
      as.numeric(bondSet$DateMaturity[nrow(bondSet)]-bondSet$DateMaturity[1]+1)
    bondSet$CCSoutYield[i] <- CCSpline(x,y,xstar)
    bondSet$CCSoutPrice[i] <- PricePU(bondSet$CCSoutYield[i],
                                      as.numeric(rev(seq(from=bondSet$DateMaturity[i]+1,to=dateReference, by="-6 month"))-dateReference),
                                      faceValue,
                                      bondSet$Coupon[i]/2)
  }
  
  ###############################################################################
  ## Run Differential Evolution for parametric curves
  ##
  
  NP <- 100 # population size
  nmaxG = 500 # number of generations
  set.seed(0)
  lower = c(Beta1=0,Beta2=-0.05,Beta3=-0.5,Beta4=-0.5,Beta5=-0.5,Lambda1=0,Lambda2=0);
  upper = c(Beta1=1,Beta2=0.05,Beta3=0.5,Beta4=0.5,Beta5=-0.5,Lambda1=1.25,Lambda2=1.25);
  
  ### Restrictions
  gx = list();
  gx[[1]] = function(params){return(-params["Beta1"]-params["Beta2"])}
  
  nsim <- nsimDE # number of simulations
  
  for (m in models) {
    cat("\n","Reference date ",format(dateReference,"%Y-%m-%d")," Calculating Model - Differential Evolution = ",m,"\n")
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

    poscol <- ncol(bondSet)
    df <- data.frame(matrix(ncol=4, nrow=nrow(bondSet)))
    names(df) <- c(paste(m,"inPrice",sep=""),paste(m,"inYield",sep=""),paste(m,"outPrice",sep=""),paste(m,"outYield",sep=""))
    bondSet <- cbind(bondSet,df)
    
    for (i in 1:nrow(bondSet)) {
      # get in-sample fit
      bondSet[i,poscol+2] <- ModelFunction[[m]](params=paramsOpt[,which.min(errors)],
                                                t=as.numeric(as.Date(bondSet$DateMaturity[i])-dateReference)/daysYear)
      bondSet[i,poscol+1] <- PricePU(bondSet[i,poscol+2],
                                     as.numeric(rev(seq(from=bondSet$DateMaturity[i],to=dateReference, by="-6 month"))-dateReference),
                                     faceValue,
                                     bondSet$Coupon[i]/2)

      # get out-of-sample fit
      bondSet[i,poscol+4] <- ModelFunction[[m]](params=paramsOpt[,which.min(errors)],
                                                      t=as.numeric(as.Date(bondSet$DateMaturity[i])-dateReference+1)/daysYear)
      bondSet[i,poscol+3] <- PricePU(bondSet[i,poscol+4],
                                     as.numeric(rev(seq(from=bondSet$DateMaturity[i]+1,to=dateReference, by="-6 month"))-dateReference),
                                     faceValue,
                                     bondSet$Coupon[i]/2)
    } # end-for
  } # end-for
  
  ############################################################################
  ## KERNEL-RIDGE MODEL

  cat("\n","Reference date ",format(dateReference,"%Y-%m-%d")," Calculating Model - Kernel Ridge\n")
  
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

  # converts objects from R to python
  priceVector <- r_to_py(priceVector)
  cashflowMatrix <- r_to_py(cashflowMatrix)
  
  # Compute in-sample fitted KR curve in python
  g_hat <- KRcurve(dateReference,priceVector,cashflowMatrix,maxYearsToMaturity)[[1]]
  fit <- FitKR(cashflowMatrix,g_hat,dateReference)
  bondSet$KRinPrice <- fit[[1]]
  bondSet$KRinYield <- fit[[2]]

  # out-of-sample (next day t+1 for model estimated on day t)
  # Fit vector of dates
  cashflowMatrix <- matrix(0,nrow=nrow(bondSet),ncol=maxYearsToMaturity*daysYear)
  ## Load bonds cashflows
  for (i in 1:nrow(bondSet)) {
    dateCoupon <- rev(seq(from=bondSet$DateMaturity[i]+1,to=dateReference, by="-6 month"))
    timeToCoupon <- as.numeric(dateCoupon-dateReference)
    for (j in 1:length(dateCoupon)) {
      cashflowMatrix[i,timeToCoupon[j]] <- (trunc(j/length(dateCoupon))+bondSet$Coupon[i]/2)*faceValue  
    } #end-for
  } # end-for

  fit <- FitKR(cashflowMatrix,g_hat,dateReference)
  bondSet$KRoutPrice <- fit[[1]]
  bondSet$KRoutYield <- fit[[2]]
  bondSet$Dur <- fit[[3]]
  bondSet$TTM <- fit[[4]]
  bondSet$InvW <- fit[[5]]
 
  # out-of-sample errors
  bondSet$NSSoutErrPrice <- (bondSet$NSSoutPrice-bondSet$PricePU)
  bondSet$KRoutErrPrice <- (bondSet$KRoutPrice-bondSet$Price)
  
  bondSet$NSSoutErrYield <- bondSet$NSSoutYield-bondSet$Yield
  bondSet$KRoutErrYield <- bondSet$KRoutYield-bondSet$Yield
  
  # in-sample errors
  bondSet$NSSinErrPrice <- (bondSet$NSSinPrice-bondSet$PricePU)
  bondSet$KRinErrPrice <- (bondSet$KRinPrice-bondSet$Price)
  
  bondSet$NSSinErrYield <- bondSet$NSSinYield-bondSet$Yield
  bondSet$KRinErrYield <- bondSet$KRinYield-bondSet$Yield
  
  
  ###############################################################################
  # Summarize by reference date
  
  posrow <- nrow(curvas)+1
  curvas[posrow,1] <- dateReference
  
  # out-of-sample errors
  curvas$NSSoutYtm[posrow] <- sqrt(sum(bondSet$NSSoutErrYield^2))*1e4
  curvas$NSSoutRelPrice[posrow] <- sqrt(sum((bondSet$NSSoutErrPrice/bondSet$PricePU)^2))*1e4
  curvas$NSSoutDurW[posrow] <- sqrt(sum(bondSet$NSSoutErrPrice^2/bondSet$InvW))*1e4
  curvas$NSSoutMatW[posrow] <- sqrt(sum(bondSet$NSSoutErrYield^2*bondSet$TTM/max(bondSet$TTM)))*1e4
  
  curvas$KRoutYtm[posrow] <- sqrt(sum(bondSet$KRoutErrYield^2))*1e4
  curvas$KRoutRelPrice[posrow] <- sqrt(sum((bondSet$KRoutErrPrice/bondSet$PricePU)^2))*1e4
  curvas$KRoutDurW[posrow] <- sqrt(sum(bondSet$KRoutErrPrice^2/bondSet$InvW))*1e4
  curvas$KRoutMatW[posrow] <- sqrt(sum(bondSet$KRoutErrYield^2*bondSet$TTM/max(bondSet$TTM)))*1e4
  
  
  # in-sample errors
  curvas$NSSinYtm[posrow] <- sqrt(sum(bondSet$NSSinErrYield^2))*1e4
  curvas$NSSinRelPrice[posrow] <- sqrt(sum((bondSet$NSSinErrPrice/bondSet$PricePU)^2))*1e4
  curvas$NSSinDurW[posrow] <- sqrt(sum(bondSet$NSSinErrPrice^2/bondSet$InvW))*1e4
  curvas$NSSinMatW[posrow] <- sqrt(sum(bondSet$NSSinErrYield^2*bondSet$TTM/max(bondSet$TTM)))*1e4
  
  curvas$KRinYtm[posrow] <- sqrt(sum(bondSet$KRinErrYield^2))*1e4
  curvas$KRinRelPrice[posrow] <- sqrt(sum((bondSet$KRinErrPrice/bondSet$PricePU)^2))*1e4
  curvas$KRinDurW[posrow] <- sqrt(sum(bondSet$KRinErrPrice^2/bondSet$InvW))*1e4
  curvas$KRinMatW[posrow] <- sqrt(sum(bondSet$KRinErrYield^2*bondSet$TTM/max(bondSet$TTM)))*1e4
  
  ###############################################################################
  # Smoothness - Tension and curvature

  h <- 1/daysYear  
  # KR
  # f(x)
  f_x <- bondSet$KRinYield
  
  # f(x+h)
  f_xplus <-bondSet$KRoutYield
  
  # f(x-h)
  cashflowMatrix <- matrix(0,nrow=nrow(bondSet),ncol=maxYearsToMaturity*daysYear)  
  for (i in 1:nrow(bondSet)) {
    dateCoupon <- rev(seq(from=bondSet$DateMaturity[i]-1,to=dateReference, by="-6 month"))
    timeToCoupon <- as.numeric(dateCoupon-dateReference)
    for (j in 1:length(dateCoupon)) {
      cashflowMatrix[i,timeToCoupon[j]] <- (trunc(j/length(dateCoupon))+bondSet$Coupon[i]/2)*faceValue  
    } #end-for
  } # end-for
 
  f_xminus <- FitKR(cashflowMatrix,g_hat,dateReference)[[2]]
  
  # KR first derivative
  bondSet$KRinD1 <- (f_xplus-f_xminus)/(2*h)
  # KR second derivative
  bondSet$KRinD2 <- (f_xplus-2*f_x+f_xminus)/(h^2)
  
  curvas$KRinTension[posrow] <- sum(bondSet$KRinD1^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
  curvas$KRinCurvature[posrow] <- sum(bondSet$KRinD2^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
    
  # NSS
  f_x <- bondSet$NSSinYield
  f_xplus <- bondSet$NSSoutYield
  for (i in 1:nrow(bondSet)) {
    f_xminus[i] <- ModelFunction[[m]](params=paramsOpt[,which.min(errors)],
                                      t=as.numeric(as.Date(bondSet$DateMaturity[i])-dateReference-1)/daysYear)
  }
  
  # NSS first derivative
  bondSet$NSSinD1 <- (f_xplus-f_xminus)/(2*h)
  # NSS second derivative
  bondSet$NSSinD2 <- (f_xplus-2*f_x+f_xminus)/(h^2)
  
  # smoothness
  curvas$NSSinTension[posrow] <- sum(bondSet$NSSinD1^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
  curvas$NSSinCurvature[posrow] <- sum(bondSet$NSSinD2^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
  
  ###############################################################################
  # Append results
  
  results <- rbind(results,bondSet)
  
  # save output errors
  write.xlsx(bondSet,file=paste(outputDirectory,"yc_rmse_",format(dateReference,"%Y-%m-%d"),".xlsx",sep=""),sheetName="BondSet",append=FALSE)
  
  ###############################################################################
  # plot results

  # Plot model curves
  p <- ggplot(data=bondSet,aes(x=DateMaturity))+
    geom_point(data=bondSet,aes(y=Yield,colour="Taxas Observadas"),shape=4,size=3)+ #geom_line(data=bondSet,aes(y=Yield,colour="Taxas Observadas"))+
    geom_point(data=bondSet,aes(y=NSSinYield,colour="NSS"),shape=15,size=2)+geom_line(aes(y=NSSinYield,colour="NSS"))+
    geom_point(data=bondSet,aes(y=KRinYield,colour="KR"),shape=19,size=2)+geom_line(aes(y=KRinYield,colour="KR"))+
    theme(legend.position = c(0.8, 0.3),
          legend.direction = "vertical") +
    scale_colour_manual(name=NULL,
                        breaks = c("Taxas Observadas","NSS","KR"),
                        values = c("darkgreen","blue","red"))+
    labs(title = paste("Data referência",format(dateReference,"%Y-%m-%d")),
         y="Taxa",
         x="Vencimento [ano]")
  
  # Save plot
  ggsave(filename = paste(outputDirectory,"yc_plot_",format(dateReference,"%Y-%m-%d"),".png",sep=""),
         units = "in",
         width = 8, height = 6,
         dpi = 100)
  
  datePlots[[k-nDateInitial+1]] <- p
  print(p)
  
  ###############################################################################
  # Filter dataset
  
  bondFilter <- subset(bondSet,DateMaturity-dateReference>=90)
  bondFilter <- subset(bondFilter,NSSinErrYield<=3*sd(NSSinErrYield))
  bondFilter <- subset(bondFilter,KRinErrYield<=3*sd(KRinErrYield))
  
  # Sumarize by reference date
  posrow <- nrow(curvasF)+1
  curvasF[posrow,1] <- dateReference
  
  # out-of-sample errors
  curvasF$NSSoutYtm[posrow] <- sqrt(sum(bondFilter$NSSoutErrYield^2))*1e4
  curvasF$NSSoutRelPrice[posrow] <- sqrt(sum((bondFilter$NSSoutErrPrice/bondFilter$PricePU)^2))*1e4
  curvasF$NSSoutDurW[posrow] <- sqrt(sum(bondFilter$NSSoutErrPrice^2/bondFilter$InvW))*1e4
  curvasF$NSSoutMatW[posrow] <- sqrt(sum(bondFilter$NSSoutErrYield^2*bondFilter$TTM/max(bondFilter$TTM)))*1e4
  
  curvasF$KRoutYtm[posrow] <- sqrt(sum(bondFilter$KRoutErrYield^2))*1e4
  curvasF$KRoutRelPrice[posrow] <- sqrt(sum((bondFilter$KRoutErrPrice/bondFilter$PricePU)^2))*1e4
  curvasF$KRoutDurW[posrow] <- sqrt(sum(bondFilter$KRoutErrPrice^2/bondFilter$InvW))*1e4
  curvasF$KRoutMatW[posrow] <- sqrt(sum(bondFilter$KRoutErrYield^2*bondFilter$TTM/max(bondFilter$TTM)))*1e4

  # smoothness out-of-sample
  curvasF$KRoutTension[posrow] <- sum(bondFilter$KRoutD1^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  curvasF$KRoutCurvature[posrow] <- sum(bondFilter$KRoutD2^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  
  curvasF$NSSoutTension[posrow] <- sum(bondFilter$NSSoutD1^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  curvasF$NSSoutCurvature[posrow] <- sum(bondFilter$NSSoutD2^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  
    
  # in-sample errors
  curvasF$NSSinYtm[posrow] <- sqrt(sum(bondFilter$NSSinErrYield^2))*1e4
  curvasF$NSSinRelPrice[posrow] <- sqrt(sum((bondFilter$NSSinErrPrice/bondFilter$PricePU)^2))*1e4
  curvasF$NSSinDurW[posrow] <- sqrt(sum(bondFilter$NSSinErrPrice^2/bondFilter$InvW))*1e4
  curvasF$NSSinMatW[posrow] <- sqrt(sum(bondFilter$NSSinErrYield^2*bondFilter$TTM/max(bondFilter$TTM)))*1e4
  
  curvasF$KRinYtm[posrow] <- sqrt(sum(bondFilter$KRinErrYield^2))*1e4
  curvasF$KRinRelPrice[posrow] <- sqrt(sum((bondFilter$KRinErrPrice/bondFilter$PricePU)^2))*1e4
  curvasF$KRinDurW[posrow] <- sqrt(sum(bondFilter$KRinErrPrice^2/bondFilter$InvW))*1e4
  curvasF$KRinMatW[posrow] <- sqrt(sum(bondFilter$KRinErrYield^2*bondFilter$TTM/max(bondFilter$TTM)))*1e4

  # smoothness in-sample
  curvasF$KRinTension[posrow] <- sum(bondFilter$KRinD1^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  curvasF$KRinCurvature[posrow] <- sum(bondFilter$KRinD2^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4

  curvasF$NSSinTension[posrow] <- sum(bondFilter$NSSinD1^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  curvasF$NSSinCurvature[posrow] <- sum(bondFilter$NSSinD2^2)*h/(as.numeric(max(bondFilter$DateMaturity)-dateReference)/daysYear)*1e4
  
  resultsF <- rbind(resultsF,bondFilter)
  
  # skip to next date
} # end-for

###############################################################################
# Calculate results

cols <- c("Model","OutYtm","OutRelPrice","OutDurW","OutMatW","InYtm","InRelPrice","InDurW","InMatW","InTension","InCurvature")
cols1 <- c("outYtm","outRelPrice","outDurW","outMatW","inYtm","inRelPrice","inDurW","inMatW","inTension","inCurvature")
rows <- c("NSS","KR")

# rmse
rmse <- data.frame(matrix(ncol=11,nrow=2))
colnames(rmse) <- cols
rownames(rmse) <- rows
rmse$Model <- rows

rmseF <- rmse

for (r in rows) {
  for (i in 1:length(cols1)) {
    rmse[r,cols[i+1]] <- mean(curvas[c(paste(r,cols1[i],sep=""))][,1])
  }
}

###############################################################################
# rmse FILTERED

for (r in rows) {
  for (i in 1:(length(cols1))) {
    rmseF[r,cols[i+1]] <- mean(curvasF[c(paste(r,cols1[i],sep=""))][,1])
  }
}

###############################################################################
# Plot comparison

f_plot <- function(df,c0,t0,d=1) {
  ggplot(data=df,aes(x=Model,y=get(c0),fill=Model)) + 
    geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
    geom_text(aes(label=round(get(c0),d)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
    scale_fill_manual(values = c("red","blue") ) +
    theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
    labs(title=t0,y=NULL,x=NULL)
}

bot <- c("Todos os Dados (sem filtros)","Todos os Dados (sem filtros)","Dados Filtrados","Dados Filtrados")
top <- c("Erros de precificação para conjunto de treinamento (dentro da amostra)","Erros de precificação para conjunto de teste (fora da amostra)","Erros de precificação para conjunto de treinamento (dentro da amostra)","Erros de precificação para conjunto de teste (fora da amostra)")
tit <- c("Taxa REQM (BPS)","Preço relativo REQM (BPS)","Duração ponderada REQM (BPS)","Vencimento ponderado REQM (BPS)")
lab <- c("Todos os Dados (sem filtros)","Dados Filtrados")
col <- data.frame(matrix(nrow=2,ncol=4))
col[1,] <- c("InYtm","InRelPrice","InDurW","InMatW")
col[2,] <- c("OutYtm","OutRelPrice","OutDurW","OutMatW")
#col[3,] <- c("InYtm","InRelPrice","InDurW","InMatW")
#col[4,] <- c("OutYtm","OutRelPrice","OutDurW","OutMatW")

for (i in 1:4) {
  if (i<=2) dat <- rmse else dat <- rmseF
  j=c(1,2,1,2)[i]
  grid.arrange(f_plot(dat,col[j,1],tit[1]),f_plot(dat,col[j,2],tit[2]),f_plot(dat,col[j,3],tit[3]),f_plot(dat,col[j,4],tit[4]),
               nrow=2,ncol=2,top=textGrob(top[i],gp=gpar(fontsize=15,font=1)),
               bottom=textGrob(bot[i],gp=gpar(fontsize=15,font=1)))
} # end-if


# Presentation format
for (i in 1:2) {
  if (i==1) dat <- rmse else dat <- rmseF
  grid.arrange(f_plot(dat,col[1,1],tit[1]),f_plot(dat,col[1,2],tit[2]),f_plot(dat,col[1,3],tit[3]),f_plot(dat,col[1,4],tit[4]),
               f_plot(dat,col[2,1],tit[1]),f_plot(dat,col[2,2],tit[2]),f_plot(dat,col[2,3],tit[3]),f_plot(dat,col[2,4],tit[4]),
               nrow=2,ncol=4,top=textGrob("Erros de precificação para conjuntos de treinamento e teste",gp=gpar(fontsize=15,font=1)),
               bottom=textGrob(lab[trunc(i/2)+1],gp=gpar(fontsize=15,font=1)))
} # end-if

# Smoothness plot
titS <- c("Tensão (BPS)","Curvatura (BPS)")
colS <- c("InTension","InCurvature")
for (i in 1:2) {
  if (i==1) dat <- rmse else dat <- rmseF
  grid.arrange(f_plot(dat,colS[1],titS[1],3),f_plot(dat,colS[2],titS[2],3),
               nrow=1,ncol=2,top=textGrob("Suavidade",gp=gpar(fontsize=15,font=1)),
               bottom=textGrob(lab[i],gp=gpar(fontsize=15,font=1)))
} # end-forf


###############################################################################
# Separate bonds in Buckets by maturity

buckets <- data.frame(matrix(ncol=9,nrow=10))

# out-of-sample
colnames(buckets) <- c("Range","NSSoutYtm","KRoutYtm","NSSoutRelPrice","KRoutRelPrice","NSSoutDurW","KRoutDurW","NSSoutMatW","KRoutMatW")
buckets$Range <- c("0 a 1a","1a a 2a","2a a 3a","3a a 4a","4a a 5a","5a a 6a","6a a 7a","7a a 8a","8a a 9a","9a<")
rownames(buckets) <- buckets$Range 

for (i in 1:10) {
  resultsB <- subset(results,(Years<=i)&(Years>i-1))

  if (nrow(resultsB)==0) buckets$NSSoutYtm[i] <- 0 else buckets$NSSoutYtm[i] <- sqrt(sum(resultsB$NSSoutErrYield^2))*1e4
  buckets$NSSoutRelPrice[i] <- sqrt(sum((resultsB$NSSoutErrPrice/resultsB$PricePU)^2))*1e4
  buckets$NSSoutDurW[i] <- sqrt(sum(resultsB$NSSoutErrPrice^2/resultsB$InvW))*1e4
  buckets$NSSoutMatW[i] <- sqrt(sum(resultsB$NSSoutErrYield^2*resultsB$TTM/max(resultsB$TTM)))*1e4
  
  if (nrow(resultsB)==0) buckets$KRoutYtm[i] <- 0 else buckets$KRoutYtm[i] <- sqrt(sum(resultsB$KRoutErrYield^2))*1e4
  buckets$KRoutRelPrice[i] <- sqrt(sum((resultsB$KRoutErrPrice/resultsB$PricePU)^2))*1e4
  buckets$KRoutDurW[i] <- sqrt(sum(resultsB$KRoutErrPrice^2/resultsB$InvW))*1e4
  buckets$KRoutMatW[i] <- sqrt(sum(resultsB$KRoutErrYield^2*resultsB$TTM/max(resultsB$TTM)))*1e4
  
} # end-for

# plot results errors for different maturities

f_plot2 <- function(df,c1,c2,t1,t2) {
    ggplot(data=df,aes(x=Range,group=1))+
    geom_point(aes(y=get(c1),colour="NSS"))+geom_line(aes(y=get(c1),colour="NSS"))+
    geom_point(aes(y=get(c2),colour="KR"))+geom_line(aes(y=get(c2),colour="KR"))+
    theme(legend.position=c(0.8, 0.8),legend.direction="vertical") +
    scale_colour_manual(name=NULL,breaks=c("NSS","KR"),values=c("blue","red"))+
    theme(plot.title=element_text(hjust=0.5,size=14,face="bold")) +
    labs(title = t1,x=NULL,y=t2)
}

tit1 <- c("Taxa REQM","Preço relativo REQM","Duração ponderada REQM","Vencimento ponderado REQM")
tit2 <- c("Taxa REQM (BPS)","Preço relativo REQM (BPS)","Duração ponderada REQM (BPS)","Vencimento ponderado REQM (BPS)")

col1 <- c("NSSoutYtm","NSSoutRelPrice","NSSoutDurW","NSSoutMatW")
col2 <- c("KRoutYtm","KRoutRelPrice","KRoutDurW","KRoutMatW")

grid.arrange(f_plot2(buckets,col1[1],col2[1],tit1[1],tit2[1]),f_plot2(buckets,col1[2],col2[2],tit1[2],tit2[2]),
             f_plot2(buckets,col1[3],col2[3],tit1[3],tit2[3]),f_plot2(buckets,col1[4],col2[4],tit1[4],tit2[4]),
             nrow=2,ncol=2,top=textGrob("Erros de precificação para conjunto de teste em diferentes vencimentos",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Todos os Dados (sem filtros)",gp=gpar(fontsize=15,font=1)))


# plot dates
if (length(datePlots)==4) { 
  grid.arrange(datePlots[[1]],datePlots[[2]],datePlots[[3]],datePlots[[4]],
               nrow=2,ncol=2,top=textGrob("Comparação entre múltiplas datas",gp=gpar(fontsize=15,font=1)),
               bottom=textGrob("Todos os Dados (sem filtros)",gp=gpar(fontsize=15,font=1)))
} # end-if

###############################################################################
# Write data

write.xlsx(rmse,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="RMSE",append=FALSE)
write.xlsx(rmseF,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="RMSE_F",append=TRUE)
write.xlsx(curvas,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Curves",append=TRUE)
write.xlsx(curvasF,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Curves_F",append=TRUE)
write.xlsx(buckets,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Buckets",append=TRUE)
write.xlsx(results,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Results",append=TRUE)
write.xlsx(resultsF,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Results_F",append=TRUE)

####################################################################################################
