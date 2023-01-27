
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

nsimDE <- 5 # number of simulations differential evolution
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


##############################################################################
## Load data

# run data from tesouro direto
bondData <- ReadTesouroDireto(inputDirectory,coupon=0.1)

bondData$YieldMid <- (bondData$RateBuy+bondData$RateSell)/2
bondData <- bondData[c(1,2,3,8,9)]
bondData$Years <- as.numeric((bondData$DateMaturity-bondData$DateReference)/daysYear)
bondData <- cbind(bondData[c(1,2,3,4)],bondData[c(6)],bondData[c(5)])
bondData <- subset(bondData,DateReference!=DateMaturity)

q=as.data.frame(table(bondData$DateReference)) # conta títulos por data

df <- bondData[c(1,5)] %>%
  group_by(DateReference) %>%
  summarize(max_value=max(Years))

# Plot model curves
bondYears <- bondData[c(1,3)]
bondYears$Years <- as.numeric((bondYears$DateMaturity-bondYears$DateReference)/daysYear)

ggplot(data=bondYears,aes(x=DateReference))+
  geom_point(data=bondYears,aes(y=Years,colour="Maturity"),size=0.3) +
  geom_line(data=df,aes(y=max_value,colour="Max Maturity"),size=0.5)+
  theme(legend.position = c(0.8, 0.3),
        legend.direction = "vertical") +
  scale_colour_manual(name=NULL,
    breaks = c("Maturity","Max Maturity"),
    values = c("blue","red"))+
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold")) +
  labs(title = "Maturity Ranges",
       y="Max Maturity [years]",
       x="Reference Date [Calendar Year]")

# temp working data frame of unique reference dates
dates <- as.data.frame(unique(bondData$DateReference))
colnames(dates) <- c("date")
dates$yearmonth <- substr(as.character(dates$date),1,7)

# vector of end-of-month reference dates
uniqueYearMonth <- unique(dates$yearmonth)

workingDates <- c()

for (i in 1:length(uniqueYearMonth)) {
  subdates <- subset(dates,yearmonth==uniqueYearMonth[i])
  n <- nrow(subdates)
  workingDates <- append(workingDates,sort(subdates$date,partial=n-1)[n-1])
}

rm(list=c("uniqueYearMonth","dates","subdates","i"))

workingDates=c(as.Date("2021-09-30"),as.Date("2021-10-29"),as.Date("2021-11-30"),as.Date("2021-12-30"));k=1

nDateInitial <- 1
nDateFinal <- length(workingDates)

# main loop to compute all reference dates
for (k in nDateInitial:nDateFinal) {
  dateReference <- workingDates[k]
  
  cat("\n","=> Working with reference date",format(dateReference,"%Y-%m-%d"),sprintf(" Date %d/%d\n",k-nDateInitial+1,nDateFinal-nDateInitial+1),"\n")
  
  # Set to generate price vector and cashfow matrix (subset of BondPrinces) 
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
      couponValue <- bondSet$Coupon[i]/2
      if (j==length(dateCoupon)) {
        cashflowMatrix[i,timeToCoupon[j]] <- (1+couponValue)*faceValue  
      } else {
        cashflowMatrix[i,timeToCoupon[j]] <-  couponValue*faceValue  
      }
    } #end-for
  } # end-for
  # remove temporary variables
  rm(list=c("i","j","dateCoupon"))

  # View reduced cashflow matrix
  print(as.vector(priceVector))
  print(ReduceSparseMatrix(cashflowMatrix))

  # converts objects from R to python
  priceVector <- r_to_py(priceVector)
  cashflowMatrix <- r_to_py(cashflowMatrix)
  
  # Compute in-sample fitted KR curve in python
  g_hat <- KRcurve(dateReference,priceVector,cashflowMatrix,maxYearsToMaturity)[[1]]
  fit <- FitKR(cashflowMatrix,g_hat,dateReference)
  bondSet$KRinPrice <- fit[[2]]
  bondSet$KRinYield <- fit[[1]]

  # out-of-sample (next day t+1 for model estimated on day t)
  # Fit vector of dates
  cashflowMatrix <- matrix(0,nrow=nrow(bondSet),ncol=maxYearsToMaturity*daysYear)
  ## Load bonds cashflows
  for (i in 1:nrow(bondSet)) {
    dateCoupon <- rev(seq(from=bondSet$DateMaturity[i]+1,to=dateReference, by="-6 month"))
    timeToCoupon <- as.numeric(dateCoupon-dateReference)
    for (j in 1:length(dateCoupon)) {
      couponValue <- bondSet$Coupon[i]/2
      if (j==length(dateCoupon)) {
        cashflowMatrix[i,timeToCoupon[j]] <- (1+couponValue)*faceValue  
      } else {
        cashflowMatrix[i,timeToCoupon[j]] <-  couponValue*faceValue  
      }
    } #end-for
  } # end-for

  fit <- FitKR(cashflowMatrix,g_hat,dateReference)
  bondSet$KRoutPrice <- fit[[2]]
  bondSet$KRoutYield <- fit[[1]]
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
  # Tension and curvature
  
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
      couponValue <- bondSet$Coupon[i]/2
      if (j==length(dateCoupon)) {
        cashflowMatrix[i,timeToCoupon[j]] <- (1+couponValue)*faceValue  
      } else {
        cashflowMatrix[i,timeToCoupon[j]] <-  couponValue*faceValue  
      }
    } #end-for
  } # end-for
 
  f_xminus <- FitKR(cashflowMatrix,g_hat,dateReference)[[1]]
  
  h <- 1/daysYear
  # KR first derivate
  bondSet$KRd1 <- (f_xplus-f_xminus)/(2*h)
  # KR second derivate
  bondSet$KRd2 <- (f_xplus-2*f_x+f_xminus)/(h^2)
  
  curvas$KRtension[posrow] <- sum(bondSet$KRd1^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
  curvas$KRcurvature[posrow] <- sum(bondSet$KRd2^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
    
  # NSS
  f_x <- bondSet$NSSinYield
  f_xplus <- bondSet$NSSoutYield
  for (i in 1:nrow(bondSet)) {
    f_xminus[i] <- ModelFunction[[m]](params=paramsOpt[,which.min(errors)],
                                      t=as.numeric(as.Date(bondSet$DateMaturity[i])-dateReference-1)/daysYear)
  }
  
  # NSS first derivate
  bondSet$NSSd1 <- (f_xplus-f_xminus)/(2*h)
  # NSS second derivate
  bondSet$NSSd2 <- (f_xplus-2*f_x+f_xminus)/(h^2)
  
  curvas$NSStension[posrow] <- sum(bondSet$NSSd1^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
  curvas$NSScurvature[posrow] <- sum(bondSet$NSSd2^2)*h/(as.numeric(max(bondSet$DateMaturity)-dateReference)/daysYear)*1e4
  
  ###############################################################################
  # Append results
  
  results <- rbind(results,bondSet)
  
  ###############################################################################
  # plot results
  # write.xlsx(bondSet,file=paste(outputDirectory,"yc_rmse_",format(dateReference,"%Y-%m-%d"),".xlsx",sep=""),sheetName="BondSet",append=FALSE)
  # Plot model curves curve
  p <- ggplot(data=bondSet,aes(x=DateMaturity))+
    geom_point(data=bondSet,aes(y=Yield,colour="Market"))+geom_line(data=bondSet,aes(y=Yield,colour="Market"))+
    geom_point(data=bondSet,aes(y=NSSinYield,colour="NSS"))+geom_line(aes(y=NSSinYield,colour="NSS"))+
    geom_point(data=bondSet,aes(y=KRinYield,colour="KR"))+geom_line(aes(y=KRinYield,colour="KR"))+
    theme(legend.position = c(0.8, 0.3),
          legend.direction = "vertical") +
    scale_colour_manual(name=NULL,
                        breaks = c("Market","NSS","KR"),
                        values = c("darkgreen","blue","red"))+
    labs(title = paste("Reference date",format(dateReference,"%Y-%m-%d")),
         y="Yield",
         x="Maturity [year]")
  
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
  
  
  # in-sample errors
  curvasF$NSSinYtm[posrow] <- sqrt(sum(bondFilter$NSSinErrYield^2))*1e4
  curvasF$NSSinRelPrice[posrow] <- sqrt(sum((bondFilter$NSSinErrPrice/bondFilter$PricePU)^2))*1e4
  curvasF$NSSinDurW[posrow] <- sqrt(sum(bondFilter$NSSinErrPrice^2/bondFilter$InvW))*1e4
  curvasF$NSSinMatW[posrow] <- sqrt(sum(bondFilter$NSSinErrYield^2*bondFilter$TTM/max(bondFilter$TTM)))*1e4
  
  curvasF$KRinYtm[posrow] <- sqrt(sum(bondFilter$KRinErrYield^2))*1e4
  curvasF$KRinRelPrice[posrow] <- sqrt(sum((bondFilter$KRinErrPrice/bondFilter$PricePU)^2))*1e4
  curvasF$KRinDurW[posrow] <- sqrt(sum(bondFilter$KRinErrPrice^2/bondFilter$InvW))*1e4
  curvasF$KRinMatW[posrow] <- sqrt(sum(bondFilter$KRinErrYield^2*bondFilter$TTM/max(bondFilter$TTM)))*1e4
  
  resultsF <- rbind(resultsF,bondFilter)
  
  # skip to next date
} # end-for

# Calculate results
###############################################################################
# rmse

rmse <- data.frame(matrix(ncol=11,nrow=2))

# out-of-sample
colnames(rmse) <- c("Model","OutYtm","OutRelPrice","OutDurW","OutMatW",
                    "InYtm","InRelPrice","InDurW","InMatW",
                    "Tension","Curvature")
rownames(rmse) <- c("NSS","KR")
rmse$Model <- c("NSS","KR")

rmse["NSS","OutYtm"] <- mean(curvas$NSSoutYtm)
rmse["KR","OutYtm"] <- mean(curvas$KRoutYtm)

rmse["NSS","OutRelPrice"] <- mean(curvas$NSSoutRelPrice)
rmse["KR","OutRelPrice"] <- mean(curvas$KRoutRelPrice)

rmse["NSS","OutDurW"] <- mean(curvas$NSSoutDurW)
rmse["KR","OutDurW"] <- mean(curvas$KRoutDurW)

rmse["NSS","OutMatW"] <- mean(curvas$NSSoutMatW)
rmse["KR","OutMatW"] <- mean(curvas$KRoutMatW)

rmse["NSS","InYtm"] <- mean(curvas$NSSinYtm)
rmse["KR","InYtm"] <- mean(curvas$KRinYtm)

rmse["NSS","InRelPrice"] <- mean(curvas$NSSinRelPrice)
rmse["KR","InRelPrice"] <- mean(curvas$KRinRelPrice)

rmse["NSS","InDurW"] <- mean(curvas$NSSinDurW)
rmse["KR","InDurW"] <- mean(curvas$KRinDurW)

rmse["NSS","InMatW"] <- mean(curvas$NSSinMatW)
rmse["KR","InMatW"] <- mean(curvas$KRinMatW)

rmse["NSS","Tension"] <- mean(curvas$NSStension)
rmse["KR","Tension"] <- mean(curvas$KRtension)

rmse["NSS","Curvature"] <- mean(curvas$NSScurvature)
rmse["KR","Curvature"] <- mean(curvas$KRcurvature)


rmsePlots[[1]] <- ggplot(data=rmse,aes(x=Model,y=InYtm,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InYtm,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Year-To-Maturity RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[2]] <- ggplot(data=rmse,aes(x=Model,y=InRelPrice,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InRelPrice,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Relative Pricing RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[3]] <- ggplot(data=rmse,aes(x=Model,y=InDurW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InDurW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Duration Weighted RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[4]] <- ggplot(data=rmse,aes(x=Model,y=InMatW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InMatW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Maturity Weighted RMSE (BPS)",y=NULL,x=NULL)

grid.arrange(rmsePlots[[1]],rmsePlots[[2]],rmsePlots[[3]],rmsePlots[[4]],nrow=2,ncol=2,
             top=textGrob("Aggregated in-sample pricing errors",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))

rmsePlots[[5]] <- ggplot(data=rmse,aes(x=Model,y=OutYtm,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutYtm,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Year-To-Maturity RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[6]] <- ggplot(data=rmse,aes(x=Model,y=OutRelPrice,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutRelPrice,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Relative Pricing RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[7]] <- ggplot(data=rmse,aes(x=Model,y=OutDurW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutDurW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Duration Weighted RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[8]] <- ggplot(data=rmse,aes(x=Model,y=OutMatW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutMatW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Maturity Weighted RMSE (BPS)",y=NULL,x=NULL)

grid.arrange(rmsePlots[[5]],rmsePlots[[6]],rmsePlots[[7]],rmsePlots[[8]],nrow=2,ncol=2,
             top=textGrob("Aggregated out-of-sample pricing errors",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))

###############################################################################
# rmse FILTERED

rmseF <- data.frame(matrix(ncol=11,nrow=2))

# out-of-sample
colnames(rmseF) <- c("Model","OutYtm","OutRelPrice","OutDurW","OutMatW",
                    "InYtm","InRelPrice","InDurW","InMatW",
                    "Tension","Curvature")
rownames(rmseF) <- c("NSS","KR")
rmseF$Model <- c("NSS","KR")

rmseF["NSS","OutYtm"] <- mean(curvasF$NSSoutYtm)
rmseF["KR","OutYtm"] <- mean(curvasF$KRoutYtm)

rmseF["NSS","OutRelPrice"] <- mean(curvasF$NSSoutRelPrice)
rmseF["KR","OutRelPrice"] <- mean(curvasF$KRoutRelPrice)

rmseF["NSS","OutDurW"] <- mean(curvasF$NSSoutDurW)
rmseF["KR","OutDurW"] <- mean(curvasF$KRoutDurW)

rmseF["NSS","OutMatW"] <- mean(curvasF$NSSoutMatW)
rmseF["KR","OutMatW"] <- mean(curvasF$KRoutMatW)

rmseF["NSS","InYtm"] <- mean(curvasF$NSSinYtm)
rmseF["KR","InYtm"] <- mean(curvasF$KRinYtm)

rmseF["NSS","InRelPrice"] <- mean(curvasF$NSSinRelPrice)
rmseF["KR","InRelPrice"] <- mean(curvasF$KRinRelPrice)

rmseF["NSS","InDurW"] <- mean(curvasF$NSSinDurW)
rmseF["KR","InDurW"] <- mean(curvasF$KRinDurW)

rmseF["NSS","InMatW"] <- mean(curvasF$NSSinMatW)
rmseF["KR","InMatW"] <- mean(curvasF$KRinMatW)

#rmseF["NSS","Tension"] <- mean(curvasF$NSStension)
#rmseF["KR","Tension"] <- mean(curvasF$KRtension)

#rmseF["NSS","Curvature"] <- mean(curvasF$NSScurvature)
#rmseF["KR","Curvature"] <- mean(curvasF$KRcurvature)


rmsePlots[[9]] <- ggplot(data=rmseF,aes(x=Model,y=InYtm,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InYtm,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Year-To-Maturity RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[10]] <- ggplot(data=rmseF,aes(x=Model,y=InRelPrice,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InRelPrice,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Relative Pricing RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[11]] <- ggplot(data=rmseF,aes(x=Model,y=InDurW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InDurW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Duration Weighted RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[12]] <- ggplot(data=rmseF,aes(x=Model,y=InMatW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(InMatW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Maturity Weighted RMSE (BPS)",y=NULL,x=NULL)

grid.arrange(rmsePlots[[9]],rmsePlots[[10]],rmsePlots[[11]],rmsePlots[[12]],nrow=2,ncol=2,
             top=textGrob("Aggregated in-sample pricing errors",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Filtered Data",gp=gpar(fontsize=15,font=1)))


rmsePlots[[13]] <- ggplot(data=rmseF,aes(x=Model,y=OutYtm,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutYtm,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Year-To-Maturity RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[14]] <- ggplot(data=rmseF,aes(x=Model,y=OutRelPrice,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutRelPrice,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Relative Pricing RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[15]] <- ggplot(data=rmseF,aes(x=Model,y=OutDurW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutDurW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Duration Weighted RMSE (BPS)",y=NULL,x=NULL)

rmsePlots[[16]] <- ggplot(data=rmseF,aes(x=Model,y=OutMatW,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(OutMatW,1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Maturity Weighted RMSE (BPS)",y=NULL,x=NULL)

grid.arrange(rmsePlots[[13]],rmsePlots[[14]],rmsePlots[[15]],rmsePlots[[16]],nrow=2,ncol=2,
             top=textGrob("Aggregated out-of-sample pricing errors",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Filtered Data",gp=gpar(fontsize=15,font=1)))


rmsePlots[[17]] <- ggplot(data=rmse,aes(x=Model,y=Tension,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(Tension,digits=3)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Tension (BPS)",y=NULL,x=NULL)


rmsePlots[[18]] <- ggplot(data=rmse,aes(x=Model,y=Curvature,fill=Model)) + 
  geom_bar(stat="identity",width=0.5,colour="black",alpha=0.7) +
  geom_text(aes(label=round(Curvature,digits=1)),colour="black",vjust=1.1,size=3.5,fontface="bold") +
  scale_fill_manual(values = c("red","blue") ) +
  theme(plot.title=element_text(hjust=0.5,size=8,face="bold"),legend.position="none") +
  labs(title = "Curvature (BPS)",y=NULL,x=NULL)

grid.arrange(rmsePlots[[17]],rmsePlots[[18]],nrow=1,ncol=2,
             top=textGrob("Aggregated smoothness",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))



###############################################################################
# Separate bonds in Buckets by maturity

buckets <- data.frame(matrix(ncol=9,nrow=10))

# out-of-sample
colnames(buckets) <- c("Range","NSSoutYtm","KRoutYtm","NSSoutRelPrice","KRoutRelPrice","NSSoutDurW","KRoutDurW","NSSoutMatW","KRoutMatW")
buckets$Range <- c("0 to 1Y","1Y to 2Y","2Y to 3Y","3Y to 4Y","4Y to 5Y","5Y to 6Y","6Y to 7Y","7Y to 8Y","8Y to 9Y","9Y<")
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

# plot results

rmsePlots[[19]] <- ggplot(data=buckets,aes(x=Range,group=1))+
  geom_point(aes(y=NSSoutYtm,colour="NSS"))+geom_line(aes(y=NSSoutYtm,colour="NSS"))+
  geom_point(aes(y=KRoutYtm,colour="KR"))+geom_line(aes(y=KRoutYtm,colour="KR"))+
  theme(legend.position=c(0.8, 0.8),legend.direction="vertical") +
  scale_colour_manual(name=NULL,breaks=c("NSS","KR"),values=c("blue","red"))+
  theme(plot.title=element_text(hjust=0.5,size=14,face="bold")) +
  labs(title = "Yield RMSE",x=NULL,y="YTM RMSE (MPS)")

rmsePlots[[20]] <- ggplot(data=buckets,aes(x=Range,group=1))+
  geom_point(aes(y=NSSoutRelPrice,colour="NSS"))+geom_line(aes(y=NSSoutRelPrice,colour="NSS"))+
  geom_point(aes(y=KRoutRelPrice,colour="KR"))+geom_line(aes(y=KRoutRelPrice,colour="KR"))+
  theme(legend.position=c(0.8, 0.8),legend.direction="vertical") +
  scale_colour_manual(name=NULL,breaks=c("NSS","KR"),values=c("blue","red"))+
  theme(plot.title=element_text(hjust=0.5,size=14,face="bold")) +
  labs(title = "Relative Pricing RMSE",x=NULL,y="Relative Pricing RMSE (MPS)")

rmsePlots[[21]] <- ggplot(data=buckets,aes(x=Range,group=1))+
  geom_point(aes(y=NSSoutDurW,colour="NSS"))+geom_line(aes(y=NSSoutDurW,colour="NSS"))+
  geom_point(aes(y=KRoutDurW,colour="KR"))+geom_line(aes(y=KRoutDurW,colour="KR"))+
  theme(legend.position=c(0.8, 0.8),legend.direction="vertical") +
  scale_colour_manual(name=NULL,breaks=c("NSS","KR"),values=c("blue","red"))+
  theme(plot.title=element_text(hjust=0.5,size=14,face="bold")) +
  labs(title = "Duration Weighted RMSE",x=NULL,y="Duration Weighted RMSE (MPS)")

rmsePlots[[22]] <- ggplot(data=buckets,aes(x=Range,group=1))+
  geom_point(aes(y=NSSoutMatW,colour="NSS"))+geom_line(aes(y=NSSoutMatW,colour="NSS"))+
  geom_point(aes(y=KRoutMatW,colour="KR"))+geom_line(aes(y=KRoutMatW,colour="KR"))+
  theme(legend.position=c(0.8, 0.8),legend.direction="vertical") +
  scale_colour_manual(name=NULL,breaks=c("NSS","KR"),values=c("blue","red"))+
  theme(plot.title=element_text(hjust=0.5,size=14,face="bold")) +
  labs(title = "Maturity Weighted RMSE",x=NULL,y="Maturity Weighted RMSE (MPS)")

grid.arrange(rmsePlots[[19]],rmsePlots[[20]],rmsePlots[[21]],rmsePlots[[22]],nrow=2,ncol=2,
             top=textGrob("Out-of-sample pricing errors for different maturities",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))


# Presentation format
grid.arrange(rmsePlots[[1]],rmsePlots[[2]],rmsePlots[[3]],rmsePlots[[4]],
             rmsePlots[[5]],rmsePlots[[6]],rmsePlots[[7]],rmsePlots[[8]],nrow=2,ncol=4,
             top=textGrob("Aggregated in-sample and out-of-sample pricing errors",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))


grid.arrange(rmsePlots[[9]],rmsePlots[[10]],rmsePlots[[11]],rmsePlots[[12]],
             rmsePlots[[13]],rmsePlots[[14]],rmsePlots[[15]],rmsePlots[[16]],nrow=2,ncol=4,
             top=textGrob("Aggregated in-sample and out-of-sample pricing errors",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Filtered Data",gp=gpar(fontsize=15,font=1)))




# plot dates
grid.arrange(datePlots[[1]],datePlots[[2]],datePlots[[3]],datePlots[[4]],
             nrow=2,ncol=2,
             top=textGrob("Multiple Dates Comparison",gp=gpar(fontsize=15,font=1)),
             bottom=textGrob("Full Data (no outlier removal)",gp=gpar(fontsize=15,font=1)))



 
###############################################################################
# Write data

write.xlsx(rmse,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="RMSE",append=FALSE)
write.xlsx(curvas,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Curves",append=TRUE)
write.xlsx(results,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Buckets",append=TRUE)
write.xlsx(results,file=paste(outputDirectory,"yc_fulldata.xlsx",sep=""),sheetName="Results",append=TRUE)







####################################################################################################
