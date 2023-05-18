
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
#library(reticulate)
#library(Matrix)
library(ggplot2)
#library(grid)
#library(gridExtra)
library(dplyr)
library(bizdays)
#library(lubridate)

source("../lib/yc_lib.r")
source("../lib/yc_loaddata.r")


nsimDE <- 10 # number of simulations differential evolution

faceValue <- 1000
daysYear <- 365


## Define vertices
vertices <- data.frame(DateReference=as.Date(character()),Vertice=as.numeric(),DateVertice=as.Date(character()),Yield=as.numeric()) # summary of curves by reference date
paramsNSS <- data.frame(Beta1=as.numeric(),Beta2=as.numeric(),Beta3=as.numeric(),Beta4=as.numeric(),Beta5=as.numeric(),Lambda1=as.numeric(),Lambda2=as.numeric())
calendar = calendars()[["Brazil/ANBIMA"]]

###############################################################################
## Define parametric curves
##

models <- c("NSS") # Models to optimize

## Define model functions as a list
ModelFunction = list();


ModelFunction[["NSS"]] = function(params,t) {
  ## Description: Compute Nelson-Siegel-Svensson (1994) function
  return(params["Beta1"]+
           params["Beta2"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"]))+
           params["Beta3"]*((1-exp(-t/params["Lambda1"]))/(t/params["Lambda1"])-exp(-t/params["Lambda1"]))+
           params["Beta4"]*((1-exp(-t/params["Lambda2"]))/(t/params["Lambda2"])-exp(-t/params["Lambda2"])))
} # end-function

##############################################################################
## Load data

# run data from tesouro direto

#bondData <- ReadTesouroDireto(inputDirectory,coupon=0.1,y1a=2003,y1b=2004,y2=2023)
bondData <- ReadTesouroDireto(inputDirectory,coupon=0.1,y1a=2002,y1b=2004,y2=2023)

bondData$YieldMid <- (bondData$RateBuy+bondData$RateSell)/2
bondData <- bondData[c(1,2,3,8,9)]
bondData$Years <- as.numeric((bondData$DateMaturity-bondData$DateReference)/daysYear)
bondData <- cbind(bondData[c(1,2,3,4)],bondData[c(6)],bondData[c(5)])
bondData <- subset(bondData,DateReference!=DateMaturity)

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
  #n <- trunc(nrow(subdates)/2)
  n=1
  #workingDates <- append(workingDates,sort(subdates$date,partial=n-1)[n-1])
  workingDates <- append(workingDates,sort(subdates$date,partial=n)[n])
}

#workingDates=c(as.Date("2021-10-15"),as.Date("2021-11-16"),as.Date("2021-12-14"),as.Date("2022-01-14"));k=1 # TESTE
#workingDates=c(as.Date("2023-05-12"));k=1 # TESTE


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
  ## Vertices
  
  par <- paramsOpt[,which.min(errors)]
  
  #par=c(Beta1=0.1258, Beta2=0.0064, Beta3=0.0940,Beta4=-0.1054, Beta5=0.0, Lambda1=2.4475, Lambda2=1.1826)
  
  mm <- nrow(paramsNSS)+1
  for (nn in 1:7) paramsNSS[mm,nn] <- par[nn]

  dateBusiness <- bizseq(from=dateReference+1, to=max(bondSet$DateMaturity),calendar)
  #numVertices <- trunc(length(dateBusiness)/63)+2
  numVertices <- trunc(length(dateBusiness)/126)
  l <- nrow(vertices)
  
  for (j in 1:numVertices) {
    #if (j==1) v <- 21 else if (j==2) v <- 42 else v <- (j-2)*63
    v <- j*126
    vertices[l+j,1] <- dateReference
    vertices[l+j,2] <- v
    vertices[l+j,3] <- dateBusiness[v]
    vertices[l+j,4] <- ModelFunction[[m]](params=par,t=as.numeric(vertices[l+j,3]-dateReference)/daysYear)
  }
  
  cat("\n","Params NSS =>",par,"\n")
  l1=l+1; l2=l+j
  print(vertices[l1:l2,])
  
  
  ###############################################################################
  # plot results

  # Plot model curves
  p = ggplot(data=bondSet,aes(x=DateMaturity))+
    geom_point(data=bondSet,aes(y=Yield,colour="Taxas Observadas"),shape=4,size=3)+ #geom_line(data=bondSet,aes(y=Yield,colour="Taxas Observadas"))+
    geom_point(data=bondSet,aes(y=NSSinYield,colour="NSS"),shape=15,size=2)+geom_line(aes(y=NSSinYield,colour="NSS"))+
    #geom_point(data=bondSet,aes(y=KRinYield,colour="KR"),shape=19,size=2)+geom_line(aes(y=KRinYield,colour="KR"))+
    theme(legend.position = c(0.8, 0.3),
          legend.direction = "vertical") +
    scale_colour_manual(name=NULL,
                        breaks = c("Taxas Observadas","NSS","KR"),
                        values = c("darkgreen","blue","red"))+
    labs(title = paste("Data referência",format(dateReference,"%Y-%m-%d")),
         y="Taxa",
         x="Vencimento [ano]")
  print(p)
  # skip to next date
} # end-for


###############################################################################
# Write data

write.xlsx(vertices,file=paste(outputDirectory,"yc_ettj.xlsx",sep=""),sheetName="ETTJ",append=FALSE)
write.xlsx(paramsNSS,file=paste(outputDirectory,"yc_ettj.xlsx",sep=""),sheetName="NSS",append=TRUE)

####################################################################################################
