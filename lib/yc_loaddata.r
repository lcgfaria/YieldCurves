# Identification ----------------------------------------------------------
# 
# FGV - Fundacao Getulio Vargas
# EESP - Escola de Economia São Paulo
# MPEF - Mestrado Profissional em Economia e Finanças
#
# YIELD CURVES COMPARISON
# Autor: Luis Giovanni Faria
#
# Date: 13/09/2022
#

## Functions Library


################################################################################
## Read historic price files downloaded from brazilian treasure site
## https://www.tesourodireto.com.br/titulos/historico-de-precos-e-taxas.htm
################################################################################


ReadBondFile <- function(path, type, year) {
  ## Function to read a given bond file with multiple maturity dates
  ## path: path file
  ## type: bond type (LTN or NT-F)
  ## year: year with reference dates and prices
  pathFileName <- paste(path, type, "_", year, ".xls", sep="")
  cat("Reading file: ", pathFileName, "\n")
  # Initialize variables
  numberOfBonds <- length(excel_sheets(pathFileName))
  bonds <- data.frame(matrix(ncol=8, nrow=0))
  # Loop to read multiple sheets
  for (j in 1:numberOfBonds) {
    d <- read_excel(pathFileName, sheet = j, range="B1:B1", col_names=c("DateMaturity"))[[1]]
    maturity <- as.Date(d, format ="%d/%m/%Y")
    if (year<=2011) {tipoData <- "date"} else {tipoData <- "text"}
    prices <- read_excel(pathFileName, sheet = j, skip=2, 
                         col_types = c(tipoData, "numeric", "numeric", "numeric", "numeric", "numeric"), 
                         col_names=c("DateReference", "RateBuy", "RateSell", "PriceBuy", "PriceSell", "PricePU"))
    prices <- cbind(type,maturity,prices)
    prices[,3] <- as.Date(prices[,3], format="%d/%m/%Y")
    bonds <- rbind(bonds, prices)
    cat("File=", type, "_", year, " Sheet=", j, " Maturity=", as.character(maturity), " Records=", nrow(prices), "\n", sep="")
  }
  return(bonds)  
} # end-function


ReadTesouroDireto <- function (path="../database/",coupon=0.1,y1a=2004,y1b=2004,y2=2022) {
  ## Read bond prices from multiple official treasure files
  ## path: path file
  ## coupon: coupon value for NTN-F (default=0.01)
  ##
  # Initialize data frame of prices
  
  bondPrices <- data.frame(matrix(ncol=8, nrow=0))
  
  # Year range
  bondTypes <- c("LTN", "NTN-F")
  bondRange <- data.frame(matrix(ncol=length(bondTypes), nrow=2))
  colnames(bondRange) <- bondTypes
  rownames(bondRange) <- c("YearBegin", "YearEnd")
  bondRange[,1] <- c(y1a, y2) # LTN range
  bondRange[,2] <- c(y1b, y2) # NTN-F range

  # Read Multiple files
  # Loop to read multiple files (bond types x years)
  for (t in bondTypes) {
    for (i in bondRange[1,t]:bondRange[2,t]) {
      # File names
      prices <- ReadBondFile(path, t, i)
      bondPrices <- rbind(bondPrices, prices)
    } # end-for
  } # end-for

  colnames(bondPrices) <- c("BondType","DateMaturity","DateReference","RateBuy","RateSell","PriceBuy","PriceSell","PricePU")  
  bondPrices$Coupon <- c(0)

  # load coupons
  for (i in 1:nrow(bondPrices)) {
    if (bondPrices$BondType[i]=="LTN") bondPrices$Coupon[i] <- 0 else bondPrices$Coupon[i] <- coupon
  } # end-for  

  bondPrices <- cbind(bondPrices[c(3)],bondPrices[c(1,2,4,5,6,7,9,8)])
  
  return(subset(bondPrices, !is.na(PricePU)))
} # end-function


ReadExampleDateKR <- function(path="../database/",dateReference) {
  ## read original example files
  ## path: path file
  ## dateReference
  ## 
  
  # read python files
  np <- import("numpy",convert=FALSE)
  priceVector <- as.matrix(np$load(paste(path,"precoR_",as.character(dateReference),".npy",sep="")))
  
  sps <- import("scipy.sparse")
  cashflowMatrix <- as.matrix(sps$load_npz(paste(path,"fluxoR_",as.character(dateReference),".npz",sep="")))
  ##fluxoSparse <- Matrix(fluxo,sparse=TRUE)

  bondSet <- data.frame(matrix(ncol=4, nrow=nrow(priceVector)))
  colnames(bondSet) <- c("BondType","DateMaturity","PricePU","Coupon")
  
  bondSet$BondType <- rep("BOND",nrow(bondSet))
  bondSet$DateMaturity <- as.Date(bondSet$DateMaturity)
  
  # get prices    
  bondSet$PricePU <- priceVector
  # get maturity dates and coupons
  for (i in 1:nrow(bondSet)) {
    x=as.Date(dateReference+which.max(cashflowMatrix[i,]))
    bondSet$DateMaturity[i] <- as.Date(dateReference+which.max(cashflowMatrix[i,]))
    bondSet$Coupon[i] <- (max(cashflowMatrix[i,])-100)*2/100
  } # end-for  
  return(list(bondSet,priceVector,cashflowMatrix))
}






