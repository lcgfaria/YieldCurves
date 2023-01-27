## Testa datas exemplo do KR

# Remove working set
rm(list=ls())

# Working Directory
sysInfo <- Sys.info()
if (sysInfo[[4]] == "GIO-YOGA") {
  workingDirectory <- "C:/Users/giovanni/OneDrive/$FGV/Computação Aplicada/Yield Curves/R/dev"
} else if (sysInfo[[4]] == "DESKTOP-G3EH8EA") {
  workingDirectory <- "C:/Users/lcgfa/OneDrive/$FGV/Computação Aplicada/Yield Curves/R/dev"
} else {
  workingDirectory <- getwd()
}
setwd(workingDirectory)
getwd()

# Load R libraries
library(readxl)
library(bizdays)

library(reticulate)
library(Matrix)

source_python("../lib/kr_model.py")
source_python("../lib/kr_utils.py")
source_python("../lib/yc_utils.py")
source_python("../lib/yc_krdate.py")

dirData <- "../database/"


## roda KR com arquivos sample originais
dateReference <- as.Date("2013-12-31")

np <- import("numpy",convert=FALSE)

##np <- import("numpy",convert=FALSE)
filePathName <- paste(dirData,"precoR_2013-12-31.npy",sep="")
priceVector <- as.matrix(np$load(filePathName))

sps <- import("scipy.sparse")
filePathName <- paste(dirData,"fluxoR_2013-12-31.npz",sep="")
cashflowMatrix <- as.matrix(sps$load_npz(filePathName))
##fluxoSparse <- Matrix(fluxo,sparse=TRUE)


# converte objetos de R para python
priceVector <- r_to_py(priceVector)
cashflowMatrix <- r_to_py(cashflowMatrix)

## chama rotina em python
KernelRidgeDate(dateReference,priceVector,cashflowMatrix)

##### FIM
#######################################################




## executa programa em python
py_run_file("krdate_test.py")


