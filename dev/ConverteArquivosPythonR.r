

### Manipula arquivos
np <- import("numpy",convert=FALSE)
dirData <- "../database/"

filePathName <- paste(dirData,"price_2013-12-31.npy",sep="")
preco <- as.matrix(np$load(filePathName))

filePathNameNew <- paste(dirData,"precoR_2013-12-31.npy",sep="")
np$save(filePathNameNew,preco)

priceVector <- as.matrix(np$load(filePathNameNew))


library(Matrix)
sps = import("scipy.sparse")
filePathName <- paste(dirData,"cashflow_2013-12-31.npz",sep="")
fluxo <- as.matrix(sps$load_npz(filePathName))
fluxoSparse <- Matrix(fluxo,sparse=TRUE)

filePathNameNew <- paste(dirData,"fluxoR_2013-12-31.npz",sep="")
sps$save_npz(filePathNameNew,fluxoSparse)

cashflowMatrix <- as.matrix(sps$load_npz(filePathNameNew))

KernelRidgeDate(dateReference,r_to_py(priceVector),r_to_py(cashflowMatrix))


## FIM