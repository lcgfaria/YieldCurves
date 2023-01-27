
## Clear environment
##for obj in dir():
##  if not obj.startswith("__"):
##    del globals()[obj]
##del obj



### Import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps 
import time
import sys
import os

import pandas as pd
import openpyxl


dir_kr = 'C:/Users/lcgfa/OneDrive/$FGV/Computação Aplicada/Yield Curves/R/'

print(dir_kr+'lib/')
sys.path.append(dir_kr+'lib/')
import kr_model
import kr_utils
import yc_utils 
import yc_krdate


### load example data
## dir_data=dir_kr+'example_data/'
dir_data='C:/Users/lcgfa/OneDrive/$FGV/Computação Aplicada/Yield Curves/R/database/'
example_dates=['1961-06-30', '2013-12-31']
date=example_dates[1]
# load price vector and cashflow matrix
B=np.load(dir_data+'preco_{}.npy'.format(date))
C=sps.load_npz(dir_data+'fluxo_{}.npz'.format(date)).toarray()
M=B.shape[0]

yc_krdate.KernelRidgeDate(date,B,C)




