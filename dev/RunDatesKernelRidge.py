# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

## Clear environment
for obj in dir():
  if not obj.startswith("__"):
    del globals()[obj]
del obj


### Import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps 
import time
import sys
import os

import pandas as pd
import openpyxl

### download KR example files from Github
### if not os.path.exists('KR_example'):
  ###  !git clone https://github.com/yye9701/KR_example.git
dir_kr = 'C:/Users/lcgfa/OneDrive/$FGV/Computação Aplicada/05 Trabalho Final/Yield Curves/R/'

print(dir_kr+'lib/')
sys.path.append(dir_kr+'lib/')
import kr_model
import kr_utils
import yc_utils 


### settings
# kernel hyper-parameters
alpha=0.05
delta=0.00

# max time to maturity in days
N=30*365

start_time=time.time()
K=kr_model.generate_kernel_matrix(alpha, delta, N, N)
end_time=time.time()

print('Time elapsed for generating a {}-by-{} kernel matrix with alpha = {} and delta = {}: {:.1f} sec'.format(N,N,alpha, delta, end_time-start_time))

### load example data
## dir_data=dir_kr+'example_data/'
dir_data='C:/Users/lcgfa/OneDrive/$FGV/Computação Aplicada/05 Trabalho Final/Yield Curves/R/database/'
example_dates=['1961-06-30', '2013-12-31']


## Gera matrizes para cada data

start_time=time.time()

dict_data={}

for date in example_dates:
    # load price vector and cashflow matrix
    B=np.load(dir_data+'preco_{}.npy'.format(date))
    C=sps.load_npz(dir_data+'fluxo_{}.npz'.format(date)).toarray()
    M=B.shape[0]

    print('Load Date: {}; Number of securities: {}'.format(date, M))

    dict_data[date]={'B':B, 'C':C}

for date in example_dates:
    # load price vector and cashflow matrix
    B=dict_data[date]["B"]
    C=dict_data[date]["C"]
    M=B.shape[0]
    
    print('Calculate Date: {}; Number of securities: {}'.format(date, M))
    
    # Calculate YTM and Duration
    dict_data[date]=yc_utils.ComputeMatrices(date,B,C,M)
    
end_time=time.time()

print('Time elapsed for loading and preparing data for {} example dates: {:.1f} sec'.format(len(example_dates),end_time-start_time))    



import yc_utils


for date in example_dates:
    # Plot Prices and Yields
    yc_utils.PlotPriceYield(dict_data,date)
    
    # Plot Cashflow
    yc_utils.PlotCashFlow(dict_data,date)
    


### fit KR model on example dates
# KR ridge penalty term 
ridge=1

start_time=time.time()

dict_fit={}
for date in example_dates:
    
    dict_fit[date]=kr_model.KR(C=dict_data[date]['C'], # cashflow matrix
                         B=dict_data[date]['B'], # price vector
                         ridge=ridge, # ridge hyper-parameter
                         inv_w=dict_data[date]['inv_w'], # inverse of the weighting vector
                         K=K # kernel matrix
                        )
                               
end_time=time.time()

print('Time elapsed for fitting KR model for {} example dates: {:.1f} sec'.format(len(example_dates),end_time-start_time))


for date in example_dates:
    # Plot fitted prices and yields
    yc_utils.PlotFitted(dict_data,dict_fit,date)



