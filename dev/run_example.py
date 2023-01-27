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

## Code to manipulate files (LGF)

## Grava preços e cashflow em excel a partir dos aruivos npy originais

gravaPlanilha = True

if gravaPlanilha:
    for date in example_dates:
        B=np.load(dir_data+'price_{}.npy'.format(date))
        print(B[:5])
        C=sps.load_npz(dir_data+'cashflow_{}.npz'.format(date)).toarray()
        print(C[:5])
        dados = np.load(dir_data+'cashflow_{}.npz'.format(date), allow_pickle=True)
                
        M=B.shape[0]
        print(M)
        
        print('Date: {}; Number of securities: {}'.format(date, M))
    
        Bdf = pd.DataFrame(B)
        Bdf.to_excel(dir_data+'preco_{}.xlsx'.format(date), sheet_name=date, index=False, header=False)
    
        Cdf = pd.DataFrame(C)
        Cdf.to_excel(dir_data+'fluxo_{}.xlsx'.format(date), sheet_name=date, index=False, header=False)
    
   

## Converte arquivos Excel para npz/npy

lePlanilha = True

if lePlanilha:
    for date in example_dates:
        Bdf = pd.read_excel(dir_data+'preco_{}.xlsx'.format(date), header=None, engine="openpyxl") 
        print(Bdf.head())
        B = np.array(Bdf[0])
        print(B.shape)
        np.save(dir_data+'preco_{}.npy'.format(date), B)
        
        Cdf = pd.read_excel(dir_data+'fluxo_{}.xlsx'.format(date), header=None, engine="openpyxl") 
        print(Cdf.head())
        C = np.asarray(Cdf)
        sps.save_npz(dir_data+'fluxo_{}.npz'.format(date), sps.csr_matrix(C))
        
        M=B.shape[0]
        print(M)


## Gera matrizes para cada data

import yc_utils 

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
    dict_data[date]=yc_utils.yc_ComputeMatrices(date,B,C,M)
    
end_time=time.time()

print('Time elapsed for loading and preparing data for {} example dates: {:.1f} sec'.format(len(example_dates),end_time-start_time))    


    
    
for date in example_dates:
    # Plot Prices and Yields
    yc_PlotPriceYield(dict_data,date)
    
    # Plot Cashflow
    yc_PlotCashflow(dict_data,date)
    
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
    yc_PlotFitted(dict_data,dict_fit,date)







######## Original code

start_time=time.time()

dict_data={}
for date in example_dates:
    # load price vector and cashflow matrix
    B=np.load(dir_data+'price_{}.npy'.format(date))
    C=sps.load_npz(dir_data+'cashflow_{}.npz'.format(date)).toarray()
    M=B.shape[0]
    
    print('Date: {}; Number of securities: {}'.format(date, M))

    # get YTM and duration
    ytm, dur=np.zeros(M), np.zeros(M) # YTM and duration
    ttm = np.zeros(M) # time to maturity in days
    for i in range(M):
        time_to_cashflow_inday=np.where(C[i]!=0)[0]+1
        ytm[i], dur[i] = kr_utils.get_ytm_and_duration(C[i][time_to_cashflow_inday-1], time_to_cashflow_inday, B[i])
        ttm[i] = max(time_to_cashflow_inday)

    # get inverse of weights for fitting
    # weights w is computed as w=1/inv_w
    inv_w=(dur*B)**2*M
    
    dict_data[date]={'B':B, 'C':C, 'ytm':ytm, 'dur':dur, 'ttm':ttm, 'inv_w':inv_w}

end_time=time.time()

print('Time elapsed for loading and preparing data for {} example dates: {:.1f} sec'.format(len(example_dates),end_time-start_time))    

for date in example_dates:

    B=dict_data[date]['B']
    ttm=dict_data[date]['ttm']
    ytm=dict_data[date]['ytm']

    fig=plt.figure(figsize=(15,6))
    ax_1,ax_2=fig.add_subplot(1,2,1), fig.add_subplot(1,2,2)
    
    ax_1.scatter(ttm/365, B)
    ax_2.scatter(ttm/365, ytm)
    
    ax_1.set_title('Observed prices on {}'.format(date))
    ax_2.set_title('YTM on {}'.format(date))
    
    for ax in [ax_1, ax_2]:
        ax.set_xlabel('Time to maturity in years');
    
    plt.show()

for date in example_dates:

    B=dict_data[date]['B']
    C=dict_data[date]['C']
    M=B.shape[0]

    time_to_cashflow_inday=[]
    for i in range(M):
        time_to_cashflow_inday.append(np.where(C[i]!=0)[0]+1)

    xs=np.concatenate(time_to_cashflow_inday)
    ys=np.concatenate([i* np.ones(time_to_cashflow_inday[i].shape) for i in range(M)])
    color=np.concatenate([C[i, time_to_cashflow_inday[i]-1] for i in range(M)])


    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(1,1,1)
    im = ax.scatter(x=xs/365, y=1+ys, c=color, cmap='Reds', alpha=0.7, s=15)
    fig.colorbar(im,ax=ax)

    ax.set_xlabel('Time to cashflow in years');
    ax.set_ylabel('Security index');
    ax.set_title('Cashflow matrix on {}'.format(date))

    plt.show()

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
    
    # load price vector and cashflow matrix
    C=dict_data[date]['C']
    M=C.shape[0]
  
    # calculate implied prices by fitted discount curve
    B_fitted=C@dict_fit[date]['g_solved'][:C.shape[1]]

    # get YTM and duration
    ytm_fitted = np.zeros(M) # YTM and duration
    ttm = np.zeros(M) # time to maturity in days
    
    for i in range(M):
        time_to_cashflow_inday=np.where(C[i]!=0)[0]+1
        ytm_fitted[i], _ = kr_utils.get_ytm_and_duration(C[i][time_to_cashflow_inday-1],
                                                         time_to_cashflow_inday,
                                                         B_fitted[i])
        ttm[i] = max(time_to_cashflow_inday)
        
    
    # plot
    fig=plt.figure(figsize=(15,6))
    ax_1,ax_2=fig.add_subplot(1,2,1), fig.add_subplot(1,2,2)
    

    ax_1.scatter(ttm/365, dict_data[date]['B'], marker='o', s=15, label='Observed Price');
    ax_1.scatter(ttm/365, B_fitted, marker='x', label='Fitted Price');
    
    ax_2.scatter(ttm/365, dict_data[date]['ytm'], marker='o', s=15, label='Observed YTM');
    ax_2.scatter(ttm/365, ytm_fitted, marker='x', label='Fitted YTM');
  

    ax_1.set_title('Observed and fitted price on {}'.format(date));
    ax_2.set_title('Observed and fitted YTM on {}'.format(date));
    
    for ax in [ax_1, ax_2]:
        ax.legend();
        ax.set_xlabel('Time to maturity in years');
    
    plt.show()


######################
