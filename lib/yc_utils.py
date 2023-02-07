# Identification ----------------------------------------------------------
# 
# FGV - Fundacao Getulio Vargas
# EESP - Escola de Economia São Paulo
# MPEF - Mestrado Profissional em Economia e Finanças
#
# YIELD CURVES COMPARISON
# Autor: Luis Giovanni Faria
#
# Date: 23/01/2023
#

## Functions Library

### Import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps 
import time
import sys
import os
import pandas as pd
import openpyxl
import pdb

dir_kr = 'C:/Users/lcgfa/OneDrive/$FGV/Yield Curves/R/'
sys.path.append(dir_kr+'lib/')
#print(os.getcwd())
#sys.path.append("../lib/")

import kr_model
import kr_utils

def KRcurve(date,B,C,years):
    ## Description: Solve KR yield curve
    ## Arguments:
    ##  date: reference date
    ##  B: price vector
    ##  C: cashflow matrix
    ##  years: maximum number of years to maturity
    ## Return
    ##  dict_out: yield curve solved (g_solved, y_solved)
    ##
    M=B.shape[0]
    B=B.flatten()

    ### settings
    # kernel hyper-parameters
    alpha=0.05
    delta=0.00
    # KR ridge penalty term (original = 1)
    ridge=1

    # max time to maturity in days (original was 30 years)
    N=years*365
    K=kr_model.generate_kernel_matrix(alpha, delta, N, N)

    # get YTM and duration
    ytm, dur=np.zeros(M), np.zeros(M) # YTM and duration
    
    ttm = np.zeros(M) # time to maturity in days
    
    for i in range(M):
        time_to_cashflow_inday=np.where(C[i]!=0)[0]+1
        ytm[i], dur[i] = kr_utils.get_ytm_and_duration(C[i][time_to_cashflow_inday-1], time_to_cashflow_inday, B[i])
        ttm[i] = max(time_to_cashflow_inday)

    # get inverse of weights for fitting. Weights w is computed as w=1/inv_w
    inv_w=(dur*B)**2*M

    dict_in={'B':B, 'C':C, 'ytm':ytm, 'dur':dur, 'ttm':ttm, 'inv_w':inv_w}
    
    ### fit KR model on example dates
    dict_out=kr_model.KR(C=dict_in['C'], # cashflow matrix
                         B=dict_in['B'], # price vector
                         ridge=ridge, # ridge hyper-parameter
                         inv_w=dict_in['inv_w'], # inverse of the weighting vector
                         K=K # kernel matrix
                         )

    return dict_out


###############################################################################

def FitKR(C,g_solved,date):
    ## Description: Get fitted values
    ## Arguments:
    ##  date: reference date
    ##  C: cashflow matrix
    ##  g_hat: yield estimator
    ## Return
    ##  output (B, ytm, dur, ttm, inv_w)
    ##

    # load price vector and cashflow matrix
    M=C.shape[0]
  
    # calculate implied prices by fitted discount curve
    B=C@g_solved

    # get YTM and duration
    ytm, dur = np.zeros(M), np.zeros(M) # YTM and duration
    ttm = np.zeros(M) # time to maturity in days
    
    for i in range(M):
        time_to_cashflow_inday=np.where(C[i]!=0)[0]+1
        ytm[i],dur[i] = kr_utils.get_ytm_and_duration(C[i][time_to_cashflow_inday-1],time_to_cashflow_inday,B[i])
        ttm[i] = max(time_to_cashflow_inday)
  
    # get inverse of weights for fitting w=1/inv_w
    inv_w=(dur*B)**2*M
        
    return B,ytm,dur,ttm,inv_w
  
###############################################################################


def GetFittedOLD(dict_data,dict_fit,date):
    # load price vector and cashflow matrix
    C=dict_data[date]['C']
    M=C.shape[0]
  
    # calculate implied prices by fitted discount curve
    B_fitted=C@dict_fit[date]['g_solved'][:C.shape[1]]

    # get YTM and duration
    ytm_fitted = np.zeros(M) # YTM and duration
    #ttm = np.zeros(M) # time to maturity in days
    
    for i in range(M):
        time_to_cashflow_inday=np.where(C[i]!=0)[0]+1
        ytm_fitted[i], _ = kr_utils.get_ytm_and_duration(C[i][time_to_cashflow_inday-1],
                                                         time_to_cashflow_inday,
                                                         B_fitted[i])
    return ytm_fitted,B_fitted
  


###############################################################################


def PlotPriceYieldOLD(dict_data,date):
    p=plt.figure()
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
    
    return
    
###############################################################################

def PlotCashFlowOLD(dict_data,date):
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
    
    return

###############################################################################

def PlotFittedOLD(dict_data,dict_fit,date):
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
    
    return
  

###############################################################################

def ComputeMatricesOLD(date,B,C,M):
    '''
    ## - Calculate annualized YTM (not in %) and duration (in years) of a security
    - YTM is estimated using Newton's method
    - Assume (1) continuous compounding and (2) each year has 365 days
    - Args:
        - cashflow (numpy array): amount of cashflow
        - time_to_cashflow_inday (numpy array): time to cashflow in days
        - B_i (float): price of the security
        - y_guess (float): initial guess for YTM in Newton's method
    - Returns:
        - ytm_solved (float): estimated YTM
        - dur_solved (float): estimated duration in years
    '''
    # print('Date: {}; Number of securities: {}'.format(date, M))
    
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

    dict={'B':B, 'C':C, 'ytm':ytm, 'dur':dur, 'ttm':ttm, 'inv_w':inv_w}
    
    return dict

###############################################################################

