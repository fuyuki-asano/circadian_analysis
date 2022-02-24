# -*- coding: utf-8 -*-
"""
Created on Fri May 14 23:06:21 2021

@author: Sleepymouse
"""

import glob
import numpy as np
import pandas as pd
from scipy.stats import chi2
import sympy
import csv
import matplotlib.pyplot as plt
import sys

#wheel managerから排出されたファイルは.csv形式（名前は.xls）なので、.xlsに変換sにて保存しなおしてから使用すること
#mac
#path = glob.glob('/Users/asanofuyuki/Desktop/python working/*.xls')
#win
path = glob.glob('.\\*.xls')
df = pd.read_excel(path[0], sheet_name = 0, header = None)
Ch = np.concatenate((np.array(['Channel']).reshape(1,1), np.array(df.iloc[10, 1:]).reshape(-1,1)), axis=0)

tau_LL = ['tau_LL']
sig_LL = ['significant_LL']


df = df.drop(0, axis=1)
df = df.drop(range(0,11,1), axis=0)
df = df.reset_index(drop=True)

basal_LD_days = 7
LL_days = 21
skip_LL_days = 7


time_points = len(df.index)
rec_days = time_points//(60*24)
bloks_number = 60*8 #period between 20:00 and 28:00 h

period_range = range(60*20, 60*28, 1)
sig_array = [chi2.ppf(q=0.99**(1/len(period_range)), df=per) for per in range(60*20-1, 60*28-1, 1)] #df=P-1
#sig_array2 = [chi2.ppf(q=0.01/len(period_range), df=per) for per in range(60*20-1, 60*28-1, 1)] #df=P-1

sig_array_for_exp = np.array(sig_array).reshape((1,-1))


for i in range(len(df.columns)):
    column = (np.array(df.iloc[:, i])).astype(np.float64)
    c = np.argwhere(np.isnan(column))
    for d in np.reshape(c, (-1)):
        if d <= time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):(d + 11)]))
        elif d > time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):]))
#    df.iloc[:,i] = column

    Qp_array_LL = []

    val_LL = column[(basal_LD_days + skip_LL_days)*24*60 : (basal_LD_days + LL_days)*24*60]
    
    for ii in period_range:
        
        P = ii

        K_folded_days_LL = len(val_LL)//ii        
        calc_values_LL = val_LL[:ii*K_folded_days_LL]
        N_LL = len(calc_values_LL)
        folded_values_LL = calc_values_LL.reshape(-1, ii)
        M_LL = np.mean(calc_values_LL)        
        Mh_LL = np.mean(folded_values_LL, axis=0)
        
        numerator_LL = np.sum(np.square(Mh_LL - M_LL))
        denom_LL = np.sum(np.square(calc_values_LL - M_LL)) / (K_folded_days_LL * P * K_folded_days_LL)
        Qp_LL = numerator_LL / denom_LL

        Qp_array_LL.append(Qp_LL)

        
    p_idx_LL = np.argmax(Qp_array_LL)
    period_LL = (period_range[p_idx_LL])
    tau_LL.append(period_LL)
    
    Qp_max_LL = Qp_array_LL[p_idx_LL]
    p_sig_LL = sig_array[p_idx_LL]
    
    if Qp_max_LL >= p_sig_LL:
        sig_LL.append('True')
    else:
        sig_LL.append('False')
        
        
    Qp = np.array(Qp_array_LL).reshape((1,-1))
    sig_array_for_exp = np.concatenate((sig_array_for_exp, Qp), axis=0)
    
    
    plt.plot(period_range, Qp_array_LL)
    plt.plot(period_range, sig_array, 'red')
    plt.text(period_LL, (Qp_max_LL + 500), ('P=' + str(period_LL)))
#    plt.plot(period_range, sig_array2, 'orange')
    plt.xlim([60*20,60*28])
    plt.ylim(0, round((Qp_max_LL + 1600), -3))
    plt.xlabel("min")
    plt.ylabel("Qp")
    name_LL = str(Ch[i+1, 0]) + "_X2_periodgram_LL.tiff"
    plt.savefig(name_LL)
    plt.show()
        

tau_LL = np.array(tau_LL).reshape(-1, 1)
sig_LL = np.array(sig_LL).reshape(-1, 1)

data_tau = np.concatenate((Ch, tau_LL, sig_LL,), axis=1)
data = np.concatenate((data_tau, sig_array_for_exp), axis=1)

with open ('tau_X2_periodgram_Qp_amplitude.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(data)


    


    