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
tau_bLD = ['tau_basal_LD']
sig_bLD = ['significant_basal_LD']
tau_DD = ['tau_DD']
sig_DD = ['significant_DD']


df = df.drop(0, axis=1)
df = df.drop(range(0,11,1), axis=0)
df = df.reset_index(drop=True)

basal_LD_days = 14
DD_days = 21
skip_DD_days = 4


time_points = len(df.index)
rec_days = time_points//(60*24)
bloks_number = 60*8 #period between 20:00 and 28:00 h

period_range = range(60*20, 60*28, 1)
sig_array = [chi2.ppf(q=0.99**(1/len(period_range)), df=per) for per in range(60*20-1, 60*28-1, 1)] #df=P-1
#sig_array2 = [chi2.ppf(q=0.01/len(period_range), df=per) for per in range(60*20-1, 60*28-1, 1)] #df=P-1

sig_array_for_exp = np.concatenate((np.array(sig_array).reshape((1,-1)), \
                                    np.array([' ']).reshape((1,-1)), \
                                        np.array(sig_array).reshape((1,-1))), axis=1)


for i in range(len(df.columns)):
    column = (np.array(df.iloc[:, i])).astype(np.float64)
    c = np.argwhere(np.isnan(column))
    for d in np.reshape(c, (-1)):
        if d <= time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):(d + 11)]))
        elif d > time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):]))
#    df.iloc[:,i] = column

            
    Qp_array_bLD = []
    Qp_array_DD = []

    val_bLD = column[:basal_LD_days*24*60]
    val_DD = column[(basal_LD_days + skip_DD_days)*24*60 : (basal_LD_days + DD_days)*24*60]
    
    for ii in period_range:
        
        P = ii
        
        K_folded_days_bLD = len(val_bLD)//ii        
        calc_values_bLD = val_bLD[:ii*K_folded_days_bLD]
        N_bLD = len(calc_values_bLD)
        folded_values_bLD = calc_values_bLD.reshape(-1, ii)
        M_bLD = np.mean(calc_values_bLD)        
        Mh_bLD = np.mean(folded_values_bLD, axis=0)
        
        numerator_bLD = np.sum(np.square(Mh_bLD - M_bLD))
        denom_bLD = np.sum(np.square(calc_values_bLD - M_bLD)) / (K_folded_days_bLD * P * K_folded_days_bLD)
        Qp_bLD = numerator_bLD / denom_bLD
        
        Qp_array_bLD.append(Qp_bLD)
        
        
        K_folded_days_DD = len(val_DD)//ii        
        calc_values_DD = val_DD[:ii*K_folded_days_DD]
        N_DD = len(calc_values_DD)
        folded_values_DD = calc_values_DD.reshape(-1, ii)
        M_DD = np.mean(calc_values_DD)        
        Mh_DD = np.mean(folded_values_DD, axis=0)
        
        numerator_DD = np.sum(np.square(Mh_DD - M_DD))
        denom_DD = np.sum(np.square(calc_values_DD - M_DD)) / (K_folded_days_DD * P * K_folded_days_DD)
        Qp_DD = numerator_DD / denom_DD

        Qp_array_DD.append(Qp_DD)

    
    p_idx_bLD = np.argmax(Qp_array_bLD)
    period_bLD = (period_range[p_idx_bLD])
    tau_bLD.append(period_bLD)
    
    Qp_max_bLD = Qp_array_bLD[p_idx_bLD]
    p_sig_bLD = sig_array[p_idx_bLD]
    
    if Qp_max_bLD >= p_sig_bLD:
        sig_bLD.append('True')
    else:
        sig_bLD.append('False')
        
        
    p_idx_DD = np.argmax(Qp_array_DD)
    period_DD = (period_range[p_idx_DD])
    tau_DD.append(period_DD)
    
    Qp_max_DD = Qp_array_DD[p_idx_DD]
    p_sig_DD = sig_array[p_idx_DD]
    
    if Qp_max_DD >= p_sig_DD:
        sig_DD.append('True')
    else:
        sig_DD.append('False')
        
        
    Qp = np.concatenate((np.array(Qp_array_bLD).reshape((1,-1)), \
                         np.array([' ']).reshape((1,-1)), \
                             np.array(Qp_array_DD).reshape((1,-1))), axis=1)
    sig_array_for_exp = np.concatenate((sig_array_for_exp, Qp), axis=0)
    
                
    plt.plot(period_range, Qp_array_bLD)
    plt.plot(period_range, sig_array, 'red')
    plt.text(period_bLD, (Qp_max_bLD + 500), ('P=' + str(period_bLD)))
#    plt.plot(period_range, sig_array2, 'orange')
    plt.xlim([60*20,60*28])
    plt.ylim(0, round((Qp_max_bLD + 1600), -3))
    plt.xlabel("min")
    plt.ylabel("Qp")
    name_bLD = str(Ch[i+1, 0]) + "_X2_periodgram_basal_LD.tiff"
    plt.savefig(name_bLD)
    plt.show()
    
    
    plt.plot(period_range, Qp_array_DD)
    plt.plot(period_range, sig_array, 'red')
    plt.text(period_DD, (Qp_max_DD + 500), ('P=' + str(period_DD)))
#    plt.plot(period_range, sig_array2, 'orange')
    plt.xlim([60*20,60*28])
    plt.ylim(0, round((Qp_max_DD + 1600), -3))
    plt.xlabel("min")
    plt.ylabel("Qp")
    name_DD = str(Ch[i+1, 0]) + "_X2_periodgram_DD.tiff"
    plt.savefig(name_DD)
    plt.show()
        

tau_bLD = np.array(tau_bLD).reshape(-1, 1)
sig_bLD = np.array(sig_bLD).reshape(-1, 1)

tau_DD = np.array(tau_DD).reshape(-1, 1)
sig_DD = np.array(sig_DD).reshape(-1, 1)

data_tau = np.concatenate((Ch, tau_bLD, sig_bLD, tau_DD, sig_DD,), axis=1)
data = np.concatenate((data_tau, sig_array_for_exp), axis=1)

with open ('tau_X2_periodgram_Qp_amplitude.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(data)


    


    