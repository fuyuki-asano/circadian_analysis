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
tau = ['tau']
sig = ['significant']


df = df.drop(0, axis=1)
df = df.drop(range(0,11,1), axis=0)
df = df.reset_index(drop=True)

basal_LD_days = 14
DD_days = 21
skip_DD_days = 4 #DDはじめの4日間を除いて、データを折りたたむ
skip_points = (basal_LD_days + skip_DD_days)*24*60

time_points = len(df.index)
rec_days = time_points//(60*24)
bloks_number = 60*8 #period between 20:00 and 28:00 h

period_range = range(60*20, 60*28, 1)
#q = sympy.root(0.99, len(period_range))
sig_array = [chi2.ppf(q=0.99**(1/len(period_range)), df=per) for per in range(60*20-1, 60*28-1, 1)] #df=P-1
#sig_array2 = [chi2.ppf(q=0.01/len(period_range), df=per) for per in range(60*20-1, 60*28-1, 1)] #df=P-1


for i in range(len(df.columns)):
    column = (np.array(df.iloc[:, i])).astype(np.float64)
    c = np.argwhere(np.isnan(column))
    for d in np.reshape(c, (-1)):
        if d <= time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):(d + 11)]))
        elif d > time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):]))
#    df.iloc[:,i] = column
            
    Qp_array = []

    values = column[skip_points: (basal_LD_days + DD_days)*24*60]
    for ii in period_range:
        
        P = ii        
        K_folded_days = len(values)//ii        
        calc_values = values[:ii*K_folded_days]
        N = len(calc_values)
        folded_values = calc_values.reshape(-1, ii)
        M = np.mean(calc_values)        
        Mh = np.mean(folded_values, axis=0)
        
        numerator = np.sum(np.square(Mh-M))
        denom = np.sum(np.square(calc_values - M)) / (K_folded_days * P * K_folded_days)
        Qp = numerator / denom

        Qp_array.append(Qp)
    
    p_idx = np.argmax(Qp_array)
    period = period_range[p_idx]
    tau.append(period)
    
    Qp_max = Qp_array[p_idx]
    p_sig = sig_array[p_idx]
    
    if Qp_max >= p_sig:
        sig.append('True')
    else:
        sig.append('False')

        
    plt.plot(period_range, Qp_array)
    plt.plot(period_range, sig_array, 'red')
    plt.text(period, (Qp_max + 500), ('P=' + str(period)))
#    plt.plot(period_range, sig_array2, 'orange')
    plt.xlim([60*20,60*28])
    plt.ylim(0, round((Qp_max + 1600), -3))
    plt.xlabel("min")
    plt.ylabel("Qp")
    name = str(Ch[i+1, 0]) + "_X2_periodgram_.tiff"
    plt.savefig(name)
    plt.show()
        

tau = np.array(tau).reshape(-1, 1)
sig = np.array(sig).reshape(-1, 1)

data = np.concatenate((Ch, tau, sig), axis=1)

with open ('tau_X2_periodgram.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(data)


    


    