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
bins = np.array(range(0,144,1)).reshape(1,-1)

df = df.drop(0, axis=1)
df = df.drop(range(0,11,1), axis=0)
df = df.reset_index(drop=True)

basal_LD_days = 14
DD_days = 21
recovery_LD_days = 7

time_points = len(df.index)
rec_days = time_points//(60*24)

if basal_LD_days + DD_days + recovery_LD_days == rec_days:
    pass
else:
    raise ValueError('Check rec schedule')


for i in range(len(df.columns)):
    column = (np.array(df.iloc[:, i])).astype(np.float64)
    c = np.argwhere(np.isnan(column))
    for d in np.reshape(c, (-1)):
        if d <= time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):(d + 11)]))
        elif d > time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):]))
            
    basal_LD_act = column[0:basal_LD_days*60*24]
    basal_LD_act_10 = np.sum(basal_LD_act.reshape(-1, 10), axis=1)
    basal_LD_act_10r = basal_LD_act_10.reshape(basal_LD_days, -1)
    act_mean = np.mean(basal_LD_act_10r, axis=0).reshape(1,-1)
    
    bins = np.concatenate((bins, act_mean), axis=0)

data = np.concatenate((Ch, bins), axis=1)

with open ('basal_LD_mean_act.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(data)


    


    