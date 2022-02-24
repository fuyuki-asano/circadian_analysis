# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 20:15:56 2021

@author: Fuyuki Asano
"""

import glob
import numpy as np
import pandas as pd
from scipy.stats import chi2
import sympy
import csv
import matplotlib.pyplot as plt
import sys


#wheel managerから排出されたファイルは.txt形式（名前は.xls）なので、.xlsに変換に保存しなおしてから使用すること
#mac
#path = glob.glob('/Users/asanofuyuki/Desktop/python working/*.xls')
#win
path = glob.glob('.\\*.xls')

df = pd.read_excel(path[0], sheet_name = 0, header = None)
Ch = np.concatenate((np.array(['Channel']).reshape(1,1), np.array(df.iloc[10, 1:]).reshape(-1,1)), axis=0)

df = df.drop(0, axis=1)
df = df.drop(range(0,11,1), axis=0)
df = df.reset_index(drop=True)


basal_LD_days = 14
DD_days = 21
recovery_LD_days = 7

header1 = np.array(['Day' + str(x) for x in range(1, basal_LD_days+1, 1)]).reshape((1,-1))
header2 = np.array(['Day' + str(x) for x in range(1, recovery_LD_days+1, 1)]).reshape((1,-1))

time_points = len(df.index)

for i in range(len(df.columns)):
    
    #nan padding(average around nan)
    column = (np.array(df.iloc[:, i])).astype(np.float64)
    c = np.argwhere(np.isnan(column))
    for d in np.reshape(c, (-1)):
        if  d < 10:
            column[d] = np.floor(np.nanmean(column[:(d+11)]))
        elif d > time_points-10:
            column[d] = np.floor(np.nanmean(column[(d-10):]))
        else:
            column[d] = np.floor(np.nanmean(column[(d-10):(d+11)]))
            
            
    basal_LD_val = np.array(column[:basal_LD_days*24*60])
    r_basal_val = basal_LD_val.reshape((basal_LD_days, -1))
    recovery_LD_val = np.array(column[(basal_LD_days + DD_days)*24*60:])
    r_recovery_val = recovery_LD_val.reshape((recovery_LD_days, -1))
    
    
    stack1 = []
    for ii in range(basal_LD_days):   
        s = 60*6
        a1 = np.array([np.sum(r_basal_val[ii,:][x:x+s]) for x in range(60*24-2*s)])
        a2 = np.array([np.sum(r_basal_val[ii,:][x+s:x+2*s]) for x in range(60*24-2*s)])
        gap1 = a2-a1
        onset = ((np.argmax(gap1)+s))/60
        stack1.append(onset)
    
    astack1 = np.array(stack1).reshape((1,-1))
    header1 = np.concatenate((header1, astack1), axis=0)
    
    
    stack2 = []
    for ii in range(recovery_LD_days):   
        s = 60*6
        a3 = np.array([np.sum(r_recovery_val[ii,:][x:x+s]) for x in range(60*24-2*s)])
        a4 = np.array([np.sum(r_recovery_val[ii,:][x+s:x+2*s]) for x in range(60*24-2*s)])
        gap2 = a4-a3
        onset = ((np.argmax(gap2)+s))/60
        stack2.append(onset)
    
    astack2 = np.array(stack2).reshape((1,-1))
    header2 = np.concatenate((header2, astack2), axis=0)
    
    
data = np.concatenate((Ch, header1, header2), axis=1)
with open ('LD_activity_onset.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(data)


