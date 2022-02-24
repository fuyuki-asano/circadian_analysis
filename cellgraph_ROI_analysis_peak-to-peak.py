# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 23:35:24 2021

@author: Fuyuki Asano
"""


import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import os
from scipy import signal
import csv
import math
import statistics
import scipy


#mac
#path = glob.glob('/Users/asanofuyuki/Desktop/python working/*.xls')
#win
path_list = glob.glob('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis\\data_folder\\*.csv')


column1 = ['index', 'time_stamp', 'time_h', 'time_min']
column2 = ['dosal' + str(x) for x in range(1,21,1)]
column3 = ['ventral' + str(x) for x in range(1,21,1)]
column = column1 + column2 + column3
name = np.array(['name']).reshape(1,1)
region = np.array(['region']).reshape(1,1)
data_per = np.array(['period', 'Amplitude', 'Damping factor', 'Phase angle', 'Non 0 Amplitude']).reshape(1,-1)
data_det = np.array([x for x in np.arange(0, 72, 0.5)]).reshape(1,-1)
data_sm = np.array([x for x in np.arange(0, 72, 0.5)]).reshape(1,-1)
emp = np.zeros(len(path_list)*40+1)
emp[:] = np.nan
emp = emp.reshape(-1,1)



for i in range(len(path_list)):
#win
    df = pd.read_csv(path_list[i], names=column, encoding='ANSI', skiprows=2, engine='python')

#mac
#    df = pd.read_csv(path_list[i], names=column, skipfooter=2, engine='python')
    
    file_name = path_list[i].split('\\')[-1]
    name_list = np.array([file_name for x in np.arange(0,40,1)]).reshape(-1,1)
    name = np.concatenate((name, name_list), axis=0)
    region_list1 = np.array(['dosal' + str(x) for x in range(1,21,1)]).reshape(-1,1)
    region_list2 = np.array(['ventral' + str(x) for x in range(1,21,1)]).reshape(-1,1)
    region = np.concatenate((region, region_list1, region_list2), axis=0)
    
    
    #find zt0
    time_points = pd.to_datetime(df.iloc[:,1])
    start_time = time_points[0]  
    next_day8 = start_time + np.timedelta64(1,'D')
    next_day8 = next_day8.replace(hour=8)
    next_day8 = next_day8.replace(minute=0)
    next_day8 = next_day8.replace(second=0)       
    td_zt8 = np.array(time_points - next_day8)
    onset = np.argwhere(td_zt8 > np.timedelta64(0, 's')).reshape(-1)[0]
    offset = onset + 144
    hr_zt = np.array(df.iloc[onset:offset, 2])
    hr_0 =  np.array(df.iloc[onset:offset, 2]) - np.array(df.iloc[onset:offset, 2])[0]

    
#data analysis 
    for ii in range(0,40,1):
        per = []
        dif = []
        sm = []
        
        luc = np.array(df.iloc[:,4+ii])
        
        #detrend 24h
        v = np.ones(48)/48
        conv = np.convolve(luc, v, mode='same')
        det_luc = (luc - conv)[24:192-24]
        narray = np.zeros(24)
        narray[:] = np.nan
        narray = narray.reshape(-1,1)
        df['det_luc'] = np.concatenate((narray, det_luc.reshape(-1,1), narray), axis=0)
        
        #Savitzky-Golay filter
        window = 25 #about 12h
        order = 3
        sm_luc = signal.savgol_filter(det_luc, window, order)
        df['sm_luc'] = np.concatenate((narray, sm_luc.reshape(-1,1), narray), axis=0)
        
        det_luc_zt = np.array(df.iloc[onset:offset, 44])
        sm_luc_zt = np.array(df.iloc[onset:offset, 45])
        
        
        #difference between zt0 and rec timestamp  
        zt0_onset = pd.to_datetime(df.iloc[onset,1])
        zt0_0sec = zt0_onset.replace(minute=0)
        zt0_0sec = zt0_0sec.replace(second=0)
        zt0_dif = zt0_onset - zt0_0sec
        zt0_dif = zt0_dif.seconds
        
        #calc amplitude and damping coefficient
        #amplitude = first peak
        #damping factor = first/third, logarithmic decrement
        #phase angle = duration to first peak
        plus_bool = np.array(sm_luc_zt >= 0)
        plus_val = sm_luc_zt[plus_bool]
        plus_idx = np.argwhere(plus_bool).reshape(-1)
        transition = (plus_idx != np.roll((plus_idx +1), 1))
        start_idx = plus_idx[transition]
        end_idx = plus_idx[np.roll(transition, -1)]
        
        
        if len(start_idx)>=3:
            if 0 in (start_idx - end_idx):
                non0amp = 'no'
                per.append(' ')
                per.append(' ')
                per.append(' ')
                per.append(' ')
                per.append(non0amp)
                
            
                #plot
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                line1 = ax.plot(hr_0, det_luc_zt, label='detrended')
                line2 = ax.plot(hr_0, sm_luc_zt, label='smoothing')
                ax.hlines(y=0, xmin=hr_0[0], xmax=hr_0[-1], color='k', linestyles='solid', linewidths=1)
                ax.set_xlim(left=hr_0[0], right=hr_0[-1])
                ax.legend()
                ax.set_xlabel('Time (hr)')
                ax.set_ylabel('RLU')
                ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
                ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
                fig.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis_fig\\' + file_name + '_' + region[ii+1,0] + '.png')
                plt.show()
            

            else:
                a1 = sm_luc_zt[start_idx[0]:end_idx[0]]
                a1_max_idx = np.round((start_idx[0] + end_idx[0])/2).astype('int64')
                a1_max = sm_luc_zt[a1_max_idx]
                
                a2 = sm_luc_zt[start_idx[1]:end_idx[1]]
                a2_max_idx = np.round((start_idx[1] + end_idx[1])/2).astype('int64')
                a2_max = sm_luc_zt[a2_max_idx]
                
                trough_idx = np.round(((end_idx[0]+1) + (start_idx[1]-1))/2).astype('int64')
                trough = sm_luc_zt[trough_idx]
                
                amp = a1_max - trough
                per_dif1 = ((hr_zt[a2_max_idx] - hr_zt[a1_max_idx]))
                ang = hr_zt[a1_max_idx] - hr_zt[0] + (zt0_dif/3600)
                
                delta1to2 = np.log(a1_max/a2_max)
                epsilon = delta1to2 / (2*np.pi)
                
                
                #When amplitude != 0, go simulation
                amp_sd_range = 6 #almost equal 3h
                amp_range_val = det_luc_zt[(a1_max_idx - amp_sd_range) : (a1_max_idx + amp_sd_range)]
                amp_mean = np.mean(amp_range_val)
                amp_sd = statistics.stdev(amp_range_val)
    
    
                if (amp_mean - 2*amp_sd) > 0:
                    non0amp = 'Yes'
                    per.append(per_dif1)
                    per.append(amp)
                    per.append(epsilon)
                    per.append(ang)
                    per.append(non0amp)
            
                else:
                    non0amp = 'no'
                    per.append(' ')
                    per.append(' ')
                    per.append(' ')
                    per.append(' ')
                    per.append(non0amp)
                    
                #plot  
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                line1 = ax.plot(hr_0, det_luc_zt, label='detrended')
                line2 = ax.plot(hr_0, sm_luc_zt, label='smoothing')
                ax.hlines(y=0, xmin=hr_0[0], xmax=hr_0[-1], color='k', linestyles='solid', linewidths=1)
                ax.set_xlim(left=hr_0[0], right=hr_0[-1])
                ax.legend()
                ax.set_xlabel('Time (hr)')
                ax.set_ylabel('RLU')
                ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
                ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
                fig.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis_fig\\' + file_name + '_' + region[ii+1,0] + '.png')
                plt.show()
            
            
        else:
            non0amp = 'no'
            per.append(' ')
            per.append(' ')
            per.append(' ')
            per.append(' ')
            per.append(non0amp)
        
            #plot
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(hr_0, det_luc_zt, label='detrended')
            line2 = ax.plot(hr_0, sm_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=hr_0[0], xmax=hr_0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=hr_0[0], right=hr_0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
            fig.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis_fig\\' + file_name + '_' + region[ii+1,0] + '.png')
            plt.show()
            
            
        data_per = np.concatenate((data_per, np.array(per).reshape(1,-1)), axis=0)
        data_det = np.concatenate((data_det, det_luc_zt.reshape(1,-1)), axis=0)    
        data_sm = np.concatenate((data_sm, sm_luc_zt.reshape(1,-1)), axis=0)


output = np.concatenate((name, region, data_per, data_det, emp, data_sm), axis=1)

with open ('cellgraph_ROI_analysis.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(output)

    

