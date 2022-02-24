# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 18:47:46 2021

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
path_list = glob.glob('F:\\circadian_experiment\\python_scripts\\kronos\\Kronos_data_file\\*.csv')
#path = 'F:\\circadian_experiment\\python_scripts\\20200923-1922\\'
#folderlist = sorted(os.listdir(path))
#foldercount = len(folderlist)


column = ['time', 'p_sec', 'p_hour', 'F0', 'temp', 'CO2']
output1 = np.array(['Sample', 'date', 'start_time', 'Period_1'\
                   , 'Period_2', 'period_3', 'period_4', 'Amplitude'\
                       , 'Damping factor', 'Phase angle', 'Non 0 Amplitude']).reshape(1,-1)
    
h1 = np.array([n for n in range(1100)]).reshape(1,-1)
#h1 = h1.reshape(1,-1)
space = np.array(['']).reshape(1,-1)
output2 = np.concatenate((h1, space, h1), axis=1) 


for i in range(len(path_list)):
#なぜか最後に2行追加される(3行の場合もある?)から、skipfooter=2
#win
    df = pd.read_csv(path_list[i], names=column, encoding='ANSI', skipfooter=2, engine='python') 
#mac
#    df = pd.read_csv(path_list[i], names=column, skipfooter=2, engine='python')

    data = []
    name = path_list[i].split('\\')[-1]
    date = path_list[i].split('\\')[-2]
    data.append(name)
    data.append(date)
    
    
    interval = int(df.iloc[5,1])
    df = df.drop(range(0,12,1), axis=0) #should be modified depending on pause times
    df = df.drop(df.columns[[4,5]], axis=1)
    df = df.reset_index(drop=True)
    
    data.append(df.iloc[0,0])

    time_point = np.array(pd.to_datetime(df.iloc[:,0]))
    lost_idx = np.argwhere(np.isnat(time_point))


#time nan padding   
    for ii in np.reshape(lost_idx, (-1)):        
        time_point[ii] = time_point[ii-1] + np.timedelta64(interval, 's')
        
    df.iloc[:, 0] = time_point


#p_sec nan padding    
    p_sec = df.iloc[:, 1].astype('float64')
    p_sec_interp = p_sec.interpolate()
#    p_sec_interp2 = p_sec.interpolate('spline', order=1)
    p_sec_interp = p_sec_interp.round()
    df.iloc[:, 1] = p_sec_interp


#p_hour nan padding    
    p_hour = df.iloc[:, 2].astype('float64')
    p_hour_interp = p_hour.interpolate()
#    p_hour_interp2 = p_hour.interpolate('spline', order=1)
    df.iloc[:, 2] = p_hour_interp
    
    
#F0 nan padding    
    df.iloc[:,3] = df.iloc[:,3].astype('float64')
    f0 = df.iloc[:, 3]
    f0_interp = f0.interpolate()
#    f0_interp2 = f0.interpolate('spline', order=1)
    df.iloc[:,3] = f0_interp


#detrend 12h
#    luc = np.array(df.iloc[:,3])    
#    first12h = df.iloc[0,2] + 12
#    sub12h = np.argwhere(np.array(df.iloc[:,2]) > first12h).reshape(-1)
#    
#    num = len(df.iloc[:,2]) - len(sub12h) +1
#    w = np.ones(num)/num
#    ma = np.convolve(luc, w, mode='valid')
#    det_luc = luc[sub12h] - ma
#    narray = np.zeros(num)
#    narray[:] = np.nan
#    df['det_luc'] = np.concatenate((narray.reshape(-1,1), det_luc.reshape(-1,1)), axis=0)
    
    
#detrend 24h
    luc = np.array(df.iloc[:,3])
    sec = np.array(df.iloc[:,1])   
    first12h = df.iloc[0,1] + 12*60*60
    last12h = df.iloc[-1,1] - 12*60*60
    sub12hrs = np.argwhere((sec > first12h) * (sec < last12h)).reshape(-1)
    ma =[]
    for j in sub12hrs:
        st =sec[j] - 12*60*60 - (interval-1)
        se = np.argwhere(sec > st).reshape(-1)[0]
        et = sec[j] + 12*60*60 + (interval-1)
        ee = np.argwhere(sec < et).reshape(-1)[-1]
        mean = np.mean(luc[se:ee])
        ma.append(mean)
        
    det_luc = luc[sub12hrs] - np.array(ma)
    narray = np.zeros(sub12hrs[0])
    narray[:] = np.nan
    narray = narray.reshape(-1,1)
    df['det_luc'] = np.concatenate((narray, det_luc.reshape(-1,1), narray), axis=0)


#Savitzky-Golay filter, moving average smoothing <- 12h Savitzky-Golay filter has enough power, skipped moving average
    window = 103 #about 12h
    order = 3
    sav = signal.savgol_filter(det_luc, window, order)
#    conv_point = 10
#    v = np.ones(conv_point)/conv_point
#    conv = np.convolve(sav, v, mode='same')
    df['filt_luc'] = np.concatenate((narray, sav.reshape(-1,1), narray), axis=0)
#    df['filt_luc'] = np.concatenate((narray, conv.reshape(-1,1), narray), axis=0)
    

#data from ZT0 (8:00am)
    time_stamp = df.iloc[:,0]
    later12h = time_stamp[0] + np.timedelta64(12,'h')    
    next_day8 = time_stamp[0] + np.timedelta64(1,'D')
    next_day8 = next_day8.replace(hour=8)
    next_day8 = next_day8.replace(minute=0)
    next_day8 = next_day8.replace(second=0)       
    td = np.array(time_stamp - next_day8)
    onset = np.argwhere(td > np.timedelta64(0, 's')).reshape(-1)[0]
    offset = onset + 1100
        

#difference between zt0 and rec timestamp    
    zt0_onset = pd.to_datetime(df.iloc[onset,0])
    zt0_0sec = zt0_onset.replace(minute=0)
    zt0_0sec = zt0_0sec.replace(second=0)
    zt0_dif = zt0_onset - zt0_0sec
    zt0_dif = zt0_dif.seconds
    
    
    det_luc_zt = np.array(df.iloc[onset:offset, 4])    
    filt_luc_zt = np.array(df.iloc[onset:offset, 5])
    sec_zt = np.array(df.iloc[onset:offset, 1])
    hr_zt = np.array(df.iloc[onset:offset, 2])
    set0 = hr_zt - hr_zt[0]
    
    luc_array = np.concatenate((det_luc_zt.reshape(1,-1), space, filt_luc_zt.reshape(1,-1)), axis=1)
    output2 = np.concatenate((output2, luc_array), axis=0)


#hilbert transform for envelope calcuration <- Not good
#    hil = np.abs(signal.hilbert(filt_luc_zt))
    

#calc amplitude and damping coefficient
#amplitude = first peak
#damping factor = first/third, logarithmic decrement
#phase angle = duration to first peak
    plus_bool = np.array(filt_luc_zt >= 0)
    plus_val = filt_luc_zt[plus_bool]
    plus_idx = np.argwhere(plus_bool).reshape(-1)
    transition = (plus_idx != np.roll((plus_idx +1), 1))
    start_idx = plus_idx[transition]
    end_idx = plus_idx[np.roll(transition, -1)]



#######################################################    
    if len(start_idx)>=5:
        if 0 in (start_idx - end_idx):
            non0amp = 'no'
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()
            

        else:
            a1 = filt_luc_zt[start_idx[0]:end_idx[0]]
            a1_max_idx = np.round((start_idx[0] + end_idx[0])/2).astype('int64')
            a1_max = filt_luc_zt[a1_max_idx]
        
            a2 = filt_luc_zt[start_idx[1]:end_idx[1]]
            a2_max_idx = np.round((start_idx[1] + end_idx[1])/2).astype('int64')
            a2_max = filt_luc_zt[a2_max_idx]
        
            a3 = filt_luc_zt[start_idx[2]:end_idx[2]]
            a3_max_idx = np.round((start_idx[2] + end_idx[2])/2).astype('int64')
            a3_max = filt_luc_zt[a3_max_idx]
        
            a4 = filt_luc_zt[start_idx[3]:end_idx[3]]
            a4_max_idx = np.round((start_idx[3] + end_idx[3])/2).astype('int64')
            a4_max = filt_luc_zt[a4_max_idx]
        
            a5 = filt_luc_zt[start_idx[4]:end_idx[4]]
            a5_max_idx = np.round((start_idx[4] + end_idx[4])/2).astype('int64')
            a5_max = filt_luc_zt[a5_max_idx]

            
            trough_idx = np.round(((end_idx[0]+1) + (start_idx[1]-1))/2).astype('int64')
            trough = filt_luc_zt[trough_idx]
        
            amp = a1_max - trough
            per_dif1 = ((sec_zt[a2_max_idx] - sec_zt[a1_max_idx])) / (60*60)
            per_dif2 = ((sec_zt[a3_max_idx] - sec_zt[a2_max_idx])) / (60*60)
            per_dif3 = ((sec_zt[a4_max_idx] - sec_zt[a3_max_idx])) / (60*60)
            per_dif4 = ((sec_zt[a5_max_idx] - sec_zt[a4_max_idx])) / (60*60)
            ang = (sec_zt[a1_max_idx] - sec_zt[0] + zt0_dif) / (60*60)
            
            delta1to3 = (1/2)*np.log(a1_max/a3_max)
            epsilon = delta1to3 / (2*np.pi)
        
        
            #When amplitude != 0, go simulation
            amp_sd_range = 26 #almost equal 3h
            amp_range_val = det_luc_zt[(a1_max_idx - amp_sd_range) : (a1_max_idx + amp_sd_range)]
            amp_mean = np.mean(amp_range_val)
            amp_sd = statistics.stdev(amp_range_val)
    
    
            if (amp_mean - 2*amp_sd) > 0:
                non0amp = 'Yes'
                data.append(per_dif1)
                data.append(per_dif2)
                data.append(per_dif3)
                data.append(per_dif4)
                data.append(amp)
                data.append(epsilon)
                data.append(ang)
                data.append(non0amp)
            
            else:
                non0amp = 'no'
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()

            

            
#######################################################    
    elif len(start_idx)>=4:
        if 0 in (start_idx - end_idx):
            non0amp = 'no'
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()
            

        else:
            a1 = filt_luc_zt[start_idx[0]:end_idx[0]]
            a1_max_idx = np.round((start_idx[0] + end_idx[0])/2).astype('int64')
            a1_max = filt_luc_zt[a1_max_idx]
        
            a2 = filt_luc_zt[start_idx[1]:end_idx[1]]
            a2_max_idx = np.round((start_idx[1] + end_idx[1])/2).astype('int64')
            a2_max = filt_luc_zt[a2_max_idx]
        
            a3 = filt_luc_zt[start_idx[2]:end_idx[2]]
            a3_max_idx = np.round((start_idx[2] + end_idx[2])/2).astype('int64')
            a3_max = filt_luc_zt[a3_max_idx]
        
            a4 = filt_luc_zt[start_idx[3]:end_idx[3]]
            a4_max_idx = np.round((start_idx[3] + end_idx[3])/2).astype('int64')
            a4_max = filt_luc_zt[a4_max_idx]
            
            
            trough_idx = np.round(((end_idx[0]+1) + (start_idx[1]-1))/2).astype('int64')
            trough = filt_luc_zt[trough_idx]
        
            amp = a1_max - trough            
            per_dif1 = ((sec_zt[a2_max_idx] - sec_zt[a1_max_idx])) / (60*60)
            per_dif2 = ((sec_zt[a3_max_idx] - sec_zt[a2_max_idx])) / (60*60)
            per_dif3 = ((sec_zt[a4_max_idx] - sec_zt[a3_max_idx])) / (60*60)
            ang = (sec_zt[a1_max_idx] - sec_zt[0] + zt0_dif) / (60*60)
            
            delta1to3 = (1/2)*np.log(a1_max/a3_max)
            epsilon = delta1to3 / (2*np.pi)
        
        
            #When amplitude != 0, go simulation
            amp_sd_range = 26 #almost equal 3h
            amp_range_val = det_luc_zt[(a1_max_idx - amp_sd_range) : (a1_max_idx + amp_sd_range)]
            amp_mean = np.mean(amp_range_val)
            amp_sd = statistics.stdev(amp_range_val)
    
    
            if (amp_mean - 2*amp_sd) > 0:
                non0amp = 'Yes'
                data.append(per_dif1)
                data.append(per_dif2)
                data.append(per_dif3)
                data.append(' ')
                data.append(amp)
                data.append(epsilon)
                data.append(ang)
                data.append(non0amp)
            
            else:
                non0amp = 'no'
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()



#######################################################    
    elif len(start_idx)>=3:
        if 0 in (start_idx - end_idx):
            non0amp = 'no'
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()
            

        else:
            a1 = filt_luc_zt[start_idx[0]:end_idx[0]]
            a1_max_idx = np.round((start_idx[0] + end_idx[0])/2).astype('int64')
            a1_max = filt_luc_zt[a1_max_idx]
        
            a2 = filt_luc_zt[start_idx[1]:end_idx[1]]
            a2_max_idx = np.round((start_idx[1] + end_idx[1])/2).astype('int64')
            a2_max = filt_luc_zt[a2_max_idx]
        
            a3 = filt_luc_zt[start_idx[2]:end_idx[2]]
            a3_max_idx = np.round((start_idx[2] + end_idx[2])/2).astype('int64')
            a3_max = filt_luc_zt[a3_max_idx]
        
        
            trough_idx = np.round(((end_idx[0]+1) + (start_idx[1]-1))/2).astype('int64')
            trough = filt_luc_zt[trough_idx]
        
            amp = a1_max - trough
            per_dif1 = ((sec_zt[a2_max_idx] - sec_zt[a1_max_idx])) / (60*60)
            per_dif2 = ((sec_zt[a3_max_idx] - sec_zt[a2_max_idx])) / (60*60)
            ang = (sec_zt[a1_max_idx] - sec_zt[0] + zt0_dif) / (60*60)
            
            delta1to3 = (1/2)*np.log(a1_max/a3_max)
            epsilon = delta1to3 / (2*np.pi)
        
        
            #When amplitude != 0, go simulation
            amp_sd_range = 26 #almost equal 3h
            amp_range_val = det_luc_zt[(a1_max_idx - amp_sd_range) : (a1_max_idx + amp_sd_range)]
            amp_mean = np.mean(amp_range_val)
            amp_sd = statistics.stdev(amp_range_val)
    
    
            if (amp_mean - 2*amp_sd) > 0:
                non0amp = 'Yes'
                data.append(per_dif1)
                data.append(per_dif2)
                data.append(' ')
                data.append(' ')
                data.append(amp)
                data.append(epsilon)
                data.append(ang)
                data.append(non0amp)
            
            else:
                non0amp = 'no'
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()
            
            
            
#######################################################    
    elif len(start_idx)>=2:
        if 0 in (start_idx - end_idx):
            non0amp = 'no'
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(' ')
            data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()
            

        else:
            a1 = filt_luc_zt[start_idx[0]:end_idx[0]]
            a1_max_idx = np.round((start_idx[0] + end_idx[0])/2).astype('int64')
            a1_max = filt_luc_zt[a1_max_idx]
        
            a2 = filt_luc_zt[start_idx[1]:end_idx[1]]
            a2_max_idx = np.round((start_idx[1] + end_idx[1])/2).astype('int64')
            a2_max = filt_luc_zt[a2_max_idx]        
        
            trough_idx = np.round(((end_idx[0]+1) + (start_idx[1]-1))/2).astype('int64')
            trough = filt_luc_zt[trough_idx]
        
            amp = a1_max - trough
            per_dif1 = ((sec_zt[a2_max_idx] - sec_zt[a1_max_idx])) / (60*60)
            ang = (sec_zt[a1_max_idx] - sec_zt[0] + zt0_dif) / (60*60)
            
            delta1to2 = np.log(a1_max/a2_max)
            epsilon = delta1to2 / (2*np.pi)
        
        
            #When amplitude != 0, go simulation
            amp_sd_range = 26 #almost equal 3h
            amp_range_val = det_luc_zt[(a1_max_idx - amp_sd_range) : (a1_max_idx + amp_sd_range)]
            amp_mean = np.mean(amp_range_val)
            amp_sd = statistics.stdev(amp_range_val)
    
    
            if (amp_mean - 2*amp_sd) > 0:
                non0amp = 'Yes'
                data.append(per_dif1)
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(amp)
                data.append(epsilon)
                data.append(ang)
                data.append(non0amp)
            
            else:
                non0amp = 'no'
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(' ')
                data.append(non0amp)
            
            #plot
#            raw = np.array(df.iloc[onset:offset, 3])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            line1 = ax.plot(set0, det_luc_zt, label='detrended')
            line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
            ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
            ax.set_xlim(left=set0[0], right=set0[-1])
            ax.legend()
            ax.set_xlabel('Time (hr)')
            ax.set_ylabel('RLU')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
            ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#            fig.savefig(str(name) + '.png')
            plt.show()
            
            
            
#######################################################    
    else:
        non0amp = 'no'
        data.append(' ')
        data.append(' ')
        data.append(' ')
        data.append(' ')
        data.append(' ')
        data.append(' ')
        data.append(' ')
        data.append(non0amp)
        
        #plot
#        raw = np.array(df.iloc[onset:offset, 3])
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        line1 = ax.plot(set0, det_luc_zt, label='detrended')
        line2 = ax.plot(set0, filt_luc_zt, label='smoothing')
        ax.hlines(y=0, xmin=set0[0], xmax=set0[-1], color='k', linestyles='solid', linewidths=1)
        ax.set_xlim(left=set0[0], right=set0[-1])
        ax.legend()
        ax.set_xlabel('Time (hr)')
        ax.set_ylabel('RLU')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
        ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
#        fig.savefig(str(name) + '.png')
        plt.show()
            

#########################################################
            
    output1 = np.concatenate((output1, np.array(data).reshape(1,-1)), axis=0)


output = np.concatenate((output1, output2), axis=1)

with open ('per2_luc_summary.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(output)

