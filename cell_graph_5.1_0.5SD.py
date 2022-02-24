# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 22:18:08 2021

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
import cv2
import tifffile
import statistics
import time
import matplotlib.colors as mc 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from scipy.ndimage import gaussian_filter
from skimage import img_as_float
from skimage.morphology import reconstruction


file = '12-22-20-1755-26'


#calc threshold
start = time.time()
P_min = 0
P_max = 48
P_interval = 0.5
ang_min = 0
ang_max = 48
e_min = -1
e_max = 1
amp_min = 0
amp_max = 10000
trough_min = 0
trough_max = 48


#export preparation
cutoff = np.zeros((512,512))
period_dif = np.zeros((512,512))
period_sim = np.zeros((512,512))
amplitude = np.zeros((512,512))
phase_angle_peak = np.zeros((512,512))
damping_factor = np.zeros((512,512))
non0amp = np.zeros((512,512))
phase_angle_trough = np.zeros((512,512))
data_stack = []
counter = 0

#damping oscillation simulation
def damping_eq_0(amp, epsilon, P, t, M):
    damp_amp = amp*(np.e**((-1*epsilon*2*np.pi*t)/(P*60)))
    osci = np.cos((np.sqrt(1-epsilon*epsilon)*2*np.pi*t)/(P*60))
    curve = damp_amp*osci + M
#    x = amp*math.pow(math.e, ((-1*epsilon*2*np.pi*t)/(P*60*60)))*np.cos((math.sqrt(1-epsilon*epsilon)*2*np.pi*t)/(P*60*60))
    return curve


#16F to 8U convertor
def conv16Fto8U(data, dmin, dmax):
    if dmin == dmax:
        return np.zeros_like(data)

    scaled = 255* (data - dmin) / (dmax - dmin)
    scaled = scaled.round()
    scaled = scaled.astype(np.uint8)

    return scaled

def normalizer0_1(data, dmin, dmax):
    norm = (data - dmin) / (dmax - dmin)
    return norm
 


#gray scale, 16bit image
#win
path_list = sorted(glob.glob('F:\\circadian_experiment\\python_scripts\\cellgraph\\' + file + '\\TIF\\*.TIF'))
csv_path = sorted(glob.glob('F:\\circadian_experiment\\python_scripts\\cellgraph\\'+ file +'\\*.csv'))

#mac
#path_list = sorted(glob.glob('/Users/asanofuyuki/Desktop/python_scripts/12-15-20-1714-31/TIF/*.TIF'))
#csv_path = sorted(glob.glob('/Users/asanofuyuki/Desktop/python_scripts/12-15-20-1714-31/*.csv'))

if len(csv_path) != 1:
    raise ValueError('Check the number of csv files')


#image stacking
i = 0
stacked_im = cv2.imread(path_list[i], cv2.IMREAD_ANYDEPTH) #gray_scale, 16bit
stacked_im_3d = stacked_im.reshape(1, 512, 512)

while True:
    i = i+1
    if i == len(path_list):
        break
    elif i < len(path_list):
        im_gray16 = cv2.imread(path_list[i], cv2.IMREAD_ANYDEPTH)
        im_gray16_3d = im_gray16.reshape(1, 512, 512)
        stacked_im_3d = np.concatenate((stacked_im_3d, im_gray16_3d), axis=0)
        
        
#threshold calc, mean_im, mean + 1SD        
mean_im = stacked_im_3d.mean(axis=0)
threshold = mean_im.mean() + np.std(mean_im)*0.5

ret, img_binary = cv2.threshold(mean_im, threshold, 2**16, cv2.THRESH_BINARY)


mean_im_U6 = conv16Fto8U(mean_im, np.min(mean_im), np.max(mean_im))
cv2.imwrite('mean.png', mean_im_U6)
cv2.imwrite('binary.png', img_binary)
 

#import dataframe       
interval = 30*60 #sec
df = pd.read_csv(csv_path[0], skiprows=[0], header=[0], encoding='SHIFT-JIS')
df = df.drop(df.columns[0], axis=1)
#minutes = np.array(df.iloc[:,2])
df2 = df.iloc[:,0:3]


#find zt0
time_points = pd.to_datetime(df.iloc[:,0])
start_time = time_points[0]  
next_day8 = start_time + np.timedelta64(1,'D')
next_day8 = next_day8.replace(hour=8)
next_day8 = next_day8.replace(minute=0)
next_day8 = next_day8.replace(second=0)       
td_zt8 = np.array(time_points - next_day8)
onset = np.argwhere(td_zt8 > np.timedelta64(0, 's')).reshape(-1)[0]
offset = onset + 144
hr_zt = np.array(df.iloc[onset:offset, 1])
hr_0 =  np.array(df.iloc[onset:offset, 1]) - np.array(df.iloc[onset:offset, 1])[0]


narray = np.zeros(24)
narray[:] = np.nan
narray = narray.reshape(-1,1)


#difference between zt0 and rec timestamp      
zt0_onset = pd.to_datetime(df.iloc[onset,0])
zt0_0sec = zt0_onset.replace(minute=0)
zt0_0sec = zt0_0sec.replace(second=0)
zt0_dif = zt0_onset - zt0_0sec
zt0_dif = zt0_dif.seconds


binary_pass = 0


#start calc
for x in range(0,512):
    for y in range(0,512):
        luc = stacked_im_3d[:, x, y]
        

#threshhold    
        if img_binary[x,y] == 0:
            binary_pass = binary_pass + 1
            pass
        
        else:
            #de-trend 24h
            v = np.ones(48)/48
            conv = np.convolve(luc, v, mode='same')
            det_luc = (luc - conv)[24:192-24]
            df2['det_luc'] = np.concatenate((narray, det_luc.reshape(-1,1), narray), axis=0)
            

            #Savitzky-Golay filter      
            window = 25 #about 24h
            order = 3
            sm_luc = signal.savgol_filter(det_luc, window, order)
            df2['sm_luc'] = np.concatenate((narray, sm_luc.reshape(-1,1), narray), axis=0)
            M = np.mean(sm_luc)
#            smoothed_subM = sm_luc - M
            
        
            #data from ZT0
            det_luc_zt = np.array(df2.iloc[onset:offset, 3])
            sm_luc_zt = np.array(df2.iloc[onset:offset, 4])
#            smoothed_subM_zt0 = smoothed_subM[onset:]
        

            #amplitude check preparation
            plus_bool = np.array(sm_luc_zt >= 0)
            plus_val = sm_luc_zt[plus_bool]
            plus_idx = np.argwhere(plus_bool).reshape(-1)
            transition = (plus_idx != np.roll((plus_idx +1), 1))
            start_idx = plus_idx[transition]
            end_idx = plus_idx[np.roll(transition, -1)]
            

#if the number of positive value is less than 2, pass following steps and go next y
            if len(start_idx)>=3:
                if 0 in (start_idx - end_idx):
                    #plot
                    #fig = plt.figure()
                    #ax = fig.add_subplot(1,1,1)
                    #line1 = ax.plot(hr_0, det_luc_zt, label='detrended')
                    #line2 = ax.plot(hr_0, sm_luc_zt, label='smoothing')
                    #ax.hlines(y=0, xmin=hr_0[0], xmax=hr_0[-1], color='k', linestyles='solid', linewidths=1)
                    #ax.set_xlim(left=hr_0[0], right=hr_0[-1])
                    #ax.legend()
                    #ax.set_xlabel('Time (hr)')
                    #ax.set_ylabel('RLU')
                    #ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
                    #ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
                    #plt.show()
                    pass
                    
                
                else:
                    a1 = sm_luc_zt[start_idx[0]:end_idx[0]]
                    a1_max_idx = np.round((start_idx[0] + end_idx[0])/2).astype('int64')
                    a1_max = sm_luc_zt[a1_max_idx]
                
                    a2 = sm_luc_zt[start_idx[1]:end_idx[1]]
                    a2_max_idx = np.round((start_idx[1] + end_idx[1])/2).astype('int64')
                    a2_max = sm_luc_zt[a2_max_idx]
                
                    trough_idx = np.round(((end_idx[0]+1) + (start_idx[1]-1))/2).astype('int64')
                    trough = sm_luc_zt[trough_idx]
                    
                    if a1_max_idx < 6:
                        pass
                        
                    else:
                        amp = a1_max - trough
                        per_dif1 = ((hr_zt[a2_max_idx] - hr_zt[a1_max_idx]))
                        ang = hr_zt[a1_max_idx] - hr_zt[0] + (zt0_dif/3600)
                
                        delta1to2 = np.log(a1_max/a2_max)
                        epsilon = delta1to2 / (2*np.pi)
                        if epsilon < -1:
                            print(x,y)
                            pass
                        else:
                            #When amplitude != 0, go simulation
                            amp_sd_range = 6 #almost equal 3h
                            amp_range_val = det_luc_zt[(a1_max_idx - amp_sd_range) : (a1_max_idx + amp_sd_range)]
                            amp_mean = np.mean(amp_range_val)
                            amp_sd = statistics.stdev(amp_range_val)
                    
                    
                            if (amp_mean - 2*amp_sd) > 0:
                                non0amp[x,y] = pow(2,16)
                                phase_angle_peak[x,y] = ang
                                period_dif[x,y] = per_dif1
                                amplitude[x,y] = amp
                                damping_factor[x,y] = epsilon
                        
            
                            else:
                                pass
                    
                            #plot  
                            #fig = plt.figure()
                            #ax = fig.add_subplot(1,1,1)
                            #line1 = ax.plot(hr_0, det_luc_zt, label='detrended')
                            #line2 = ax.plot(hr_0, sm_luc_zt, label='smoothing')
                            #ax.hlines(y=0, xmin=hr_0[0], xmax=hr_0[-1], color='k', linestyles='solid', linewidths=1)
                            #ax.set_xlim(left=hr_0[0], right=hr_0[-1])
                            #ax.legend()
                            #ax.set_xlabel('Time (hr)')
                            #ax.set_ylabel('RLU')
                            #ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
                            #ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
                            #fig.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis_fig\\' + file_name + '_' + region[ii+1,0] + '.png')
                            #plt.show()
    
    
            else:
                #plot  
                #fig = plt.figure()
                #ax = fig.add_subplot(1,1,1)
                #line1 = ax.plot(hr_0, det_luc_zt, label='detrended')
                #line2 = ax.plot(hr_0, sm_luc_zt, label='smoothing')
                #ax.hlines(y=0, xmin=hr_0[0], xmax=hr_0[-1], color='k', linestyles='solid', linewidths=1)
                #ax.set_xlim(left=hr_0[0], right=hr_0[-1])
                #ax.legend()
                #ax.set_xlabel('Time (hr)')
                #ax.set_ylabel('RLU')
                #ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
                #ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=1)
                #fig.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis_fig\\' + file_name + '_' + region[ii+1,0] + '.png')
                #plt.show()
                pass
            
                    
                    
#save psuedo color fig, per_sim
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.axis('off')
#norm_per = Normalize(0, P_max)
#cmap = get_cmap('jet')
#mappable = ScalarMappable(cmap=cmap, norm=norm_per)
#mappable._A = []
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="3%", pad="3%")
#ax.imshow(norm_per(period_sim), cmap=cmap)
#cbar = plt.colorbar(mappable, cax=cax)
#ticks = np.linspace(0, P_max, 5)
#cbar.set_ticks(ticks)
#plt.savefig('F:\\circadian_experiment\\python_scripts\\'+ file +'\\TIF\\pc_period_sim.png')
#plt.show()


#save psuedo color fig, per_dif
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axis('off')
norm_per = norm_e = mc.DivergingNorm(vcenter=24, vmin=16, vmax=32)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm_per)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(norm_per(period_dif), cmap=cmap)
cbar = plt.colorbar(mappable, cax=cax)
ticks = np.linspace(16, 32, 5)
cbar.set_ticks(ticks)
plt.savefig('pc_period_dif_0.5SD.png')
plt.show()


#save psuedo color fig, ang
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axis('off')
norm_ang = norm_e = mc.DivergingNorm(vcenter=12, vmin=0, vmax=24)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm_ang)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(norm_ang(phase_angle_peak), cmap=cmap)
cbar = plt.colorbar(mappable, cax=cax)
ticks = np.linspace(0, 24, 7)
cbar.set_ticks(ticks)
plt.savefig('pc_ang_0.5SD.png')
plt.show()


#save psuedo color fig, amp
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axis('off')
norm_amp = norm_e = mc.DivergingNorm(vcenter=5000, vmin=0, vmax=10000)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm_amp)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(norm_amp(amplitude), cmap=cmap)
cbar = plt.colorbar(mappable, cax=cax)
ticks = np.linspace(0, 10000, 5)
cbar.set_ticks(ticks)
plt.savefig('pc_amp_0.5SD.png')
plt.show()


#save psuedo color fig, damp
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axis('off')
norm_e = mc.DivergingNorm(vcenter=0.05, vmin=0, vmax=0.1)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm_e)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(norm_e(damping_factor), cmap=cmap)
cbar = plt.colorbar(mappable, cax=cax)
ticks = np.linspace(0, 0.1, 5)
cbar.set_ticks(ticks)
plt.savefig('pc_damp_0.5SD.png')
plt.show()




#data export
#if onset >= 24:
#    data_stack = np.array(data_stack)
#    data_stack = data_stack.reshape(counter,-1)
#    with open ('F:\\circadian_experiment\\python_scripts\\'+ file +'\\pixel_amp.csv', 'w', newline = '') as f:
#        writer = csv.writer(f)
#        writer.writerows(data_stack)


end = time.time()
print(end - start)  