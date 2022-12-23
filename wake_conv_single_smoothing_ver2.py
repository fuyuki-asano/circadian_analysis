# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 16:28:13 2021

@author: Fuyuki Asano
"""


import glob
import numpy as np
import pandas as pd
import sympy
import csv
import matplotlib.pyplot as plt
import sys
import os
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from scipy.stats import gaussian_kde


path = 'C:\\Users\\Sleepymouse\\Desktop\\experiment'
folderlist = sorted(os.listdir(path))
foldercount = len(folderlist)

p1 = 30
h = 540
p2 = int(h/p1)

exp_W = np.concatenate((np.array(['name', 'stage', 'peak_latency', 'peak_len', 'peak_Amp']).reshape((1,-1))\
                        , np.array([x for x in range(int(4320/p1))]).reshape((1,-1))), axis=1)
exp_NR = np.concatenate((np.array(['name', 'stage', 'peak_latency', 'peak_len', 'peak_Amp']).reshape((1,-1))\
                        , np.array([x for x in range(int(4320/p1))]).reshape((1,-1))), axis=1)
exp_R = np.concatenate((np.array(['name', 'stage', 'peak_latency', 'peak_len', 'peak_Amp']).reshape((1,-1))\
                        , np.array([x for x in range(int(4320/p1))]).reshape((1,-1))), axis=1)


plot_W = np.array([x for x in range(int(4320/p1))]).reshape(1, -1)
plot_NR = np.array([x for x in range(int(4320/p1))]).reshape(1, -1)
plot_R = np.array([x for x in range(int(4320/p1))]).reshape(1, -1)


for i in range(foldercount):
    filepath = sorted(glob.glob(path+'\\'+folderlist[i]+'\\*.xls'))
    filecount = len(filepath)

#    for ii in range(filecount-3):
    for ii in range(filecount-1): #for sleepy
#        ii = ii+1 #DD1
        
        name = filepath[ii].split('\\')[-1]
        data_W = np.array([name, 'W']).reshape(1, -1)
        data_NR = np.array([name, 'NR']).reshape(1, -1)
        data_R = np.array([name, 'R']).reshape(1, -1)
        

        df = pd.read_excel(filepath[ii], sheet_name=0)
        epoch = np.array(df.iloc[:,0])
        stage = np.array(df.iloc[:,1])
        
        
        wake = np.where(stage=='W', 1, 0)
        NR = np.where(stage=='NR', 1, 0)
        R = np.where(stage=='R', 1, 0)
        
        x = np.arange(len(wake))
        wake_idx = np.where(wake==1)
        wake_kde = gaussian_kde(wake_idx)
        y = wake_kde(x)
        

##################Wake####################################        
        wake_p1 = np.sum(wake.reshape(-1, p1), axis=1)
        wake_p1_norm = wake_p1/np.max(wake_p1)
#        wake_p1_th = np.where(wake_p1>=np.mean(wake_p1), np.max(wake_p1), 0)
        
        W_p2 = wake_p1[len(wake_p1)-int(p2/2):].tolist() + wake_p1.tolist() + wake[:int(p2/2)].tolist()
        W_p2 = np.array(W_p2)
        a_W1 = np.array([np.sum(W_p2[x:x+p2]) for x in range(int(4320/p1))])
#        a_W1_norm = a_W1/h
        a_W1_norm = a_W1/np.max(a_W1)
        a_W1_norm_th = np.where(a_W1_norm>=np.mean(a_W1_norm), 1, 0)
#        a_W1_norm_th = np.where(a_W1_norm>0.5, 1, 0)

        plot_W = np.concatenate((plot_W, a_W1_norm.reshape(1, -1)), axis=0)           

#peak latency detection
        temp_W = np.array([0] + a_W1_norm_th.tolist())   
        ar1_W = np.arange(len(temp_W))
        ar1_W[temp_W>0] = 0
        left_W = np.maximum.accumulate(ar1_W)
        
        temp_W_r = temp_W[::-1]
        ar2_W = np.arange(len(temp_W_r))
        ar2_W[temp_W_r>0] = 0
        acc_W_r = np.maximum.accumulate(ar2_W)
        right_W = len(temp_W_r) -1 -acc_W_r[::-1]
        
        W_ntbout = right_W - (left_W + 1)
        W_ntbout[temp_W==0] = 0
        W_ntbout = W_ntbout[1:]
        
        W_latency = np.argmax(W_ntbout)*30*20/3600
        W_len = np.sum(np.where(W_ntbout==np.max(W_ntbout), 1, 0))*p1*20/60
        W_amp = np.sum(wake_p1[W_ntbout==np.max(W_ntbout)]*20/60)/(W_len/60)  
        
        

##################NREMS####################################        
        NR_p1 = np.sum(NR.reshape(-1, p1), axis=1)
#        NR_p1_norm = NR_p1/np.max(NR_p1)
#        NR_p1_th = np.where(NR_p1>=np.mean(NR_p1), np.max(NR_p1), 0)
        
        NR_p2 = NR_p1[int(len(NR_p1)/2-p2/2):].tolist() + NR_p1[:int(len(NR_p1)/2+p2/2)].tolist()
        NR_p2 = np.array(NR_p2)
        a_NR1 = np.array([np.sum(NR_p2[x:x+p2]) for x in range(int(4320/p1))])
        a_NR1_norm = a_NR1/np.max(a_NR1)
        a_NR1_norm_th = np.where(a_NR1_norm>=np.mean(a_NR1_norm), 1, 0)
        
        plot_NR = np.concatenate((plot_NR, a_NR1_norm.reshape(1, -1)), axis=0)

#peak latency detection
        temp_NR = np.array([0] + a_NR1_norm_th.tolist())   
        ar1_NR = np.arange(len(temp_NR))
        ar1_NR[temp_NR>0] = 0
        left_NR = np.maximum.accumulate(ar1_NR)
        
        temp_NR_r = temp_NR[::-1]
        ar2_NR = np.arange(len(temp_NR_r))
        ar2_NR[temp_NR_r>0] = 0
        acc_NR_r = np.maximum.accumulate(ar2_NR)
        right_NR = len(temp_NR_r) -1 -acc_NR_r[::-1]
        
        NR_ntbout = right_NR - (left_NR + 1)
        NR_ntbout[temp_NR==0] = 0
        NR_ntbout = NR_ntbout[1:]
        
        NR_latency = (np.argmax(NR_ntbout) - (4320/p1)/2)*30*20/3600
        NR_len = np.sum(np.where(NR_ntbout==np.max(NR_ntbout), 1, 0))*p1*20/60
        NR_amp = np.sum(NR_p1[NR_ntbout==np.max(NR_ntbout)]*20/60)/(NR_len/60)      
        
        
##################REMS####################################        
        R_p1 = np.sum(R.reshape(-1, p1), axis=1)
#        R_p1_norm = R_p1/np.max(R_p1)
#        R_p1_th = np.where(R_p1>=np.mean(R_p1), np.max(R_p1), 0)
        
        R_p2 = R_p1[int(len(R_p1)/2-p2/2):].tolist() + R_p1[:int(len(R_p1)/2+p2/2)].tolist()
        R_p2 = np.array(R_p2)
        a_R1 = np.array([np.sum(R_p2[x:x+p2]) for x in range(int(4320/p1))])
        a_R1_norm = a_R1/np.max(a_R1)
        a_R1_norm_th = np.where(a_R1_norm>=np.mean(a_R1_norm), 1, 0)
        
        plot_R = np.concatenate((plot_R, a_R1_norm.reshape(1, -1)), axis=0)

#peak latency detection
        temp_R = np.array([0] + a_R1_norm_th.tolist())   
        ar1_R = np.arange(len(temp_R))
        ar1_R[temp_R>0] = 0
        left_R = np.maximum.accumulate(ar1_R)
        
        temp_R_r = temp_R[::-1]
        ar2_R = np.arange(len(temp_R_r))
        ar2_R[temp_R_r>0] = 0
        acc_R_r = np.maximum.accumulate(ar2_R)
        right_R = len(temp_R_r) -1 -acc_R_r[::-1]
        
        R_ntbout = right_R - (left_R + 1)
        R_ntbout[temp_R==0] = 0
        R_ntbout = R_ntbout[1:]
        
        R_latency = (np.argmax(R_ntbout) - (4320/p1)/2)*30*20/3600
        R_len = np.sum(np.where(R_ntbout==np.max(R_ntbout), 1, 0))*p1*20/60
        R_amp = np.sum(R_p1[R_ntbout==np.max(R_ntbout)]*20/60)/(R_len/60)     
        
        

        
        data_W = np.concatenate((data_W, W_latency.reshape(1,-1), W_len.reshape(1,-1), W_amp.reshape(1,-1)\
                                 , a_W1_norm.reshape(1,-1)), axis=1)
        data_NR = np.concatenate((data_NR, NR_latency.reshape(1,-1), NR_len.reshape(1,-1), NR_amp.reshape(1,-1)\
                                 , a_NR1_norm.reshape(1,-1)), axis=1)
        data_R = np.concatenate((data_R, R_latency.reshape(1,-1), R_len.reshape(1,-1), R_amp.reshape(1,-1)\
                                 , a_R1_norm.reshape(1,-1)), axis=1)
            
            
        exp_W = np.concatenate((exp_W, data_W), axis=0)
        exp_NR = np.concatenate((exp_NR, data_NR), axis=0)
        exp_R = np.concatenate((exp_R, data_R), axis=0)
        
        
        
        xx = [x for x in range(len(a_W1_norm))]
#        xx = np.linspace(0, 24, int(4320/p1)+1)[1:]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
#        ax.plot(xx, wake_p1_norm, label='10')
#        ax.plot(x, wake, label='w')
#        ax.plot(x, y, label='kde')
        ax.plot(xx, a_W1_norm, label='W1_norm')
        ax.plot(xx, a_W1_norm_th, label='W1_norm_th')
        ax.set_xlim(left=0, right=np.max(xx))
#        ax.set_ylim(bottom=-0.1, top=np.max(a_W1_norm)+0.1)
#        ax.set_ylim(bottom=0, top=np.max(y))
        ax.legend()
        plt.show()
        


exp = np.concatenate((exp_W, exp_NR, exp_R) ,axis=0)
exp_W2 = exp_W[1:, 2:,].astype('float64')
exp_W2_sort = exp_W2[np.argsort(exp_W2[:, 0])]
plot_W2 = exp_W2_sort[:, 3:]
               
fig = plt.figure()
fig.set_figheight(10)
fig.set_figwidth(10)
ax = fig.add_subplot(1,1,1)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
norm = Normalize(0, 1)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
#ax.imshow(plot_W[1:, :], cmap=cmap, vmin=0, vmax=1)
ax.imshow(plot_W2, cmap=cmap, vmin=0, vmax=1)
cbar = plt.colorbar(mappable, cax=cax)
#ax.scatter((exp_W2_sort[:, 0]*6), np.arange(exp_W2_sort.shape[0]), marker=".", color="k")
ax.set_xlim(left=0, right=143)
ax.set_ylim(bottom=25, top=0)
ax.xaxis.set_major_locator(ticker.MultipleLocator(12))
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
ax.set_xlabel('x10 min')
ax.set_ylabel('Mouse#')
plt.savefig('C:\\Users\\Sleepymouse\\Desktop\\experiment\\wake_peak_hm.pdf')
plt.show()

      

with open ('C:\\Users\\Sleepymouse\\Desktop\\experiment\\wake_10m_conv_3h_smoothing.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(exp)