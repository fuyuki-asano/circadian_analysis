# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 12:39:31 2021

@author: Fuyuki Asano
"""

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
import sys


path_list = glob.glob('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis\\for_heatmap\\for_heatmap.xlsx')
df = pd.read_excel(path_list[0],  engine='openpyxl')
n = 60

data = np.array(df.iloc[:, 3:])
sm_max = np.nanmax(data)
sm_min = np.nanmin(data)
data_n = (data - sm_min) / (sm_max - sm_min)


fd_n = np.array(data_n[0:60, :])
fv_n = np.array(data_n[60:120, :])
vd_n = np.array(data_n[120:180, :])
vv_n = np.array(data_n[180:240, :])


fd_n_max = np.nanmax(fd_n, axis=1)
fd_n_idx = np.argsort(fd_n_max)
fd_n_s = np.zeros(144)
fd_n_s[:] = np.nan
fd_n_s = fd_n_s.reshape(1,-1)

for i in range(60):
    fd_n_s = np.concatenate((fd_n_s, np.array(fd_n[fd_n_idx[59-i],:]).reshape(1,-1)), axis=0)
    
    
fv_n_max = np.nanmax(fv_n, axis=1)
fv_n_idx = np.argsort(fv_n_max)
fv_n_s = np.zeros(144)
fv_n_s[:] = np.nan
fv_n_s = fv_n_s.reshape(1,-1)

for i in range(60):
    fv_n_s = np.concatenate((fv_n_s, np.array(fv_n[fv_n_idx[59-i],:]).reshape(1,-1)), axis=0)


vd_n_max = np.nanmax(vd_n, axis=1)
vd_n_idx = np.argsort(vd_n_max)
vd_n_s = np.zeros(144)
vd_n_s[:] = np.nan
vd_n_s = vd_n_s.reshape(1,-1)

for i in range(60):
    vd_n_s = np.concatenate((vd_n_s, np.array(vd_n[vd_n_idx[59-i],:]).reshape(1,-1)), axis=0)


vv_n_max = np.nanmax(vv_n, axis=1)
vv_n_idx = np.argsort(vv_n_max)
vv_n_s = np.zeros(144)
vv_n_s[:] = np.nan
vv_n_s = vv_n_s.reshape(1,-1)

for i in range(60):
    vv_n_s = np.concatenate((vv_n_s, np.array(vv_n[vv_n_idx[59-i],:]).reshape(1,-1)), axis=0)

#save psuedo color fig
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
norm = Normalize(0, 1)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(fd_n_s, cmap=cmap, vmin=0, vmax=1)
cbar = plt.colorbar(mappable, cax=cax)
ax.set_xlim(left=0, right=144)
ax.set_ylim(bottom=60, top=0)
ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.set_xlabel('Time (hr)')
ax.set_ylabel('ROI')
plt.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis\\for_heatmap\\Sik3-ex3-flox_ROI_dosal_heatmap_sorted_nm.png')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
norm = Normalize(0, 1)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(fv_n_s, cmap=cmap, vmin=0, vmax=1)
cbar = plt.colorbar(mappable, cax=cax)
ax.set_xlim(left=0, right=144)
ax.set_ylim(bottom=60, top=0)
ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.set_xlabel('Time (hr)')
ax.set_ylabel('ROI')
plt.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis\\for_heatmap\\Sik3-ex3-flox_ROI_ventral_heatmap_sorted_nm.png')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
norm = Normalize(0, 1)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(vd_n_s, cmap=cmap, vmin=0, vmax=1)
cbar = plt.colorbar(mappable, cax=cax)
ax.set_xlim(left=0, right=144)
ax.set_ylim(bottom=60, top=0)
ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.set_xlabel('Time (hr)')
ax.set_ylabel('ROI')
plt.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis\\for_heatmap\\Vgat-Cre;Sik3-ex3-flox_ROI_dosal_heatmap_sorted_nm.png')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
norm = Normalize(0, 1)
cmap = get_cmap('jet')
mappable = ScalarMappable(cmap=cmap, norm=norm)
mappable._A = []
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad="3%")
ax.imshow(vv_n_s, cmap=cmap, vmin=0, vmax=1)
cbar = plt.colorbar(mappable, cax=cax)
ax.set_xlim(left=0, right=144)
ax.set_ylim(bottom=60, top=0)
ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.set_xlabel('Time (hr)')
ax.set_ylabel('ROI')
plt.savefig('F:\\circadian_experiment\\python_scripts\\cellgraph\\ROI_analysis\\for_heatmap\\Vgat-Cre;Sik3-ex3-flox_ROI_ventral_heatmap_sorted_nm.png')
plt.show()






