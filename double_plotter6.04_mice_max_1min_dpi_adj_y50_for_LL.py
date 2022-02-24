
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:07:49 2018

@author: asanofuyuki

.tiff
Normalized_by_each_mouse_max/ 1 min bin

出来上がったデータを.xlsで再保存して、doubleplotを描く。
データは9:00(ZT0)スタート、9:00(ZT0)終了でエクスポートすること。 (最終日の8:59までのデータがある。)
"""

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import sys
from PIL import Image
import time

t1 = time.time()

#好きな場所に好きな名前のフォルダを作り、その中にdouble_plot.pyと解析したい.xlsファイル1つをおく
path = glob.glob("./*.xls") #xlsファイルを見つける。
df = pd.read_excel(path[0], sheet_name = 0, header = None)
Ch = np.array(df.iloc[10, 1:])


df = df.drop(0, axis=1)
df = df.drop(range(0,11,1), axis=0)
df = df.reset_index(drop=True)


num = len(df.columns)
time_points = len(df.index)
recdays = int(len(df) / (24 * 60))


#マウス1匹ごとに処理していく
for ID in range(num): #エクセル上で一番左のマウスから処理していく    

    column = (np.array(df.iloc[:, ID])).astype(np.float64)
    c = np.argwhere(np.isnan(column))
    for d in np.reshape(c, (-1)):
        if d <= time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):(d + 11)]))
        elif d > time_points-10:
            column[d] = np.floor(np.nanmean(column[(d - 10):]))
#    df.iloc[:,i] = column
    

    maxi = np.max(column[7*1440 : (len(column)-7*1440)]) #LL中の最大の回転数を見つける
#    mean = np.mean(column[7*1440 : (len(column)-7*1440)]) #LL中の最大の回転数を見つける
    s = np.concatenate((np.zeros(60 * 24), column[:], np.zeros(60 * 48)), axis = 0) #データの前後に0を適切な数埋める
    t = np.reshape(s, (recdays + 3, 60 * 24)) #行にちょうど1日分のデータが並ぶ形に並び替え、
    u = np.concatenate((t[:recdays + 1, 60 * 18:], t[1:recdays + 2, :], t[2:recdays + 3, :60 * 18]), axis = 1) #tをずらしたり削ったりして1行分のデータが横に並ぶようにする
    
    v = u / (maxi*0.25) #最大値で割って、ratioにする
#    v = u / mean #最大値で割って、ratioにする
    v = np.round(v, 2)
    v = np.where(v>1, 1, v)

    
    ref = np.concatenate((np.zeros((50, 60*6+1))\
                          , np.ones((50, 60*12))\
                          , np.zeros((50, 60*12))\
                          , np.ones((50, 60*12))\
                          , np.zeros((50, 60*6+1))), axis=1)
    
    ref = np.concatenate((np.zeros((1, 60*24*2+2))\
                          , ref, np.zeros((1, 60*24*2+2))), axis=0)
    
    actogram = ref
    
    for i in range(recdays + 1):
        w = np.array(v[i,:])
        w = np.round(w*100/2).astype(np.int64)
        
        row = np.zeros((51, 1))
        
        for ii in range(len(w)):
            stac = np.concatenate((np.ones((50-w[ii], 1)), np.zeros((w[ii]+1, 1))), axis=0)
            row = np.concatenate((row, stac), axis=1)
            
        row = np.concatenate((row,  np.zeros((51, 1))), axis=1)
        actogram = np.concatenate((actogram, row), axis=0)
     
    name = str(Ch[ID]) + '_y50_LL.tiff'
    pil_img = Image.fromarray(actogram).save(name)

    
t2 = time.time()
print(t2-t1)
        
        
